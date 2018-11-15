#include "stdafx.h"
#include "gwmodel.h"
#include "iboundarycondition.h"
#include "element.h"

using namespace std;

void GWModel::update()
{
  if(m_currentDateTime < m_endDateTime)
  {
    applyBoundaryConditions(m_currentDateTime);

    if(m_component)
      m_component->applyInputValues();

    m_prevTimeStep = m_timeStep;

    m_timeStep = computeTimeStep();

    computeDerivedHydraulics();

    solveHydHead(m_timeStep);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(int j = 0; j < 2; j++)
    {
      switch (j)
      {
        case 0:
          {
            if(m_solveHeatTransport)
              solveHeatTransport(m_timeStep);
          }
          break;
        case 1:
          {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
            for(int i = 0 ; i < (int)m_solutes.size(); i++)
            {
              solveSoluteTransport(i, m_timeStep);
            }
          }
          break;
      }
    }

    m_prevDateTime = m_currentDateTime;
    m_currentDateTime = m_currentDateTime + m_timeStep / 86400.0;

    prepareForNextTimeStep();

    if(m_currentDateTime >= m_nextOutputTime)
    {
      writeOutput();
      m_nextOutputTime = std::min(m_nextOutputTime + m_outputInterval / 86400.0 , m_endDateTime);
    }

    if(m_verbose)
    {
      printStatus();
    }
  }
}

void GWModel::prepareForNextTimeStep()
{
  m_minTemp = m_minHead = std::numeric_limits<double>::max();
  m_maxTemp = m_maxHead = std::numeric_limits<double>::lowest();

  std::fill(m_maxSolute.begin(), m_maxSolute.end(), m_maxTemp);
  std::fill(m_minSolute.begin(), m_minSolute.end(), m_minTemp);

  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];

    for(size_t j = 0; j < element->elementCells.size(); j++)
    {
      ElementCell *elementCell = element->elementCells[j];

      elementCell->computeDVolumeDt();

      elementCell->computeMassBalance(m_timeStep);
      m_totalMassBalance = elementCell->totalMassBalance;

      elementCell->computeHeatBalance(m_timeStep);
      m_totalHeatBalance += elementCell->totalHeatBalance;

      m_totalExternalHeatBalance += elementCell->totalExternalHeatBalance;


      elementCell->prevHydHead.copy(elementCell->hydHead);
      elementCell->prevTemperature.copy(elementCell->temperature);

      m_minTemp = min(m_minTemp , elementCell->temperature.value);
      m_maxTemp = max(m_maxTemp , elementCell->temperature.value);

      m_minHead = min(m_minHead , elementCell->hydHead.value);
      m_maxHead = max(m_maxHead , elementCell->hydHead.value);

      for(size_t j = 0; j < m_solutes.size(); j++)
      {
        elementCell->computeSoluteBalance(m_timeStep, j);
        m_totalSoluteMassBalance[j] += elementCell->totalSoluteMassBalance[j];

        elementCell->prevSoluteConcs[j].copy(elementCell->soluteConcs[j]);

        m_minSolute[j] = min(m_minSolute[j] , elementCell->soluteConcs[j].value);
        m_maxSolute[j] = max(m_maxSolute[j] , elementCell->soluteConcs[j].value);
      }
    }
  }
}

void GWModel::applyInitialConditions()
{

  //Initialize heat and solute balance trackers
  m_totalMassBalance = 0.0;
  m_totalHeatBalance = 0.0;
  m_totalExternalHeatBalance = 0.0;

  std::fill(m_totalSoluteMassBalance.begin(), m_totalSoluteMassBalance.end(), 0.0);
  std::fill(m_totalExternalSoluteMassBalance.begin(), m_totalExternalSoluteMassBalance.end(), 0.0);

  applyBoundaryConditions(m_currentDateTime);

  //Write initial output
  writeOutput();

  //Set next output time
  m_nextOutputTime += m_outputInterval / 86400.0;

  for(size_t i = 0; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    ElementCell *leftCell = element->elementCells[99];
    ElementCell *rightCell = element->elementCells[100];
    leftCell->edgeHydHead[0].isBC = true;
    rightCell->edgeHydHead[2].isBC = true;

    double amplitude = 1.0;
    double phase = 0.0;
    double elapsed = (m_currentDateTime - m_startDateTime) * 86400.0;
    double head = amplitude * sin(0.000189804556154383 * elapsed + phase);
    leftCell->edgeHydHead[0].value = head;
    rightCell->edgeHydHead[2].value = head;
  }
}

void GWModel::applyBoundaryConditions(double dateTime)
{
  //reset external fluxes
#ifdef USE_OPENMMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];

    for(size_t j = 0; j < element->elementCells.size(); j++)
    {
      ElementCell *elementCell = element->elementCells[j];
      elementCell->externalInflow = 0;
      elementCell->externalHeatFluxes = 0;

      for(size_t j = 0; j < m_solutes.size(); j++)
      {
        elementCell->externalSoluteFluxes[j] = 0.0;
      }
    }
  }

#ifdef USE_OPENMMP
#pragma omp parallel for
#endif
  for(size_t i = 0; i < m_boundaryConditions.size() ; i++)
  {
    IBoundaryCondition *boundaryCondition = m_boundaryConditions[i];
    boundaryCondition->applyBoundaryConditions(dateTime);
  }


  for(size_t i = 0; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    ElementCell *leftCell = element->elementCells[99];
    ElementCell *rightCell = element->elementCells[100];
    leftCell->edgeHydHead[0].isBC = true;
    rightCell->edgeHydHead[2].isBC = true;


    double amplitude = 1.0;
    double phase = 0.0;
    double elapsed = (m_currentDateTime - m_startDateTime) * 86400.0;
    double head = amplitude * sin(0.000189804556154383 * elapsed + phase);
    leftCell->edgeHydHead[0].value = head;
    rightCell->edgeHydHead[2].value = head;
  }
}

double GWModel::computeTimeStep()
{
  double timeStep = m_maxTimeStep;
  double maxCourantFactor = 0.0;

  if(m_numCurrentInitFixedTimeSteps < m_numInitFixedTimeSteps)
  {
    timeStep = m_minTimeStep;
    m_numCurrentInitFixedTimeSteps++;
  }
  else if(m_useAdaptiveTimeStep)
  {
#ifdef _WIN32
    {

      for(int i = 0 ; i < (int)m_elements.size()  ; i++)
      {
        Element *element = m_elements[i];

        for(int j = 0; j < m_totalCellsPerElement; j++)
        {
          double courantFactor = element->elementCells[j]->computeCourantFactor();

          if(!std::isinf(courantFactor) && courantFactor > maxCourantFactor)
          {
            maxCourantFactor = courantFactor;
          }
        }
      }
    }
#else
    {
      for(int i = 0 ; i < (int)m_elements.size()  ; i++)
      {
        Element *element = m_elements[i];

        for(int j = 0; j < m_totalCellsPerElement; j++)
        {
          double courantFactor = element->elementCells[j]->computeCourantFactor();

          if(!(std::isinf(courantFactor) || std::isnan(courantFactor)) && courantFactor > maxCourantFactor)
          {
            maxCourantFactor = courantFactor;
          }
        }
      }
    }
#endif

    timeStep = maxCourantFactor ? m_timeStepRelaxationFactor / maxCourantFactor : m_maxTimeStep;
  }

  double nextTime = m_currentDateTime + timeStep / 86400.0;

  if(nextTime > m_nextOutputTime)
  {
    timeStep = std::max(m_minTimeStep,  (m_nextOutputTime - m_currentDateTime) *  86400.0);
  }

  timeStep = std::min(std::max(timeStep, m_minTimeStep), m_maxTimeStep);

  return timeStep;
}

void GWModel::computeDerivedHydraulics()
{
  for(int i = 0 ; i < (int)m_elements.size()  ; i++)
  {
    Element *element = m_elements[i];

    for(int j = 0; j < m_totalCellsPerElement; j++)
    {
      ElementCell *elementCell = element->elementCells[j];
      elementCell->computeDerivedHydraulics();
    }
  }
}

void GWModel::solveHydHead(double timeStep)
{
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0 ; i < (int)m_elements.size(); i++)
  {
    Element *element = m_elements[i];

    for(int j = 0; j < m_totalCellsPerElement; j++)
    {
      ElementCell *elementCell = element->elementCells[j];
      m_currHydHead[elementCell->index] = elementCell->hydHead.value;
      m_outHydHead[elementCell->index] = elementCell->hydHead.value;
    }
  }

  //Solve using ODE solver
  SolverUserData solverUserData; solverUserData.model = this; solverUserData.variableIndex = -2;

  if(m_hydHeadSolver->solve(m_currHydHead.data(), m_totalCellsPerElement * m_elements.size() , m_currentDateTime * 86400.0, timeStep,
                            m_outHydHead.data(), &GWModel::computeDHydHeadDt, &solverUserData))
  {
    printf("GWModel HydHead solver failed \n");
  }

  //Apply computed values;
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0 ; i < (int)m_elements.size(); i++)
  {
    Element *element = m_elements[i];

    for(int j = 0; j < m_totalCellsPerElement; j++)
    {
      ElementCell *elementCell = element->elementCells[j];
      double outHydHead = m_outHydHead[elementCell->index];
      elementCell->hydHead.value = outHydHead;
    }
  }
}

void GWModel::solveHeatTransport(double timeStep)
{
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0 ; i < (int)m_elements.size(); i++)
  {
    Element *element = m_elements[i];

    for(int j = 0; j < m_totalCellsPerElement; j++)
    {
      ElementCell *elementCell = element->elementCells[j];
      m_currTemps[elementCell->index] = elementCell->temperature.value;
      m_outTemps[elementCell->index] = elementCell->temperature.value;
    }
  }

  //Solve using ODE solver
  SolverUserData solverUserData; solverUserData.model = this; solverUserData.variableIndex = -1;

  if(m_heatSolver->solve(m_currTemps.data(), m_totalCellsPerElement * m_elements.size(), m_currentDateTime * 86400.0, timeStep,
                         m_outTemps.data(), &GWModel::computeDTDt, &solverUserData))
  {
    printf("GWModel heat solver failed \n");
  }

  //Apply computed values;
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0 ; i < (int)m_elements.size(); i++)
  {
    Element *element = m_elements[i];

    for(int j = 0; j < m_totalCellsPerElement; j++)
    {
      ElementCell *elementCell = element->elementCells[j];
      double outTemp = m_outTemps[elementCell->index];
      elementCell->temperature.value = outTemp;
    }
  }
}

void GWModel::solveSoluteTransport(int soluteIndex, double timeStep)
{
  std::vector<double> &currentSoluteConcs = m_currSoluteConcs[soluteIndex];
  std::vector<double> &outputSoluteConcs = m_outSoluteConcs[soluteIndex];

  //Set initial values.
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];

    for(int j = 0; j < m_totalCellsPerElement; j++)
    {
      ElementCell *elementCell = element->elementCells[j];
      currentSoluteConcs[elementCell->index] = elementCell->soluteConcs[soluteIndex].value;
      outputSoluteConcs[elementCell->index] = elementCell->soluteConcs[soluteIndex].value;
    }
  }

  //Solve using ODE solver
  SolverUserData solverUserData; solverUserData.model = this; solverUserData.variableIndex = soluteIndex;

  if(m_soluteSolvers[soluteIndex]->solve(outputSoluteConcs.data(), m_totalCellsPerElement * m_elements.size(), m_currentDateTime * 86400.0, timeStep,
                                         outputSoluteConcs.data(), &GWModel::computeDSoluteDt, &solverUserData))
  {
    printf("GWModel Solute solver failed \n");
  }

  //Apply computed values;
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(size_t i = 0 ; i < m_elements.size(); i++)
  {
    Element *element = m_elements[i];

    for(int j = 0; j < m_totalCellsPerElement; j++)
    {
      ElementCell *elementCell = element->elementCells[j];
      elementCell->soluteConcs[soluteIndex].value = outputSoluteConcs[elementCell->index];
    }
  }
}

void GWModel::computeDHydHeadDt(double t, double y[], double dydt[], void *userData)
{
  SolverUserData *solverUserData = (SolverUserData*) userData;
  GWModel *modelInstance = solverUserData->model;
  double dt = t - (modelInstance->m_currentDateTime *  86400.0);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0 ; i < (int)modelInstance->m_elements.size(); i++)
  {
    Element *element = modelInstance->m_elements[i];

    for(int j = 0 ; j < modelInstance->m_totalCellsPerElement; j++)
    {
      ElementCell *elementCell = element->elementCells[j];
      double DHeadDt = elementCell->computeDHydHeadDt(dt,y);
      dydt[elementCell->index] = DHeadDt;
    }
  }
}

void GWModel::computeDTDt(double t, double y[], double dydt[], void *userData)
{
  SolverUserData *solverUserData = (SolverUserData*) userData;
  GWModel *modelInstance = solverUserData->model;
  double dt = t - (modelInstance->m_currentDateTime *  86400.0);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0 ; i < (int)modelInstance->m_elements.size(); i++)
  {
    Element *element = modelInstance->m_elements[i];

    for(int j = 0 ; j < modelInstance->m_totalCellsPerElement; j++)
    {
      ElementCell  *elementCell = element->elementCells[j];
      double DTDt = elementCell->computeDHydHeadDt(dt,y);
      dydt[elementCell->index] = DTDt;
    }
  }
}

void GWModel::computeDSoluteDt(double t, double y[], double dydt[], void *userData)
{
  SolverUserData *solverUserData = (SolverUserData*) userData;
  GWModel *modelInstance = solverUserData->model;
  double dt = t - (modelInstance->m_currentDateTime *  86400.0);

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0 ; i < (int)modelInstance->m_elements.size(); i++)
  {
    Element *element = modelInstance->m_elements[i];

    for(int j = 0 ; j < modelInstance->m_totalCellsPerElement; j++)
    {
      ElementCell  *elementCell = element->elementCells[j];
      dydt[elementCell->index] = elementCell->computeDSoluteDt(dt,y,solverUserData->variableIndex);
    }
  }
}

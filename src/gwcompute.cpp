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

    calculatePreComputedHydHeadVariables();
    calculatePreComputedTempVariables();
    calculatePreComputedSoluteVariables();

    solve(m_timeStep);

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

      elementCell->computeMassBalance(m_timeStep);
      m_totalMassBalance += elementCell->totalMassBalance;

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

  calculatePreComputedHydHeadVariables();

  calculatePreComputedTempVariables();

  calculatePreComputedSoluteVariables();

  //Write initial output
  writeOutput();

  //Set next output time
  m_nextOutputTime += m_outputInterval / 86400.0;

}

void GWModel::applyBoundaryConditions(double dateTime)
{
  //reset external fluxes

  if(m_simulateWaterAge)
  {
#ifdef USE_OPENMMP
#pragma omp parallel for
#endif
    for(int i = 0 ; i < (int)m_elements.size(); i++)
    {
      Element *element = m_elements[i];

#ifdef USE_OPENMMP
#pragma omp parallel for
#endif
      for(int j = 0; j < (int)element->elementCells.size(); j++)
      {
        ElementCell *elementCell = element->elementCells[j];
        elementCell->externalInflow = 0.0;
        elementCell->externalHeatFluxes = 0.0;
        elementCell->externalSoluteFluxes[m_numSolutes] = elementCell->volume / 86400.0;

        for(int k = 0; k < m_numSolutes; k++)
        {
          elementCell->externalSoluteFluxes[k] = 0.0;
        }
      }
    }
  }
  else
  {
#ifdef USE_OPENMMP
#pragma omp parallel for
#endif
    for(int i = 0 ; i < (int)m_elements.size(); i++)
    {
      Element *element = m_elements[i];

#ifdef USE_OPENMMP
#pragma omp parallel for
#endif
      for(int j = 0; j < (int)element->elementCells.size(); j++)
      {
        ElementCell *elementCell = element->elementCells[j];
        elementCell->externalInflow = 0.0;
        elementCell->externalHeatFluxes = 0.0;

        for(int k = 0; k < m_numSolutes; k++)
        {
          elementCell->externalSoluteFluxes[k] = 0.0;
        }
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
          ElementCell *elementCell = element->elementCells[j];
          double courantFactor = max(elementCell->computeCourantFactor(), elementCell->computeDiffusionFactor());

          if(/*!(std::isinf(courantFactor) || std::isnan(courantFactor)) &&*/ courantFactor > maxCourantFactor)
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

void GWModel::calculatePreComputedHydHeadVariables()
{
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0 ; i < (int)m_elements.size()  ; i++)
  {
    Element *element = m_elements[i];
    element->channelInflow = 0.0;
    element->channelInflowFlux = 0.0;

    for(int j = 0; j < m_totalCellsPerElement; j++)
    {
      ElementCell *elementCell = element->elementCells[j];
      elementCell->calculatePreComputedHydHeadVariables();
    }
  }
}

void GWModel::calculatePreComputedTempVariables()
{
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0 ; i < (int)m_elements.size()  ; i++)
  {
    Element *element = m_elements[i];
    element->channelHeatFlux = 0.0;
    element->channelHeatRate = 0.0;

    for(int j = 0; j < m_totalCellsPerElement; j++)
    {
      ElementCell *elementCell = element->elementCells[j];
      elementCell->calculatePreComputedTempVariables();
    }
  }
}

void GWModel::calculatePreComputedSoluteVariables()
{
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0 ; i < (int)m_elements.size()  ; i++)
  {
    Element *element = m_elements[i];

    for(size_t j = 0; j < m_solutes.size(); j++)
    {
      element->channelSoluteFlux[j] = 0.0;
      element->channelSoluteRate[j] = 0.0;
    }

    for(int j = 0; j < m_totalCellsPerElement; j++)
    {
      ElementCell *elementCell = element->elementCells[j];
      elementCell->calculatePreComputedSoluteVariables();
    }
  }
}

void GWModel::solve(double timeStep)
{
  if(m_solveHeatTransport)
  {

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(int i = 0 ; i < (int)m_elements.size(); i++)
    {
      Element *element = m_elements[i];

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      for(int j = 0; j < m_totalCellsPerElement; j++)
      {
        ElementCell *elementCell = element->elementCells[j];
        m_solverCurrentValues[elementCell->index] = elementCell->hydHead.value;
        m_solverOutValues[elementCell->index] = elementCell->hydHead.value;
        m_solverCurrentValues[elementCell->index + m_tempIndex] = elementCell->temperature.value;
        m_solverOutValues[elementCell->index + m_tempIndex] = elementCell->temperature.value;

        for(size_t s = 0 ; s < m_solutes.size(); s++)
        {
          m_solverCurrentValues[elementCell->index + m_soluteIndexes[s]] = elementCell->soluteConcs[s].value;
          m_solverOutValues[elementCell->index + m_soluteIndexes[s]] = elementCell->soluteConcs[s].value;
        }
      }
    }

    //Solve using ODE solver
    SolverUserData solverUserData; solverUserData.model = this;

    if(m_odeSolver->solve(m_solverCurrentValues.data(), m_solverCurrentValues.size() , 0, timeStep,
                          m_solverOutValues.data(), &GWModel::computeDYDt, &solverUserData))
    {
      printf("Solver failed \n");
    }

    //Apply computed values;
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(int i = 0 ; i < (int)m_elements.size(); i++)
    {
      Element *element = m_elements[i];

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      for(int j = 0; j < m_totalCellsPerElement; j++)
      {
        ElementCell *elementCell = element->elementCells[j];

        elementCell->hydHead.value = m_solverOutValues[elementCell->index];
        elementCell->temperature.value = m_solverOutValues[elementCell->index + m_tempIndex];

        for(size_t s = 0 ; s < m_solutes.size(); s++)
        {
          elementCell->soluteConcs[s].value = m_solverOutValues[elementCell->index + m_soluteIndexes[s]];
        }

        elementCell->computeVolumeDerivative();
      }
    }
  }
  else
  {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(int i = 0 ; i < (int)m_elements.size(); i++)
    {
      Element *element = m_elements[i];

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      for(int j = 0; j < m_totalCellsPerElement; j++)
      {
        ElementCell *elementCell = element->elementCells[j];
        m_solverCurrentValues[elementCell->index] = elementCell->hydHead.value;
        m_solverOutValues[elementCell->index] = elementCell->hydHead.value;

        for(size_t s = 0 ; s < m_solutes.size(); s++)
        {
          m_solverCurrentValues[elementCell->index + m_soluteIndexes[s]] = elementCell->soluteConcs[s].value;
          m_solverOutValues[elementCell->index + m_soluteIndexes[s]] = elementCell->soluteConcs[s].value;
        }
      }
    }

    //Solve using ODE solver
    SolverUserData solverUserData; solverUserData.model = this;

    if(m_odeSolver->solve(m_solverCurrentValues.data(), m_solverCurrentValues.size() , m_currentDateTime * 86400.0, timeStep,
                          m_solverOutValues.data(), &GWModel::computeDYDt, &solverUserData))
    {
      printf("Solver failed \n");
    }

    //Apply computed values;
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(int i = 0 ; i < (int)m_elements.size(); i++)
    {
      Element *element = m_elements[i];

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      for(int j = 0; j < m_totalCellsPerElement; j++)
      {
        ElementCell *elementCell = element->elementCells[j];

        elementCell->hydHead.value = m_solverOutValues[elementCell->index];
        elementCell->temperature.value = m_solverOutValues[elementCell->index + m_tempIndex];

        for(size_t s = 0 ; s < m_solutes.size(); s++)
        {
          elementCell->soluteConcs[s].value = m_solverOutValues[elementCell->index + m_soluteIndexes[s]];
        }

        elementCell->computeVolumeDerivative();
      }
    }
  }
}

void GWModel::computeDYDt(double t, double y[], double dydt[], void *userData)
{
  SolverUserData *solverUserData = (SolverUserData*) userData;
  GWModel *modelInstance = solverUserData->model;

  if(modelInstance->m_solveHeatTransport)
  {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(int i = 0 ; i < (int)modelInstance->m_elements.size(); i++)
    {
      Element *element = modelInstance->m_elements[i];

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      for(int j = 0 ; j < modelInstance->m_totalCellsPerElement; j++)
      {
        ElementCell  *elementCell = element->elementCells[j];
        dydt[elementCell->index] = elementCell->computeDHydHeadDt(t, y);
      }

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      for(int j = 0 ; j < modelInstance->m_totalCellsPerElement; j++)
      {
        ElementCell  *elementCell = element->elementCells[j];
        dydt[elementCell->index + modelInstance->m_tempIndex] = elementCell->computeDTDt(t, &y[modelInstance->m_tempIndex]);

        for(size_t s = 0 ; s < modelInstance->m_solutes.size(); s++)
        {
          dydt[elementCell->index + modelInstance->m_soluteIndexes[s]] = elementCell->computeDSoluteDt(t, &y[modelInstance->m_soluteIndexes[s]], s);
        }
      }
    }
  }
  else
  {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for(int i = 0 ; i < (int)modelInstance->m_elements.size(); i++)
    {
      Element *element = modelInstance->m_elements[i];

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      for(int j = 0 ; j < modelInstance->m_totalCellsPerElement; j++)
      {
        ElementCell  *elementCell = element->elementCells[j];
        dydt[elementCell->index] = elementCell->computeDHydHeadDt(t, y);
      }

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      for(int j = 0 ; j < modelInstance->m_totalCellsPerElement; j++)
      {
        ElementCell  *elementCell = element->elementCells[j];
        for(size_t s = 0 ; s < modelInstance->m_solutes.size(); s++)
        {
          dydt[elementCell->index + modelInstance->m_soluteIndexes[s]] = elementCell->computeDSoluteDt(t, &y[modelInstance->m_soluteIndexes[s]], s);
        }
      }
    }
  }
}

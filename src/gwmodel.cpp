/*!
*  \file    gwmodel.cpp
*  \author  Caleb Amoa Buahin <caleb.buahin@gmail.com>
*  \version 1.0.0
*  \section Description
*  This file and its associated files and libraries are free software;
*  you can redistribute it and/or modify it under the terms of the
*  Lesser GNU Lesser General Public License as published by the Free Software Foundation;
*  either version 3 of the License, or (at your option) any later version.
*  fvhmcompopnent.h its associated files is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.(see <http://www.gnu.org/licenses/> for details)
*  \date 2018
*  \pre
*  \bug
*  \todo
*  \warning
*/

#include "gwmodel.h"
#include "elementjunction.h"
#include "element.h"
#include "iboundarycondition.h"

using namespace  std;

GWModel::GWModel(GWComponent *component)
  : QObject(component),
    m_timeStep(0.01),
    m_maxTimeStep(10),
    m_minTimeStep(0.01),
    m_outputInterval(14400),
    m_timeStepRelaxationFactor(0.8),
    m_numInitFixedTimeSteps(1),
    m_numCurrentInitFixedTimeSteps(0),
    m_printFrequency(10),
    m_currentPrintCount(0),
    m_flushToDiskFrequency(10),
    m_currentflushToDiskCount(0),
    m_addedSoluteCount(0),
    m_numLeftCellsPerElement(0),
    m_numRightCellsPerElement(0),
    m_totalCellsPerElement(0),
    m_useAdaptiveTimeStep(true),
    m_verbose(true),
    m_solveHeatTransport(false),
    m_hydHeadSolver(nullptr),
    m_heatSolver(nullptr),
    m_waterDensity(1000.0),
    m_cp(4184),
    m_sedDensity(1970),
    m_sedCp(2758),
    m_hydConX(0.0),
    m_hydConY(0.0),
    m_specificStorage(0.3),
    m_porosity(0.3),
    m_defaultCellWidth(1.0),
    #ifdef   USE_NETCDF
    m_outputNetCDF(nullptr),
    #endif
    m_retrieveCouplingDataFunction(nullptr),
    m_component(component)
{
  m_hydHeadSolver = new ODESolver(1, ODESolver::CVODE_ADAMS);
  m_heatSolver = new ODESolver(1, ODESolver::CVODE_ADAMS);
}

GWModel::~GWModel()
{
  for(Element *element : m_elements)
    delete element;

  m_elements.clear();
  m_elementsById.clear();


  for(ElementJunction *elementJunction : m_elementJunctions)
    delete elementJunction;

  m_elementJunctions.clear();
  m_elementJunctionsById.clear();


  delete m_hydHeadSolver;
  delete m_heatSolver;


  for(size_t i = 0; i < m_soluteSolvers.size(); i++)
  {
    delete m_soluteSolvers[i];
  }

  m_soluteSolvers.clear();

  closeOutputFiles();

  for(IBoundaryCondition *boundaryCondition : m_boundaryConditions)
    delete boundaryCondition;

  m_boundaryConditions.clear();
}

double GWModel::minTimeStep() const
{
  return m_minTimeStep;
}

void GWModel::setMinTimeStep(double timeStep)
{
  m_minTimeStep = timeStep;
}

double GWModel::maxTimeStep() const
{
  return m_maxTimeStep;
}

void GWModel::setMaxTimeStep(double timeStep)
{
  m_maxTimeStep = timeStep;
}

bool GWModel::useAdaptiveTimeStep() const
{
  return m_useAdaptiveTimeStep;
}

void GWModel::setUseAdaptiveTimeStep(bool use)
{
  m_useAdaptiveTimeStep = use;
}

double GWModel::timeStepRelaxationFactor() const
{
  return m_timeStepRelaxationFactor;
}

void GWModel::setTimeStepRelaxationFactor(double tStepRelaxFactor)
{
  m_timeStepRelaxationFactor = tStepRelaxFactor;
}

double GWModel::currentTimeStep() const
{
  return m_timeStep;
}

double GWModel::startDateTime() const
{
  return m_startDateTime;
}

void GWModel::setStartDateTime(double dateTime)
{
  m_startDateTime = dateTime;
}

double GWModel::endDateTime() const
{
  return m_endDateTime;
}

void GWModel::setEndDateTime(double dateTime)
{
  m_endDateTime = dateTime;
}

double GWModel::outputInterval() const
{
  return m_outputInterval;
}

void GWModel::setOutputInterval(double interval)
{
  m_outputInterval = interval;
}

double GWModel::currentDateTime() const
{
  return m_currentDateTime;
}

ODESolver *GWModel::hydHeadSolver() const
{
  return m_hydHeadSolver;
}

ODESolver *GWModel::heatSolver() const
{
  return m_heatSolver;
}

std::vector<ODESolver*> GWModel::soluteSolvers() const
{
  return m_soluteSolvers;
}

double GWModel::waterDensity() const
{
  return m_waterDensity;
}

void GWModel::setWaterDensity(double value)
{
  m_waterDensity = value;
}

double GWModel::defaultSedimentDensity() const
{
  return m_sedDensity;
}

void GWModel::setDefaultSedimentDensity(double value)
{
  m_sedDensity = value;
}

double GWModel::specificHeatCapacityWater() const
{
  return m_cp;
}

void GWModel::setSpecificHeatCapacityWater(double value)
{
  m_cp = value;
}

double GWModel::defaultSedimentSpecificHeatCapacity() const
{
  return m_sedCp;
}

void GWModel::setDefaultSedimentSpecificHeatCapacity(double value)
{
  m_sedCp = value;
}

double GWModel::defaultPorosity() const
{
  return m_porosity;
}

void GWModel::setDefaultPorosity(double porosity)
{
  m_porosity = porosity;
}

double GWModel::defaultHydraulicConductivityX() const
{
  return m_hydConX;
}

void GWModel::setDefaultHydraulicConductivityX(double hydConX)
{
  m_hydConX = hydConX;
}

double GWModel::defaultHydraulicConductivityY() const
{
  return m_hydConY;
}

void GWModel::setDefaultHydraulicConductivityY(double hydConY)
{
  m_hydConY = hydConY;
}

int GWModel::numSolutes() const
{
  return (int) m_solutes.size();
}

void GWModel::setNumSolutes(int numSolutes)
{
  for(size_t i = 0; i < m_soluteSolvers.size(); i++)
  {
    delete m_soluteSolvers[i];
  }

  m_soluteSolvers.clear();

  if(numSolutes >= 0)
  {
    m_solutes.resize(numSolutes);
    m_maxSolute.resize(numSolutes, 0.0);
    m_minSolute.resize(numSolutes, 0.0);
    m_totalSoluteMassBalance.resize(numSolutes, 0.0);
    m_totalExternalSoluteMassBalance.resize(numSolutes, 0.0);

    for(size_t i = 0 ; i < m_solutes.size(); i++)
    {
      m_soluteSolvers.push_back(new ODESolver(m_elements.size(), ODESolver::CVODE_ADAMS));
      m_solutes[i] = "Solute_" + std::to_string(i + 1);
    }
  }
}

void GWModel::setSoluteName(int soluteIndex, const std::string &soluteName)
{
  m_solutes[soluteIndex] = soluteName;
}

std::string GWModel::solute(int soluteIndex) const
{
  return m_solutes[soluteIndex];
}

int GWModel::numElementJunctions() const
{
  return (int)m_elements.size();
}

ElementJunction *GWModel::addElementJunction(const std::string &id, double x, double y, double z)
{
  if(m_elementJunctionsById.find(id) == m_elementJunctionsById.end())
  {
    ElementJunction *eJunction = new ElementJunction(id, x, y, z, this);
    eJunction->index = m_elementJunctions.size();
    m_elementJunctions.push_back(eJunction);
    m_elementJunctionsById[id] = eJunction;
    return eJunction;
  }

  return nullptr;
}

void GWModel::deleteElementJunction(const std::string &id)
{
  std::unordered_map<string,ElementJunction*>::iterator eJIter =  m_elementJunctionsById.find(id) ;

  if(eJIter != m_elementJunctionsById.end())
  {
    ElementJunction *eJunction = eJIter->second;
    m_elementJunctionsById.erase(eJIter);

    std::vector<ElementJunction*>::iterator it = std::find(m_elementJunctions.begin(), m_elementJunctions.end(), eJunction);
    if(it != m_elementJunctions.end())
    {
      m_elementJunctions.erase(it);
    }

    delete eJunction;
  }
}

void GWModel::deleteElementJunction(int index)
{
  ElementJunction *eJunction = m_elementJunctions[index];

  m_elementJunctionsById.erase(eJunction->id);

  std::vector<ElementJunction*>::iterator it = std::find(m_elementJunctions.begin(), m_elementJunctions.end(), eJunction);
  if(it != m_elementJunctions.end())
    m_elementJunctions.erase(it);

  delete eJunction;
}

ElementJunction *GWModel::getElementJunction(const std::string &id)
{
  return m_elementJunctionsById[id];
}

ElementJunction *GWModel::getElementJunction(int index)
{
  return m_elementJunctions[index];
}

int GWModel::numLeftElementCells() const
{
  return m_numLeftCellsPerElement;
}

void GWModel::setNumLeftElementCells(int value)
{
  m_numLeftCellsPerElement = value;
}

int GWModel::numRightElementCells() const
{
  return m_numRightCellsPerElement;
}

void GWModel::setNumRightElementCells(int value)
{
  m_numRightCellsPerElement = value;
}

int GWModel::numElements() const
{
  return (int)m_elements.size();
}

Element *GWModel::addElement(const std::string &id, ElementJunction *upStream, ElementJunction *downStream)
{
  if(upStream && downStream)
  {
    Element *element = new Element(id, upStream, downStream, this);
    element->index = m_elements.size();
    m_elements.push_back(element);
    m_elementsById[id] = element;
    return element;
  }

  return nullptr;
}

void GWModel::deleteElement(const std::string &id)
{
  unordered_map<string,Element*>::iterator eIter = m_elementsById.find(id);

  if(eIter != m_elementsById.end())
  {
    Element *element = eIter->second;
    m_elementsById.erase(eIter);

    vector<Element*>::iterator it = std::find(m_elements.begin() , m_elements.end(), element);
    if(it != m_elements.end())
      m_elements.erase(it);

    delete element;
  }
}

void GWModel::deleteElement(int index)
{
  Element *element = m_elements[index];
  m_elementJunctionsById.erase(element->id);

  vector<Element*>::iterator it = std::find(m_elements.begin() , m_elements.end(), element);

  if(it != m_elements.end())
    m_elements.erase(it);

  delete element;
}

Element *GWModel::getElement(const std::string &id)
{
  return m_elementsById[id];
}

Element *GWModel::getElement(int index)
{
  return m_elements[index];
}

RetrieveCouplingData GWModel::retrieveCouplingDataFunction() const
{

}

void GWModel::setRetrieveCouplingDataFunction(RetrieveCouplingData retrieveCouplingDataFunction)
{
  m_retrieveCouplingDataFunction = retrieveCouplingDataFunction;
}

bool GWModel::initialize(std::list<std::string> &errors)
{
  bool initialized = initializeInputFiles(errors) &&
                     initializeTimeVariables(errors) &&
                     initializeElements(errors) &&
                     initializeSolvers(errors) &&
                     initializeBoundaryConditions(errors) &&
                     initializeOutputFiles(errors);


  if(initialized)
  {
    applyInitialConditions();
  }

  return initialized;
}

bool GWModel::finalize(std::list<std::string> &errors)
{
  closeOutputFiles();

  for(IBoundaryCondition *boundaryCondition : m_boundaryConditions)
    delete boundaryCondition;

  m_boundaryConditions.clear();

  return true;
}

bool GWModel::initializeTimeVariables(std::list<string> &errors)
{
  if(m_startDateTime >= m_endDateTime)
  {
    errors.push_back("End datetime must be greater than startdatetime");
    return false;
  }

  if( (m_endDateTime - m_startDateTime) *  86400.0 < m_minTimeStep )
  {
    errors.push_back("Make sure timestep is less than the simulation interval");
    return false;
  }

  if(m_minTimeStep <=  0 || m_maxTimeStep <= 0)
  {
    errors.push_back("Make sure time steps are greater 0");
    return false;
  }

  if(m_minTimeStep > m_maxTimeStep)
  {
    errors.push_back("");
    return false;
  }

  m_numCurrentInitFixedTimeSteps = 0;

  m_currentDateTime = m_startDateTime;
  m_nextOutputTime = m_currentDateTime;

  m_currentPrintCount = 0;
  m_currentflushToDiskCount = 0;

  return true;
}

bool GWModel::initializeElements(std::list<string> &errors)
{

  m_totalCellsPerElement = m_numLeftCellsPerElement + m_numRightCellsPerElement;

  for(int i = 0 ; i < (int)m_elementJunctions.size()  ; i++)
  {
    ElementJunction *elementJunction = m_elementJunctions[i];
    elementJunction->index = i;
  }

  for(int i = 0 ; i < (int)m_elements.size(); i++)
  {
    Element *element = m_elements[i];
    element->index = i;
    double runningY = 0.0;

    for(int j = 0; j < m_totalCellsPerElement; j++)
    {
      ElementCell *elementCell = element->elementCells[j];
      elementCell->index  = i * m_totalCellsPerElement + j;
    }

    for(int f = m_numLeftCellsPerElement - 1; f > -1; f--)
    {
      runningY -= element->elementCells[f]->width;
      element->elementCells[f]->centerY  = runningY + element->elementCells[f]->width / 2.0;
    }

    runningY = 0.0;

    for(int f = m_numLeftCellsPerElement; f < m_totalCellsPerElement; f++)
    {
      runningY += element->elementCells[f]->width;
      element->elementCells[f]->centerY  = runningY - element->elementCells[f]->width / 2.0;
    }
  }

  for(int i = 0 ; i < (int)m_elements.size()  ; i++)
  {
    Element *element = m_elements[i];
    element->initialize();
  }

  m_currHydHead.resize(m_totalCellsPerElement * m_elements.size(), 0.0);
  m_outHydHead.resize(m_totalCellsPerElement * m_elements.size(), 0.0);

  m_currTemps.resize(m_totalCellsPerElement * m_elements.size(), 0.0);
  m_outTemps.resize(m_totalCellsPerElement * m_elements.size(), 0.0);

  m_currSoluteConcs.resize(m_solutes.size(), std::vector<double>(m_totalCellsPerElement * m_elements.size(), 0.0));
  m_outSoluteConcs.resize(m_solutes.size(), std::vector<double>(m_totalCellsPerElement * m_elements.size(), 0.0));


  return true;
}

bool GWModel::initializeSolvers(std::list<string> &errors)
{
  m_hydHeadSolver->setSize(m_elements.size() * m_totalCellsPerElement);
  m_hydHeadSolver->initialize();

  if(m_solveHeatTransport)
  {
    m_heatSolver->setSize(m_elements.size() * m_totalCellsPerElement);
    m_heatSolver->initialize();
  }

  for(size_t i = 0; i < m_soluteSolvers.size(); i++)
  {
    ODESolver *solver =  m_soluteSolvers[i];
    solver->setSize(m_elements.size() * m_totalCellsPerElement);
    solver->initialize();
  }

  return true;
}

bool GWModel::initializeBoundaryConditions(std::list<string> &errors)
{
  for(size_t i = 0; i < m_boundaryConditions.size() ; i++)
  {
    IBoundaryCondition *boundaryCondition = m_boundaryConditions[i];
    boundaryCondition->clear();
    boundaryCondition->findAssociatedGeometries();
    boundaryCondition->prepare();
  }

  return true;
}

bool GWModel::findProfile(Element *from, Element *to, std::list<Element *> &profile)
{
  for(Element *outgoing : from->downstreamJunction->outgoingElements)
  {
    if(outgoing == to)
    {
      profile.push_back(from);
      profile.push_back(outgoing);
      return true;
    }
    else if(findProfile(outgoing, to, profile))
    {
      profile.push_front(from);
      return true;
    }
  }

  return false;
}

void GWModel::calculateDistanceFromUpstreamJunction(Element *element)
{
  //  if(element->distanceFromUpStreamJunction == 0)
  //  {
  //    if(element->upstreamElement != nullptr)
  //    {
  //      if(element->upstreamElement->distanceFromUpStreamJunction == 0)
  //      {
  //        calculateDistanceFromUpstreamJunction(element->upstreamElement);
  //      }

  //      element->distanceFromUpStreamJunction = element->upstreamElement->distanceFromUpStreamJunction + element->length / 2.0;
  //    }
  //    else
  //    {
  //      element->distanceFromUpStreamJunction = element->length / 2.0;
  //    }
  //  }
}








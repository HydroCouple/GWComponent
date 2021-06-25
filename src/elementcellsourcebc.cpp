#include "stdafx.h"
#include "elementcellsourcebc.h"
#include "gwmodel.h"
#include "temporal/timeseries.h"
#include "element.h"
#include "elementcell.h"
#include "core/datacursor.h"

ElementCellSourceBC::ElementCellSourceBC(Variable variable,
                                         GeometryMultiplier geometryMultiplier,
                                         Element *startElement,
                                         Element *endElement,
                                         int startElementCell,
                                         int endElementCell,
                                         GWModel *model):
  QObject(model),
  m_variable(variable),
  m_geometryMultiplier(geometryMultiplier),
  m_startElement(startElement),
  m_endElement(endElement),
  m_startElementCell(startElementCell),
  m_endElementCell(endElementCell),
  m_timeSeries(nullptr),
  m_applyBCFunction(nullptr),
  m_parentModel(model)
{
  m_dataCursor = new DataCursor();
}

ElementCellSourceBC::~ElementCellSourceBC()
{
  delete m_dataCursor;
}

void ElementCellSourceBC::findAssociatedGeometries()
{
  m_profile.clear();
  m_parentModel->findProfile(m_startElement, m_endElement, m_profile);
}

void ElementCellSourceBC::prepare()
{
  m_match = false;
  m_cellLength = std::max(0, m_endElementCell - m_startElementCell + 1);


  switch (m_geometryMultiplier)
  {
    case Area:
      m_getGeometryMultiplier =  &ElementCellSourceBC::getArea;
      break;
    case Volume:
      m_getGeometryMultiplier =  &ElementCellSourceBC::getVolume;
      break;
    case Length:
      m_getGeometryMultiplier =  &ElementCellSourceBC::getLength;
      break;
    case Width:
      m_getGeometryMultiplier =  &ElementCellSourceBC::getWidth;
      break;
    case SaturatedThickness:
      m_getGeometryMultiplier =  &ElementCellSourceBC::getSaturatedThickness;
      break;
    default:
      m_getGeometryMultiplier =  &ElementCellSourceBC::getZero;
      break;
  }


  switch (m_variable)
  {
    case FLOW:
      {
        m_applyBCFunction = &ElementCellSourceBC::applyFlowBC;
      }
      break;
    case HEAT:
      {
        m_applyBCFunction = &ElementCellSourceBC::applyHeatBC;
      }
      break;
    case SOLUTE:
      {
        m_applyBCFunction = &ElementCellSourceBC::applySoluteBC;
      }
      break;
  }

  if(m_timeSeries->numColumns() == m_cellLength * m_profile.size())
  {
    m_match = true;
  }
}

void ElementCellSourceBC::applyBoundaryConditions(double dateTime)
{
  (this->*m_applyBCFunction)(dateTime);
}

void ElementCellSourceBC::clear()
{
  m_profile.clear();

}

ElementCellSourceBC::Variable ElementCellSourceBC::variable() const
{
  return m_variable;
}

void ElementCellSourceBC::setVariable(Variable variable)
{
  m_variable = variable;
}

ElementCellSourceBC::GeometryMultiplier ElementCellSourceBC::geometryMultiplier() const
{
  return m_geometryMultiplier;
}

void ElementCellSourceBC::setGeometryMultiplier(GeometryMultiplier multiplier)
{
  m_geometryMultiplier = multiplier;
}

int ElementCellSourceBC::soluteIndex() const
{
  return m_soluteIndex;
}

void ElementCellSourceBC::setSoluteIndex(int soluteIndex)
{
  m_soluteIndex = soluteIndex;
}

Element *ElementCellSourceBC::startElement() const
{
  return m_startElement;
}

void ElementCellSourceBC::setStartElement(Element *element)
{
  m_startElement = element;
}

Element *ElementCellSourceBC::endElement() const
{
  return m_endElement;
}

void ElementCellSourceBC::setEndElement(Element *element)
{
  m_endElement = element;
}

std::vector<Element*> ElementCellSourceBC::profiles() const
{
  return m_profile;
}

QSharedPointer<TimeSeries> ElementCellSourceBC::timeSeries() const
{
  return m_timeSeries;
}

void ElementCellSourceBC::setTimeSeries(const QSharedPointer<TimeSeries> &timeseries)
{
  m_timeSeries = timeseries;
  m_dataCursor->setMin(0);
  m_dataCursor->setMax(timeseries->numRows() - 1);
}

void ElementCellSourceBC::applyFlowBC(double dateTime)
{
  double value = 0;

  if(m_match)
  {
    for(size_t j = 0; j < m_profile.size(); j++)
    {
      for(int i = 0; i < m_cellLength; i++)
      {
        ElementCell *elementCell = m_profile[j]->elementCells[m_startElementCell + i];

        if(m_timeSeries->interpolate(dateTime, j, m_dataCursor, value))
        {
          value = value * (this->*m_getGeometryMultiplier)(elementCell);
          elementCell->externalInflow += value;

//          if(value < 0.0)
          {
            elementCell->externalHeatFluxes += m_parentModel->m_waterDensity * m_parentModel->m_cp * value * elementCell->temperature.value;
            elementCell->externalSoluteFluxes[m_soluteIndex] += value * elementCell->soluteConcs[m_soluteIndex].value;
          }
        }
      }
    }
  }
  else
  {
    if(m_timeSeries->interpolate(dateTime, 0, m_dataCursor, value))
    {
      for(size_t j = 0; j < m_profile.size(); j++)
      {
        for(int i = 0; i < m_cellLength; i++)
        {
          ElementCell *elementCell = m_profile[j]->elementCells[m_startElementCell + i];

          double tempV = value * (this->*m_getGeometryMultiplier)(elementCell);
          elementCell->externalInflow += tempV;

//          if(value < 0.0)
          {
            elementCell->externalHeatFluxes += m_parentModel->m_waterDensity * m_parentModel->m_cp * tempV * elementCell->temperature.value;
            elementCell->externalSoluteFluxes[m_soluteIndex] += tempV * elementCell->soluteConcs[m_soluteIndex].value;
          }
        }
      }
    }
  }
}

void ElementCellSourceBC::applyHeatBC(double dateTime)
{
  double value = 0;

  if(m_match)
  {
    for(size_t j = 0; j < m_profile.size(); j++)
    {
      for(int i = 0; i < m_cellLength; i++)
      {
        ElementCell *elementCell = m_profile[j]->elementCells[m_startElementCell + i];

        if(m_timeSeries->interpolate(dateTime, j, m_dataCursor, value))
        {
          value = value * (this->*m_getGeometryMultiplier)(elementCell);
          elementCell->externalHeatFluxes += value;
        }
      }
    }
  }
  else
  {
    if(m_timeSeries->interpolate(dateTime, 0, m_dataCursor, value))
    {
      for(size_t j = 0; j < m_profile.size(); j++)
      {
        for(int i = 0; i < m_cellLength; i++)
        {
          ElementCell *elementCell = m_profile[j]->elementCells[m_startElementCell + i];

          double tempV = value * (this->*m_getGeometryMultiplier)(elementCell);
          elementCell->externalHeatFluxes += tempV;
        }
      }
    }
  }
}

void ElementCellSourceBC::applySoluteBC(double dateTime)
{

  double value = 0;

  if(m_match)
  {
    for(size_t j = 0; j < m_profile.size(); j++)
    {
      for(int i = 0; i < m_cellLength; i++)
      {
        ElementCell *elementCell = m_profile[j]->elementCells[m_startElementCell + i];

        if(m_timeSeries->interpolate(dateTime, j, m_dataCursor, value))
        {
          value = value * (this->*m_getGeometryMultiplier)(elementCell);
          elementCell->externalSoluteFluxes[m_soluteIndex] += value;
        }
      }
    }
  }
  else
  {
    if(m_timeSeries->interpolate(dateTime, 0, m_dataCursor, value))
    {
      for(size_t j = 0; j < m_profile.size(); j++)
      {
        for(int i = 0; i < m_cellLength; i++)
        {
          ElementCell *elementCell = m_profile[j]->elementCells[m_startElementCell + i];
          double tempV = value * (this->*m_getGeometryMultiplier)(elementCell);
          elementCell->externalSoluteFluxes[m_soluteIndex] += tempV;
        }
      }
    }
  }
}

double ElementCellSourceBC::getZero(ElementCell *elementCell)
{
  return 0.0;
}

double ElementCellSourceBC::getArea(ElementCell *elementCell)
{
  return elementCell->parentElement->length * elementCell->width;
}

double ElementCellSourceBC::getVolume(ElementCell *elementCell)
{
  return elementCell->volume;
}

double ElementCellSourceBC::getLength(ElementCell *elementCell)
{
  return elementCell->parentElement->length;
}

double ElementCellSourceBC::getWidth(ElementCell *elementCell)
{
  return elementCell->width;
}

double ElementCellSourceBC::getSaturatedThickness(ElementCell *elementCell)
{
  return elementCell->depth;
}

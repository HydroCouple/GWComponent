#include "stdafx.h"
#include "edgebc.h"
#include "gwmodel.h"
#include "temporal/timeseries.h"
#include "element.h"
#include "core/datacursor.h"

EdgeBC::EdgeBC(Variable variable,
               Edge edge,
               Element *startElement,
               Element *endElement,
               int startElementCell,
               int endElementCell,
               GWModel *model):
  QObject(model),
  m_variable(variable),
  m_edge(edge),
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

EdgeBC::~EdgeBC()
{
  delete m_dataCursor;
}

void EdgeBC::findAssociatedGeometries()
{
  m_profile.clear();
  m_neighbourCells.clear();
  m_parentModel->findProfile(m_startElement, m_endElement, m_profile);

  switch (m_edge)
  {
    case RIGHT:
      {
        m_nEdge = LEFT;
      }
      break;
    case LEFT:
      {
        m_nEdge = RIGHT;
      }
      break;
    case UP:
      {
        m_nEdge = DOWN;
      }
      break;
    case DOWN:
      {
        m_nEdge = UP;
      }
      break;
  }
}

void EdgeBC::prepare()
{
  m_match = false;
  m_cellLength = std::max(0, m_endElementCell - m_startElementCell + 1);

  switch (m_variable)
  {
    case HYD_HEAD:
      {
        m_applyBCFunction = &EdgeBC::applyHydHeadBC;

        for(size_t i = 0 ; i < m_profile.size(); i++)
        {
          for(int j = 0; j < m_cellLength; j++)
          {
            int cellIndex = m_startElementCell + j;
            ElementCell *elementCell = m_profile[i]->elementCells[cellIndex];
            elementCell->edgeHydHead[m_edge].isBC = true;

            ElementCell *nCell = nullptr;

            if((nCell = elementCell->neighbors[m_nEdge]))
            {
              nCell->edgeHydHead[m_nEdge].isBC = true;
            }
          }
        }
      }
      break;
    case GRAD_HYD_HEAD:
      {
        m_applyBCFunction = &EdgeBC::applyGradHydHeadBC;

        for(size_t i = 0 ; i < m_profile.size(); i++)
        {
          for(int j = 0; j < m_cellLength; j++)
          {
            int cellIndex = m_startElementCell + j;
            ElementCell *elementCell = m_profile[i]->elementCells[cellIndex];
            elementCell->edgeGradHydHead[m_edge].isBC = true;

            ElementCell *nCell = nullptr;

            if((nCell = elementCell->neighbors[m_nEdge]))
            {
              nCell->edgeGradHydHead[m_nEdge].isBC = true;
            }
          }
        }
      }
      break;
    case TEMPERATURE:
      {
        m_applyBCFunction = &EdgeBC::applyTempBC;

        for(size_t i = 0 ; i < m_profile.size(); i++)
        {
          for(int j = 0; j < m_cellLength; j++)
          {
            int cellIndex = m_startElementCell + j;
            ElementCell *elementCell = m_profile[i]->elementCells[cellIndex];
            elementCell->edgeTemperatures[m_edge].isBC = true;

            ElementCell *nCell = nullptr;

            if((nCell = elementCell->neighbors[m_nEdge]))
            {
              nCell->edgeTemperatures[m_nEdge].isBC = true;
            }
          }
        }
      }
      break;
    case GRAD_TEMPERATURE:
      {
        m_applyBCFunction = &EdgeBC::applyGradTempBC;

        for(size_t i = 0 ; i < m_profile.size(); i++)
        {
          for(int j = 0; j < m_cellLength; j++)
          {
            int cellIndex = m_startElementCell + j;
            ElementCell *elementCell = m_profile[i]->elementCells[cellIndex];
            elementCell->edgeGradTemperatures[m_edge].isBC = true;

            ElementCell *nCell = nullptr;

            if((nCell = elementCell->neighbors[m_nEdge]))
            {
              nCell->edgeGradTemperatures[m_nEdge].isBC = true;
            }
          }
        }
      }
      break;
    case SOLUTE:
      {
        m_applyBCFunction = &EdgeBC::applySoluteBC;

        for(size_t i = 0 ; i < m_profile.size(); i++)
        {
          for(int j = 0; j < m_cellLength; j++)
          {
            int cellIndex = m_startElementCell + j;
            ElementCell *elementCell = m_profile[i]->elementCells[cellIndex];
            elementCell->edgeSoluteConcs[m_soluteIndex][m_edge].isBC = true;

            ElementCell *nCell = nullptr;

            if((nCell = elementCell->neighbors[m_nEdge]))
            {
              nCell->edgeSoluteConcs[m_soluteIndex][m_nEdge].isBC = true;
            }
          }
        }
      }
      break;
    case GRAD_SOLUTE:
      {
        m_applyBCFunction = &EdgeBC::applyGradSoluteBC;

        for(size_t i = 0 ; i < m_profile.size(); i++)
        {
          for(int j = 0; j < m_cellLength; j++)
          {
            int cellIndex = m_startElementCell + j;
            ElementCell *elementCell = m_profile[i]->elementCells[cellIndex];
            elementCell->edgeGradSoluteConcs[m_soluteIndex][m_edge].isBC = true;

            ElementCell *nCell = nullptr;

            if((nCell = elementCell->neighbors[m_nEdge]))
            {
              nCell->edgeGradSoluteConcs[m_soluteIndex][m_nEdge].isBC = true;
            }
          }
        }
      }
      break;
  }

  if(m_timeSeries->numColumns() == m_cellLength * m_profile.size())
  {
    m_match = true;
  }
}

void EdgeBC::applyBoundaryConditions(double dateTime)
{
  (this->*m_applyBCFunction)(dateTime);
}

void EdgeBC::clear()
{
  m_profile.clear();

}

EdgeBC::Variable EdgeBC::variable() const
{
  return m_variable;
}

void EdgeBC::setVariable(Variable variable)
{
  m_variable = variable;
}

EdgeBC::Edge EdgeBC::edge() const
{
  return m_edge;
}

void EdgeBC::setEdge(EdgeBC::Edge edge)
{
  m_edge = edge;
}

int EdgeBC::soluteIndex() const
{
  return m_soluteIndex;
}

void EdgeBC::setSoluteIndex(int soluteIndex)
{
  m_soluteIndex = soluteIndex;
}

Element *EdgeBC::startElement() const
{
  return m_startElement;
}

void EdgeBC::setStartElement(Element *element)
{
  m_startElement = element;
}

Element *EdgeBC::endElement() const
{
  return m_endElement;
}

void EdgeBC::setEndElement(Element *element)
{
  m_endElement = element;
}

std::vector<Element*> EdgeBC::profiles() const
{
  return m_profile;
}

QSharedPointer<TimeSeries> EdgeBC::timeSeries() const
{
  return m_timeSeries;
}

void EdgeBC::setTimeSeries(const QSharedPointer<TimeSeries> &timeseries)
{
  m_timeSeries = timeseries;
  m_dataCursor->setMin(0);
  m_dataCursor->setMax(timeseries->numRows() - 1);
}

void EdgeBC::applyHydHeadBC(double dateTime)
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
          elementCell->edgeHydHead[m_edge].value = value;

          ElementCell *nCell = nullptr;

          if((nCell = elementCell->neighbors[m_nEdge]))
          {
            nCell->edgeHydHead[m_nEdge].value = value;
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
          elementCell->edgeHydHead[m_edge].value = value;

          ElementCell *nCell = nullptr;

          if((nCell = elementCell->neighbors[m_nEdge]))
          {
            nCell->edgeHydHead[m_nEdge].value = value;
          }
        }
      }
    }
  }
}

void EdgeBC::applyGradHydHeadBC(double dateTime)
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
          elementCell->edgeGradHydHead[m_edge].value = value;

          ElementCell *nCell = nullptr;

          if((nCell = elementCell->neighbors[m_nEdge]))
          {
            nCell->edgeGradHydHead[m_nEdge].value = value;
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
          elementCell->edgeGradHydHead[m_edge].value = value;

          ElementCell *nCell = nullptr;

          if((nCell = elementCell->neighbors[m_nEdge]))
          {
            nCell->edgeGradHydHead[m_nEdge].value = value;
          }
        }
      }
    }
  }
}

void EdgeBC::applyTempBC(double dateTime)
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
          elementCell->edgeTemperatures[m_edge].value = value;

          ElementCell *nCell = nullptr;

          if((nCell = elementCell->neighbors[m_nEdge]))
          {
            nCell->edgeTemperatures[m_nEdge].value = value;
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
          elementCell->edgeTemperatures[m_edge].value = value;

          ElementCell *nCell = nullptr;

          if((nCell = elementCell->neighbors[m_nEdge]))
          {
            nCell->edgeTemperatures[m_nEdge].value = value;
          }
        }
      }
    }
  }
}

void EdgeBC::applyGradTempBC(double dateTime)
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
          elementCell->edgeGradTemperatures[m_edge].value = value;

          ElementCell *nCell = nullptr;

          if((nCell = elementCell->neighbors[m_nEdge]))
          {
            nCell->edgeGradTemperatures[m_nEdge].value = value;
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
          elementCell->edgeGradTemperatures[m_edge].value = value;

          ElementCell *nCell = nullptr;

          if((nCell = elementCell->neighbors[m_nEdge]))
          {
            nCell->edgeGradTemperatures[m_nEdge].value = value;
          }
        }
      }
    }
  }
}

void EdgeBC::applySoluteBC(double dateTime)
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
          elementCell->edgeSoluteConcs[m_soluteIndex][m_edge].value = value;

          ElementCell *nCell = nullptr;

          if((nCell = elementCell->neighbors[m_nEdge]))
          {
            nCell->edgeSoluteConcs[m_soluteIndex][m_nEdge].value = value;
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
          elementCell->edgeSoluteConcs[m_soluteIndex][m_edge].value = value;

          ElementCell *nCell = nullptr;

          if((nCell = elementCell->neighbors[m_nEdge]))
          {
            nCell->edgeSoluteConcs[m_soluteIndex][m_nEdge].value = value;
          }
        }
      }
    }
  }
}

void EdgeBC::applyGradSoluteBC(double dateTime)
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
          elementCell->edgeGradSoluteConcs[m_soluteIndex][m_edge].value = value;

          ElementCell *nCell = nullptr;

          if((nCell = elementCell->neighbors[m_nEdge]))
          {
            nCell->edgeGradSoluteConcs[m_soluteIndex][m_nEdge].value = value;
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
          elementCell->edgeGradSoluteConcs[m_soluteIndex][m_edge].value = value;

          ElementCell *nCell = nullptr;

          if((nCell = elementCell->neighbors[m_nEdge]))
          {
            nCell->edgeGradSoluteConcs[m_soluteIndex][m_nEdge].value = value;
          }
        }
      }
    }
  }
}

/*!
 * \file elementoutput.cpp
 * \author  Caleb Amoa Buahin <caleb.buahin@gmail.com>
 * \version 1.0.0
 * \description
 * \license
 * This file and its associated files, and libraries are free software.
 * You can redistribute it and/or modify it under the terms of the
 * Lesser GNU Lesser General Public License as published by the Free Software Foundation;
 * either version 3 of the License, or (at your option) any later version.
 * This file and its associated files is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.(see <http://www.gnu.org/licenses/> for details)
 * \copyright Copyright 2014-2018, Caleb Buahin, All rights reserved.
 * \date 2014-2018
 * \pre
 * \bug
 * \warning
 * \todo
 */

#include "stdafx.h"
#include "elementoutput.h"
#include "gwcomponent.h"
#include "gwmodel.h"
#include "temporal/timedata.h"
#include "hydrocouplespatial.h"
#include "element.h"

using namespace HydroCouple;
using namespace HydroCouple::Spatial;
using namespace SDKTemporal;

ElementOutput::ElementOutput(const QString &id,
                             Dimension *timeDimension,
                             Dimension *geometryDimension,
                             ValueDefinition *valueDefinition,
                             VariableType variableType,
                             GWComponent *modelComponent)
  : TimeGeometryOutputDouble(id, IGeometry::LineString,
                             timeDimension, geometryDimension,
                             valueDefinition, modelComponent),
    m_component(modelComponent),
    m_variableType(variableType),
    m_soluteIndex(0)
{

}

ElementOutput::~ElementOutput()
{

}

void ElementOutput::updateValues(HydroCouple::IInput *querySpecifier)
{
  if(!m_component->workflow())
  {
    ITimeComponentDataItem* timeExchangeItem = dynamic_cast<ITimeComponentDataItem*>(querySpecifier);
    QList<IOutput*>updateList;

    if(timeExchangeItem)
    {
      double queryTime = timeExchangeItem->time(timeExchangeItem->timeCount() - 1)->julianDay();

      while (m_component->modelInstance()->currentDateTime() < queryTime &&
             m_component->status() == IModelComponent::Updated)
      {
        m_component->update(updateList);
      }
    }
    else
    {
      if(m_component->status() == IModelComponent::Updated)
      {
        m_component->update(updateList);
      }
    }
  }

  refreshAdaptedOutputs();
}

void ElementOutput::updateValues()
{
  moveDataToPrevTime();

  int currentTimeIndex = timeCount() - 1;
  DateTime *currentDateTime = m_times[currentTimeIndex];

  currentDateTime->setJulianDay(m_component->modelInstance()->currentDateTime());
  resetTimeSpan();

  switch (m_variableType)
  {
    case ChannelOutflow:
      {
        for(int i = 0 ; i < (int)m_geometries.size() ; i++)
        {
          Element *element = m_component->modelInstance()->getElement(i);
          double value = -element->channelInflow;
          setValue(currentTimeIndex,i,&value);
        }
      }
      break;
    case ChannelOuflowHeatFlux:
      {
        for(int i = 0 ; i < (int)m_geometries.size() ; i++)
        {
          Element *element = m_component->modelInstance()->getElement(i);
          double value = -element->channelHeatFlux;
          setValue(currentTimeIndex,i,&value);
        }
      }
      break;
    case ChannelOutflowSoluteFlux:
      {
        for(int i = 0 ; i < (int)m_geometries.size() ; i++)
        {
          Element *element = m_component->modelInstance()->getElement(i);
          double value = -element->channelSoluteFlux[m_soluteIndex];
          setValue(currentTimeIndex,i,&value);
        }
      }
      break;
  }

  refreshAdaptedOutputs();
}

ElementOutput::VariableType ElementOutput::variableType() const
{
  return m_variableType;
}

void ElementOutput::setVariableType(VariableType variableType)
{
  m_variableType = variableType;
}

int ElementOutput::soluteIndex() const
{
  return m_soluteIndex;
}

void ElementOutput::setSoluteIndex(int soluteIndex)
{
  m_soluteIndex = soluteIndex;
}

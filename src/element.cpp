/*!
*  \file    element.cpp
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

#include "stdafx.h"
#include "element.h"
#include "elementjunction.h"
#include "gwmodel.h"

Element::Element(const std::string &id, ElementJunction *upstream, ElementJunction *downstream,  GWModel *model)
  : id(id),
    upstreamJunction(upstream),
    downstreamJunction(downstream),
    channelSoluteConcs(nullptr),
    channelSoluteRate(nullptr),
    channelSoluteFlux(nullptr),
    model(model)
{

  channelTemperature = 0;
  hydConZ = model->m_hydConZ;

  upstream->outgoingElements.insert(this);
  downstream->incomingElements.insert(this);

  x = (upstream->x +  downstream->x) / 2.0;
  y = (upstream->y +  downstream->y) / 2.0;
  z = (upstream->z +  downstream->z) / 2.0;

  initializeElementCells();
  initializeSolutes();
}

Element::~Element()
{

  for(ElementCell *cell : elementCells)
    delete cell;

  elementCells.clear();

  if(channelSoluteConcs != nullptr)
  {
    delete[] channelSoluteConcs;
    delete[] channelSoluteRate;
    delete[] channelSoluteFlux;
  }

  upstreamJunction->outgoingElements.erase(this);
  downstreamJunction->incomingElements.erase(this);
}

void Element::initialize()
{
  channelTemperature = 0;

  setUpstreamElement();
  setDownStreamElement();

  channelInflow = 0.0;
  channelHeatRate = 0.0;
  channelHeatFlux = 0.0;

  for(int i = 0; i < model->m_totalCellsPerElement; i++)
  {
    elementCells[i]->initialize();
  }
}

void Element::initializeElementCells()
{
  for(ElementCell *cell : elementCells)
    delete cell;

  elementCells.clear();

  for(int i = 0; i < model->m_totalCellsPerElement; i++)
  {
    ElementCell *cell = new ElementCell(i, this);
    elementCells.push_back(cell);
  }
}

void Element::initializeSolutes()
{
  if(channelSoluteConcs)
  {
    delete[] channelSoluteConcs; channelSoluteConcs = nullptr;
    delete[] channelSoluteFlux; channelSoluteFlux = nullptr;
    delete[] channelSoluteRate; channelSoluteRate = nullptr;
  }

  if(model->m_solutes.size() > 0)
  {
    channelSoluteConcs = new double[model->m_solutes.size()]();
    channelSoluteFlux = new double[model->m_solutes.size()]();
    channelSoluteRate = new double[model->m_solutes.size()]();
  }

  for(int i = 0; i < model->m_totalCellsPerElement; i++)
  {
    elementCells[i]->initializeSolutes();
  }
}

void Element::setUpstreamElement()
{
  upstreamElement = nullptr;

  if(upstreamJunction->incomingElements.size() + upstreamJunction->outgoingElements.size() == 2)
  {
    for(Element *element : upstreamJunction->incomingElements)
    {
      if(element != this)
      {
        upstreamElement = element;
        return;
      }
    }

    for(Element *element : upstreamJunction->outgoingElements)
    {
      if(element != this)
      {
        upstreamElement = element;
        return;
      }
    }
  }
}

void Element::setDownStreamElement()
{
  downstreamElement = nullptr;

  if(downstreamJunction->incomingElements.size() + downstreamJunction->outgoingElements.size() == 2)
  {
    for(Element *element : downstreamJunction->outgoingElements)
    {
      if(element != this)
      {
        downstreamElement = element;
        return;
      }
    }

    for(Element *element : upstreamJunction->incomingElements)
    {
      if(element != this)
      {
        downstreamElement = element;
        return;
      }
    }
  }
}

/*!
*  \file    elementcell.cpp
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
#include "gwmodel.h"

using namespace std;

ElementCell::ElementCell(int cindex, Element *cparent)
  : elementCellIndex(cindex),
    edgeHydHead(nullptr),
    gradHydHead(nullptr),
    soluteConcs(nullptr),
    prevSoluteConcs(nullptr),
    externalSoluteFluxes(nullptr),
    totalSoluteMassBalance(nullptr),
    edgeFlows(nullptr),
    edgeDepths(nullptr),
    edgeBottomElevs(nullptr),
    edgeHydCons(nullptr),
    volume(0.0),
    prevVolume(0.0),
    dVolumedt(0.0),
    parentElement(cparent)
{
  hydConX = parentElement->model->m_hydConX;
  hydConY = parentElement->model->m_hydConY;
  sedDensity = parentElement->model->m_sedDensity;
  sedCp = parentElement->model->m_sedCp;
  porosity = parentElement->model->m_porosity;
  specificStorage = parentElement->model->m_specificStorage;
  width = parentElement->model->m_defaultCellWidth;

  edgeHydHead = new Variable[4];
  gradHydHead = new Variable[4];
  edgeFlows = new double[4]();
  edgeDepths = new double[4]();
  edgeHydCons = new double[4]();
  edgeBottomElevs = new double[4]();
}

ElementCell::~ElementCell()
{
  delete[] edgeHydHead;
  delete[] gradHydHead;
  delete[] edgeDepths;
  delete[] edgeHydCons;
  delete[] edgeBottomElevs;

  if(soluteConcs)
  {
    delete[] soluteConcs;
    delete[] prevSoluteConcs;
    delete[] externalSoluteFluxes;
    delete[] totalSoluteMassBalance;
  }
}

void ElementCell::initialize()
{
  start = true;
  dVolumedt = 0.0;

  //left
  if(elementCellIndex > 0)
  {
    ElementCell *left = parentElement->elementCells[elementCellIndex - 1];
    edgeBottomElevs[2] = bedRockElev + (width / 2.0) * (left->bedRockElev - bedRockElev) / (width / 2.0 + left->width / 2.0);
  }
  else
  {
    edgeBottomElevs[2] = bedRockElev;
  }

  //right
  if(elementCellIndex < parentElement->model->m_totalCellsPerElement - 1)
  {
    ElementCell *right = parentElement->elementCells[elementCellIndex + 1];
    edgeBottomElevs[0] = bedRockElev + (width / 2.0) * (right->bedRockElev - bedRockElev) / (width / 2.0 + right->width / 2.0);
  }
  else
  {
    edgeBottomElevs[0] = bedRockElev;
  }

  //upstream
  if(parentElement->index > 0)
  {
    Element *upElement = parentElement->model->m_elements[parentElement->index - 1];
    ElementCell *up = upElement->elementCells[elementCellIndex];
    edgeBottomElevs[3] = bedRockElev + (parentElement->length / 2.0) *
                         (up->bedRockElev - bedRockElev) / (parentElement->length / 2.0 + up->parentElement->length / 2.0);
  }
  else
  {
    edgeBottomElevs[3] = bedRockElev;
  }

  //downstream
  if(parentElement->index < parentElement->model->m_elements.size() - 1)
  {
    Element *downElement = parentElement->model->m_elements[parentElement->index + 1];
    ElementCell *down = downElement->elementCells[elementCellIndex];
    edgeBottomElevs[1] = bedRockElev + (parentElement->length / 2.0) *
                         (down->bedRockElev - bedRockElev) / (parentElement->length / 2.0 + down->parentElement->length / 2.0);
  }
  else
  {
    edgeBottomElevs[1] = bedRockElev;
  }
}

void ElementCell::initializeSolutes()
{
  if(soluteConcs)
  {
    delete[] soluteConcs; soluteConcs = nullptr;
    delete[] prevSoluteConcs; prevSoluteConcs = nullptr;
    delete[] externalSoluteFluxes; externalSoluteFluxes = nullptr;
    delete[] totalSoluteMassBalance; totalSoluteMassBalance = nullptr;
  }


}

double ElementCell::computeDHydHeadDt(double dt, double H[])
{
  double DHydHeadDt = 0.0;
  computeEdgeDepths(H);

  //Left
  if(edgeHydHead[2].isBC)
  {
    double edgeFlow = edgeHydCons[2] *  (edgeHydHead[2].value - H[index]) * edgeDepths[2] * parentElement->length / (width / 2.0);
    DHydHeadDt += edgeFlow;
    edgeFlows[2] = edgeFlow;
  }
  else if(gradHydHead[2].isBC)
  {
    double edgeFlow = edgeHydCons[2] * gradHydHead[2].value * depth * parentElement->length;
    edgeFlows[2] = edgeFlow;
  }
  else
  {
    if(elementCellIndex > 0)
    {
      ElementCell *left = parentElement->elementCells[elementCellIndex - 1];
      double edgeFlow = edgeHydCons[2] * (H[left->index] - H[index]) * depth * parentElement->length /
                        ((width / 2.0) + (left->width / 2.0));
      DHydHeadDt += edgeFlow;
      edgeFlows[2] = edgeFlow;
    }
    else
    {
      edgeFlows[2] = 0.0;
    }
  }


  //Right
  if(edgeHydHead[0].isBC)
  {
    double edgeFlow = edgeHydCons[0] *  (edgeHydHead[0].value - H[index]) * edgeDepths[0] * parentElement->length / (width / 2.0);
    DHydHeadDt += edgeFlow;
    edgeFlows[0] = edgeFlow;
  }
  else if(gradHydHead[0].isBC)
  {
    double edgeFlow = edgeHydCons[0] * gradHydHead[0].value * depth * parentElement->length;
    DHydHeadDt += edgeFlow;
    edgeFlows[0] = edgeFlow;
  }
  else
  {
    if(elementCellIndex < parentElement->model->m_totalCellsPerElement - 1)
    {
      ElementCell *right = parentElement->elementCells[elementCellIndex + 1];
      double edgeFlow = edgeHydCons[0] *  (H[right->index] - H[index]) * depth * parentElement->length /
                        ((width / 2.0) + (right->width / 2.0));

      DHydHeadDt += edgeFlow;
      edgeFlows[0] = edgeFlow;
    }
    else
    {
      edgeFlows[0] = 0.0;
    }
  }


  //Upstream
  if(edgeHydHead[3].isBC)
  {
    double edgeFlow = edgeHydCons[3] *  (edgeHydHead[3].value - H[index]) * edgeDepths[3] * width /
                      (parentElement->length / 2.0);

    DHydHeadDt += edgeFlow;
    edgeFlows[3] = edgeFlow;
  }
  else if(gradHydHead[3].isBC)
  {
    double edgeFlow = edgeHydCons[3] * gradHydHead[3].value * depth * width;
    DHydHeadDt += edgeFlow;
    edgeFlows[3] = edgeFlow;
  }
  else
  {
    if(parentElement->index > 0)
    {
      Element *upElement = parentElement->model->m_elements[parentElement->index - 1];
      ElementCell *up = upElement->elementCells[elementCellIndex];
      double edgeFlow = edgeHydCons[3] *  (H[up->index] - H[index]) * depth * width /
                        ((parentElement->length / 2.0) + (upElement->length / 2.0));

      DHydHeadDt += edgeFlow;
      edgeFlows[3] = edgeFlow;
    }
    else
    {
      edgeFlows[3] = 0.0;
    }
  }


  //Downstream
  if(edgeHydHead[1].isBC)
  {
    double edgeFlow = edgeHydCons[1] *  (edgeHydHead[1].value - H[index]) * edgeDepths[1] * width /
                      (parentElement->length / 2.0);

    DHydHeadDt += edgeFlow;
    edgeFlows[1] = edgeFlow;
  }
  else if(gradHydHead[1].isBC)
  {
    double edgeFlow = edgeHydCons[1] * gradHydHead->value * depth * width ;

    DHydHeadDt += edgeFlow;
    edgeFlows[1] = edgeFlow;
  }
  else
  {

    if(parentElement->index < parentElement->model->m_elements.size() - 1)
    {
      Element *downElement = parentElement->model->m_elements[parentElement->index + 1];
      ElementCell *down = downElement->elementCells[elementCellIndex];
      double edgeFlow = edgeHydCons[1] *  (H[down->index] - H[index]) * depth * width /
                        ((parentElement->length / 2.0) + (downElement->length / 2.0));

      DHydHeadDt += edgeFlow;
      edgeFlows[1] = edgeFlow;
    }
    else
    {
      edgeFlows[1] = 0.0;
    }
  }

  //compute river inflow
  DHydHeadDt = DHydHeadDt /*- specificStorage * H[index] * dVolumedt / depth*/  + externalInflow;
  DHydHeadDt = DHydHeadDt / (specificStorage * parentElement->length * width);

  return DHydHeadDt;
}

double ElementCell::computeDTDt(double dt, double T[])
{

}

double ElementCell::computeDSoluteDt(double dt, double S[], int soluteIndex)
{

}

double ElementCell::computeCourantFactor() const
{
  return max({
               hydConY * edgeDepths[0] / (specificStorage * width * width),
               hydConX * edgeDepths[1] / (specificStorage * parentElement->length * parentElement->length),
               hydConY * edgeDepths[2] / (specificStorage * width * width),
               hydConX * edgeDepths[3] / (specificStorage * parentElement->length * parentElement->length)
             });
}

void ElementCell::computeDVolumeDt()
{
  prevVolume = max(prevHydHead.value - bedRockElev, 0.000000001) * width * parentElement->length;
  volume = max(hydHead.value - bedRockElev, 0.00000000001) * width * parentElement->length;
  dVolumedt = (volume - prevVolume)/ parentElement->model->m_timeStep;
}

void ElementCell::computeDerivedHydraulics()
{
  //  if(start)
  //  {
  //    prevVolume = volume = (hydHead.value - bedRockElev) * width * parentElement->length;
  //    start = false;
  //  }
  //  else
  //  {
  //    prevVolume = volume;
  //    volume = (hydHead.value - bedRockElev) * width * parentElement->length;
  //  }

  //  dVolumedt = (volume - prevVolume) / parentElement->model->m_timeStep;

  computeEdgeHydCons();

  computeEdgeDepths();
}

void ElementCell::computeMassBalance(double timeStep)
{

}

void ElementCell::computeHeatBalance(double timeStep)
{

}

void ElementCell::computeSoluteBalance(double timeStep, int soluteIndex)
{

}

void ElementCell::computeEdgeHydCons()
{
  //calculate edge hydraulic conductivities
  if(hydConY > 0)
  {
    //left
    if(elementCellIndex > 0)
    {
      ElementCell *left = parentElement->elementCells[elementCellIndex - 1];

      if(left->hydConY > 0)
      {
        edgeHydCons[2] = ((width / 2.0) + (left->width / 2.0)) /
                         ((width / 2.0 / hydConY)+(left->width / 2.0 / left->hydConY));
      }
      else
      {
        edgeHydCons[2] = 0.0;
      }
    }
    else
    {
      edgeHydCons[2] = hydConY;
    }

    //right
    if(elementCellIndex < parentElement->model->m_totalCellsPerElement - 1)
    {
      ElementCell *right = parentElement->elementCells[elementCellIndex + 1];

      if(right->hydConY > 0)
      {
        edgeHydCons[0] = ((width / 2.0) + (right->width / 2.0)) /
                         ((width / 2.0 / hydConY)+(right->width / 2.0 / right->hydConY));
      }
      else
      {
        edgeHydCons[0] = 0.0;
      }
    }
    else
    {
      edgeHydCons[0] = hydConY;
    }
  }
  else
  {
    edgeHydCons[0] = edgeHydCons[2] = 0.0;
  }


  if(hydConX > 0.0)
  {

    //upstream
    if(parentElement->index > 0)
    {
      Element *upElement = parentElement->model->m_elements[parentElement->index - 1];
      ElementCell *up = upElement->elementCells[elementCellIndex];

      if(up->hydConX > 0)
      {
        edgeHydCons[3] = ((parentElement->length / 2.0) + (upElement->length / 2.0)) /
                         ((parentElement->length / 2.0 / hydConX)+(upElement->length / 2.0 / up->hydConX));
      }
      else
      {
        edgeHydCons[3] = 0.0;
      }
    }
    else
    {
      edgeHydCons[3] = hydConX;
    }

    //downstream
    if(parentElement->index < parentElement->model->m_elements.size() - 1)
    {
      Element *downElement = parentElement->model->m_elements[parentElement->index + 1];
      ElementCell *down = downElement->elementCells[elementCellIndex];

      if(down->hydConX > 0)
      {
        edgeHydCons[1] = ((parentElement->length / 2.0) + (downElement->length / 2.0)) /
                         ((parentElement->length / 2.0 / hydConX)+(downElement->length / 2.0 / down->hydConX));
      }
      else
      {
        edgeHydCons[1] = 0.0;
      }
    }
    else
    {
      edgeHydCons[1] = hydConX;
    }
  }
  else
  {
    edgeHydCons[1] = edgeHydCons[3] = 0.0;
  }
}

void ElementCell::computeEdgeDepths()
{

  depth = max(hydHead.value - bedRockElev, 0.0);

  //Left
  if(!edgeHydHead[2].isBC)
  {
    if(elementCellIndex > 0)
    {
      ElementCell *left = parentElement->elementCells[elementCellIndex - 1];
      edgeHydHead[2].value = hydHead.value + (width / 2.0) * (left->hydHead.value - hydHead.value) / (width / 2.0 + left->width / 2.0);
    }
    else
    {
      edgeHydHead[2].value = hydHead.value;
    }
  }

  edgeDepths[2] = max(edgeHydHead[2].value - edgeBottomElevs[2], 0.0);

  //Right
  if(!edgeHydHead[0].isBC)
  {
    if(elementCellIndex < parentElement->model->m_totalCellsPerElement - 1)
    {
      ElementCell *right = parentElement->elementCells[elementCellIndex + 1];
      edgeHydHead[0].value = hydHead.value + (width / 2.0) * (right->hydHead.value - hydHead.value) / (width / 2.0 + right->width / 2.0);
    }
    else
    {
      edgeHydHead[0].value = hydHead.value;
    }
  }

  edgeDepths[0] = max(edgeHydHead[0].value - edgeBottomElevs[0], 0.0);

  //Upstream
  if(!edgeHydHead[3].isBC)
  {
    if(parentElement->index > 0)
    {
      Element *upElement = parentElement->model->m_elements[parentElement->index - 1];
      ElementCell *up = upElement->elementCells[elementCellIndex];
      edgeHydHead[3].value = hydHead.value + (parentElement->length / 2.0) *
                             (up->hydHead.value - hydHead.value) / (parentElement->length / 2.0 + up->parentElement->length / 2.0);
    }
    else
    {
      edgeHydHead[3].value = hydHead.value;
    }
  }

  edgeDepths[3] = max(edgeHydHead[3].value - edgeBottomElevs[3], 0.0);


  //Downstream
  if(!edgeHydHead[1].isBC)
  {
    if(parentElement->index < parentElement->model->m_elements.size() - 1)
    {
      Element *downElement = parentElement->model->m_elements[parentElement->index + 1];
      ElementCell *down = downElement->elementCells[elementCellIndex];
      edgeHydHead[1].value = hydHead.value + (parentElement->length / 2.0) *
                             (down->hydHead.value - hydHead.value) / (parentElement->length / 2.0 + down->parentElement->length / 2.0);
    }
    else
    {
      edgeHydHead[1].value = hydHead.value;
    }
  }

  edgeDepths[1] = max(edgeHydHead[1].value - edgeBottomElevs[1], 0.0);

}

void ElementCell::computeEdgeDepths(double H[])
{

  depth = max(H[index] - bedRockElev, 0.000001);

  //  //Left
  //  if(!edgeHydHead[2].isBC)
  //  {
  //    if(elementCellIndex > 0)
  //    {
  //      ElementCell *left = parentElement->elementCells[elementCellIndex - 1];
  //      edgeHydHead[2].value = H[index] + (width / 2.0) * (left->hydHead.value - H[index]) / (width / 2.0 + left->width / 2.0);
  //    }
  //    else
  //    {
  //      edgeHydHead[2].value = H[index];
  //    }
  //  }

  //  edgeDepths[2] = max(edgeHydHead[2].value - edgeBottomElevs[2], 0.0);

  //  //Right
  //  if(!edgeHydHead[0].isBC)
  //  {
  //    if(elementCellIndex < parentElement->model->m_totalCellsPerElement - 1)
  //    {
  //      ElementCell *right = parentElement->elementCells[elementCellIndex + 1];
  //      edgeHydHead[0].value = H[index] + (width / 2.0) * (right->hydHead.value - H[index]) / (width / 2.0 + right->width / 2.0);
  //    }
  //    else
  //    {
  //      edgeHydHead[0].value = H[index];
  //    }
  //  }

  //  edgeDepths[0] = max(edgeHydHead[0].value - edgeBottomElevs[0], 0.0);

  //  //Upstream
  //  if(!edgeHydHead[3].isBC)
  //  {
  //    if(parentElement->index > 0)
  //    {
  //      Element *upElement = parentElement->model->m_elements[parentElement->index - 1];
  //      ElementCell *up = upElement->elementCells[elementCellIndex];
  //      edgeHydHead[3].value = H[index] + (parentElement->length / 2.0) *
  //                             (up->hydHead.value - H[index]) / (parentElement->length / 2.0 + up->parentElement->length / 2.0);
  //    }
  //    else
  //    {
  //      edgeHydHead[3].value = H[index];
  //    }
  //  }

  //  edgeDepths[3] = max(edgeHydHead[3].value - edgeBottomElevs[3], 0.0);


  //  //Downstream
  //  if(!edgeHydHead[1].isBC)
  //  {
  //    if(parentElement->index < parentElement->model->m_elements.size() - 1)
  //    {
  //      Element *downElement = parentElement->model->m_elements[parentElement->index + 1];
  //      ElementCell *down = downElement->elementCells[elementCellIndex];
  //      edgeHydHead[1].value = H[index] + (parentElement->length / 2.0) *
  //                             (down->hydHead.value - H[index]) / (parentElement->length / 2.0 + down->parentElement->length / 2.0);
  //    }
  //    else
  //    {
  //      edgeHydHead[1].value = H[index];
  //    }
  //  }

  //  edgeDepths[1] = max(edgeHydHead[1].value - edgeBottomElevs[1], 0.0);

}

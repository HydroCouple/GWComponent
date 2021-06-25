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
#include "elementcell.h"
#include "gwmodel.h"

using namespace std;

const int ElementCell::nEdgeIndex[] = {2, 3, 0, 1, 5, 4};
const double ElementCell::ef1[] = {0.0, 1.0, 0.0, 1.0, 2.0 , 2.0};
const double ElementCell::ef2[] = {1.0, 0.0, 1.0, 0.0};

ElementCell::ElementCell(int cindex, Element *cparent)
  : elementCellIndex(cindex),
    edgeHydHead(nullptr),
    edgeGradHydHead(nullptr),
    soluteConcs(nullptr),
    prevSoluteConcs(nullptr),
    externalSoluteFluxes(nullptr),
    channelSoluteFlux(nullptr),
    channelSoluteRate(nullptr),
    totalSoluteMassBalance(nullptr),
    edgeFlows(nullptr),
    edgeHydCons(nullptr),
    parentElement(cparent),
    retardationFactor(nullptr),
    topBedCell(nullptr),
    bedCells(nullptr),
    isBedCell(false)
{
  edgeEffectiveKe = new double[6]();
  edgePorosity = new double[6]();
  edgeFlowsComp = nullptr;

  hydCon = new double[3]();
  hydCon[0] = parentElement->model->m_hydConY;
  hydCon[1] = parentElement->model->m_hydConX;
  hydCon[2] = parentElement->model->m_hydConZ;

  sedDensity = parentElement->model->m_sedDensity;
  sedCp = parentElement->model->m_sedCp;
  porosity = parentElement->model->m_porosity;
  specificYield = parentElement->model->m_specificYield;
  specificStorage = parentElement->model->m_specificStorage;
  width = parentElement->model->m_defaultCellWidth;

  dispersivity = new double[3]();
  dispersivity[0] = parentElement->model->m_dispersivityY;
  dispersivity[1] = parentElement->model->m_dispersivityX;
  dispersivity[2] = parentElement->model->m_dispersivityZ;

  sedThermalConductivity = parentElement->model->m_sedThermalConductivity;

  bedRockElevs = new double[4]();
  edgeHydHead = new Variable[6];
  edgeGradHydHead = new Variable[6];
  edgeFlows = new double[6]();
  edgeHydCons = new double[6]();
  edgeMechDispersionCoeff = new double[6]();
  edgeHeatPecletNumbers = new double[4]();
  edgeSolutePecletNumbers = new double[4]();
  edgeTemperatures = new Variable[6];
  edgeGradTemperatures = new Variable[6];
  edgeHeatFluxes = new double[6]();
  computeEdgeHeadDerivs = new ComputeEdgeDeriv[6];
  computeEdgeTempDerivs = new ComputeEdgeDeriv[6];

}

ElementCell::~ElementCell()
{
  delete[] edgeHydHead;
  delete[] edgeEffectiveKe;
  delete[] edgePorosity;
  delete[] edgeGradHydHead;
  delete[] edgeHydCons;
  delete[] edgeTemperatures;
  delete[] edgeGradTemperatures;
  delete[] edgeHeatFluxes;
  delete[] edgeMechDispersionCoeff;
  delete[] edgeSolutePecletNumbers;
  delete[] computeEdgeHeadDerivs;
  delete[] computeEdgeTempDerivs;
  delete[] bedRockElevs;

  if(edgeFlowsComp != nullptr)
  {
    for(int i = 0 ; i < 4; i++)
    {
      if(edgeFlowsComp[i] != nullptr)
      {
        delete[]  edgeFlowsComp[i];
      }
    }

    delete[] edgeFlowsComp;
  }

  deleteSoluteVariables();
  deleteBedCells();
}

ElementCell *ElementCell::neighbour(int edge, bool top)
{
  switch (edge)
  {
    case 0:
      {
        int nindex = elementCellIndex + 1;

        if(nindex >= 0 && nindex < parentElement->model->m_totalCellsPerElement)
        {
          ElementCell *cell = parentElement->elementCells[nindex];

          if(top)
          {
            return cell;
          }
          else
          {
            if(isBedCell && cell->isBedCell)
            {
              return cell->topBedCell->bedCells[zIndex];
            }
            else /*if(!isBedCell && !cell->isBedCell)*/
            {
              return  cell;
            }
          }
        }

      }
      break;
    case 1:
      {
        int nindex = parentElement->index + 1;

        if(nindex >= 0 && nindex < (int)parentElement->model->m_elements.size())
        {
          ElementCell *cell = parentElement->model->m_elements[nindex]->elementCells[elementCellIndex];

          if(top)
          {
            return cell;
          }
          else
          {
            if(isBedCell && cell->isBedCell)
            {
              return cell->topBedCell->bedCells[zIndex];
            }
            else /*if(!isBedCell && !cell->isBedCell)*/
            {
              return  cell;
            }
          }
        }

      }
      break;
    case 2:
      {
        int nindex = elementCellIndex - 1;

        if(nindex >= 0 && nindex < parentElement->model->m_totalCellsPerElement)
        {
          ElementCell *cell = parentElement->elementCells[nindex];

          if(top)
          {
            return cell;
          }
          else
          {
            if(isBedCell && cell->isBedCell)
            {
              return cell->topBedCell->bedCells[zIndex];
            }
            else /*if(!isBedCell && !cell->isBedCell)*/
            {
              return  cell;
            }
          }
        }

      }
      break;
    case 3:
      {
        int nindex = parentElement->index - 1;

        if(nindex >= 0 && nindex < (int)parentElement->model->m_elements.size())
        {
          ElementCell *cell = parentElement->model->m_elements[nindex]->elementCells[elementCellIndex];

          if(top)
          {
            return cell;
          }
          else
          {
            if(isBedCell && cell->isBedCell)
            {
              return cell->topBedCell->bedCells[zIndex];
            }
            else /*if(!isBedCell && !cell->isBedCell)*/
            {
              return  cell;
            }
          }
        }
      }
      break;
    case 4:
      {
        int ind = zIndex - 1;

        if(isBedCell && ind >= 0 && ind < parentElement->model->m_numBedZCells)
        {
          return topBedCell->bedCells[ind];
        }
      }
      break;
    case 5:
      {
        int ind = zIndex + 1;

        if(isBedCell && ind >= 0 && ind < parentElement->model->m_numBedZCells)
        {
          return topBedCell->bedCells[ind];
        }
      }
      break;
  }

  return nullptr;
}

double ElementCell::flowLengthP(int edge)
{
  switch (edge)
  {
    case 0:
    case 2:
      {
        return  width / 2.0;
      }
    case 1:
    case 3:
      {
        return  parentElement->length / 2.0;
      }
    case 4:
      {
        if(isTopCell())
        {
          return max(0.0, hydHead.value - bottomElev) / 2.0;
        }
        else
        {
          return (topElev - bottomElev) / 2.0;
        }
      }
    case 5:
      {
        if(isTopCell())
        {
          return max(0.0, hydHead.value - bottomElev) / 2.0;
        }
        else
        {
          return (topElev - bottomElev) / 2.0;
        }
      }
  }

  return 0.0;
}

double ElementCell::flowLengthN(int edge, ElementCell *n)
{
  if(n)
  {
    switch (edge)
    {
      case 0:
      case 2:
        {
          return  n->width / 2.0;
        }
      case 1:
      case 3:
        {
          return  n->parentElement->length / 2.0;
        }
      case 4:
        {
          if(n)
          {
            if(n->isTopCell())
            {
              return 0.0;
            }
            else
            {
              return (n->topElev - n->bottomElev) / 2.0;
            }
          }
          else
          {
            return 0.0;
          }
        }
      case 5:
        {
          if(n)
          {
            if (n->isTopCell())
            {
              return max(0.0, n->hydHead.value - n->bottomElev)/2.0;
            }
            else
            {
              return (n->topElev - n->bottomElev) / 2.0;
            }
          }
          else
          {
            return  0.0;
          }
        }
    }

    return 0.0;
  }

  return 0.0;
}

double ElementCell::flowArea(int edge, int layer)
{
  switch (edge)
  {
    case 0:
    case 2:
      {
        return computeEdgeDepth(edge, layer) * parentElement->length;
      }
    case 1:
    case 3:
      {
        return computeEdgeDepth(edge, layer) * width;
      }
    case 4:
    case 5:
      {
        return width * parentElement->length;
      }
  }

  return 0.0;
}

ElementCell *ElementCell::layer(int layer)
{
  if(isBedCell)
  {
    return bedCells[layer];
  }

  return this;
}

void ElementCell::initializeBedCells()
{
  deleteBedCells();

  if(isBedCell && isTopCell())
  {
    bedCells = new ElementCell*[parentElement->model->m_numBedZCells];
    topBedCell = this;
    this->zIndex = 0;
    bedCells[0] = this;

    double totalDz = topElev - bedRockElev;

    double currentTop = topElev - totalDz * parentElement->model->m_zcellFactors[0];
    bottomElev = currentTop;

    for (int i = 1 ; i < parentElement->model->m_numBedZCells; i++)
    {
      ElementCell *elementCell = new ElementCell(elementCellIndex, parentElement);
      elementCell->width = width;
      elementCell->elementCellIndex = elementCellIndex;
      elementCell->topBedCell = this;
      elementCell->hydHead.value = elementCell->prevHydHead.value = hydHead.value;
      elementCell->temperature.value = elementCell->prevTemperature.value = temperature.value;
      elementCell->zIndex = i;
      elementCell->topElev = currentTop;
      currentTop = currentTop - totalDz * parentElement->model->m_zcellFactors[i];
      elementCell->bottomElev = currentTop;
      elementCell->bedRockElev = bedRockElev;
      elementCell->specificYield = 0.0;
      elementCell->isBedCell = true;
      elementCell->sedDensity = sedDensity;
      elementCell->sedCp = elementCell->sedCp;
      elementCell->porosity = porosity;
      bedCells[i] = elementCell;

      elementCell->initializeSolutes();

      int numSolutes = parentElement->model->m_solutes.size();

      for(int l = 0; l < numSolutes; l++ )
      {
        elementCell->soluteConcs[l].value = soluteConcs[l].value;
      }
    }
  }
}

void ElementCell::deleteBedCells()
{
  if(bedCells != nullptr)
  {
    for(int i = 1; i < parentElement->model->m_numBedZCells; i++)
    {
      delete bedCells[i];
    }

    delete[] bedCells; bedCells = nullptr;
  }
}

bool ElementCell::isTopCell()
{
  return  this == topBedCell || topBedCell == nullptr;
}

void ElementCell::initialize()
{

  start = true;
  dvolume_dt = 0.0;

  totalMassBalance = 0.0;
  totalHeatBalance = 0.0;
  channelInflow = 0.0;
  channelHeatFlux = 0.0;
  channelHeatRate = 0.0;
  externalHeatFluxes = 0.0;
  externalInflow = 0.0;
  depth = 0;

  if(isBedCell)
  {
    if(zIndex == 0)
    {
      specificStorage = 0.0;
    }
    else
    {
      specificYield = 0.0;
    }
  }
  else
  {
    specificStorage = 0.0;
  }

  for(int i = 0; i < 4; i++)
  {
    ElementCell *n = neighbour(i);

    if(n)
    {
      double flp = flowLengthP(i);
      double fln = flowLengthN(i, n);
      double slp = (n->bedRockElev - bedRockElev) / (flp + fln);
      bedRockElevs[i] = bedRockElev + flp * slp;
    }
    else if((n = neighbour(nEdgeIndex[i])))
    {
      int ind = nEdgeIndex[i];
      double flp = flowLengthP(ind);
      double fln = flowLengthN(ind, n);
      double slp = (n->bedRockElev - bedRockElev) / (flp + fln);
      bedRockElevs[i] = bedRockElev - flp * slp;
    }
    else
    {
      bedRockElevs[i] = bedRockElev;
    }
  }

  switch (parentElement->model->m_advectionMode)
  {
    case GWModel::AdvectionDiscretizationMode::Central:
      {
        computeTempAdv = &ElementCell::computeDTDtCentral;
        computeSoluteAdv = &ElementCell::computeDSoluteDtCentral;
      }
      break;
    case GWModel::AdvectionDiscretizationMode::Hybrid:
      {
        computeTempAdv = &ElementCell::computeDTDtHybrid;
        computeSoluteAdv = &ElementCell::computeDSoluteDtHybrid;
      }
      break;
    case GWModel::AdvectionDiscretizationMode::TVD:
      {
        computeTempAdv = &ElementCell::computeDTDtTVD;
        computeSoluteAdv = &ElementCell::computeDSoluteDtTVD;
      }
      break;
    default:
      {
        computeTempAdv = &ElementCell::computeDTDtUpwind;
        computeSoluteAdv = &ElementCell::computeDSoluteDtUpwind;
      }
      break;
  }
}

void ElementCell::initializeSolutes()
{
  deleteSoluteVariables();

  if(parentElement->model->m_solutes.size() > 0)
  {
    int numSolutes = parentElement->model->m_solutes.size();

    soluteConcs = new Variable[numSolutes];
    prevSoluteConcs = new Variable[numSolutes];
    externalSoluteFluxes = new double[numSolutes]();
    channelSoluteFlux = new double[numSolutes]();
    channelSoluteRate = new double[numSolutes]();
    totalSoluteMassBalance = new double[numSolutes]();
    computeEdgeSoluteDerivs = new ComputeEdgeVariableDeriv *[numSolutes];

    edgeSoluteConcs = new Variable*[numSolutes];
    edgeGradSoluteConcs = new Variable*[numSolutes];
    edgeSoluteConcFluxes = new double*[numSolutes];
    retardationFactor = new double[numSolutes]();

    for(int i = 0; i < numSolutes; i++)
    {
      edgeSoluteConcs[i] = new Variable[6];
      edgeSoluteConcFluxes[i] = new double[6]();
      edgeGradSoluteConcs[i] = new Variable[6];
      computeEdgeSoluteDerivs[i] = new ComputeEdgeVariableDeriv[6];

    }
  }
}

double ElementCell::computeDHydHeadDt(double dt, double H[])
{
  double DHydHeadDt = 0.0;

  for(int i = 0; i < 6; i++)
  {
    edgeFlows[i] = 0.0;
    DHydHeadDt += (this->*computeEdgeHeadDerivs[i])(i, dt, H);
  }

  //compute river inflow
  DHydHeadDt += channelInflow;
  DHydHeadDt += externalInflow;
  //  DHydHeadDt /= ((specificYield * volume)  + (specificStorage * volume * depth));
  DHydHeadDt /= ((specificYield * parentElement->length * width)  + (specificStorage * volume));
  return DHydHeadDt;
}

double ElementCell::computeDTDt(double dt, double T[])
{
  double DTDt = 0.0;

  //Compute advection
  DTDt += (this->*computeTempAdv)(dt, T);

  //Compute dispersion
  DTDt += computeDTDtDispersion(dt, T);

  //Add external sources
  {
    DTDt += externalHeatFluxes / (rhom_Cm * volume);
  }

  //Channel Heat
  {
    DTDt += channelHeatRate / (rhom_Cm * volume);
  }

//  DTDt -= T[index] * dvolume_dt / volume;

  return DTDt;
}

double ElementCell::computeDTDtDispersion(double dt, double T[])
{
  double DTDt = 0.0;

  for(int i = 0; i < 6; i++)
  {
    DTDt += (this->*computeEdgeTempDerivs[i])(i, dt, T);
  }

  DTDt = DTDt / (rhom_Cm * volume);

  return DTDt;
}

double ElementCell::computeDTDtUpwind(double dt, double T[])
{
  double DTDt = 0.0;
  double heat = 0.0;

  for(int i = 0; i < 6; i++)
  {
    double flow = edgeFlows[i];
    edgeHeatFluxes[i] = 0.0;

    ElementCell *neigh = nullptr;

    if(flow >= 0.0)
    {
      if(edgeTemperatures[i].isBC)
      {
        heat = parentElement->model->m_waterDensity * parentElement->model->m_cp * flow * edgeTemperatures[i].value;

        DTDt += heat;
        edgeHeatFluxes[i] += heat;
      }
      else if((neigh = neighbour(i)))
      {
        heat = 0.0;

        if(i < 4)
        {
          if(!isBedCell && neigh->isBedCell)
          {

            for(int z = 0 ; z  < parentElement->model->m_numBedZCells; z++)
            {
              ElementCell *lay = neigh->topBedCell->bedCells[z];
              double tFlow = edgeFlowsComp[i][z];
              heat += parentElement->model->m_waterDensity * parentElement->model->m_cp * tFlow * T[lay->index];
            }
          }
          else if(isBedCell && neigh->isBedCell)
          {
            ElementCell *lay = neigh->topBedCell->bedCells[zIndex];
            heat = parentElement->model->m_waterDensity * parentElement->model->m_cp * flow * T[lay->index];
          }
          else
          {
            heat = parentElement->model->m_waterDensity * parentElement->model->m_cp * flow * T[neigh->index];
          }
        }
        else
        {
          heat = parentElement->model->m_waterDensity * parentElement->model->m_cp * flow * T[neigh->index];
        }

        DTDt += heat;
        edgeHeatFluxes[i] += heat;

      }
    }
    else
    {
      heat =  parentElement->model->m_waterDensity * parentElement->model->m_cp * flow * T[index];

      DTDt += heat;
      edgeHeatFluxes[i] += heat;
    }
  }

  DTDt = DTDt / (rhom_Cm * volume);

  return DTDt;
}

double ElementCell::computeDTDtCentral(double dt, double T[])
{
  double DTDt = 0.0;
  double heat = 0.0;

  for(int i = 0; i < 6; i++)
  {
    double edgeFlow = edgeFlows[i];
    edgeHeatFluxes[i] = 0.0;
    ElementCell *n = nullptr;

    if(edgeTemperatures[i].isBC)
    {
      heat =  parentElement->model->m_waterDensity * parentElement->model->m_cp *
              edgeFlow * edgeTemperatures[i].value;

      DTDt += heat;
      edgeHeatFluxes[i] += heat;
    }
    else if((n = neighbour(i)))
    {
      heat = 0.0;

      double fp = flowLengthP(i);
      double fn = flowLengthN(i,n);
      double ft = fp + fn;

      if(i < 4)
      {
        if (!isBedCell && n->isBedCell)
        {
          for(int k = 0; k < parentElement->model->m_numBedZCells; k++)
          {
            ElementCell *clayer = n->bedCells[k];
            double tFlow = edgeFlowsComp[i][k];
            double interpTemp =  T[index] * fn / ft + T[clayer->index] * fp / ft;
            heat += parentElement->model->m_waterDensity * parentElement->model->m_cp * tFlow * interpTemp;
          }
        }
        else
        {
          heat = parentElement->model->m_waterDensity * parentElement->model->m_cp * edgeFlow * (T[index] / fp   +  T[n->index] / fn) /
              (1.0 / fp + 1.0 / fn);
        }
      }
      else
      {
        heat = parentElement->model->m_waterDensity * parentElement->model->m_cp * edgeFlow * (T[index] / fp   +  T[n->index] / fn) /
            (1.0 / fp + 1.0 / fn);
      }


      DTDt += heat;
      edgeHeatFluxes[i] += heat;
    }
    else if(edgeFlow < 0)
    {
      heat =  parentElement->model->m_waterDensity * parentElement->model->m_cp * edgeFlow * T[index];
      DTDt += heat;
      edgeHeatFluxes[i] += heat;
    }
  }

  DTDt = DTDt / (rhom_Cm * volume);

  return DTDt;
}

double ElementCell::computeDTDtHybrid(double dt, double T[])
{
  double DTDt = 0.0;

  return DTDt;
}

double ElementCell::computeDTDtTVD(double dt, double T[])
{

  double DTDt = 0.0;

  for(int i = 0; i < 6; i++)
  {
    double flow = edgeFlows[i];

    if(flow > 0.0)
    {
      if(edgeTemperatures[i].isBC)
      {
        DTDt +=  flow * edgeTemperatures[i].value;
      }
      else if(neighbour(i))
      {
        DTDt +=  flow * T[neighbour(i)->index];
      }
    }
    else
    {
      DTDt += flow * T[index];
    }
  }

  DTDt = DTDt * parentElement->model->m_waterDensity * parentElement->model->m_cp / (rhom_Cm * volume);

  return DTDt;

}

double ElementCell::computeDSoluteDt(double dt, double S[], int soluteIndex)
{
  double DSoluteDt = 0;

  {
    //Compute advection
    DSoluteDt += (this->*computeSoluteAdv)(dt, S, soluteIndex);

    //Compute dispersion
    DSoluteDt += computeDSoluteDtDispersion(dt, S, soluteIndex);

    //First order reaction
    {
      DSoluteDt -= (soluteConcs[soluteIndex].value * parentElement->model->m_solute_first_order_k[soluteIndex]) / (retardationFactor[soluteIndex]);
    }

    double rf = retardationFactor[soluteIndex];

    //Add channel solute
    {
      DSoluteDt += channelSoluteRate[soluteIndex] / (rf * volume);
    }

    //Add external sources
    {
      DSoluteDt += externalSoluteFluxes[soluteIndex] / (rf * volume);
    }

  }


  dsolute_dt = DSoluteDt;

  //  DSoluteDt = 0.0;
  return DSoluteDt;
}

double ElementCell::computeDSoluteDtDispersion(double dt, double S[], int soluteIndex)
{
  double DSoluteDt = 0.0;

  for(int i = 0; i < 6; i++)
  {
    DSoluteDt += (this->*computeEdgeSoluteDerivs[soluteIndex][i])(i, dt, S, soluteIndex);
  }

  DSoluteDt = DSoluteDt / (retardationFactor[soluteIndex] * volume);

  return DSoluteDt;
}

double ElementCell::computeDSoluteDtUpwind(double dt, double S[], int soluteIndex)
{
  double DSoluteDt = 0.0;

  for(int i = 0; i < 6; i++)
  {
    ElementCell *neigh = nullptr;
    edgeSoluteConcFluxes[soluteIndex][i] = 0.0;
    double flow = edgeFlows[i];
    double edgeSoluteFlux = 0.0;

    if(flow >= 0.0)
    {
      if(edgeSoluteConcs[soluteIndex][i].isBC)
      {
        edgeSoluteFlux = flow * edgeSoluteConcs[soluteIndex][i].value;

        DSoluteDt += edgeSoluteFlux;
        edgeSoluteConcFluxes[soluteIndex][i] += edgeSoluteFlux;
      }
      else if((neigh = neighbour(i)))
      {
        double flp = flowLengthP(i) / 2.0;
        double fln = flowLengthN(i,neigh) / 2.0;

        if (i < 4)
        {
          if(!isBedCell && neigh->isBedCell)
          {

            double numer = S[index]/flp;
            double inFactor = 1.0 / fln;
            double denom = 1.0 / flp;

            for(int k = 0; k < parentElement->model->m_numBedZCells; k++)
            {
              ElementCell *clayer = neigh->bedCells[k];
              numer += S[clayer->index] * inFactor;
              denom += inFactor;
            }

            double edgeSoluteConc = numer/denom;
            edgeSoluteFlux = flow * edgeSoluteConc;
            edgeSoluteConcFluxes[soluteIndex][i] += edgeSoluteFlux;
          }
          else if(isBedCell && neigh->isBedCell)
          {
            ElementCell *lay = neigh->topBedCell->bedCells[zIndex];
            edgeSoluteFlux = flow * S[lay->index];
            edgeSoluteConcFluxes[soluteIndex][i] += edgeSoluteFlux;
          }
          else
          {
            edgeSoluteFlux = flow * S[neigh->index];
            edgeSoluteConcFluxes[soluteIndex][i] += edgeSoluteFlux;
          }
        }
        else
        {
          edgeSoluteFlux = flow * S[neigh->index];
          edgeSoluteConcFluxes[soluteIndex][i] += edgeSoluteFlux;
        }

        DSoluteDt += edgeSoluteFlux;
      }
    }
    else
    {
      edgeSoluteFlux =  flow * S[index];

      DSoluteDt += edgeSoluteFlux;
      edgeSoluteConcFluxes[soluteIndex][i] += edgeSoluteFlux;
    }
  }

  DSoluteDt = DSoluteDt / (retardationFactor[soluteIndex] * volume * porosity);

  return DSoluteDt;


}

double ElementCell::computeDSoluteDtCentral(double dt, double S[], int soluteIndex)
{
  double DSoluteDt = 0.0;
  double edgeSoluteFlux = 0.0;

  for(int i = 0; i < 6; i++)
  {
    double edgeFlow = edgeFlows[i];
    edgeSoluteConcFluxes[soluteIndex][i] = 0.0;

    if(edgeSoluteConcs[soluteIndex][i].isBC)
    {
      edgeSoluteFlux =  edgeFlow * edgeSoluteConcs[soluteIndex][i].value;
      DSoluteDt += edgeSoluteFlux;
      edgeSoluteConcFluxes[soluteIndex][i] = edgeSoluteFlux;
    }
    else if(neighbour(i))
    {
      ElementCell *n = neighbour(i);
      double flp = flowLengthP(i) / 2.0;
      double fln = flowLengthN(i,n) / 2.0;

      if (i < 4 && !isBedCell && n->isBedCell)
      {

        double numer = S[index]/flp;
        double inFactor = 1.0 / fln;
        double denom = 1.0 / flp;

        for(int k = 0; k < parentElement->model->m_numBedZCells; k++)
        {
          ElementCell *clayer = n->bedCells[k];
          numer += S[clayer->index] * inFactor;
          denom += inFactor;
        }

        double edgeSoluteConc = numer/denom;
        edgeSoluteFlux = edgeFlow * edgeSoluteConc;
      }
      else
      {
        edgeSoluteFlux = edgeFlow * (S[index] / flp   +  S[n->index] / fln) / (1.0 / flp + 1.0 / fln);
      }

      DSoluteDt += edgeSoluteFlux;
      edgeSoluteConcFluxes[soluteIndex][i] = edgeSoluteFlux;
    }
    else if(edgeFlow < 0.0)
    {
      edgeSoluteFlux =  edgeFlow * S[index];
      DSoluteDt += edgeSoluteFlux;
      edgeSoluteConcFluxes[soluteIndex][i] = edgeSoluteFlux;
    }
  }

  DSoluteDt = DSoluteDt / (retardationFactor[soluteIndex] * volume * porosity);

  return DSoluteDt;

}

double ElementCell::computeDSoluteDtHybrid(double dt, double S[], int soluteIndex)
{

  return 0.0;
}

double ElementCell::computeDSoluteDtTVD(double dt, double S[], int soluteIndex)
{

  return 0.0;
}

double ElementCell::computeCourantFactor() const
{
  double alpha = hydCon[0] * depth / (specificYield * width * width);
  double beta  = hydCon[1] * depth / (specificYield * parentElement->length * parentElement->length);
  double factor = max(alpha, beta);

  return factor ;
}

double ElementCell::computeDiffusionFactor()
{
  //Change to solute and heat dispersion based on (Woods, Teubner, Simmons, & Narayan, 2003)

  double diffFactor = 0;

  for(int i = 0; i < 6 ; i++)
  {
    double flp = flowLengthP(i);
    diffFactor = max(diffFactor, ((edgeEffectiveKe[i] / rhom_Cm) + (edgeMechDispersionCoeff[i] / (parentElement->model->m_waterDensity * parentElement->model->m_cp)))
                     / (flp * flp));
  }

  return diffFactor;
}

void ElementCell::calculatePreComputedHydHeadVariables()
{
  ElementCell *n = nullptr;

  for(int i = 0; i < 6; i++)
  {
    if(edgeHydHead[i].isBC)
    {
      computeEdgeHeadDerivs[i] = &ElementCell::computeEdgeHydHeadBC;
    }
    else if(edgeGradHydHead[i].isBC)
    {
      computeEdgeHeadDerivs[i] = &ElementCell::computeEdgeGradHydHeadBC;
    }
    else if((n = neighbour(i)))
    {
      if(i < 4)
      {
        if(!isBedCell && n->isBedCell)
        {
          computeEdgeHeadDerivs[i] = &ElementCell::computeNeighborLayersHydHeadBC;
          //          computeEdgeHeadDerivs[i] = &ElementCell::computeZeroHydHeadBC;

          if(edgeFlowsComp == nullptr)
          {
            edgeFlowsComp = new double*[4];

            for(int f = 0; f < 4; f++)
            {
              edgeFlowsComp[f] = new double[parentElement->model->m_numBedZCells]();
            }
          }
        }
        //        else if(isBedCell && !n->isBedCell)
        //        {

        //          computeEdgeHeadDerivs[i] = &ElementCell::computeZeroHydHeadBC;
        //        }
        else
        {
          computeEdgeHeadDerivs[i] = &ElementCell::computeNeighborHydHeadBC;
        }
      }
      else
      {
        computeEdgeHeadDerivs[i] = &ElementCell::computeNeighborHydHeadBC;
      }
    }
    else
    {
      computeEdgeHeadDerivs[i] = &ElementCell::computeZeroHydHeadBC;
    }
  }

  computeDepth();

  computeEdgeHydCons();

  computeEdgeDispersionCoefficients();

}

void ElementCell::calculatePreComputedTempVariables()
{
  ElementCell *n = nullptr;

  for(int i = 0; i < 6; i++)
  {
    if(edgeTemperatures[i].isBC)
    {
      computeEdgeTempDerivs[i] = &ElementCell::computeEdgeTempBC;
    }
    else if(edgeGradHydHead[i].isBC)
    {
      computeEdgeTempDerivs[i] = &ElementCell::computeEdgeGradTempBC;
    }
    else if((n = neighbour(i)))
    {
      if(i < 4 && !isBedCell && n->isBedCell)
      {
        computeEdgeTempDerivs[i] = &ElementCell::computeNeighborLayersTempBC;
      }
      else
      {
        computeEdgeTempDerivs[i] = &ElementCell::computeNeighborTempBC;
      }
    }
    else
    {
      computeEdgeTempDerivs[i] = &ElementCell::computeZeroTempBC;
    }
  }

}

void ElementCell::calculatePreComputedSoluteVariables()
{

  ElementCell *n = nullptr;
  double dryBulkDensity = (1.0 - porosity) * sedDensity;

  for(int i = 0; i < (int)parentElement->model->m_solutes.size(); i++)
  {
    for(int j = 0; j < 6; j++)
    {
      if(edgeSoluteConcs[i][j].isBC)
      {
        computeEdgeSoluteDerivs[i][j] = &ElementCell::computeEdgeSoluteBC;
      }
      else if(edgeGradSoluteConcs[i][i].isBC)
      {
        computeEdgeSoluteDerivs[i][j] = &ElementCell::computeEdgeGradSoluteBC;
      }
      else if((n = neighbour(j)))
      {
        if(i < 4)
        {
          if((!isBedCell && n->isBedCell) || (isBedCell && !n->isBedCell))
          {
            computeEdgeSoluteDerivs[i][j] = &ElementCell::computeZeroSoluteBC;
          }
          else
          {
            computeEdgeSoluteDerivs[i][j] = &ElementCell::computeNeighborSoluteBC;
          }
        }
        else
        {
          computeEdgeSoluteDerivs[i][j] = &ElementCell::computeNeighborSoluteBC;
        }
      }
      else
      {
        computeEdgeSoluteDerivs[i][j] = &ElementCell::computeZeroSoluteBC;
      }
    }


    retardationFactor[i] =  1.0 + dryBulkDensity * parentElement->model->m_solute_kd[i] / porosity;
  }
}

void ElementCell::computeMassBalance(double timeStep)
{
  totalMassBalance += (hydHead.value - prevHydHead.value) * width * parentElement->length / timeStep;
}

void ElementCell::computeHeatBalance(double timeStep)
{
  totalHeatBalance += rhom_Cm * (prevTemperature.value - temperature.value)  * volume / 1000.0;
}

void ElementCell::computeSoluteBalance(double timeStep, int soluteIndex)
{
  totalSoluteMassBalance[soluteIndex] =  (prevSoluteConcs[soluteIndex].value - soluteConcs[soluteIndex].value)  * volume * retardationFactor[soluteIndex];
}

void ElementCell::deleteSoluteVariables()
{
  if(soluteConcs)
  {
    delete[] soluteConcs; soluteConcs = nullptr;
    delete[] prevSoluteConcs; prevSoluteConcs = nullptr;
    delete[] externalSoluteFluxes; externalSoluteFluxes = nullptr;
    delete[] channelSoluteFlux; channelSoluteFlux = nullptr;
    delete[] channelSoluteRate; channelSoluteRate = nullptr;
    delete[] totalSoluteMassBalance; totalSoluteMassBalance = nullptr;

    for(int i = 0; i < parentElement->model->m_solutes.size(); i++)
    {
      delete[] edgeSoluteConcs[i];
      delete[] edgeSoluteConcFluxes[i];
      delete[] edgeGradSoluteConcs[i];
      delete[] computeEdgeSoluteDerivs[i];
    }

    delete[] computeEdgeSoluteDerivs; computeEdgeSoluteDerivs = nullptr;
    delete[] edgeSoluteConcs; edgeSoluteConcs = nullptr;
    delete[] edgeSoluteConcFluxes; edgeSoluteConcFluxes = nullptr;
    delete[] edgeGradSoluteConcs; edgeGradSoluteConcs = nullptr;
    delete[] retardationFactor; retardationFactor = nullptr;
  }
}

void ElementCell::computeEdgeHydCons()
{
  for(int i = 0; i < 6; i++)
  {
    int f = ef1[i];
    ElementCell *neigh = nullptr;

    if((neigh = neighbour(i)))
    {
      ElementCell *neighbor = neighbour(i);

      double flp = flowLengthP(i) * 2;
      double fln = flowLengthN(i, neigh) * 2;

      double hydConP = hydCon[f];
      double hydConN = neighbor->hydCon[f];
      double hydConLocal = hydConP > 0 && hydConN > 0 ? (flp +  fln) / ((flp / hydConP) + (fln / hydConN)) : hydCon[f];
      edgeHydCons[i] = hydConLocal;
    }
    else
    {
      edgeHydCons[i] = max(0.0, hydCon[f]);
    }
  }
}

void ElementCell::computeDepth()
{
  if(isTopCell())
    depth = max(hydHead.value - bottomElev, 0.0);
  else
    depth = topElev - bottomElev;

  prev_volume = volume;
  volume = depth * width * parentElement->length;

  dvolume_dt = (volume - prev_volume) / parentElement->model->m_timeStep;
}

double ElementCell::computeEdgeDepth(int edge, int layerIndex)
{
  double edgeDepth = 0.0;

  if(edgeHydHead->isBC)
  {
    if(isTopCell())
      edgeDepth = max(0.0, edgeHydHead->value - bottomElev);
    else
      edgeDepth = topElev - bottomElev;
  }
  else if(edgeGradHydHead->isBC)
  {
    if(isTopCell())
    {
      double edgeHead = edgeGradHydHead->value * flowLengthP(edge) + hydHead.value;
      edgeDepth = max(0.0, edgeHead - bottomElev);
    }
    else
    {
      edgeDepth = topElev - bottomElev;
    }
  }
  else
  {
    ElementCell *neigh = neighbour(edge, true);

    if(neigh)
    {
      neigh = neigh->layer(layerIndex) ;

      if(isBedCell && neigh->isBedCell)
      {
        if(isTopCell())
        {
          neigh = neigh->topBedCell;
          double fp = flowLengthP(edge);
          double fn = flowLengthN(edge, neigh);
          double ft = fp + fn;
          double interpHead =  hydHead.value ; // * fn / ft + neigh->hydHead.value * fp / ft;
          double interpBottom = bottomElev; // * fn / ft + neigh->bottomElev * fp / ft;
          edgeDepth = max(0.0, interpHead - interpBottom);
        }
        else
        {
          double fp = flowLengthP(edge);
          double fn = flowLengthN(edge, neigh);
          double ft = fp + fn;
          double interpHead =  topElev ; //* fn / ft + neigh->topElev * fp / ft;
          double interpBottom = bottomElev; // * fn / ft + neigh->bottomElev * fp / ft;
          edgeDepth = max(0.0, interpHead - interpBottom);
        }
      }
      else if(isBedCell && !neigh->isBedCell)
      {
        if(isTopCell())
        {
          double fp = flowLengthP(edge);
          double fn = flowLengthN(edge, neigh);
          double ft = fp + fn;
          double interpHead = hydHead.value; // * fn / ft + neigh->hydHead.value * fp / ft;
          double interpBottom = bottomElev;
          edgeDepth = max(0.0, interpHead - interpBottom);
        }
        else
        {
          double fp = flowLengthP(edge);
          double fn = flowLengthN(edge, neigh);
          double ft = fp + fn;
          double interpHead =  topElev;
          double interpBottom = bottomElev;

          //          if (zIndex ==  parentElement->model->m_numBedZCells - 1)
          //          {
          //            interpBottom = bottomElev * fn / ft + neigh->bottomElev * fp / ft;
          //          }

          edgeDepth = max(0.0, interpHead - interpBottom);
        }
      }
      else if(!isBedCell && neigh->isBedCell)
      {
        if(neigh->isTopCell())
        {
          double fp = flowLengthP(edge);
          double fn = flowLengthN(edge, neigh);
          double ft = fp + fn;
          double interpHead =  hydHead.value; // * fn/ft + neigh->hydHead.value * fp / ft;
          double interpBottom = neigh->bottomElev;
          edgeDepth = max(0.0, interpHead - interpBottom);
        }
        else
        {
          double fp = flowLengthP(edge);
          double fn = flowLengthN(edge, neigh);
          double ft = fp + fn;
          double interpHead = neigh->topElev;
          double interpBottom = neigh->bottomElev;
          edgeDepth = max(0.0, interpHead - interpBottom);
        }
      }
      else if(!isBedCell && !neigh->isBedCell)
      {
        double fp = flowLengthP(edge);
        double fn = flowLengthN(edge, neigh);
        double ft = fp + fn;
        double interpHead =  hydHead.value ; //* fn / ft + neigh->hydHead.value * fp / ft;
        double interpBottom = bottomElev; // * fn / ft + neigh->bottomElev * fp / ft;
        edgeDepth = max(0.0, interpHead - interpBottom);
      }
    }
    else
    {
      edgeDepth = depth;
    }
  }

  return edgeDepth;

}

void ElementCell::computeVolumeDerivative()
{
  //  prev_volume = volume;
  //  depth = max(hydHead.value - bottomElev, 0.0);
  //  volume = depth * width * parentElement->length;
  //  dvolume_dt = (volume - prev_volume) / parentElement->model->m_timeStep;
}

void ElementCell::computeEdgeDispersionCoefficients()
{
  rhom_Cm = porosity * parentElement->model->m_waterDensity * parentElement->model->m_cp +
            (1.0 - porosity) * sedDensity * sedCp;

  double kep = porosity * parentElement->model->m_waterThermalConductivity + (1.0 - porosity) * sedThermalConductivity;


  for(int i = 0; i < 6; i++)
  {
    int f = ef1[i];

    double disp = dispersivity[f];
    ElementCell *neighbor = neighbour(i);
    double farea = flowArea(i, zIndex);

    if(neighbor)
    {
      double ken = neighbor->porosity * parentElement->model->m_waterThermalConductivity + (1.0 - porosity) * neighbor->sedThermalConductivity;

      double flp = flowLengthP(i);
      double fln = flowLengthN(i, neighbor);

      edgeEffectiveKe[i] = kep > 0 && ken > 0 ?
                             (flp +  fln) / ((flp / kep) + (fln / ken)) : 0.0;

      double dise = neighbor->dispersivity[f];

      double combDisp = disp > 0 && dise > 0 ?
                          (flp +  fln) / ((flp / disp) + (fln / dise)) : 0.0;

      double vlp = volume / 2.0;
      double vln = neighbor->volume / 2.0;

      edgePorosity[i] = (porosity * vlp + neighbor->porosity * vln) /(vlp + vln);
      double vel = fabs(edgeFlows[i] / farea);
      edgeMechDispersionCoeff[i] = combDisp * fabs(vel / edgePorosity[i]);

    }
    else
    {
      edgeEffectiveKe[i] = kep;
      edgePorosity[i] = porosity;
      double vel = fabs(edgeFlows[i] / farea);
      edgeMechDispersionCoeff[i] = dispersivity[f] * fabs(vel / porosity);
    }
  }
}

void ElementCell::computeChannelMassFlux()
{
  channelInflow = 0.0;
  channelInflowFlux = 0.0;

  if(isTopCell())
  {
    double lower = fabs(centerY) - width / 2.0;
    double upper = fabs(centerY) + width / 2.0;

    wettedWidth = 0.0;

    if(parentElement->channelWidth / 2.0 > lower)
    {
      wettedWidth = min(upper, parentElement->channelWidth / 2.0) - lower;

      if(parentElement->channelBedHydCond > 0.0 && parentElement->hydConZ > 0.0)
      {
        double hydConZCom = parentElement->hydConZ;
        double dz = parentElement->channelWSE - hydHead.value;
        double gradH = dz / ((topElev - bottomElev) / 10.0);
        channelInflow = hydConZCom * gradH * wettedWidth * parentElement->length;
      }

      channelInflowFlux = channelInflow / (wettedWidth * parentElement->length);

    }
  }
}

void ElementCell::computeChannelHeatFlux()
{
  channelHeatFlux = 0.0;
  channelHeatRate = 0.0;

  if(isTopCell() && wettedWidth > 0)
  {

    if(channelInflow >= 0.0)
    {
      channelHeatRate = parentElement->model->m_waterDensity * parentElement->model->m_cp *
                        channelInflow * parentElement->channelTemperature;
    }
    else
    {
      channelHeatRate = parentElement->model->m_waterDensity * parentElement->model->m_cp *
                        channelInflow * temperature.value ;
    }

    if(parentElement->channelWSE > topElev)
    {

      double t = (topElev - bottomElev)/2.0;
      double gradT = (parentElement->channelTemperature - temperature.value) / t;
      double mech = edgeMechDispersionCoeff[4];
      double kep = porosity * parentElement->model->m_waterThermalConductivity + (1.0 - porosity) * sedThermalConductivity;
      double heatDisCoeff = mech *  parentElement->model->m_waterDensity * parentElement->model->m_cp * porosity + kep;
      channelHeatRate += heatDisCoeff * gradT * wettedWidth * parentElement->length;
    }

    channelHeatFlux = channelHeatRate /  (wettedWidth * parentElement->length);
  }
}

void ElementCell::computeChannelSoluteFlux(int soluteIndex)
{
  channelSoluteFlux[soluteIndex] = 0;
  channelSoluteRate[soluteIndex] = 0;

  if(isTopCell() && wettedWidth > 0)
  {

    if(channelInflow >= 0.0)
    {
      channelSoluteRate[soluteIndex] = channelInflow * parentElement->channelSoluteConcs[soluteIndex];
    }
    else
    {
      channelSoluteRate[soluteIndex] = channelInflow * soluteConcs[soluteIndex].value;
    }

    if(parentElement->channelWSE > topElev)
    {
      double t = (topElev - bottomElev)/4.0;
      double gradSolute = (parentElement->channelSoluteConcs[soluteIndex] - soluteConcs[soluteIndex].value) / t;
      double mech = edgeMechDispersionCoeff[4];
      double soluteDisCoeff = mech + parentElement->model->m_solute_molecular_diff[soluteIndex];
      channelSoluteRate[soluteIndex] += soluteDisCoeff *  gradSolute * wettedWidth * parentElement->length;
    }

    channelSoluteFlux[soluteIndex] = channelSoluteRate[soluteIndex] / (wettedWidth * parentElement->length);
  }
}

double ElementCell::computeEdgeHydHeadBC(int edgeIndex, double dt, double H[])
{
  double flp = flowLengthP(edgeIndex);
  double farea = flowArea(edgeIndex, zIndex);

  double gradH = (edgeHydHead[edgeIndex].value - H[index]) / (flp);
  double edgeFlow = edgeHydCons[edgeIndex] *  gradH * farea;

  edgeFlows[edgeIndex] = edgeFlow;
  edgeGradHydHead[edgeIndex].value = gradH;

  return edgeFlow;
}

double ElementCell::computeEdgeGradHydHeadBC(int edgeIndex, double dt, double H[])
{
  double edgeFlow = edgeHydCons[edgeIndex] * edgeGradHydHead[edgeIndex].value * flowArea(edgeIndex, zIndex);
  edgeFlows[edgeIndex] = edgeFlow;
  return edgeFlow;
}

double ElementCell::computeNeighborHydHeadBC(int edgeIndex, double dt, double H[])
{

  ElementCell *neigh = neighbour(edgeIndex);
  int nIndex = neigh->index;

  double flp = flowLengthP(edgeIndex);
  double fln = flowLengthN(edgeIndex, neigh);
  double farea = flowArea(edgeIndex, zIndex);

  double gradH = (H[nIndex] - H[index]) / (flp + fln);
  double edgeFlow = edgeHydCons[edgeIndex] * gradH *  farea;

  edgeFlows[edgeIndex] = edgeFlow;
  edgeGradHydHead[edgeIndex].value = gradH;

  return edgeFlow;
}

double ElementCell:: computeNeighborLayersHydHeadBC(int edgeIndex, double dt, double H[])
{
  ElementCell *neigh = neighbour(edgeIndex);

  double gradH = 0.0;
  double edgeFlow = 0.0;

  for(int i = 0; i < parentElement->model->m_numBedZCells; i++)
  {
    ElementCell *lay = neigh->topBedCell->bedCells[i];
    int nIndex = lay->index;

    double flp = flowLengthP(edgeIndex);
    double fln = flowLengthN(edgeIndex, lay);
    double farea = flowArea(edgeIndex, zIndex);

    double gradHT = (H[nIndex] - H[index]) / (flp + fln);
    double flow = edgeHydCons[edgeIndex] * gradHT *  farea;
    edgeFlow += flow;
    edgeFlowsComp[edgeIndex][i] = flow;
    gradH += gradHT;
  }

  edgeFlows[edgeIndex] = edgeFlow;
  edgeGradHydHead[edgeIndex].value = gradH / ( parentElement->model->m_numBedZCells * 1.0);

  return edgeFlow;
}

double ElementCell::computeZeroHydHeadBC(int edgeIndex, double dt, double H[])
{
  edgeFlows[edgeIndex] = 0.0;
  edgeGradHydHead[edgeIndex].value = 0.0;
  return 0.0;
}

double ElementCell::computeEdgeTempBC(int edgeIndex, double dt, double T[])
{
  double gradT = (edgeTemperatures[edgeIndex].value - T[index]) / flowLengthP(edgeIndex);
  double heatDisCoeff = edgeMechDispersionCoeff[edgeIndex] * parentElement->model->m_waterDensity *
                        parentElement->model->m_cp * edgePorosity[edgeIndex] + edgeEffectiveKe[edgeIndex];

  edgeGradTemperatures[edgeIndex].value = gradT;

  double edgeFlow = heatDisCoeff *  gradT * flowArea(edgeIndex, zIndex);
  edgeHeatFluxes[edgeIndex] += edgeFlow;

  return edgeFlow;
}

double ElementCell::computeEdgeGradTempBC(int edgeIndex, double dt, double T[])
{
  double heatDisCoeff = edgeMechDispersionCoeff[edgeIndex] * parentElement->model->m_waterDensity *
                        parentElement->model->m_cp * edgePorosity[edgeIndex] + edgeEffectiveKe[edgeIndex];

  double edgeFlow = heatDisCoeff *  edgeGradTemperatures[edgeIndex].value * flowArea(edgeIndex, zIndex);
  edgeHeatFluxes[edgeIndex] += edgeFlow;
  return edgeFlow;
}

double ElementCell::computeNeighborTempBC(int edgeIndex, double dt, double T[])
{
  ElementCell *n = neighbour(edgeIndex);

  double flp = flowLengthP(edgeIndex);
  double fln = flowLengthN(edgeIndex, n);

  double gradT = (T[n->index] - T[index]) / (flp + fln);

  double heatDisCoeff = edgeMechDispersionCoeff[edgeIndex] * parentElement->model->m_waterDensity *
                        parentElement->model->m_cp * edgePorosity[edgeIndex] + edgeEffectiveKe[edgeIndex];

  edgeGradTemperatures[edgeIndex].value = gradT;
  double edgeFlow = heatDisCoeff * gradT *  flowArea(edgeIndex, zIndex);
  edgeHeatFluxes[edgeIndex] += edgeFlow;

  return edgeFlow;
}

double ElementCell::computeNeighborLayersTempBC(int edgeIndex, double dt, double T[])
{
  ElementCell *n = neighbour(edgeIndex);

  double edgeFlow = 0;
  double gradT = 0.0;

  for(int i = 0; i < parentElement->model->m_numBedZCells; i++)
  {
    ElementCell *lay = n->topBedCell->bedCells[i];

    double flp = flowLengthP(edgeIndex);
    double fln = flowLengthN(edgeIndex, lay);
    double farea = flowArea(edgeIndex, i);

    double gradTT = (T[lay->index] - T[index]) / (flp + fln);

    double heatDisCoeff = edgeMechDispersionCoeff[edgeIndex] * parentElement->model->m_waterDensity *
                          parentElement->model->m_cp * edgePorosity[edgeIndex] + edgeEffectiveKe[edgeIndex];

    edgeFlow += heatDisCoeff * gradTT *  farea;
    gradT += gradTT;
  }

  edgeGradTemperatures[edgeIndex].value = gradT / (parentElement->model->m_numBedZCells * 1.0);
  edgeHeatFluxes[edgeIndex] += edgeFlow;

  return edgeFlow;
}

double ElementCell::computeZeroTempBC(int edgeIndex, double dt, double T[])
{
  edgeGradTemperatures[edgeIndex].value = 0.0;

  return 0.0;
}

double ElementCell::computeEdgeSoluteBC(int edgeIndex, double dt, double S[], int soluteIndex)
{
  double gradSolute = (edgeSoluteConcs[soluteIndex][edgeIndex].value - S[index]) / flowLengthP(edgeIndex);
  double heatDisCoeff = edgeMechDispersionCoeff[edgeIndex] + parentElement->model->m_solute_molecular_diff[soluteIndex];

  edgeGradSoluteConcs[soluteIndex][edgeIndex].value = gradSolute;
  double edgeFlow = heatDisCoeff *  gradSolute * flowArea(edgeIndex, zIndex);

  edgeSoluteConcFluxes[soluteIndex][edgeIndex] += edgeFlow;

  return edgeFlow;
}

double ElementCell::computeEdgeGradSoluteBC(int edgeIndex, double dt, double S[], int soluteIndex)
{
  double heatDisCoeff = edgeMechDispersionCoeff[edgeIndex] + parentElement->model->m_solute_molecular_diff[soluteIndex];
  double edgeFlow = heatDisCoeff *  edgeGradSoluteConcs[soluteIndex][edgeIndex].value * flowArea(edgeIndex, zIndex);

  edgeSoluteConcFluxes[soluteIndex][edgeIndex] += edgeFlow;

  return edgeFlow;
}

double ElementCell::computeNeighborSoluteBC(int edgeIndex, double dt, double S[], int soluteIndex)
{
  ElementCell *n = neighbour(edgeIndex);

  double flp = flowLengthP(edgeIndex);
  double fln = flowLengthN(edgeIndex, n);

  double gradSolute = (S[n->index] - S[index]) / (flp + fln);

  double solDisCoeff = edgeMechDispersionCoeff[edgeIndex] +
                       parentElement->model->m_solute_molecular_diff[soluteIndex];


  edgeGradSoluteConcs[soluteIndex][edgeIndex].value = gradSolute;
  double edgeFlow = solDisCoeff * gradSolute * flowArea(edgeIndex,zIndex);

  edgeSoluteConcFluxes[soluteIndex][edgeIndex] += edgeFlow;

  return edgeFlow;
}

double ElementCell::computeNeighborLayersSoluteBC(int edgeIndex, double dt, double S[], int soluteIndex)
{

  ElementCell *n = neighbour(edgeIndex);

  double edgeFlow = 0;
  double gradSolute = 0.0;

    for(int i = 0; i < parentElement->model->m_numBedZCells; i++)
    {
      ElementCell *lay = n->topBedCell->bedCells[i];

      double flp = flowLengthP(edgeIndex);
      double fln = flowLengthN(edgeIndex, lay);
      double farea = flowArea(edgeIndex, i);

      double gradSoluteSS = (S[lay->index] - S[index]) / (flp + fln);

      double solDisCoeff = edgeMechDispersionCoeff[edgeIndex] + parentElement->model->m_solute_molecular_diff[soluteIndex];

      edgeFlow += solDisCoeff * gradSoluteSS *  farea;
      gradSolute += gradSoluteSS;
    }

  edgeSoluteConcFluxes[soluteIndex][edgeIndex] += edgeFlow;
  edgeGradSoluteConcs[soluteIndex][edgeIndex].value = gradSolute / (parentElement->model->m_numBedZCells * 1.0);

  return edgeFlow;
}

double ElementCell::computeZeroSoluteBC(int edgeIndex, double dt, double S[], int soluteIndex)
{
  edgeGradSoluteConcs[soluteIndex][edgeIndex].value = 0.0;
  return 0.0;
}



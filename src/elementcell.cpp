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

const int ElementCell::nEdgeIndex[] = {2, 3, 0, 1};
const double ElementCell::ef1[] = {0.0, 1.0, 0.0, 1.0};
const double ElementCell::ef2[] = {1.0, 0.0, 1.0, 0.0};

ElementCell::ElementCell(int cindex, Element *cparent)
  : elementCellIndex(cindex),
    edgeHydHead(nullptr),
    edgeGradHydHead(nullptr),
    soluteConcs(nullptr),
    prevSoluteConcs(nullptr),
    externalSoluteFluxes(nullptr),
    channelSoluteFlux(nullptr),
    totalSoluteMassBalance(nullptr),
    edgeFlows(nullptr),
    edgeHydCons(nullptr),
    parentElement(cparent),
    retardationFactor(nullptr)
{

  edgeBedRockElevs = new double[4]();
  edgeEffectiveKe = new double[4]();
  edgePorosity = new double[4]();
  hydCon = new double[2]();
  hydCon[0] = parentElement->model->m_hydConY;
  hydCon[1] = parentElement->model->m_hydConX;
  sedDensity = parentElement->model->m_sedDensity;
  sedCp = parentElement->model->m_sedCp;
  porosity = parentElement->model->m_porosity;
  specificStorage = parentElement->model->m_specificStorage;
  width = parentElement->model->m_defaultCellWidth;
  dispersivity = new double[2]();
  dispersivity[0] = parentElement->model->m_dispersivityY;
  dispersivity[1] = parentElement->model->m_dispersivityX;
  sedThermalConductivity = parentElement->model->m_sedThermalConductivity;

  edgeHydHead = new Variable[4];
  edgeGradHydHead = new Variable[4];
  edgeFlows = new double[4]();
  edgeHydCons = new double[4]();
  edgeMechDispersionCoeff = new double[4]();
  edgeHeatPecletNumbers = new double[4]();
  edgeSolutePecletNumbers = new double[4]();
  edgeTemperatures = new Variable[4];
  edgeGradTemperatures = new Variable[4];
  neighbors = new ElementCell*[4];
  computeEdgeHeadDerivs = new ComputeEdgeDeriv[4];
  computeEdgeTempDerivs = new ComputeEdgeDeriv[4];
  edgeDepths = new double[4]();
  flowLengthP = new double[4]();
  flowLengthN = new double[4]();
  flowWidth = new double[4]();
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
  delete[] edgeMechDispersionCoeff;
  delete[] edgeSolutePecletNumbers;
  delete[] neighbors;
  delete[] computeEdgeHeadDerivs;
  delete[] computeEdgeTempDerivs;
  delete[] edgeDepths;
  delete[] flowLengthP;
  delete[] flowLengthN;
  delete[] flowWidth;
  delete[] edgeBedRockElevs;

  deleteSoluteVariables();
}

void ElementCell::initialize()
{

  start = true;
  dvolume_dt = 0.0;

  totalMassBalance = 0.0;
  totalHeatBalance = 0.0;
  channelInflow = 0.0;
  channelHeatFlux = 0.0;
  externalHeatFluxes = 0.0;
  externalInflow = 0.0;
  depth = 0;


  int add[4] = {1,1,-1,-1};

  //set neighbours
  for(int i = 0; i < 4 ; i++)
  {
    neighbors[i] = nullptr;
    edgeBedRockElevs[i] = bedRockElev;

    if(ef1[i] == 0)
    {
      flowLengthP[i] = width;
      flowWidth[i] = parentElement->length;

      int nindex = elementCellIndex + add[i];

      if(nindex >= 0 && nindex < parentElement->model->m_totalCellsPerElement)
      {
        ElementCell *neighbor = parentElement->elementCells[nindex];
        flowLengthN[i] = neighbor->width;
        neighbors[i] = neighbor;

        double flp = width / 2.0;
        double fln = neighbor->width / 2.0;
        edgeBedRockElevs[i] = (bedRockElev / flp +  neighbor->bedRockElev / fln) /
                              ((1.0 / flp) + (1.0 / fln));
      }
    }
    else
    {

      flowLengthP[i] = parentElement->length;
      flowWidth[i] = width;

      int nindex = parentElement->index + add[i];

      if(nindex >= 0 && nindex < (int)parentElement->model->m_elements.size())
      {
        Element *element = parentElement->model->m_elements[nindex];
        ElementCell *neighbor = element->elementCells[elementCellIndex];
        flowLengthN[i] = neighbor->parentElement->length;
        neighbors[i] = neighbor;

        double flp = parentElement->length / 2.0;
        double fln = neighbor->parentElement->length / 2.0;

        edgeBedRockElevs[i] = (bedRockElev / flp +  neighbor->bedRockElev / fln) /
                              ((1.0 / flp) + (1.0 / fln));
      }
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
    totalSoluteMassBalance = new double[numSolutes]();
    computeEdgeSoluteDerivs = new ComputeEdgeVariableDeriv *[numSolutes];

    edgeSoluteConcs = new Variable*[numSolutes];
    edgeGradSoluteConcs = new Variable*[numSolutes];
    retardationFactor = new double[numSolutes]();

    for(int i = 0; i < numSolutes; i++)
    {
      edgeSoluteConcs[i] = new Variable[4];
      edgeGradSoluteConcs[i] = new Variable[4];
      computeEdgeSoluteDerivs[i] = new ComputeEdgeVariableDeriv[4];

    }
  }
}

double ElementCell::computeDHydHeadDt(double dt, double H[])
{
  double DHydHeadDt = 0.0;

  for(int i = 0; i < 4; i++)
  {
    DHydHeadDt += (this->*computeEdgeHeadDerivs[i])(i, dt, H);
  }

  qY = (-edgeFlows[0] + edgeFlows[2]) / 2.0;
  qX = (-edgeFlows[1] + edgeFlows[3]) / 2.0;

  //  if(dt > 0)
  {
    dvolume_dt = ((max(H[index] - bedRockElev, 0.0) * parentElement->length * width)
                  - volume) / parentElement->model->m_timeStep;
  }

  //compute river inflow
  DHydHeadDt += channelInflow;
  DHydHeadDt += externalInflow;
  DHydHeadDt -= specificStorage * hydHead.value * dvolume_dt / depth;
  DHydHeadDt /= (specificStorage * parentElement->length * width);

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
    DTDt += channelHeatFlux / (rhom_Cm * volume);
  }

  //Product rule subtract volume derivative
  {
    DTDt -= T[index] * dvolume_dt / volume;
  }

  return DTDt;
}

double ElementCell::computeDTDtDispersion(double dt, double T[])
{
  double DTDt = 0.0;

  for(int i = 0; i < 4; i++)
  {
    DTDt += (this->*computeEdgeTempDerivs[i])(i, dt, T);
  }

  DTDt = DTDt / (rhom_Cm * volume);

  return DTDt;
}

double ElementCell::computeDTDtUpwind(double dt, double T[])
{
  double DTDt = 0.0;

  for(int i = 0; i < 4; i++)
  {
    double flow = edgeFlows[i];

    if(flow >= 0.0)
    {
      if(edgeTemperatures[i].isBC)
      {
        DTDt += parentElement->model->m_waterDensity * parentElement->model->m_cp * flow * edgeTemperatures[i].value;
      }
      else if(neighbors[i])
      {
        DTDt += parentElement->model->m_waterDensity * parentElement->model->m_cp * flow * T[neighbors[i]->index];
      }
    }
    else
    {
      DTDt +=  parentElement->model->m_waterDensity * parentElement->model->m_cp * flow * T[index];
    }
  }

  DTDt = DTDt / (rhom_Cm * volume);

  return DTDt;
}

double ElementCell::computeDTDtCentral(double dt, double T[])
{
  double DTDt = 0.0;

  for(int i = 0; i < 4; i++)
  {
    double edgeFlow = edgeFlows[i];

    if(edgeTemperatures[i].isBC)
    {
      DTDt +=  parentElement->model->m_waterDensity * parentElement->model->m_cp *
               edgeFlow * edgeTemperatures[i].value;
    }
    else if(neighbors[i])
    {
      double flp = flowLengthP[i] / 2.0;
      double fln = flowLengthN[i] / 2.0;

      double edgeT = (T[index] / flp   +  T[neighbors[i]->index] / fln) / (1.0 / flp + 1.0 / fln);

      DTDt += parentElement->model->m_waterDensity * parentElement->model->m_cp * edgeFlow * edgeT;
    }
    else if(edgeFlow < 0)
    {
      DTDt +=  parentElement->model->m_waterDensity * parentElement->model->m_cp * edgeFlow * T[index];
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

  for(int i = 0; i < 4; i++)
  {
    double flow = edgeFlows[i];

    if(flow > 0.0)
    {
      if(edgeTemperatures[i].isBC)
      {
        DTDt +=  flow * edgeTemperatures[i].value;
      }
      else if(neighbors[i])
      {
        DTDt +=  flow * T[neighbors[i]->index];
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

    //subtract chain rule volume derivative
    {
      DSoluteDt -= (S[index]  * dvolume_dt) / (retardationFactor[soluteIndex] * volume);
    }

    //First order reaction
    {
      DSoluteDt -= (soluteConcs[soluteIndex].value * parentElement->model->m_solute_first_order_k[soluteIndex]) / (retardationFactor[soluteIndex]);
    }

    //Add channel solute
    {
      DSoluteDt += channelSoluteFlux[soluteIndex] / (retardationFactor[soluteIndex] * volume);
    }

    //Add external sources
    {
      DSoluteDt += externalSoluteFluxes[soluteIndex] / (retardationFactor[soluteIndex] * volume);
    }
  }

  return DSoluteDt;
}

double ElementCell::computeDSoluteDtDispersion(double dt, double S[], int soluteIndex)
{
  double DSoluteDt = 0.0;

  for(int i = 0; i < 4; i++)
  {
    DSoluteDt += (this->*computeEdgeSoluteDerivs[soluteIndex][i])(i, dt, S, soluteIndex);
  }

  DSoluteDt = DSoluteDt / (retardationFactor[soluteIndex] * volume);

  return DSoluteDt;
}

double ElementCell::computeDSoluteDtUpwind(double dt, double S[], int soluteIndex)
{
  double DSoluteDt = 0.0;

  for(int i = 0; i < 4; i++)
  {
    double flow = edgeFlows[i];

    if(flow > 0.0)
    {
      if(edgeSoluteConcs[soluteIndex][i].isBC)
      {
        DSoluteDt +=  flow * edgeSoluteConcs[soluteIndex][i].value;
      }
      else if(neighbors[i])
      {
        DSoluteDt +=  flow * S[neighbors[i]->index];
      }
    }
    else
    {
      DSoluteDt += flow * S[index];
    }
  }

  DSoluteDt = DSoluteDt / (retardationFactor[soluteIndex] * volume * porosity);

  return DSoluteDt;
}

double ElementCell::computeDSoluteDtCentral(double dt, double S[], int soluteIndex)
{
  double DSoluteDt = 0.0;

  for(int i = 0; i < 4; i++)
  {
    double edgeFlow = edgeFlows[i];

    if(edgeSoluteConcs[soluteIndex][i].isBC)
    {
      DSoluteDt +=  edgeFlow * edgeSoluteConcs[soluteIndex][i].value;
    }
    else if(neighbors[i])
    {
      double flp = flowLengthP[i] / 2.0;
      double fln = flowLengthN[i] / 2.0;

      DSoluteDt += edgeFlow * (S[index] / flp   +  S[neighbors[i]->index] / fln) /
          (1.0 / flp + 1.0 / fln);
    }
    else if(edgeFlow < 0.0)
    {
      DSoluteDt +=  edgeFlow * S[index];
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
  double alpha = hydCon[0] * depth / (specificStorage * width * width);
  double beta  = hydCon[1] * depth / (specificStorage * parentElement->length * parentElement->length);
  double factor = max(alpha, beta);

  return factor ;
}

double ElementCell::computeDiffusionFactor() const
{
  //Change to solute and heat dispersion based on (Woods, Teubner, Simmons, & Narayan, 2003)

  double diffFactor = 0;

  for(int i = 0; i < 4 ; i++)
  {
    double flp = flowLengthP[i];
    diffFactor = max(diffFactor, ((edgeEffectiveKe[i] / rhom_Cm) + (edgeMechDispersionCoeff[i] / (parentElement->model->m_waterDensity * parentElement->model->m_cp)))
                     / (flp * flp));
  }

  return diffFactor;
}

void ElementCell::calculatePreComputedHydHeadVariables()
{
  for(int i = 0; i < 4; i++)
  {
    if(edgeHydHead[i].isBC)
    {
      computeEdgeHeadDerivs[i] = &ElementCell::computeEdgeHydHeadHeadBC;
    }
    else if(edgeGradHydHead[i].isBC)
    {
      computeEdgeHeadDerivs[i] = &ElementCell::computeEdgeGradHydHeadHeadBC;
    }
    else if(neighbors[i])
    {
      computeEdgeHeadDerivs[i] = &ElementCell::computeNeighborHydHeadBC;
    }
    else
    {
      computeEdgeHeadDerivs[i] = &ElementCell::computeZeroHydHeadBC;
    }
  }

  computeEdgeHydCons();

  computeEdgeDepths();

  computeEdgeDispersionCoefficients();

  computeChannelMassFlux();

}

void ElementCell::calculatePreComputedTempVariables()
{
  for(int i = 0; i < 4; i++)
  {
    if(edgeTemperatures[i].isBC)
    {
      computeEdgeTempDerivs[i] = &ElementCell::computeEdgeTempBC;
    }
    else if(edgeGradTemperatures[i].isBC)
    {
      computeEdgeTempDerivs[i] = &ElementCell::computeEdgeGradTempBC;
    }
    else if(neighbors[i])
    {
      computeEdgeTempDerivs[i] = &ElementCell::computeNeighborTempBC;
    }
    else
    {
      computeEdgeTempDerivs[i] = &ElementCell::computeZeroTempBC;
    }
  }

  computeChannelHeatFlux();
}

void ElementCell::calculatePreComputedSoluteVariables()
{

  double dryBulkDensity = (1.0 - porosity) * sedDensity;

  for(int i = 0; i < parentElement->model->m_solutes.size(); i++)
  {
    for(int j = 0; j < 4; j++)
    {
      if(edgeSoluteConcs[i][j].isBC)
      {
        computeEdgeSoluteDerivs[i][j] = &ElementCell::computeEdgeSoluteBC;
      }
      else if(edgeGradSoluteConcs[i][j].isBC)
      {
        computeEdgeSoluteDerivs[i][j] = &ElementCell::computeEdgeGradSoluteBC;
      }
      else if(neighbors[j])
      {
        computeEdgeSoluteDerivs[i][j] = &ElementCell::computeNeighborSoluteBC;
      }
      else
      {
        computeEdgeSoluteDerivs[i][j] = &ElementCell::computeZeroSoluteBC;
      }
    }

    retardationFactor[i] =  1.0 + dryBulkDensity * parentElement->model->m_solute_kd[i] / porosity;

    computeChannelSoluteFlux(i);
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

}

void ElementCell::deleteSoluteVariables()
{
  if(soluteConcs)
  {
    delete[] soluteConcs; soluteConcs = nullptr;
    delete[] prevSoluteConcs; prevSoluteConcs = nullptr;
    delete[] externalSoluteFluxes; externalSoluteFluxes = nullptr;
    delete[] channelSoluteFlux; channelSoluteFlux = nullptr;
    delete[] totalSoluteMassBalance; totalSoluteMassBalance = nullptr;

    for(int i = 0; i < parentElement->model->m_solutes.size(); i++)
    {
      delete[] edgeSoluteConcs[i];
      delete[] edgeGradSoluteConcs[i];
      delete[] computeEdgeSoluteDerivs[i];
    }

    delete[] computeEdgeSoluteDerivs; computeEdgeSoluteDerivs = nullptr;
    delete[] edgeSoluteConcs; edgeSoluteConcs = nullptr;
    delete[] edgeGradSoluteConcs; edgeGradSoluteConcs = nullptr;
    delete[] retardationFactor; retardationFactor = nullptr;
  }
}

void ElementCell::computeEdgeHydCons()
{
  for(int i = 0; i < 4; i++)
  {
    int f = ef1[i];

    if(neighbors[i])
    {
      ElementCell *neighbor = neighbors[i];

      double flp = flowLengthP[i] / 2.0;
      double fln = flowLengthN[i] / 2.0;

      double hydConP = hydCon[f];
      double hydConN = neighbor->hydCon[f];

      edgeHydCons[i] = hydConP > 0 && hydConN > 0 ?
                         (flp +  fln) / ((flp / hydConP) + (fln / hydConN)) : 0.0;
    }
    else
    {
      edgeHydCons[i] = max(0.0, hydCon[f]);
    }
  }
}

void ElementCell::computeEdgeDepths()
{
  depth = max(hydHead.value - bedRockElev, 0.0);
  volume = depth * width * parentElement->length;

  //calculate edge bottom elevations
  for(int i = 0; i < 4; i++)
  {
    if(!edgeHydHead[i].isBC)
    {
      edgeHydHead[i].value = hydHead.value;
    }

    edgeDepths[i] = max(0.0, edgeHydHead[i].value - bedRockElev);
  }
}

void ElementCell::computeVolumeDerivative()
{
  prev_volume = volume;
  depth = max(hydHead.value - bedRockElev, 0.0);
  volume = depth * width * parentElement->length;
  dvolume_dt = (volume - prev_volume) / parentElement->model->m_timeStep;
}

void ElementCell::computeEdgeDispersionCoefficients()
{
  rhom_Cm = porosity * parentElement->model->m_waterDensity * parentElement->model->m_cp +
            (1.0 - porosity) * sedDensity * sedCp;

  double kep = porosity * parentElement->model->m_waterThermalConductivity + (1.0 - porosity) * sedThermalConductivity;


  for(int i = 0; i < 4; i++)
  {
    int f = ef1[i];

    double disp = dispersivity[f];

    if(neighbors[i])
    {
      ElementCell *neighbor = neighbors[i];

      double ken = neighbor->porosity * parentElement->model->m_waterThermalConductivity + (1.0 - porosity) * neighbor->sedThermalConductivity;

      double flp = flowLengthP[i] / 2.0;
      double fln = flowLengthN[i] / 2.0;

      edgeEffectiveKe[i] = kep > 0 && ken > 0 ?
                             (flp +  fln) / ((flp / kep) + (fln / ken)) : 0.0;

      double dise = neighbor->dispersivity[f];

      double combDisp = disp > 0 && dise > 0 ?
                          (flp +  fln) / ((flp / disp) + (fln / dise)) : 0.0;

      double vlp = volume / 2.0;
      double vln = neighbor->volume / 2.0;

      edgePorosity[i] = (porosity * vlp + neighbor->porosity * vln) /(vlp + vln);
      edgeMechDispersionCoeff[i] = combDisp * fabs(edgeGradHydHead[i].value * edgeHydCons[i] / edgePorosity[i]);

    }
    else
    {
      edgeEffectiveKe[i] = kep;
      edgePorosity[i] = porosity;
      edgeMechDispersionCoeff[i] = dispersivity[f] * fabs(edgeGradHydHead[i].value * edgeHydCons[i] / porosity);
    }
  }
}

void ElementCell::computeEdgeDepths(double H[])
{
  depth = max(H[index] - bedRockElev, 0.0);
  volume = depth * width * parentElement->length;

  for(int i = 0; i < 4; i++)
  {
    if(!edgeHydHead[i].isBC)
    {
      edgeHydHead[i].value = H[index];
    }

    edgeDepths[i] = max(0.0, edgeHydHead[i].value - bedRockElev);
  }
}

void ElementCell::computeChannelMassFlux()
{
  channelInflow = 0.0;
  channelInflowFlux = 0.0;

  double lower = fabs(centerY) - width / 2.0;
  double upper = fabs(centerY) + width / 2.0;
  double wettedWidth = 0.0;

  if(parentElement->channelWidth / 2.0 > lower)
  {
    wettedWidth = min(upper, parentElement->channelWidth / 2.0) - lower;

    if(parentElement->channelBedHydCond > 0.0 && parentElement->hydConZ > 0.0)
    {

      double hydConZCom = ((parentElement->channelBedThickness) / 2.0 + (topElev - bedRockElev) / 2.0) /
                          ((parentElement->channelBedThickness) / 2.0 / parentElement->channelBedHydCond +
                           (topElev - bedRockElev) / 2.0 / parentElement->hydConZ);

      if(hydHead.value >= topElev)
      {
        double dz = parentElement->channelWSE - hydHead.value;

        if((hydHead.value > parentElement->channelWSE) ||
           (parentElement->channelWSE > hydHead.value && parentElement->channelWSE - topElev > 0))
        {
          double cond = hydConZCom / (parentElement->channelBedThickness / 2.0 + max(0.0, topElev - bedRockElev) / 2.0);
          channelInflow = cond * wettedWidth * parentElement->length * dz;
        }

        double maxFlow = dz * wettedWidth * parentElement->length / parentElement->model->m_timeStep;

        if(dz >= 0)
        {
          channelInflow = min(channelInflow, maxFlow * 0.99);
        }
        else
        {
          channelInflow = max(channelInflow, maxFlow * 0.99);
        }
      }
      else if(parentElement->channelWSE - topElev > 0.0)
      {
        double dz = parentElement->channelWSE - topElev ;
        double maxFlow = dz * wettedWidth * parentElement->length / parentElement->model->m_timeStep;

        double cond = hydConZCom / (parentElement->channelBedThickness / 2.0 + max(0.0, topElev - bedRockElev) / 2.0);
        channelInflow = min({hydConZCom * wettedWidth * parentElement->length,
                             cond * wettedWidth * parentElement->length * dz,
                             maxFlow * 0.99});
      }
    }

    channelInflowFlux = channelInflow / (wettedWidth * parentElement->length);
  }

  parentElement->channelInflow += channelInflow;
  parentElement->channelInflowFlux += parentElement->channelWidth > 0 ? channelInflowFlux * wettedWidth / parentElement->channelWidth : 0.0;
}

void ElementCell::computeChannelHeatFlux()
{
  channelHeatFlux = 0.0;

  if(channelInflow >= 0.0)
  {
    channelHeatFlux = parentElement->model->m_waterDensity * parentElement->model->m_cp *
                      channelInflow * parentElement->channelTemperature;
  }
  else
  {
    channelHeatFlux = parentElement->model->m_waterDensity * parentElement->model->m_cp *
                      channelInflow * temperature.value;
  }

  parentElement->channelHeatFlux += channelHeatFlux;
}

void ElementCell::computeChannelSoluteFlux(int soluteIndex)
{
  channelSoluteFlux[soluteIndex] = 0.0;

  if(channelInflow >= 0.0)
  {
    channelSoluteFlux[soluteIndex] = channelInflow * parentElement->channelSoluteConcs[soluteIndex];
  }
  else
  {
    channelSoluteFlux[soluteIndex] = channelInflow * soluteConcs[soluteIndex].value;
  }

  parentElement->channelSoluteFlux[soluteIndex] += channelSoluteFlux[soluteIndex];

}

double ElementCell::computeEdgeHydHeadHeadBC(int edgeIndex, double dt, double H[])
{
  double flp = flowLengthP[edgeIndex];
  double fwp = flowWidth[edgeIndex];

  double gradH = (edgeHydHead[edgeIndex].value - H[index]) / (flp / 2.0);
  double edgeFlow = edgeHydCons[edgeIndex] *  gradH * edgeDepths[edgeIndex] * fwp;

  edgeFlows[edgeIndex] = edgeFlow;
  edgeGradHydHead[edgeIndex].value = gradH;

  return edgeFlow;
}

double ElementCell::computeEdgeGradHydHeadHeadBC(int edgeIndex, double dt, double H[])
{
  double edgeFlow = edgeHydCons[edgeIndex] * edgeGradHydHead[edgeIndex].value * edgeDepths[edgeIndex] * flowWidth[edgeIndex];
  edgeFlows[edgeIndex] = edgeFlow;
  return edgeFlow;
}

double ElementCell::computeNeighborHydHeadBC(int edgeIndex, double dt, double H[])
{
  int nIndex = neighbors[edgeIndex]->index;

  double flp = flowLengthP[edgeIndex] / 2.0;
  double fln = flowLengthN[edgeIndex] / 2.0;
  double fwp = flowWidth[edgeIndex];

  double gradH = (H[nIndex] - H[index]) / (flp + fln);
  double edgeFlow = edgeHydCons[edgeIndex] * gradH *  edgeDepths[edgeIndex] * fwp;
  //  double edgeFlow = edgeHydCons[edgeIndex] * gradH *  depth * fwp;

  edgeFlows[edgeIndex] = edgeFlow;
  edgeGradHydHead[edgeIndex].value = gradH;

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
  double gradT = (edgeTemperatures[edgeIndex].value - T[index]) / (flowLengthP[edgeIndex] / 2.0);
  double heatDisCoeff = edgeMechDispersionCoeff[edgeIndex] * parentElement->model->m_waterDensity *
                        parentElement->model->m_cp * edgePorosity[edgeIndex] + edgeEffectiveKe[edgeIndex];

  edgeGradTemperatures[edgeIndex].value = gradT;
  double edgeFlow = heatDisCoeff *  gradT * edgeDepths[edgeIndex] * flowWidth[edgeIndex] ;
  return edgeFlow;
}

double ElementCell::computeEdgeGradTempBC(int edgeIndex, double dt, double T[])
{
  double heatDisCoeff = edgeMechDispersionCoeff[edgeIndex] * parentElement->model->m_waterDensity *
                        parentElement->model->m_cp * edgePorosity[edgeIndex] + edgeEffectiveKe[edgeIndex];

  double edgeFlow = heatDisCoeff *  edgeGradTemperatures[edgeIndex].value * edgeDepths[edgeIndex] * flowWidth[edgeIndex] ;
  return edgeFlow;
}

double ElementCell::computeNeighborTempBC(int edgeIndex, double dt, double T[])
{
  ElementCell *n = neighbors[edgeIndex];

  double flp = flowLengthP[edgeIndex] / 2.0;
  double fln = flowLengthN[edgeIndex] / 2.0;

  double gradT = (T[n->index] - T[index]) / (flp + fln);

  double heatDisCoeff = edgeMechDispersionCoeff[edgeIndex] * parentElement->model->m_waterDensity *
                        parentElement->model->m_cp * edgePorosity[edgeIndex] + edgeEffectiveKe[edgeIndex];

  edgeGradTemperatures[edgeIndex].value = gradT;
  double edgeFlow = heatDisCoeff * gradT *  edgeDepths[edgeIndex] * flowWidth[edgeIndex];

  return edgeFlow;
}

double ElementCell::computeZeroTempBC(int edgeIndex, double dt, double T[])
{
  edgeGradTemperatures[edgeIndex].value = 0.0;
  return 0.0;
}

double ElementCell::computeEdgeSoluteBC(int edgeIndex, double dt, double S[], int soluteIndex)
{
  double gradSolute = (edgeSoluteConcs[soluteIndex][edgeIndex].value - S[index]) / (flowLengthP[edgeIndex] / 2.0);
  double heatDisCoeff = edgeMechDispersionCoeff[edgeIndex] + parentElement->model->m_solute_molecular_diff[soluteIndex];

  edgeGradSoluteConcs[soluteIndex][edgeIndex].value = gradSolute;
  double edgeFlow = heatDisCoeff *  gradSolute * edgeDepths[edgeIndex] * flowWidth[edgeIndex];

  return edgeFlow;
}

double ElementCell::computeEdgeGradSoluteBC(int edgeIndex, double dt, double S[], int soluteIndex)
{
  double heatDisCoeff = edgeMechDispersionCoeff[edgeIndex] + parentElement->model->m_solute_molecular_diff[soluteIndex];
  double edgeFlow = heatDisCoeff *  edgeGradSoluteConcs[soluteIndex][edgeIndex].value * edgeDepths[edgeIndex] * flowWidth[edgeIndex];

  return edgeFlow;
}

double ElementCell::computeNeighborSoluteBC(int edgeIndex, double dt, double S[], int soluteIndex)
{
  ElementCell *n = neighbors[edgeIndex];

  double flp = flowLengthP[edgeIndex] / 2.0;
  double fln = flowLengthN[edgeIndex] / 2.0;

  double gradSolute = (S[n->index] - S[index]) / (flp + fln);

  double heatDisCoeff = edgeMechDispersionCoeff[edgeIndex] +
                        parentElement->model->m_solute_molecular_diff[soluteIndex];


  edgeGradSoluteConcs[soluteIndex][edgeIndex].value = gradSolute;
  double edgeFlow = heatDisCoeff * gradSolute *  edgeDepths[edgeIndex] * flowWidth[edgeIndex];

  return edgeFlow;
}

double ElementCell::computeZeroSoluteBC(int edgeIndex, double dt, double S[], int soluteIndex)
{
  edgeGradSoluteConcs[soluteIndex][edgeIndex].value = 0.0;
  return 0.0;
}



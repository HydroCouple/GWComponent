/*!
*  \file    gwmodelio.cpp
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
#include <QTextStream>
#include <QDateTime>

#include "gwmodel.h"
#include "temporal/timedata.h"
#include "element.h"
#include "elementjunction.h"
#include "temporal/timedata.h"
#include "threadsafenetcdf/threadsafencfile.h"
#include "threadsafenetcdf/threadsafencdim.h"
#include "threadsafenetcdf/threadsafencatt.h"
#include "temporal/timeseries.h"
#include "edgebc.h"
#include "elementcellsourcebc.h"
#include "channelbc.h"

#include <QDir>
#include <QDate>

using namespace std;

#ifdef USE_NETCDF

using namespace netCDF;
using namespace netCDF::exceptions;

#endif

const float GWModel::dir[] = {-1.0,-1.0,1.0,1.0};

bool GWModel::verbose() const
{
  return m_verbose;
}

void GWModel::setVerbose(bool verbose)
{
  m_verbose = verbose;
}

int GWModel::printFrequency() const
{
  return m_printFrequency;
}

void GWModel::setPrintFrequency(int printFreq)
{
  m_printFrequency = printFreq;
}

int GWModel::flushToDiskFrequency() const
{
  return m_flushToDiskFrequency;
}

void GWModel::setFlushToDiskFrequency(int diskFlushFrequency)
{
  m_flushToDiskFrequency = diskFlushFrequency;
}

QFileInfo GWModel::inputFile() const
{
  return m_inputFile;
}

void GWModel::setInputFile(const QFileInfo &inputFile)
{
  m_inputFile = inputFile;
}

QFileInfo GWModel::outputNetCDFFile() const
{
  return m_outputNetCDFFileInfo;
}

void GWModel::setOutputNetCDFFile(const QFileInfo &outputNetCDFFile)
{
  m_outputNetCDFFileInfo = outputNetCDFFile;
}

void GWModel::printStatus()
{
  m_currentPrintCount++;

  if (m_currentPrintCount >= m_printFrequency)
  {

    printf("GWModel TimeStep (s): %f\tDateTime: %f\tIters: %i/%i\t Head (m) { Min: %f\tMax: %f\tTotalMassBalance: %g (m^3/s)}", m_timeStep, m_currentDateTime,
           m_odeSolver->getIterations(), m_odeSolver->maxIterations(), m_minHead, m_maxHead, m_totalMassBalance);

    if(m_solveHeatTransport)
    {
      printf("\tTemperature (°C) { Min: %f\tMax: %f\tTotalMassBalance: %g (KJ)}", m_minTemp, m_maxTemp, m_totalHeatBalance);
    }


    for (size_t j = 0; j < m_solutes.size(); j++)
    {
      std::string &solute = m_solutes[j];

      if(solute == "WATER_AGE")
      {
        printf("\t%s (days) { Min: %f\tMax: %f}", solute.c_str(),
               m_minSolute[j], m_maxSolute[j]);
      }
      else
      {
        printf("\t%s (kg/m) { Min: %f\tMax: %f\tTotalMassBalance: %g (kg)}", solute.c_str(),
               m_minSolute[j], m_maxSolute[j], m_totalSoluteMassBalance[j]);
      }
    }

    printf("\n");

    m_currentPrintCount = 0;
  }
}

bool GWModel::initializeInputFiles(list<string> &errors)
{
  if (QFile::exists(m_inputFile.absoluteFilePath()))
  {
    QFile file(m_inputFile.absoluteFilePath());

    if (file.open(QIODevice::ReadOnly))
    {
      m_timeSeries.clear();

      m_delimiters = QRegExp("(\\,|\\t|\\;|\\s+)");
      int currentFlag = -1;
      m_addedSoluteCount = 0;

      QTextStream streamReader(&file);
      int lineCount = 0;

      while (!streamReader.atEnd())
      {
        QString line = streamReader.readLine().trimmed();
        lineCount++;

        if (!line.isEmpty() && !line.isNull())
        {
          bool readSuccess = true;
          QString error = "";

          auto it = m_inputFileFlags.find(line.toStdString());

          if (it != m_inputFileFlags.cend())
          {
            currentFlag = it->second;
          }
          else if (!QStringRef::compare(QStringRef(&line, 0, 2), ";;"))
          {
            //commment do nothing
          }
          else
          {
            switch (currentFlag)
            {
              case 1:
                readSuccess = readInputFileOptionTag(line, error);
                break;
              case 2:
                readSuccess = readInputFileOutputTag(line, error);
                break;
              case 3:
                readSuccess = readInputFileSolutesTag(line, error);
                break;
              case 4:
                readSuccess = readInputFileElementJunctionsTag(line, error);
                break;
              case 5:
                readSuccess = readInputFileElementsTag(line, error);
                break;
              case 6:
                readSuccess = readInputFileElementCellWidths(line, error);
                break;
              case 7:
                readSuccess = readInputFileElementCellHydraulics(line, error);
                break;
              case 8:
                readSuccess = readInputFileElementCellInitConditions(line, error);
                break;
              case 9:
                readSuccess = readInputFileElementCellBoundaryConditions(line, error);
                break;
              case 10:
                readSuccess = readInputFileTimeSeriesTag(line, error);
                break;
              case 11:
                readSuccess = readInputFileElementCellSources(line, error);
                break;
              case 12:
                readSuccess = readInputFileChannelBoundaryConditions(line, error);
                break;
            }
          }

          if (!readSuccess)
          {
            errors.push_back("Line " + std::to_string(lineCount) + " : " + error.toStdString());
            file.close();
            return false;
            break;
          }
        }
      }

      file.close();
    }
  }

  return true;
}

bool GWModel::initializeOutputFiles(std::list<std::string> &errors)
{
  return initializeNetCDFOutputFile(errors);
}

bool GWModel::initializeNetCDFOutputFile(std::list<std::string> &errors)
{

#ifdef USE_NETCDF


  if (m_outputNetCDFFileInfo.isRelative())
  {
    m_outputNetCDFFileInfo = relativePathToAbsolute(m_outputNetCDFFileInfo);
  }

  if (m_outputNetCDFFileInfo.absoluteFilePath().isEmpty() || m_outputNetCDFFileInfo.isDir())
  {
    return true;
  }
  else if (!m_outputNetCDFFileInfo.absoluteFilePath().isEmpty() &&
           !m_outputNetCDFFileInfo.absoluteFilePath().isNull() &&
           !m_outputNetCDFFileInfo.absoluteDir().exists())
  {
    std::string message = "NetCDF output file directory does not exist: " + m_outputNetCDFFileInfo.absoluteFilePath().toStdString();
    errors.push_back(message);
    return false;
  }

  bool returnValue = false;

  closeOutputNetCDFFile();

  try
  {
    m_outNetCDFVariables.clear();

    m_outputNetCDF = new ThreadSafeNcFile(m_outputNetCDFFileInfo.absoluteFilePath().toStdString(), NcFile::replace);

    //time variable
    ThreadSafeNcDim timeDim =  m_outputNetCDF->addDim("time");
    ThreadSafeNcVar timeVar =  m_outputNetCDF->addVar("time", NcType::nc_DOUBLE, timeDim);
    timeVar.putAtt("long_name", "Time");
    timeVar.putAtt("standard_name", "time");
    timeVar.putAtt("calendar", "julian");
    m_outNetCDFVariables["time"] = timeVar;

    //Add Solutes
    ThreadSafeNcDim solutesDim =  m_outputNetCDF->addDim("solutes",m_numSolutes);
    ThreadSafeNcVar solutes =  m_outputNetCDF->addVar("solute_names", NcType::nc_STRING, solutesDim);
    solutes.putAtt("long_name", "Solutes");
    m_outNetCDFVariables["solutes"] = solutes;

    if(m_numSolutes > 0)
    {
      char **soluteNames = new char *[m_numSolutes];

      for (int i = 0; i < m_numSolutes; i++)
      {
        string soluteName = m_solutes[i];
        soluteNames[i] = new char[soluteName.size() + 1];
        strcpy(soluteNames[i], soluteName.c_str());
      }

      solutes.putVar(soluteNames);

      for (int i = 0; i < m_numSolutes; i++)
      {
        delete[] soluteNames[i];
      }

      delete[] soluteNames;
    }

    //Add element junctions
    ThreadSafeNcDim junctionDim =  m_outputNetCDF->addDim("element_junctions", m_elementJunctions.size());

    ThreadSafeNcVar junctionIdentifiers =  m_outputNetCDF->addVar("element_junction_id", NcType::nc_STRING, junctionDim);
    junctionIdentifiers.putAtt("long_name", "Element Junction Identifiers");
    m_outNetCDFVariables["element_junction_id"] = junctionIdentifiers;


    ThreadSafeNcVar junctionX =  m_outputNetCDF->addVar("x", NcType::nc_FLOAT, junctionDim);
    junctionX.putAtt("long_name", "Junction X Coordinate");
    junctionX.putAtt("units", "m");
    m_outNetCDFVariables["x"] = junctionX;

    ThreadSafeNcVar junctionY =  m_outputNetCDF->addVar("y", NcType::nc_FLOAT, junctionDim);
    junctionY.putAtt("long_name", "Junction Y Coordinate");
    junctionY.putAtt("units", "m");
    m_outNetCDFVariables["y"] = junctionY;

    ThreadSafeNcVar junctionZ =  m_outputNetCDF->addVar("z", NcType::nc_FLOAT, junctionDim);
    junctionZ.putAtt("long_name", "Junction Z Coordinate");
    junctionZ.putAtt("units", "m");
    m_outNetCDFVariables["z"] = junctionZ;

    float *vertx = new float[m_elementJunctions.size()];
    float *verty = new float[m_elementJunctions.size()];
    float *vertz = new float[m_elementJunctions.size()];
    char **junctionIds = new char *[m_elementJunctions.size()];

    //write other relevant junction attributes here.
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < (int)m_elementJunctions.size(); i++)
    {
      ElementJunction *junction = m_elementJunctions[i];

      junctionIds[i] = new char[junction->id.size() + 1];
      strcpy(junctionIds[i], junction->id.c_str());

      vertx[i] = junction->x;
      verty[i] = junction->y;
      vertz[i] = junction->z;
    }

    junctionX.putVar(vertx);
    junctionY.putVar(verty);
    junctionZ.putVar(vertz);
    junctionIdentifiers.putVar(junctionIds);

    delete[] vertx;
    delete[] verty;
    delete[] vertz;

    for (size_t i = 0; i < m_elementJunctions.size(); i++)
    {
      delete[] junctionIds[i];
    }

    delete[] junctionIds;

    //Add Elements
    ThreadSafeNcDim elementsDim =  m_outputNetCDF->addDim("elements", m_elements.size());

    ThreadSafeNcVar elementIdentifiers =  m_outputNetCDF->addVar("element_id", NcType::nc_STRING, elementsDim);
    elementIdentifiers.putAtt("long_name", "Element Identifier");
    m_outNetCDFVariables["element_id"] = elementIdentifiers;

    ThreadSafeNcVar elementFromJunction =  m_outputNetCDF->addVar("from_junction", NcType::nc_INT64, elementsDim);
    elementFromJunction.putAtt("long_name", "Upstream Junction");
    m_outNetCDFVariables["from_junction"] = elementFromJunction;

    ThreadSafeNcVar elementToJunction =  m_outputNetCDF->addVar("to_junction", NcType::nc_INT64, elementsDim);
    elementToJunction.putAtt("long_name", "Downstream Junction");
    m_outNetCDFVariables["to_junction"] = elementToJunction;

    ThreadSafeNcVar element_x =  m_outputNetCDF->addVar("element_x", NcType::nc_FLOAT, elementsDim);
    element_x.putAtt("standard_name", "projection_x_coordinate");
    element_x.putAtt("long_name", "Element X Coordinate");
    element_x.putAtt("units", "m");
    m_outNetCDFVariables["element_x"] = element_x;

    ThreadSafeNcVar element_y =  m_outputNetCDF->addVar("element_y", NcType::nc_FLOAT, elementsDim);
    element_y.putAtt("standard_name", "projection_y_coordinate");
    element_y.putAtt("long_name", "Element Y Coordinate");
    element_y.putAtt("units", "m");
    m_outNetCDFVariables["element_y"] = element_y;


    int *fromJunctions = new int[m_elements.size()];
    int *toJunctions = new int[m_elements.size()];
    char **elementIds = new char *[m_elements.size()];
    float *elX = new float[m_elements.size()];
    float *elY = new float[m_elements.size()];

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < (int)m_elements.size(); i++)
    {
      Element *element = m_elements[i];

      elementIds[i] = new char[element->id.size() + 1];
      strcpy(elementIds[i], element->id.c_str());

      fromJunctions[i] = element->upstreamJunction->index;
      toJunctions[i] = element->downstreamJunction->index;

      elX[i] = element->x;
      elY[i] = element->y;
    }

    elementIdentifiers.putVar(elementIds);
    elementFromJunction.putVar(fromJunctions);
    elementToJunction.putVar(toJunctions);
    element_x.putVar(elX);
    element_y.putVar(elY);

    delete[] fromJunctions;
    delete[] toJunctions;
    delete[] elX;
    delete[] elY;

    for (size_t i = 0; i < m_elements.size(); i++)
    {
      delete[] elementIds[i];
    }

    delete[] elementIds;

    //ElementCells
    ThreadSafeNcDim elementCellsDim =  m_outputNetCDF->addDim("element_cells", m_totalCellsPerElement);

    //Edges
    ThreadSafeNcDim elementCellFaceDim = m_outputNetCDF->addDim("element_cell_face", 4);

    //hydraulics variables
    ThreadSafeNcVar lengthVar =  m_outputNetCDF->addVar("length", "float",
                                                        std::vector<std::string>({"elements"}));
    lengthVar.putAtt("long_name", "Element Length");
    lengthVar.putAtt("units", "m");
    m_outNetCDFVariables["length"] = lengthVar;

    std::vector<float> lengths(m_elements.size());

    for(size_t i = 0; i < m_elements.size(); i++)
    {
      lengths[i] = m_elements[i]->length;
    }

    lengthVar.putVar(lengths.data());

    ThreadSafeNcVar cellWidths = m_outputNetCDF->addVar("cell_width", NcType::nc_FLOAT, elementCellsDim);
    cellWidths.putAtt("long_name","Element Cell Width");
    cellWidths.putAtt("units","m");
    m_outNetCDFVariables["cell_width"] = cellWidths;

    std::vector<float> cell_widths(m_totalCellsPerElement, 0.0);

    for(size_t i = 0; i < m_elements.size(); i++)
    {
      for(int j  = 0; j < m_totalCellsPerElement; j++)
      {
        cell_widths[j] = m_elements[i]->elementCells[j]->width;
      }
      break;
    }

    cellWidths.putVar(cell_widths.data());

    //hydraulics variables
    ThreadSafeNcVar hydHeadVar =  m_outputNetCDF->addVar("hydraulic_head", "float",
                                                         std::vector<std::string>({"time", "elements", "element_cells"}));
    hydHeadVar.putAtt("long_name", "Hydraulic Head");
    hydHeadVar.putAtt("units", "m");
    m_outNetCDFVariables["hydraulic_head"] = hydHeadVar;


    ThreadSafeNcVar dVolumeDtVar =  m_outputNetCDF->addVar("dvolume_dt", "float",
                                                           std::vector<std::string>({"time", "elements", "element_cells"}));
    dVolumeDtVar.putAtt("long_name", "Volume Time Derivative");
    dVolumeDtVar.putAtt("units", "m");
    m_outNetCDFVariables["dvolume_dt"] = dVolumeDtVar;


    ThreadSafeNcVar satDepthVar =  m_outputNetCDF->addVar("saturated_depth", "float",
                                                          std::vector<std::string>({"time", "elements", "element_cells"}));
    satDepthVar.putAtt("long_name", "Saturated Depth");
    satDepthVar.putAtt("units", "m");
    m_outNetCDFVariables["saturated_depth"] = satDepthVar;


    ThreadSafeNcVar totalElementCellMassBalanceVar =  m_outputNetCDF->addVar("total_element_cell_mass_balance", "float",
                                                                             std::vector<std::string>({"time", "elements","element_cells"}));
    totalElementCellMassBalanceVar.putAtt("long_name", "Total Element Cell Mass Balance");
    totalElementCellMassBalanceVar.putAtt("units", "m^3");
    m_outNetCDFVariables["total_element_cell_mass_balance"] = totalElementCellMassBalanceVar;

    ThreadSafeNcVar elementExternalInflowVar =  m_outputNetCDF->addVar("element_external_inflow", "float",
                                                                       std::vector<std::string>({"time", "elements","element_cells"}));
    elementExternalInflowVar.putAtt("long_name", "Element Cell External Inflow");
    elementExternalInflowVar.putAtt("units", "m^3");
    m_outNetCDFVariables["element_cell_external_inflow"] = elementExternalInflowVar;

    ThreadSafeNcVar elementCellExternalInflowVar =  m_outputNetCDF->addVar("element_cell_external_inflow", "float",
                                                                           std::vector<std::string>({"time", "elements"}));
    elementCellExternalInflowVar.putAtt("long_name", "Element External Inflow");
    elementCellExternalInflowVar.putAtt("units", "m^3");
    m_outNetCDFVariables["element_external_inflow"] = elementCellExternalInflowVar;


    ThreadSafeNcVar elementChannelInflowVar =  m_outputNetCDF->addVar("element_channel_inflow", "float",
                                                                      std::vector<std::string>({"time", "elements"}));
    elementChannelInflowVar.putAtt("long_name", "Element Channel Inflow");
    elementChannelInflowVar.putAtt("units", "m^3/s");
    m_outNetCDFVariables["element_channel_inflow"] = elementChannelInflowVar;

    ThreadSafeNcVar elementChannelInflowFluxVar =  m_outputNetCDF->addVar("element_channel_flux", "float",
                                                                          std::vector<std::string>({"time", "elements"}));
    elementChannelInflowFluxVar.putAtt("long_name", "Element Channel Inflow Flux");
    elementChannelInflowFluxVar.putAtt("units", "m^3/m^2/s");
    m_outNetCDFVariables["element_channel_flux"] = elementChannelInflowFluxVar;


    ThreadSafeNcVar elementCellChannelInflowVar =  m_outputNetCDF->addVar("element_cell_channel_inflow", "float",
                                                                          std::vector<std::string>({"time", "elements","element_cells"}));
    elementCellChannelInflowVar.putAtt("long_name", "Element Cell Channel Inflow");
    elementCellChannelInflowVar.putAtt("units", "m^3/s");
    m_outNetCDFVariables["element_cell_channel_inflow"] = elementCellChannelInflowVar;


    ThreadSafeNcVar elementCellChannelInflowFluxVar =  m_outputNetCDF->addVar("element_cell_channel_flux", "float",
                                                                              std::vector<std::string>({"time", "elements","element_cells"}));
    elementCellChannelInflowFluxVar.putAtt("long_name", "Element Cell Channel Inflow Flux");
    elementCellChannelInflowFluxVar.putAtt("units", "m^3/m^2/s");
    m_outNetCDFVariables["element_cell_channel_flux"] = elementCellChannelInflowFluxVar;


    ThreadSafeNcVar elementCellFaceFlowVar =  m_outputNetCDF->addVar("element_cell_face_flow", "float",
                                                                     std::vector<std::string>({"time", "elements","element_cells","element_cell_face"}));
    elementCellFaceFlowVar.putAtt("long_name", "Element Cell Face Flow");
    elementCellFaceFlowVar.putAtt("units", "m^3/s");
    m_outNetCDFVariables["element_cell_face_flow"] = elementCellFaceFlowVar;


    ThreadSafeNcVar elementCellFaceHeatFlowVar =  m_outputNetCDF->addVar("element_cell_face_heat_flow", "float",
                                                                     std::vector<std::string>({"time", "elements","element_cells","element_cell_face"}));
    elementCellFaceHeatFlowVar.putAtt("long_name", "Element Cell Face Heat Flow");
    elementCellFaceHeatFlowVar.putAtt("units", "J/s");
    m_outNetCDFVariables["element_cell_face_heat_flow"] = elementCellFaceHeatFlowVar;


    //    ThreadSafeNcVar elementCellFaceFluxVar =  m_outputNetCDF->addVar("element_cell_face_flux", "float",
    //                                                                     std::vector<std::string>({"time", "elements","element_cells","element_cell_face"}));
    //    elementCellFaceFluxVar.putAtt("long_name", "Element Cell Face Flux");
    //    elementCellFaceFluxVar.putAtt("units", "m^3/m^2/s");
    //    m_outNetCDFVariables["element_cell_face_flux"] = elementCellFaceFluxVar;

    //    ThreadSafeNcVar elementCellFaceSupVelVar =  m_outputNetCDF->addVar("element_cell_face_sup_velocity", "float",
    //                                                                       std::vector<std::string>({"time", "elements","element_cells","element_cell_face"}));
    //    elementCellFaceSupVelVar.putAtt("long_name", "Element Cell Face Superficial Velocity");
    //    elementCellFaceSupVelVar.putAtt("units", "m/s");
    //    m_outNetCDFVariables["element_cell_face_sup_velocity"] = elementCellFaceSupVelVar;


    //    ThreadSafeNcVar elementCellSupVelXVar =  m_outputNetCDF->addVar("element_cell_sup_velocity_x", "float",
    //                                                                    std::vector<std::string>({"time", "elements","element_cells"}));
    //    elementCellSupVelXVar.putAtt("long_name", "Element Cell Superficial X Velocity");
    //    elementCellSupVelXVar.putAtt("units", "m/s");
    //    m_outNetCDFVariables["element_cell_sup_velocity_x"] = elementCellSupVelXVar;


    //    ThreadSafeNcVar elementCellSupVelYVar =  m_outputNetCDF->addVar("element_cell_sup_velocity_y", "float",
    //                                                                    std::vector<std::string>({"time", "elements","element_cells"}));
    //    elementCellSupVelYVar.putAtt("long_name", "Element Cell Superficial Y Velocity");
    //    elementCellSupVelYVar.putAtt("units", "m/s");
    //    m_outNetCDFVariables["element_cell_sup_velocity_y"] = elementCellSupVelYVar;


    ThreadSafeNcVar totalMassBalanceVar =  m_outputNetCDF->addVar("total_mass_balance", "float",
                                                                  std::vector<std::string>({"time"}));
    totalMassBalanceVar.putAtt("long_name", "Total Mass Balance");
    totalMassBalanceVar.putAtt("units", "m^3");
    m_outNetCDFVariables["total_mass_balance"] = totalMassBalanceVar;

    ThreadSafeNcVar tempVar =  m_outputNetCDF->addVar("temperature", "float",
                                                      std::vector<std::string>({"time", "elements", "element_cells"}));
    tempVar.putAtt("long_name", "Temperature");
    tempVar.putAtt("units", "°C");
    m_outNetCDFVariables["temperature"] = tempVar;

    ThreadSafeNcVar waterAgeVar =  m_outputNetCDF->addVar("water_age", "float",
                                                          std::vector<std::string>({"time", "elements", "element_cells"}));
    waterAgeVar.putAtt("long_name", "Water Age");
    waterAgeVar.putAtt("units", "days");
    m_outNetCDFVariables["water_age"] = waterAgeVar;


    ThreadSafeNcVar elementExternalHeatFluxVar =  m_outputNetCDF->addVar("element_external_heat_flux", "float",
                                                                         std::vector<std::string>({"time", "elements"}));
    elementExternalHeatFluxVar.putAtt("long_name", "Element External Heat Flux");
    elementExternalHeatFluxVar.putAtt("units", "J/s");
    m_outNetCDFVariables["element_external_heat_flux"] = elementExternalHeatFluxVar;


    ThreadSafeNcVar elementCellExternalHeatFluxVar =  m_outputNetCDF->addVar("element_cell_external_heat_flux", "float",
                                                                             std::vector<std::string>({"time", "elements","element_cells"}));
    elementCellExternalHeatFluxVar.putAtt("long_name", "Element Cell Channel Heat Flux");
    elementCellExternalHeatFluxVar.putAtt("units", "J/s");
    m_outNetCDFVariables["element_cell_external_heat_flux"] = elementCellExternalHeatFluxVar;


    ThreadSafeNcVar elementChannelHeatFluxVar =  m_outputNetCDF->addVar("element_channel_heat_flux", "float",
                                                                        std::vector<std::string>({"time", "elements"}));
    elementChannelHeatFluxVar.putAtt("long_name", "Element Channel Heat Flux");
    elementChannelHeatFluxVar.putAtt("units", "J/s");
    m_outNetCDFVariables["element_channel_heat_flux"] = elementChannelHeatFluxVar;


    ThreadSafeNcVar elementCellChannelHeatFluxVar =  m_outputNetCDF->addVar("element_cell_channel_heat_flux", "float",
                                                                            std::vector<std::string>({"time", "elements","element_cells"}));
    elementCellChannelHeatFluxVar.putAtt("long_name", "Element Cell Channel Heat Flux");
    elementCellChannelHeatFluxVar.putAtt("units", "J/s");
    m_outNetCDFVariables["element_cell_channel_heat_flux"] = elementCellChannelHeatFluxVar;


    ThreadSafeNcVar elementChannelWSEVar =  m_outputNetCDF->addVar("element_channel_wse", "float",
                                                                   std::vector<std::string>({"time", "elements"}));
    elementChannelWSEVar.putAtt("long_name", "Element Channel Water Surface Elevation");
    elementChannelWSEVar.putAtt("units", "m");
    m_outNetCDFVariables["element_channel_wse"] = elementChannelWSEVar;


    ThreadSafeNcVar elementChannelFlowWidthVar =  m_outputNetCDF->addVar("element_channel_flow_top_width", "float",
                                                                         std::vector<std::string>({"time", "elements"}));
    elementChannelFlowWidthVar.putAtt("long_name", "Element Channel Flow Top Width");
    elementChannelFlowWidthVar.putAtt("units", "m");
    m_outNetCDFVariables["element_channel_flow_top_width"] = elementChannelFlowWidthVar;


    ThreadSafeNcVar solutesVar =  m_outputNetCDF->addVar("solute_concentration", "float",
                                                         std::vector<std::string>({"time", "solutes", "elements", "element_cells"}));
    solutesVar.putAtt("long_name", "Solute Concentration");
    solutesVar.putAtt("units", "kg/m^3");
    m_outNetCDFVariables["solute_concentration"] = solutesVar;

    ThreadSafeNcVar elementCellFaceSoluteFlowVar =  m_outputNetCDF->addVar("element_cell_face_solute_flow", "float",
                                                                     std::vector<std::string>({"time","solutes", "elements","element_cells","element_cell_face"}));
    elementCellFaceSoluteFlowVar.putAtt("long_name", "Element Cell Face Solute Flow");
    elementCellFaceSoluteFlowVar.putAtt("units", "kg/s");
    m_outNetCDFVariables["element_cell_face_solute_flow"] = elementCellFaceSoluteFlowVar;


    ThreadSafeNcVar elementCellChannelSoluteFluxVar =  m_outputNetCDF->addVar("element_cell_channel_solute_flux", "float",
                                                                              std::vector<std::string>({"time", "solutes", "elements", "element_cells"}));
    elementCellChannelSoluteFluxVar.putAtt("long_name", "Cell Channel Solute Flux");
    elementCellChannelSoluteFluxVar.putAtt("units", "kg/s");
    m_outNetCDFVariables["element_cell_channel_solute_flux"] = elementCellChannelSoluteFluxVar;

    ThreadSafeNcVar elementChannelSoluteFluxVar =  m_outputNetCDF->addVar("element_channel_solute_flux", "float",
                                                                          std::vector<std::string>({"time", "solutes", "elements"}));
    elementChannelSoluteFluxVar.putAtt("long_name", "Channel Solute Flux");
    elementChannelSoluteFluxVar.putAtt("units", "kg/s");
    m_outNetCDFVariables["element_channel_solute_flux"] = elementChannelSoluteFluxVar;

    m_outputNetCDF->sync();

    returnValue = true;

  }
  catch (NcException &e)
  {
    std::string message = std::string(e.what());
    printf("%s\n", e.what());
    errors.push_back(message);
    returnValue = false;
  }


#endif

  return returnValue;
}

bool GWModel::readInputFileOptionTag(const QString &line, QString &errorMessage)
{
  QStringList options = line.split(m_delimiters, QString::SkipEmptyParts);
  std::string optionsFlag = options[0].toStdString();
  auto it = m_optionsFlags.find(optionsFlag);

  if (it != m_optionsFlags.end())
  {
    int optionsIndex = it->second;

    switch (optionsIndex)
    {
      case 1:
        {
          bool foundError = false;

          if (options.size() == 3)
          {
            QDateTime dateTime;
            if (SDKTemporal::DateTime::tryParse(options[1] + " " + options[2], dateTime))
            {
              m_startDateTime = SDKTemporal::DateTime::toJulianDays(dateTime);
            }
            else
            {
              foundError = true;
            }
          }

          if (foundError)
          {
            errorMessage = "Start datetime format error";
            return false;
          }
        }
        break;
      case 2:
        {
          bool foundError = false;

          if (options.size() == 3)
          {
            QDateTime dateTime;
            if (SDKTemporal::DateTime::tryParse(options[1] + " " + options[2], dateTime))
            {
              m_endDateTime = SDKTemporal::DateTime::toJulianDays(dateTime);
            }
            else
            {
              foundError = true;
            }
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Start datetime format error";
            return false;
          }
        }
        break;
      case 3:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_outputInterval = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Report interval error";
            return false;
          }
        }
        break;
      case 4:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_maxTimeStep = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Max time step error";
            return false;
          }
        }
        break;
      case 5:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_minTimeStep = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Min time step error";
            return false;
          }
        }
        break;
      case 6:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_numInitFixedTimeSteps = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Number of initial time step error";
            return false;
          }
        }
        break;
      case 7:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            m_useAdaptiveTimeStep = QString::compare(options[1], "No", Qt::CaseInsensitive);
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Use adaptive time step error";
            return false;
          }
        }
        break;
      case 8:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_timeStepRelaxationFactor = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Time step relaxation factor error";
            return false;
          }
        }
        break;
      case 9:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            std::string code = options[1].toUpper().toStdString();
            auto it = m_solverTypeFlags.find(code);

            int solverMode = -1;

            if (it != m_solverTypeFlags.end())
              solverMode = it->second;

            switch (solverMode)
            {
              case 1:
                m_odeSolver->setSolverType(ODESolver::RK4);
                break;
              case 2:
                m_odeSolver->setSolverType(ODESolver::RKQS);
                break;
              case 3:
                {
                  m_odeSolver->setSolverType(ODESolver::CVODE_ADAMS);
                  m_odeSolver->setSolverIterationMethod(ODESolver::IterationMethod::FUNCTIONAL);
                }
                break;
              case 4:
                {
                  m_odeSolver->setSolverType(ODESolver::CVODE_BDF);
                  m_odeSolver->setSolverIterationMethod(ODESolver::IterationMethod::NEWTON);
                  m_odeSolver->setLinearSolverType(ODESolver::LinearSolverType::GMRES);
                }
                break;
              case 5:
                m_odeSolver->setSolverType(ODESolver::EULER);
                break;
              default:
                foundError = true;
                break;
            }
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Hydraulic head solver type error";
            return false;
          }
        }
        break;
      case 10:
        {
          bool foundError = false;

          if (options.size() == 2)
          {

            bool ok;
            double abs_tol = options[1].toDouble(&ok);

            if (ok)
              m_odeSolver->setAbsoluteTolerance(abs_tol);

            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Hydraulic head solver absolute tolerance error";
            return false;
          }
        }
        break;
      case 11:
        {
          bool foundError = false;

          if (options.size() == 2)
          {

            bool ok;
            double rel_tol = options[1].toDouble(&ok);

            if (ok)
              m_odeSolver->setRelativeTolerance(rel_tol);

            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Hydraulic head solver relative tolerance error";
            return false;
          }
        }
        break;
      case 12:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            m_solveHeatTransport = QString::compare(options[1], "No", Qt::CaseInsensitive);
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Solve heat transport";
            return false;
          }
        }
        break;
      case 13:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_waterDensity = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Water density error";
            return false;
          }
        }
        break;
      case 14:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_cp = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Specific heat capacity of water error";
            return false;
          }
        }
        break;
      case 15:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_sedDensity = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Default sediment density";
            return false;
          }
        }
        break;
      case 16:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_sedCp = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Default sediment sepecific heat capacity";
            return false;
          }
        }
        break;
      case 17:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_porosity = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Default sediment porosity";
            return false;
          }
        }
        break;
      case 18:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_hydConX = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Hydraulic conductivity in the x-direction";
            return false;
          }
        }
        break;
      case 19:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_hydConY = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Hydraulic conductivity in the y-direction";
            return false;
          }
        }
        break;
      case 20:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_specificStorage = options[1].toDouble(&ok);
            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Default specific storage";
            return false;
          }
        }
        break;
      case 21:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_numLeftCellsPerElement = options[1].toDouble(&ok);
            foundError = !ok || m_numLeftCellsPerElement < 0;

            m_totalCellsPerElement = m_numLeftCellsPerElement + m_numRightCellsPerElement;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Number of cells per element must be even";
            return false;
          }
        }
        break;
      case 22:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_numRightCellsPerElement = options[1].toDouble(&ok);
            foundError = !ok || m_numRightCellsPerElement < 0;

            m_totalCellsPerElement = m_numLeftCellsPerElement + m_numRightCellsPerElement;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Cell width";
            return false;
          }
        }
        break;
      case 23:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool parsed = false;
            int numSolutes = options[1].toDouble(&parsed);

            if (parsed)
              setNumSolutes(numSolutes);

            foundError = !parsed;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Number of solutes error";
            return false;
          }
        }
        break;
      case 24:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            m_verbose = QString::compare(options[1], "No", Qt::CaseInsensitive);
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Verbose error";
            return false;
          }
        }
        break;
      case 25:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool parsed = false;
            m_flushToDiskFrequency = options[1].toDouble(&parsed);
            foundError = !parsed;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Flush to disk frequency error";
            return false;
          }
        }
        break;
      case 26:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool parsed = false;
            m_printFrequency = options[1].toDouble(&parsed);
            foundError = !parsed;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Print frequency error";
            return false;
          }
        }
        break;
      case 27:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            m_defaultCellWidth = options[1].toDouble(&ok);
            foundError = !ok || m_defaultCellWidth <= 0.0;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Default cell width invalid";
            return false;
          }
        }
        break;
      case 28:
        {
          bool foundError = false;

          if (options.size() == 2)
          {

            bool ok;
            m_dispersivityX = options[1].toDouble(&ok);

            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "dispersivity in x-direction error";
            return false;
          }
        }
        break;
      case 29:
        {
          bool foundError = false;

          if (options.size() == 2)
          {

            bool ok;
            m_dispersivityY = options[1].toDouble(&ok);

            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "dispersivity in Y-direction error";
            return false;
          }
        }
        break;
      case 30:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            std::string code = options[1].toUpper().toStdString();
            auto it = m_advectionFlags.find(code);

            int advectionMode = -1;

            if (it != m_advectionFlags.end())
              advectionMode = it->second;

            switch (advectionMode)
            {
              case 1:
                m_advectionMode = AdvectionDiscretizationMode::Upwind;
                break;
              case 2:
                m_advectionMode = AdvectionDiscretizationMode::Central;
                break;
              case 3:
                m_advectionMode = AdvectionDiscretizationMode::Hybrid;
                break;
              case 4:
                m_advectionMode = AdvectionDiscretizationMode::TVD;
                break;
              default:
                foundError = true;
                break;
            }
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Advection mode error";
            return false;
          }
        }
        break;
      case 31:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            bool ok;
            int fluxLimiter = options[1].toInt(&ok);
            foundError = fluxLimiter > 7 || !ok;

            m_TVDFluxLimiter = fluxLimiter;
          }
          else
          {
            foundError = true;
          }

          if(foundError)
          {
            errorMessage = "TVD flux limiter error";
            return false;
          }
        }
        break;
      case 32:
        {
          bool foundError = false;

          if (options.size() == 2)
          {

            bool ok;
            m_waterThermalConductivity = options[1].toDouble(&ok);

            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Water thermal conductivity error";
            return false;
          }
        }
        break;
      case 33:
        {
          bool foundError = false;

          if (options.size() == 2)
          {

            bool ok;
            m_sedThermalConductivity = options[1].toDouble(&ok);

            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Sediment thermal conductivity error";
            return false;
          }
        }
        break;
      case 34:
        {
          bool foundError = false;

          if (options.size() == 2)
          {

            bool ok;
            m_hydConZ = options[1].toDouble(&ok);

            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Default hydraulic conductivity Z direction";
            return false;
          }
        }
        break;
      case 35:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            std::string code = options[1].toUpper().toStdString();
            auto it = m_linearSolverTypeFlags.find(code);

            int linearSolver = -1;

            if (it != m_linearSolverTypeFlags.end())
              linearSolver = it->second;

            switch (linearSolver)
            {
              case 1:
                m_odeSolver->setLinearSolverType(ODESolver::LinearSolverType::GMRES);
                break;
              case 2:
                m_odeSolver->setLinearSolverType(ODESolver::LinearSolverType::FGMRES);
                break;
              case 3:
                m_odeSolver->setLinearSolverType(ODESolver::LinearSolverType::Bi_CGStab);
                break;
              case 4:
                m_odeSolver->setLinearSolverType(ODESolver::LinearSolverType::TFQMR);
                break;
              case 5:
                m_odeSolver->setLinearSolverType(ODESolver::LinearSolverType::PCG);
                break;
              default:
                foundError = true;
                break;
            }
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Linear solver type error";
            return false;
          }
        }
        break;
      case 36:
        {
          bool foundError = false;

          if (options.size() == 2)
          {
            m_simulateWaterAge = QString::compare(options[1], "No", Qt::CaseInsensitive);

            if(m_simulateWaterAge)
              setNumSolutes(m_numSolutes);
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Simulate water age error";
            return false;
          }
        }
        break;
      case 37:
        {
          bool foundError = false;

          if (options.size() == 2)
          {

            bool ok;
            m_dispersivityZ = options[1].toDouble(&ok);

            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "dispersivity in z-direction error";
            return false;
          }
        }
        break;
    }
  }

  return true;
}

bool GWModel::readInputFileOutputTag(const QString &line, QString &errorMessage)
{
  QStringList options = line.split(m_delimiters, QString::SkipEmptyParts);
  QString optionsFlag = options[0];

  if (options.size() == 2)
  {
    if (!QString::compare(optionsFlag, "netcdf", Qt::CaseInsensitive))
    {
      m_outputNetCDFFileInfo = QFileInfo(options[1]);
    }
  }
  else
  {
    errorMessage = "Output file error";
    return false;
  }

  return true;
}

bool GWModel::readInputFileSolutesTag(const QString &line, QString &errorMessage)
{
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() >= 4)
  {
    bool foundError = false;

    if (m_addedSoluteCount < (int)m_solutes.size())
    {
      m_solutes[m_addedSoluteCount] = columns[0].toStdString();

      bool parsed;

      double first_order_k = columns[1].toDouble(&parsed);

      if(parsed)
      {
        m_solute_first_order_k[m_addedSoluteCount] = first_order_k;
      }
      else
      {
        errorMessage = "Invalid solute first order reaction rate";
        return false;
      }


      double kd = columns[2].toDouble(&parsed);

      if(parsed)
      {
        m_solute_kd[m_addedSoluteCount] = kd;
      }
      else
      {
        errorMessage = "Invalid solute first order reaction rate";
        return false;
      }

      double molDiff = columns[3].toDouble(&parsed);

      if(parsed)
      {
        m_solute_molecular_diff[m_addedSoluteCount] = molDiff;
      }
      else
      {
        errorMessage = "Invalid solute molecular diffision";
        return false;
      }

      m_addedSoluteCount++;
    }
    else
    {
      errorMessage = "Solute error count";
      return false;
    }
  }
  else
  {
    errorMessage = "Solute error";
    return false;
  }

  return true;
}

bool GWModel::readInputFileElementJunctionsTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";

  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() == 4)
  {
    QString id = columns[0];

    bool workedX;
    bool workedY;
    bool workedZ;

    double x = columns[1].toDouble(&workedX);

    double y = columns[2].toDouble(&workedY);

    double z = columns[3].toDouble(&workedZ);

    if (workedX && workedY && workedZ)
    {
      addElementJunction(id.toStdString(), x, y, z);
    }
    else
    {
      errorMessage = "Junctions error";
      return false;
    }
  }
  else
  {
    errorMessage = "Junctions error";
    return false;
  }

  return true;
}

bool GWModel::readInputFileElementsTag(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if (columns.size() > 11)
  {
    QString id = columns[0];
    QString fromId = columns[1];
    QString toId = columns[2];

    auto fromIt = m_elementJunctionsById.find(fromId.toStdString());
    auto toIt = m_elementJunctionsById.find(toId.toStdString());

    if (fromIt != m_elementJunctionsById.end() &&
        toIt != m_elementJunctionsById.end())
    {
      ElementJunction *ej1 = m_elementJunctionsById[fromId.toStdString()];
      ElementJunction *ej2 = m_elementJunctionsById[toId.toStdString()];

      bool lengthOk;
      double length = columns[3].toDouble(&lengthOk);

      bool bottomElOk;
      double bottomEl = columns[4].toDouble(&bottomElOk);

      bool topElOk;
      double topEl = columns[5].toDouble(&topElOk);

      bool hydHeadOk ;
      double hydHead = columns[6].toDouble(&hydHeadOk);

      bool hydConZOk;
      double hydConZ = columns[7].toDouble(&hydConZOk);

      bool topWidthOk ;
      double topWidth = columns[8].toDouble(&topWidthOk);

      bool wseOk ;
      double wse = columns[9].toDouble(&wseOk);

      bool bedThicknessOk ;
      double bedThickness = columns[10].toDouble(&bedThicknessOk);

      bool channelKOk ;
      double channelK = columns[11].toDouble(&channelKOk);

      if (lengthOk && hydHeadOk && topWidthOk && bottomElOk && topElOk &&
          wseOk && bedThicknessOk && channelKOk && hydConZOk)
      {
        Element *element = addElement(id.toStdString(), ej1, ej2);
        element->length = length;
        element->channelWSE = wse;
        element->channelBedThickness = bedThickness;
        element->channelWidth = topWidth;
        element->channelBedHydCond = channelK;
        element->hydConZ = hydConZ;

        for(int j = 0 ; j < m_totalCellsPerElement; j++)
        {
          ElementCell *elementCell = element->elementCells[j];
          elementCell->hydHead.value = elementCell->prevHydHead.value = hydHead;
          elementCell->bedRockElev = bottomEl;
          elementCell->topElev = topEl;
        }

        int currentColumn = 12;

        if(m_solveHeatTransport)
        {
          if(columns.size() > 12)
          {
            bool tempOk;
            double temp = columns[12].toDouble(&tempOk);

            bool chanTempOk;
            double chanTemp = columns[13].toDouble(&chanTempOk);

            if(tempOk && chanTempOk)
            {
              element->channelTemperature = chanTemp;
              currentColumn = 14;

              for(int j = 0 ; j < m_totalCellsPerElement; j++)
              {
                ElementCell *elementCell = element->elementCells[j];
                elementCell->temperature.value = elementCell->prevTemperature.value = temp;
              }
            }
            else
            {
              errorMessage = "Temperature specified is invalid";
              return false;
            }
          }
          else
          {
            errorMessage = "";
            return false;
          }
        }

        if(columns.size() - currentColumn >= m_numSolutes * 2)
        {
          int soluteIndex = 0;

          for (int i = currentColumn; i < columns.size(); i = i + 2)
          {
            bool soluteOk ;
            double solute = columns[i].toDouble(&soluteOk);

            bool chanSoluteOk;
            double chanSolute = columns[i+1].toDouble(&chanSoluteOk);

            if (soluteOk && chanSoluteOk)
            {
              element->channelSoluteConcs[soluteIndex] = chanSolute;

              for(int j = 0 ; j < m_totalCellsPerElement; j++)
              {
                ElementCell *elementCell = element->elementCells[j];
                elementCell->soluteConcs[soluteIndex].value = elementCell->prevSoluteConcs[soluteIndex].value = solute;
              }
            }
            else
            {
              errorMessage = "Wrong initial solute concetration";
              return false;
            }
          }
        }
      }
      else
      {
        errorMessage = "Must be valid number";
        return false;
      }
    }
    else
    {
      errorMessage = "Wrong upstream or downstream junction";
      return false;
    }
  }

  return true;
}

bool GWModel::readInputFileElementCellWidths(const QString &line, QString &errorMessage)
{
  errorMessage = "";
  QStringList columns = line.split(m_delimiters, QString::SkipEmptyParts);

  if(columns.size() == m_totalCellsPerElement)
  {
    for(int i = 0; i < columns.size(); i++)
    {
      bool widthOk;
      double width = columns[i].toDouble(&widthOk);

      if(widthOk && width > 0)
      {
        for(size_t j = 0; j < m_elements.size(); j++)
        {
          Element *element = m_elements[j];
          element->elementCells[i]->width = width;
        }
      }
      else
      {
        errorMessage = "Width must be positive number";
        return false;
      }
    }
  }
  else
  {
    errorMessage = "Number of width values specified must match left cells plus right cells";
    return false;
  }

  return true;
}

bool GWModel::readInputFileElementCellHydraulics(const QString &line, QString &errorMessage)
{
  QStringList options = line.split(m_delimiters, QString::SkipEmptyParts);

  if(options.size() ==  2 + m_numLeftCellsPerElement + m_numRightCellsPerElement)
  {
    QString elementId = options[0].trimmed();

    std::vector<double> values;

    for(int i = 2; i < options.size(); i++)
    {
      bool valueOk;
      double value = options[i].toDouble(&valueOk);

      if(valueOk)
      {
        values.push_back(value);
      }
      else
      {
        errorMessage = "Specified variable is not number";
        return false;
      }
    }

    if(m_elementsById.find(elementId.toStdString()) != m_elementsById.end())
    {
      Element *element = m_elementsById[elementId.toStdString()];
      std::string optionsFlag = options[1].toStdString();
      auto it = m_hydraulicVariableFlags.find(optionsFlag);

      if (it != m_hydraulicVariableFlags.end())
      {
        int hydIndex = it->second;

        switch (hydIndex)
        {
          case 1:
            {
              for(int j = 0; j < m_numLeftCellsPerElement + m_numRightCellsPerElement; j++)
              {
                element->elementCells[j]->hydHead.value = values[j];
                element->elementCells[j]->prevHydHead.value = values[j];
              }
            }
            break;
          case 2:
            {
              for(int j = 0; j < m_numLeftCellsPerElement + m_numRightCellsPerElement; j++)
              {
                element->elementCells[j]->bedRockElev = values[j];
              }
            }
            break;
          case 3:
            {
              for(int j = 0; j < m_numLeftCellsPerElement + m_numRightCellsPerElement; j++)
              {
                element->elementCells[j]->topElev = values[j];
              }
            }
            break;
          case 4:
            {
              for(int j = 0; j < m_numLeftCellsPerElement + m_numRightCellsPerElement; j++)
              {
                ElementCell *n = element->elementCells[j];
                n->hydCon[1] = values[j];
              }
            }
            break;
          case 5:
            {
              for(int j = 0; j < m_numLeftCellsPerElement + m_numRightCellsPerElement; j++)
              {
                ElementCell *n = element->elementCells[j];
                n->hydCon[0] = values[j];
              }
            }
            break;
          case 6:
            {
              for(int j = 0; j < m_numLeftCellsPerElement + m_numRightCellsPerElement; j++)
              {
                element->elementCells[j]->specificStorage = values[j];
              }
            }
            break;
          case 7:
            {
              for(int j = 0; j < m_numLeftCellsPerElement + m_numRightCellsPerElement; j++)
              {
                element->elementCells[j]->porosity = values[j];
              }
            }
            break;
          case 8:
            {
              for(int j = 0; j < m_numLeftCellsPerElement + m_numRightCellsPerElement; j++)
              {
                element->elementCells[j]->sedDensity = values[j];
              }
            }
            break;
          case 9:
            {
              for(int j = 0; j < m_numLeftCellsPerElement + m_numRightCellsPerElement; j++)
              {
                element->elementCells[j]->sedCp = values[j];
              }
            }
            break;
        }
      }
    }
    else
    {
      errorMessage = "Element specified not found: " + elementId;
      return false;
    }
  }
  else
  {
    errorMessage = "Number of columes must be of size " + QString::number(2 + m_numLeftCellsPerElement + m_numRightCellsPerElement);
    return false;
  }


  return true;
}

bool GWModel::readInputFileElementCellInitConditions(const QString &line, QString &errorMessage)
{
  QStringList options = line.split(m_delimiters, QString::SkipEmptyParts);

  if(options.size() ==  2 + m_numLeftCellsPerElement + m_numRightCellsPerElement)
  {
    QString elementId = options[0].trimmed();

    std::vector<double> values;

    for(int i = 2; i < options.size(); i++)
    {
      bool valueOk;
      double value = options[i].toDouble(&valueOk);

      if(valueOk)
      {
        values.push_back(value);
      }
      else
      {
        errorMessage = "Specified variable is not number";
        return false;
      }
    }

    if(m_elementsById.find(elementId.toStdString()) != m_elementsById.end())
    {
      Element *element = m_elementsById[elementId.toStdString()];
      QString variable = options[1];

      if(!variable.compare("Temperature", Qt::CaseInsensitive))
      {
        for(int j = 0; j < m_numLeftCellsPerElement + m_numRightCellsPerElement; j++)
        {
          element->elementCells[j]->hydHead.value = values[j];
          element->elementCells[j]->prevHydHead.value = values[j];
        }
      }
      else
      {
        for(size_t k = 0; k < m_solutes.size(); k++)
        {
          QString solute = QString::fromStdString(m_solutes[k]);

          if(!variable.compare(solute, Qt::CaseInsensitive))
          {
            for(int j = 0; j < m_numLeftCellsPerElement + m_numRightCellsPerElement; j++)
            {
              element->elementCells[j]->soluteConcs[k].value = values[j];
              element->elementCells[j]->prevSoluteConcs[k].value = values[j];
            }
            break;
          }
        }
      }
    }
    else
    {
      errorMessage = "Element specified not found: " + elementId;
      return false;
    }
  }
  else
  {
    errorMessage = "Number of columes must be of size " + QString::number(2 + m_numLeftCellsPerElement + m_numRightCellsPerElement);
    return false;
  }

  return true;
}

bool GWModel::readInputFileElementCellBoundaryConditions(const QString &line, QString &errorMessage)
{
  QStringList options = line.split(m_delimiters, QString::SkipEmptyParts);

  auto bcflag = m_BCFlags.find(options[0].toUpper().toStdString());

  if(bcflag != m_BCFlags.end())
  {
    int flag = bcflag->second;

    switch (flag)
    {
      case 1:
      case 2:
      case 3:
      case 4:
        {
          if(options.size() == 8)
          {
            QString edge = options[1].toUpper();
            std::string stdEdge = edge.toStdString();

            auto edgeIt = m_edgeFlags.find(stdEdge);

            if(edgeIt != m_edgeFlags.end())
            {
              int eFlag = edgeIt->second;
              EdgeBC::Edge edgeFlag = (EdgeBC::Edge)eFlag;

              QString fromId = options[2];
              QString toId = options[3];

              auto itFrom = m_elementsById.find(fromId.toStdString());
              auto itTo = m_elementsById.find(toId.toStdString());

              if(itFrom != m_elementsById.end() && itTo != m_elementsById.end())
              {
                Element *fromElement = itFrom->second;
                Element *toElement = itTo->second;

                bool oks;
                bool oke;

                int startCell = options[4].toInt(&oks);
                int endCell = options[5].toInt(&oke);

                if(oks && oke && endCell >= startCell &&
                   startCell >= 0 && endCell < m_numLeftCellsPerElement + m_numRightCellsPerElement)
                {
                  QString type = options[6];

                  if(!type.compare("Value", Qt::CaseInsensitive))
                  {
                    double value = options[7].toDouble(&oks);

                    if(oks)
                    {
                      EdgeBC *boundaryCondition = new EdgeBC((EdgeBC::Variable)flag,
                                                             edgeFlag,
                                                             fromElement,
                                                             toElement,
                                                             startCell,
                                                             endCell,
                                                             this);

                      QUuid uid = QUuid::createUuid();
                      QSharedPointer<TimeSeries>ts(new TimeSeries(uid.toString(), 1, this));

                      ts->addRow(m_startDateTime, value);
                      ts->addRow(m_endDateTime, value);

                      m_timeSeries[ts->id().toStdString()] = ts;

                      boundaryCondition->setTimeSeries(ts);

                      m_boundaryConditions.push_back(boundaryCondition);
                    }
                    else
                    {
                      errorMessage = "Specified value is invalid";
                      return false;
                    }
                  }
                  else if(!type.compare("TimeSeries", Qt::CaseInsensitive))
                  {
                    std::string tsId = options[7].toStdString();
                    auto tsIt = m_timeSeries.find(tsId);

                    if(tsIt != m_timeSeries.end())
                    {
                      EdgeBC *boundaryCondition = new EdgeBC((EdgeBC::Variable)flag,
                                                             edgeFlag,
                                                             fromElement,
                                                             toElement,
                                                             startCell,
                                                             endCell,
                                                             this);

                      boundaryCondition->setTimeSeries(tsIt->second);
                      m_boundaryConditions.push_back(boundaryCondition);
                    }
                    else
                    {
                      errorMessage = "Specified timeseries was not found";
                      return false;
                    }
                  }
                  else
                  {
                    errorMessage = "Specified file type is invalid";
                    return false;
                  }
                }
                else
                {
                  errorMessage = "Start/end cell was not found or is invalid";
                  return false;
                }
              }
              else
              {
                errorMessage = "Start/end element was not found";
                return false;
              }
            }
            else
            {
              errorMessage = "Specified edge flag is invalid";
              return false;
            }
          }
          else
          {
            errorMessage = "Solute boundary condition must have 9 columns";
            return false;
          }
        }
        break;
      case 5:
      case 6:
        {
          return readInputFileElementCellSoluteBoundaryConditions(flag, options, errorMessage);
        }
        break;
    }
  }

  return true;
}

bool GWModel::readInputFileElementCellSoluteBoundaryConditions(int variable, const QStringList &options, QString &errorMessage)
{
  if(options.size() == 9)
  {
    std::string solute = options[1].toStdString();

    bool soluteFound = false;

    for(size_t m = 0; m < m_solutes.size(); m++)
    {
      if(solute == m_solutes[m])
      {
        soluteFound = true;

        QString edge = options[2].toUpper();
        std::string stdEdge = edge.toStdString();

        auto edgeIt = m_edgeFlags.find(stdEdge);

        if(edgeIt != m_edgeFlags.end())
        {
          EdgeBC::Edge edgeFlag = (EdgeBC::Edge)edgeIt->second;
          QString fromId = options[3];
          QString toId = options[4];

          auto itFrom = m_elementsById.find(fromId.toStdString());
          auto itTo = m_elementsById.find(toId.toStdString());

          if(itFrom != m_elementsById.end() && itTo != m_elementsById.end())
          {
            Element *fromElement = itFrom->second;
            Element *toElement = itTo->second;

            bool oks;
            bool oke;

            int startCell = options[5].toInt(&oks);
            int endCell = options[6].toInt(&oke);

            if(oks && oke && endCell >= startCell &&
               startCell >= 0 && endCell < m_numLeftCellsPerElement + m_numRightCellsPerElement)
            {
              QString type = options[7];

              if(!type.compare("Value", Qt::CaseInsensitive))
              {
                double value = options[8].toDouble(&oks);

                if(oks)
                {
                  EdgeBC *boundaryCondition = new EdgeBC((EdgeBC::Variable)variable,
                                                         edgeFlag,
                                                         fromElement,
                                                         toElement,
                                                         startCell,
                                                         endCell,
                                                         this);
                  boundaryCondition->setSoluteIndex(m);

                  QUuid uid = QUuid::createUuid();
                  QSharedPointer<TimeSeries>ts(new TimeSeries(uid.toString(), 1, boundaryCondition));

                  ts->addRow(m_startDateTime, value);
                  ts->addRow(m_endDateTime, value);

                  m_timeSeries[ts->id().toStdString()] = ts;

                  boundaryCondition->setTimeSeries(ts);
                  m_boundaryConditions.push_back(boundaryCondition);
                }
                else
                {
                  errorMessage = "Specified value is invalid";
                  return false;
                }
              }
              else if(!type.compare("TimeSeries", Qt::CaseInsensitive))
              {
                std::string tsId = options[8].toStdString();
                auto tsIt = m_timeSeries.find(tsId);

                if(tsIt != m_timeSeries.end())
                {
                  EdgeBC *boundaryCondition = new EdgeBC((EdgeBC::Variable)variable,
                                                         edgeFlag,
                                                         fromElement,
                                                         toElement,
                                                         startCell,
                                                         endCell,
                                                         this);

                  boundaryCondition->setSoluteIndex(m);

                  boundaryCondition->setTimeSeries(tsIt->second);
                  m_boundaryConditions.push_back(boundaryCondition);
                }
                else
                {
                  errorMessage = "Specified timeseries was not found";
                  return false;
                }
              }
              else
              {
                errorMessage = "Specified file type is invalid";
                return false;
              }
            }
            else
            {
              errorMessage = "Start/end cell was not found or is invalid";
              return false;
            }
          }
          else
          {
            errorMessage = "Start/end element was not found";
            return false;
          }
        }
        else
        {
          errorMessage = "Specified edge flag is invalid";
          return false;
        }

        break;
      }
    }

    if(!soluteFound)
    {
      errorMessage = "Specified solute not found";
      return false;
    }
  }
  else
  {
    errorMessage = "Solute boundary condition must have 9 columns";
    return false;
  }

  return true;
}

bool GWModel::readInputFileElementCellSources(const QString &line, QString &errorMessage)
{
  QStringList options = line.split(m_delimiters, QString::SkipEmptyParts);

  auto srcflag = m_srcFlags.find(options[0].toUpper().toStdString());

  if(srcflag != m_BCFlags.end())
  {
    int flag = srcflag->second;

    switch (flag)
    {
      case 1:
      case 2:
        {
          if(options.size() == 8)
          {
            QString fromId = options[1];
            QString toId = options[2];

            auto itFrom = m_elementsById.find(fromId.toStdString());
            auto itTo = m_elementsById.find(toId.toStdString());

            if(itFrom != m_elementsById.end() && itTo != m_elementsById.end())
            {
              Element *fromElement = itFrom->second;
              Element *toElement = itTo->second;

              bool oks;
              bool oke;

              int startCell = options[3].toInt(&oks);
              int endCell = options[4].toInt(&oke);

              if(oks && oke && endCell >= startCell &&
                 startCell >= 0 && endCell < m_numLeftCellsPerElement + m_numRightCellsPerElement)
              {

                auto geomMult = m_geomMultFlags.find(options[5].toStdString());

                if(geomMult != m_geomMultFlags.end())
                {
                  ElementCellSourceBC::GeometryMultiplier geomMultiplier = (ElementCellSourceBC::GeometryMultiplier)geomMult->second;

                  QString type = options[6];

                  if(!type.compare("Value", Qt::CaseInsensitive))
                  {
                    double value = options[7].toDouble(&oks);

                    if(oks)
                    {
                      ElementCellSourceBC *sourceCondition = new ElementCellSourceBC((ElementCellSourceBC::Variable)flag,
                                                                                     geomMultiplier,
                                                                                     fromElement,
                                                                                     toElement,
                                                                                     startCell,
                                                                                     endCell,
                                                                                     this);

                      QUuid uid = QUuid::createUuid();
                      QSharedPointer<TimeSeries>ts(new TimeSeries(uid.toString(), 1, sourceCondition));

                      ts->addRow(m_startDateTime, value);
                      ts->addRow(m_endDateTime, value);

                      m_timeSeries[ts->id().toStdString()] = ts;

                      sourceCondition->setTimeSeries(ts);
                      m_boundaryConditions.push_back(sourceCondition);
                    }
                    else
                    {
                      errorMessage = "Specified value is invalid";
                      return false;
                    }
                  }
                  else if(!type.compare("TimeSeries", Qt::CaseInsensitive))
                  {
                    std::string tsId = options[7].toStdString();
                    auto tsIt = m_timeSeries.find(tsId);

                    if(tsIt != m_timeSeries.end())
                    {
                      ElementCellSourceBC *sourceCondition = new ElementCellSourceBC((ElementCellSourceBC::Variable)flag,
                                                                                     geomMultiplier,
                                                                                     fromElement,
                                                                                     toElement,
                                                                                     startCell,
                                                                                     endCell,
                                                                                     this);

                      sourceCondition->setTimeSeries(tsIt->second);
                      m_boundaryConditions.push_back(sourceCondition);
                    }
                    else
                    {
                      errorMessage = "Specified timeseries was not found";
                      return false;
                    }
                  }
                  else
                  {
                    errorMessage = "Specified file type is invalid";
                    return false;
                  }
                }
                else
                {
                  errorMessage = "Specified geometry multiplier is invalid";
                  return false;
                }
              }
              else
              {
                errorMessage = "Start/end cell was not found or is invalid";
                return false;
              }
            }
            else
            {
              errorMessage = "Start/end element was not found";
              return false;
            }

          }
          else
          {
            errorMessage = "Solute boundary condition must have 9 columns";
            return false;
          }
        }
        break;
      case 3:
        {

          if(options.size() == 9)
          {
            std::string solute = options[1].toStdString();
            bool soluteFound = false;

            for(size_t m = 0; m < m_solutes.size(); m++)
            {
              if(solute == m_solutes[m])
              {
                soluteFound = true;
                QString fromId = options[2];
                QString toId = options[3];

                auto itFrom = m_elementsById.find(fromId.toStdString());
                auto itTo = m_elementsById.find(toId.toStdString());

                if(itFrom != m_elementsById.end() && itTo != m_elementsById.end())
                {
                  Element *fromElement = itFrom->second;
                  Element *toElement = itTo->second;

                  bool oks;
                  bool oke;

                  int startCell = options[4].toInt(&oks);
                  int endCell = options[5].toInt(&oke);

                  if(oks && oke && endCell >= startCell &&
                     startCell >= 0 && endCell < m_numLeftCellsPerElement + m_numRightCellsPerElement)
                  {

                    auto geomMult = m_geomMultFlags.find(options[6].toStdString());

                    if(geomMult != m_geomMultFlags.end())
                    {
                      ElementCellSourceBC::GeometryMultiplier geomMultiplier = (ElementCellSourceBC::GeometryMultiplier)geomMult->second;

                      QString type = options[7];

                      if(!type.compare("Value", Qt::CaseInsensitive))
                      {
                        double value = options[8].toDouble(&oks);

                        if(oks)
                        {
                          ElementCellSourceBC *sourceCondition = new ElementCellSourceBC((ElementCellSourceBC::Variable)flag,
                                                                                         geomMultiplier,
                                                                                         fromElement,
                                                                                         toElement,
                                                                                         startCell,
                                                                                         endCell,
                                                                                         this);

                          QUuid uid = QUuid::createUuid();
                          QSharedPointer<TimeSeries>ts(new TimeSeries(uid.toString(), 1, sourceCondition));

                          ts->addRow(m_startDateTime, value);
                          ts->addRow(m_endDateTime, value);

                          m_timeSeries[ts->id().toStdString()] = ts;

                          sourceCondition->setTimeSeries(ts);
                          m_boundaryConditions.push_back(sourceCondition);
                        }
                      }
                      else if(!type.compare("TimeSeries", Qt::CaseInsensitive))
                      {
                        std::string tsId = options[8].toStdString();
                        auto tsIt = m_timeSeries.find(tsId);

                        if(tsIt != m_timeSeries.end())
                        {
                          ElementCellSourceBC *sourceCondition = new ElementCellSourceBC((ElementCellSourceBC::Variable)flag,
                                                                                         geomMultiplier,
                                                                                         fromElement,
                                                                                         toElement,
                                                                                         startCell,
                                                                                         endCell,
                                                                                         this);

                          sourceCondition->setTimeSeries(tsIt->second);
                          m_boundaryConditions.push_back(sourceCondition);
                        }
                        else
                        {
                          errorMessage = "Specified timeseries was not found";
                          return false;
                        }
                      }
                      else
                      {
                        errorMessage = "Specified file type is invalid";
                        return false;
                      }
                    }
                    else
                    {
                      errorMessage = "Specified geometry multiplier is invalid";
                      return false;
                    }
                  }
                  else
                  {
                    errorMessage = "Start/end cell was not found or is invalid";
                    return false;
                  }
                }
                else
                {
                  errorMessage = "Start/end element was not found";
                  return false;
                }

                break;
              }
            }
            if(!soluteFound)
            {
              errorMessage = "Specified solute not found";
              return false;
            }
          }
          else
          {
            errorMessage = "Solute boundary condition must have 9 columns";
            return false;
          }

        }
        break;
    }
  }

  return true;
}

bool GWModel::readInputFileTimeSeriesTag(const QString &line, QString &errorMessage)
{
  QStringList options = line.split(m_delimiters, QString::SkipEmptyParts);

  if(options.size() ==  2)
  {
    QFileInfo fileInfo(options[1].trimmed());

    if (fileInfo.isRelative())
      fileInfo = relativePathToAbsolute(fileInfo);

    if(QFile::exists(fileInfo.absoluteFilePath()))
    {
      QSharedPointer<TimeSeries> timeSeries(TimeSeries::createTimeSeries(options[0], fileInfo, this));

      if(!timeSeries.isNull())
      {
        m_timeSeries[timeSeries->id().toStdString()] = timeSeries;
      }
      else
      {
        errorMessage = "Timeseries specified is invalid";
        return false;
      }
    }
    else
    {
      errorMessage = "Specified filepath does not exist";
      return false;
    }
  }
  else
  {
    errorMessage = "TimeSeries must have two columns";
    return false;
  }

  return true;
}

bool GWModel::readInputFileChannelBoundaryConditions(const QString &line, QString &errorMessage)
{
  QStringList options = line.split(m_delimiters, QString::SkipEmptyParts);

  auto varflag = m_chanVarFlags.find(options[0].toUpper().toStdString());

  if(varflag != m_chanVarFlags.end())
  {
    int flag = varflag->second;

    switch (flag)
    {
      case 1:
      case 2:
      case 3:
        {
          if(options.size() == 5)
          {
            QString fromId = options[1];
            QString toId = options[2];

            auto itFrom = m_elementsById.find(fromId.toStdString());
            auto itTo = m_elementsById.find(toId.toStdString());

            if(itFrom != m_elementsById.end() && itTo != m_elementsById.end())
            {
              Element *fromElement = itFrom->second;
              Element *toElement = itTo->second;

              QString type = options[3];

              if(!type.compare("Value", Qt::CaseInsensitive))
              {
                bool oks;
                double value = options[4].toDouble(&oks);

                if(oks)
                {
                  ChannelBC *channelBC = new ChannelBC((ChannelBC::Variable)flag,
                                                       fromElement,
                                                       toElement,
                                                       this);

                  QUuid uid = QUuid::createUuid();
                  QSharedPointer<TimeSeries>ts(new TimeSeries(uid.toString(), 1, channelBC));

                  ts->addRow(m_startDateTime, value);
                  ts->addRow(m_endDateTime, value);

                  m_timeSeries[ts->id().toStdString()] = ts;

                  channelBC->setTimeSeries(ts);
                  m_boundaryConditions.push_back(channelBC);
                }
                else
                {
                  errorMessage = "Specified value is invalid";
                  return false;
                }
              }
              else if(!type.compare("TimeSeries", Qt::CaseInsensitive))
              {
                std::string tsId = options[4].toStdString();
                auto tsIt = m_timeSeries.find(tsId);

                if(tsIt != m_timeSeries.end())
                {
                  ChannelBC *channelBC = new ChannelBC((ChannelBC::Variable)flag,
                                                       fromElement,
                                                       toElement,
                                                       this);

                  channelBC->setTimeSeries(tsIt->second);
                  m_boundaryConditions.push_back(channelBC);
                }
                else
                {
                  errorMessage = "Specified timeseries was not found";
                  return false;
                }
              }
              else
              {
                errorMessage = "Specified file type is invalid";
                return false;
              }
            }
            else
            {
              errorMessage = "Specified geometry multiplier is invalid";
              return false;
            }
          }
          else
          {
            errorMessage = "Solute boundary condition must have 9 columns";
            return false;
          }
        }
        break;
      case 4:
        {

          if(options.size() == 6)
          {
            std::string solute = options[1].toStdString();
            bool soluteFound = false;

            for(size_t m = 0; m < m_solutes.size(); m++)
            {
              if(solute == m_solutes[m])
              {
                soluteFound = true;
                QString fromId = options[2];
                QString toId = options[3];

                auto itFrom = m_elementsById.find(fromId.toStdString());
                auto itTo = m_elementsById.find(toId.toStdString());

                if(itFrom != m_elementsById.end() && itTo != m_elementsById.end())
                {
                  Element *fromElement = itFrom->second;
                  Element *toElement = itTo->second;

                  QString type = options[4];

                  if(!type.compare("Value", Qt::CaseInsensitive))
                  {
                    bool oks;
                    double value = options[5].toDouble(&oks);

                    if(oks)
                    {
                      ChannelBC *channelBC = new ChannelBC((ChannelBC::Variable)flag,
                                                           fromElement,
                                                           toElement,
                                                           this);

                      QUuid uid = QUuid::createUuid();
                      QSharedPointer<TimeSeries>ts(new TimeSeries(uid.toString(), 1, channelBC));
                      ts->addRow(m_startDateTime, value);
                      ts->addRow(m_endDateTime, value);

                      m_timeSeries[ts->id().toStdString()] = ts;

                      channelBC->setTimeSeries(ts);
                      m_boundaryConditions.push_back(channelBC);

                    }
                    else
                    {
                      errorMessage = "Value specified is invalid";
                      return false;
                    }
                  }
                  else if(!type.compare("TimeSeries", Qt::CaseInsensitive))
                  {
                    std::string tsId = options[5].toStdString();
                    auto tsIt = m_timeSeries.find(tsId);

                    if(tsIt != m_timeSeries.end())
                    {
                      ChannelBC *channelBC = new ChannelBC((ChannelBC::Variable)flag,
                                                           fromElement,
                                                           toElement,
                                                           this);

                      channelBC->setTimeSeries(tsIt->second);
                      m_boundaryConditions.push_back(channelBC);
                    }
                    else
                    {
                      errorMessage = "Specified timeseries was not found";
                      return false;
                    }
                  }
                  else
                  {
                    errorMessage = "Specified file type is invalid";
                    return false;
                  }

                }
                else
                {
                  errorMessage = "Start/end element was not found";
                  return false;
                }

                break;
              }
            }

            if(!soluteFound)
            {
              errorMessage = "Specified solute not found";
              return false;
            }
          }
          else
          {
            errorMessage = "Solute boundary condition must have 9 columns";
            return false;
          }

        }
        break;
    }
  }

  return true;
}

void GWModel::writeOutput()
{
  m_currentflushToDiskCount++;

  if (m_currentflushToDiskCount >= m_flushToDiskFrequency)
  {
    m_flushToDisk = true;
    m_currentflushToDiskCount = 0;
  }
  else
  {
    m_flushToDisk = false;
  }

  writeNetCDFOutput();
}

void GWModel::writeNetCDFOutput()
{
#ifdef USE_NETCDF

  if (m_outputNetCDF)
  {

    size_t currentTime = m_outNetCDFVariables["time"].getDim(0).getSize();

    //Set current dateTime
    m_outNetCDFVariables["time"].putVar(std::vector<size_t>({currentTime}), m_currentDateTime);

    float *hydHead = new float[m_elements.size() * m_totalCellsPerElement];
    float *dvolumedt = new float[m_elements.size() * m_totalCellsPerElement];
    float *satDepth = new float[m_elements.size() * m_totalCellsPerElement];
    float *elementInflow = new float[m_elements.size()]();
    float *elementCellInflow = new float[m_elements.size() * m_totalCellsPerElement];
    float *elementChannelInflow = new float[m_elements.size()]();
    float *elementChannelInflowFlux = new float[m_elements.size()]();
    float *elementCellChannelInflow = new float[m_elements.size() * m_totalCellsPerElement];
    float *elementCellChannelInflowFlux = new float[m_elements.size() * m_totalCellsPerElement];
    float *totalElementCellMassBal = new float[m_elements.size() * m_totalCellsPerElement];
    float *edgeFlow = new float[m_elements.size() * m_totalCellsPerElement * 4];
    float *edgeHeatFlow = new float[m_elements.size() * m_totalCellsPerElement * 4];
    //    float *edgeFlux = new float[m_elements.size() * m_totalCellsPerElement * 4];
    //    float *edgeSupVel = new float[m_elements.size() * m_totalCellsPerElement * 4];
    //    float *cellSupVelX = new float[m_elements.size() * m_totalCellsPerElement];
    //    float *cellSupVelY = new float[m_elements.size() * m_totalCellsPerElement];
    float *temperature = new float[m_elements.size() * m_totalCellsPerElement];
    float *waterAge = new float[m_elements.size() * m_totalCellsPerElement];
    float *elementHeatFlux = new float[m_elements.size()];
    float *elementCellHeatFlux = new float[m_elements.size() * m_totalCellsPerElement];
    float *elementChannelHeatFlux = new float[m_elements.size()]();
    float *elementCellChannelHeatFlux = new float[m_elements.size() * m_totalCellsPerElement];
    float *elementChannelWSE = new float[m_elements.size()];
    float *elementChannelWidth = new float[m_elements.size()];
    float *solutes = new float[m_elements.size() * m_numSolutes * m_totalCellsPerElement];
    float *edgeSolutesFlow = new float[m_elements.size() * m_numSolutes * m_totalCellsPerElement * 4];
    float *cellChannelSoluteFlux = new float[m_elements.size() * m_totalCellsPerElement * m_numSolutes];
    float *channelSoluteFlux = new float[m_elements.size() * m_numSolutes];

    //element_cell_face_flux

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < (int)m_elements.size(); i++)
    {
      Element *element = m_elements[i];
      elementChannelWSE[i] = element->channelWSE;
      elementChannelWidth[i] = element->channelWidth;
      elementChannelInflow[i] = element->channelInflow;
      elementChannelInflowFlux[i] = element->channelInflowFlux;
      elementChannelHeatFlux[i] = element->channelHeatFlux;

      for(int j = 0; j < m_totalCellsPerElement; j++)
      {
        ElementCell *elementCell = element->elementCells[j];

        hydHead[j + i * m_totalCellsPerElement] = elementCell->hydHead.value;
        dvolumedt[j + i * m_totalCellsPerElement] = elementCell->dvolume_dt;
        satDepth[j + i * m_totalCellsPerElement] = elementCell->depth;
        elementInflow[i] += elementCell->externalInflow;
        elementCellInflow[j + i * m_totalCellsPerElement] = elementCell->externalInflow;
        elementCellChannelInflow[j + i * m_totalCellsPerElement] = elementCell->channelInflow;
        elementCellChannelInflowFlux[j + i * m_totalCellsPerElement] = elementCell->channelInflowFlux;
        totalElementCellMassBal[j + i * m_totalCellsPerElement] = elementCell->totalMassBalance;
        temperature[j + i * m_totalCellsPerElement] = elementCell->temperature.value;
        elementHeatFlux[i] += elementCell->externalHeatFluxes;
        elementCellHeatFlux[j + i * m_totalCellsPerElement] = elementCell->externalHeatFluxes;
        elementCellChannelHeatFlux[j + i * m_totalCellsPerElement] = elementCell->channelHeatRate;

        for(int k = 0; k < 4; k++)
        {
          edgeFlow[k + j * 4 + i * 4 * m_totalCellsPerElement] = elementCell->edgeFlows[k] * dir[k];
          edgeHeatFlow[k + j * 4 + i * 4 * m_totalCellsPerElement] = elementCell->edgeHeatFluxes[k] * dir[k];
        }

        for (int k = 0; k < m_numSolutes; k++)
        {
          solutes[j + i * m_totalCellsPerElement + k * m_totalCellsPerElement * m_elements.size()] = elementCell->soluteConcs[k].value;
          cellChannelSoluteFlux[j + i * m_totalCellsPerElement + k * m_elements.size() * m_totalCellsPerElement] = elementCell->channelSoluteRate[k];

          for(int l = 0; l < 4; l++)
          {
            edgeSolutesFlow[ l + j * 4 + i * 4 * m_totalCellsPerElement + k * 4 * m_totalCellsPerElement * m_elements.size()] = elementCell->edgeSoluteConcFluxes[k][l] * dir[k];
          }
        }

//        ThreadSafeNcVar elementCellFaceSoluteFlowVar =  m_outputNetCDF->addVar("element_cell_face_solute_flow", "float",
//                                                                         std::vector<std::string>({"time","solutes", "elements","element_cells","element_cell_face"}));
      }

      for (int k = 0; k < m_numSolutes; k++)
      {
        channelSoluteFlux[i  + k * m_elements.size()] = element->channelSoluteRate[k];
      }
    }

    //    size_t size = sizeof(float) * m_totalCellsPerElement;
    m_outNetCDFVariables["hydraulic_head"].putVar(std::vector<size_t>({currentTime, 0, 0}), std::vector<size_t>({1, m_elements.size(), (size_t)m_totalCellsPerElement}), hydHead);

    m_outNetCDFVariables["dvolume_dt"].putVar(std::vector<size_t>({currentTime, 0, 0}), std::vector<size_t>({1, m_elements.size(), (size_t)m_totalCellsPerElement}), dvolumedt);

    m_outNetCDFVariables["saturated_depth"].putVar(std::vector<size_t>({currentTime, 0, 0}), std::vector<size_t>({1, m_elements.size(), (size_t)m_totalCellsPerElement}), satDepth);

    m_outNetCDFVariables["element_external_inflow"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), elementInflow);

    m_outNetCDFVariables["element_cell_external_inflow"].putVar(std::vector<size_t>({currentTime, 0, 0}), std::vector<size_t>({1, m_elements.size(), (size_t)m_totalCellsPerElement}), elementCellInflow);

    m_outNetCDFVariables["element_channel_inflow"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), elementChannelInflow);

    m_outNetCDFVariables["element_channel_flux"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), elementChannelInflowFlux);

    m_outNetCDFVariables["element_cell_channel_inflow"].putVar(std::vector<size_t>({currentTime, 0, 0}), std::vector<size_t>({1, m_elements.size(), (size_t)m_totalCellsPerElement}), elementCellChannelInflow);

    m_outNetCDFVariables["element_cell_channel_flux"].putVar(std::vector<size_t>({currentTime, 0, 0}), std::vector<size_t>({1, m_elements.size(), (size_t)m_totalCellsPerElement}), elementCellChannelInflowFlux);

    m_outNetCDFVariables["total_element_cell_mass_balance"].putVar(std::vector<size_t>({currentTime, 0, 0}), std::vector<size_t>({1, m_elements.size(), (size_t)m_totalCellsPerElement}), totalElementCellMassBal);

    m_outNetCDFVariables["total_mass_balance"].putVar(std::vector<size_t>({currentTime}), std::vector<size_t>({1}), &m_totalMassBalance);

    m_outNetCDFVariables["element_cell_face_flow"].putVar(std::vector<size_t>({currentTime, 0, 0, 0}), std::vector<size_t>({1, m_elements.size(), (size_t)m_totalCellsPerElement, 4}), edgeFlow);

    m_outNetCDFVariables["temperature"].putVar(std::vector<size_t>({currentTime, 0, 0}), std::vector<size_t>({1, m_elements.size(), (size_t)m_totalCellsPerElement}), temperature);

    m_outNetCDFVariables["element_cell_face_heat_flow"].putVar(std::vector<size_t>({currentTime, 0, 0, 0}), std::vector<size_t>({1, m_elements.size(), (size_t)m_totalCellsPerElement, 4}), edgeHeatFlow);

    m_outNetCDFVariables["element_external_heat_flux"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), elementHeatFlux);

    m_outNetCDFVariables["element_cell_external_heat_flux"].putVar(std::vector<size_t>({currentTime, 0, 0}), std::vector<size_t>({1, m_elements.size(), (size_t)m_totalCellsPerElement}), elementCellHeatFlux);

    m_outNetCDFVariables["element_channel_heat_flux"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), elementChannelHeatFlux);

    m_outNetCDFVariables["element_cell_channel_heat_flux"].putVar(std::vector<size_t>({currentTime, 0, 0}), std::vector<size_t>({1, m_elements.size(), (size_t)m_totalCellsPerElement}), elementCellChannelHeatFlux);

    m_outNetCDFVariables["element_channel_wse"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), elementChannelWSE);

    m_outNetCDFVariables["element_channel_flow_top_width"].putVar(std::vector<size_t>({currentTime, 0}), std::vector<size_t>({1, m_elements.size()}), elementChannelWidth);


    if(m_numSolutes)
    {
      m_outNetCDFVariables["solute_concentration"].putVar(std::vector<size_t>({currentTime, 0, 0, 0}), std::vector<size_t>({1, (size_t)m_numSolutes, m_elements.size(), (size_t)m_totalCellsPerElement}), solutes);

      m_outNetCDFVariables["element_cell_channel_solute_flux"].putVar(std::vector<size_t>({currentTime, 0, 0, 0}), std::vector<size_t>({1, (size_t)m_numSolutes, m_elements.size(), (size_t)m_totalCellsPerElement}), cellChannelSoluteFlux);

      m_outNetCDFVariables["element_cell_face_solute_flow"].putVar(std::vector<size_t>({currentTime, 0, 0, 0, 0}), std::vector<size_t>({1, (size_t)m_numSolutes, m_elements.size(), (size_t)m_totalCellsPerElement, 4}), edgeSolutesFlow);

      m_outNetCDFVariables["element_channel_solute_flux"].putVar(std::vector<size_t>({currentTime, 0, 0}), std::vector<size_t>({1, (size_t)m_numSolutes, m_elements.size()}), channelSoluteFlux);
    }

    if(m_simulateWaterAge)
    {

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      for (int i = 0; i < (int)m_elements.size(); i++)
      {
        Element *element = m_elements[i];

        for(int j = 0; j < m_totalCellsPerElement; j++)
        {
          ElementCell *elementCell = element->elementCells[j];
          waterAge[j + i * m_totalCellsPerElement] = elementCell->soluteConcs[m_numSolutes].value;
        }
      }

      m_outNetCDFVariables["water_age"].putVar(std::vector<size_t>({currentTime, 0, 0}), std::vector<size_t>({1, m_elements.size(), (size_t)m_totalCellsPerElement}), waterAge);
    }

    if(m_flushToDisk)
    {
      m_outputNetCDF->sync();
    }


    delete[] hydHead;
    delete[] dvolumedt;
    delete[] satDepth;
    delete[] elementInflow;
    delete[] elementCellInflow;
    delete[] elementChannelInflow;
    delete[] elementChannelInflowFlux;
    delete[] elementCellChannelInflow;
    delete[] elementCellChannelInflowFlux;
    delete[] totalElementCellMassBal;
    delete[] edgeFlow;
    delete[] temperature;
    delete[] edgeHeatFlow;
    delete[] waterAge;
    delete[] elementHeatFlux;
    delete[] elementCellHeatFlux;
    delete[] elementChannelHeatFlux;
    delete[] elementCellChannelHeatFlux;
    delete[] elementChannelWSE;
    delete[] elementChannelWidth;
    delete[] solutes;
    delete[] edgeSolutesFlow;
    delete[] cellChannelSoluteFlux;
    delete[] channelSoluteFlux;
  }

#endif
}

void GWModel::closeOutputFiles()
{
  closeOutputNetCDFFile();
}

void GWModel::closeOutputNetCDFFile()
{
#ifdef USE_NETCDF

  if(m_outputNetCDF)
  {
    m_outputNetCDF->sync();
    delete m_outputNetCDF;
    m_outputNetCDF = nullptr;
  }

#endif
}

QFileInfo GWModel::relativePathToAbsolute(const QFileInfo &fileInfo)
{
  if (fileInfo.isRelative())
  {
    if (!m_inputFile.filePath().isEmpty() &&
        !m_inputFile.filePath().isNull() &&
        QFile::exists(m_inputFile.absoluteFilePath()))
    {
      QFileInfo absoluteFilePath = m_inputFile.absoluteDir().absoluteFilePath(fileInfo.filePath());

      if (absoluteFilePath.absoluteDir().exists())
      {
        return absoluteFilePath;
      }
    }
  }

  return fileInfo;
}

const unordered_map<string, int> GWModel::m_inputFileFlags({
                                                             {"[OPTIONS]", 1},
                                                             {"[OUTPUTS]", 2},
                                                             {"[SOLUTES]", 3},
                                                             {"[ELEMENTJUNCTIONS]", 4},
                                                             {"[ELEMENTS]", 5},
                                                             {"[ELEMENT_CELL_WIDTHS]", 6},
                                                             {"[INIT_HYDRAULICS]", 7},
                                                             {"[INIT_CONDITIONS]", 8},
                                                             {"[BOUNDARY_CONDITIONS]", 9},
                                                             {"[TIMESERIES]", 10},
                                                             {"[SOURCES]", 11},
                                                             {"[CHANNEL_BOUNDARY_CONDITIONS]", 12}
                                                           });

const unordered_map<std::string, int> GWModel::m_hydraulicVariableFlags({
                                                                          {"HYD_HEAD", 1},
                                                                          {"BOTTOM_ELEV", 2},
                                                                          {"TOP_ELEV", 3},
                                                                          {"HYD_COND_X", 4},
                                                                          {"HYD_COND_Y", 5},
                                                                          {"SPECIFIC_STORAGE", 6},
                                                                          {"POROSITY", 7},
                                                                          {"SED_DENSITY", 8},
                                                                          {"SED_SPECIFIC_HEAT_CAPACITY", 9},
                                                                          {"DISPERSIVITY_X", 9},
                                                                          {"DISPERSIVITY_Y", 10},
                                                                          {"WATER_THERMAL_CONDUCTIVITY", 11},
                                                                          {"SED_THERMAL_CONDUCTIVITY", 12},
                                                                        });

const unordered_map<string, int> GWModel::m_optionsFlags({
                                                           {"START_DATETIME", 1},
                                                           {"END_DATETIME", 2},
                                                           {"REPORT_INTERVAL", 3},
                                                           {"MAX_TIME_STEP", 4},
                                                           {"MIN_TIME_STEP", 5},
                                                           {"NUM_INITIAL_FIXED_STEPS", 6},
                                                           {"USE_ADAPTIVE_TIME_STEP", 7},
                                                           {"TIME_STEP_RELAXATION_FACTOR", 8},
                                                           {"SOLVER", 9},
                                                           {"SOLVER_ABS_TOL", 10},
                                                           {"SOLVER_REL_TOL", 11},
                                                           {"SOLVE_TEMP", 12},
                                                           {"WATER_DENSITY", 13},
                                                           {"WATER_SPECIFIC_HEAT_CAPACITY", 14},
                                                           {"DEFAULT_SED_DENSITY", 15},
                                                           {"DEFAULT_SED_SPECIFIC_HEAT_CAPACITY", 16},
                                                           {"DEFAULT_POROSITY", 17},
                                                           {"DEFAULT_HYD_COND_X", 18},
                                                           {"DEFAULT_HYD_COND_Y", 19},
                                                           {"DEFAULT_SPECIFIC_STORAGE", 20},
                                                           {"NUM_LEFT_CELLS", 21},
                                                           {"NUM_RIGHT_CELLS", 22},
                                                           {"NUM_SOLUTES", 23},
                                                           {"VERBOSE", 24},
                                                           {"FLUSH_TO_DISK_FREQ", 25},
                                                           {"PRINT_FREQ", 26},
                                                           {"DEFAULT_CELL_WIDTH", 27},
                                                           {"DISPERSIVITY_X", 28},
                                                           {"DISPERSIVITY_Y", 29},
                                                           {"ADVECTION_MODE", 30},
                                                           {"TVD_FLUX_LIMITER", 31},
                                                           {"WATER_THERMAL_CONDUCTIVITY", 32},
                                                           {"SED_THERMAL_CONDUCTIVITY", 33},
                                                           {"DEFAULT_HYD_COND_Z", 34},
                                                           {"LINEAR_SOLVER", 35},
                                                           {"SIMULATE_WATER_AGE", 36},
                                                           {"DISPERSIVITY_Z", 37},
                                                         });

const unordered_map<string, int> GWModel::m_advectionFlags({
                                                             {"UPWIND", 1},
                                                             {"CENTRAL", 2},
                                                             {"HYBRID", 3},
                                                             {"TVD", 4},
                                                           });

const unordered_map<string, int> GWModel::m_BCFlags({
                                                      {"HYD_HEAD", 1},
                                                      {"GRAD_HYD_HEAD", 2},
                                                      {"TEMPERATURE", 3},
                                                      {"GRAD_TEMPERATURE", 4},
                                                      {"SOLUTE", 5},
                                                      {"GRAD_SOLUTE", 6},
                                                    });

const unordered_map<string, int> GWModel::m_srcFlags({
                                                       {"FLOW", 1},
                                                       {"HEAT", 2},
                                                       {"SOLUTE", 3},
                                                     });

const unordered_map<string, int> GWModel::m_geomMultFlags({
                                                            {"NONE", 0},
                                                            {"AREA", 1},
                                                            {"VOLUME", 2},
                                                            {"LENGTH", 3},
                                                            {"WIDTH", 4},
                                                            {"SATURATED_THICKNESS", 5},
                                                          });

const unordered_map<string, int> GWModel::m_chanVarFlags({
                                                           {"WIDTH", 1},
                                                           {"WSE", 2},
                                                           {"TEMPERATURE", 3},
                                                           {"SOLUTE", 4},
                                                         });


const unordered_map<string, int> GWModel::m_edgeFlags({
                                                        {"RIGHT", 0},
                                                        {"DOWN", 1},
                                                        {"LEFT", 2},
                                                        {"UP", 3}
                                                      });

const unordered_map<string, int> GWModel::m_solverTypeFlags({{"RK4", 1},
                                                             {"RKQS", 2},
                                                             {"ADAMS", 3},
                                                             {"BDF", 4},
                                                             {"EULER", 5}
                                                            });

const unordered_map<string, int> GWModel::m_linearSolverTypeFlags({{"GMRES", 1},
                                                                   {"FGMRES", 2},
                                                                   {"Bi_CGStab", 3},
                                                                   {"TFQMR", 4},
                                                                   {"PCG", 5}
                                                                  });

const QRegExp GWModel::m_dateTimeDelim("(\\,|\\t|\\\n|\\/|\\s+|\\:)");

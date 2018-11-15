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
#include "timeseries.h"


#include <QDir>
#include <QDate>

using namespace std;

#ifdef USE_NETCDF

using namespace netCDF;
using namespace netCDF::exceptions;

#endif

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

    printf("GWModel TimeStep (s): %f\tDateTime: %f\tHead (m) { Iters: %i/%i\tMin: %f\tMax: %f\tTotalMassBalance: %g (m^3/s)}", m_timeStep, m_currentDateTime,
           m_hydHeadSolver->getIterations(), m_hydHeadSolver->maxIterations(), m_minHead, m_maxHead, m_totalMassBalance);

    if(m_solveHeatTransport)
    {
      printf("\tTemperature (°C) { Iters: %i/%i\tMin: %f\tMax: %f\tTotalMassBalance: %g (KJ)}",
             m_heatSolver->getIterations(), m_heatSolver->maxIterations(), m_minTemp, m_maxTemp, m_totalHeatBalance);
    }


    for (size_t j = 0; j < m_solutes.size(); j++)
    {
      std::string &solute = m_solutes[j];
      ODESolver *solver = m_soluteSolvers[j];
      printf("\t%s (kg/m) { Iters: %i/%i\tMin: %f\tMax: %f\tTotalMassBalance: %g (kg)}", solute.c_str(), solver->getIterations(), solver->maxIterations(),
             m_minSolute[j], m_maxSolute[j], m_totalSoluteMassBalance[j]);
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
                readSuccess = readInputFileElementCellBCHydHeadTag(line, error);
                break;
              case 9:
                readSuccess = readInputFileElementCellBCHydHeadDerivTag(line, error);
                break;
              case 10:
                readSuccess = readInputFileElementCellBCFluxTag(line, error);
                break;
              case 11:
                readSuccess = readInputFileElementCellHeadDepFluxTag(line, error);
                break;
              case 12:
                readSuccess = readInputFileElementCellTempInitTag(line, error);
                break;
              case 13:
                readSuccess = readInputFileElementCellBCTempTag(line, error);
                break;
              case 14:
                readSuccess = readInputFileElementCellSoluteInitTag(line, error);
                break;
              case 15:
                readSuccess = readInputFileElementCellBCSoluteTag(line, error);
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
    ThreadSafeNcDim solutesDim =  m_outputNetCDF->addDim("solutes", m_solutes.size());
    ThreadSafeNcVar solutes =  m_outputNetCDF->addVar("solute_names", NcType::nc_STRING, solutesDim);
    solutes.putAtt("long_name", "Solutes");
    m_outNetCDFVariables["solutes"] = solutes;

    if (m_solutes.size())
    {
      char **soluteNames = new char *[m_solutes.size()];

      for (size_t i = 0; i < m_solutes.size(); i++)
      {
        string soluteName = m_solutes[i];
        soluteNames[i] = new char[soluteName.size() + 1];
        strcpy(soluteNames[i], soluteName.c_str());
      }

      solutes.putVar(soluteNames);

      for (size_t i = 0; i < m_solutes.size(); i++)
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

    for(int i = 0; i < m_elements.size(); i++)
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

    ThreadSafeNcVar tempVar =  m_outputNetCDF->addVar("temperature", "float",
                                                      std::vector<std::string>({"time", "elements", "element_cells"}));
    tempVar.putAtt("long_name", "Temperature");
    tempVar.putAtt("units", "°C");
    m_outNetCDFVariables["temperature"] = tempVar;


    ThreadSafeNcVar solutesVar =  m_outputNetCDF->addVar("solute_concentration", "float",
                                                         std::vector<std::string>({"time", "solutes", "elements", "element_cells"}));
    solutesVar.putAtt("long_name", "Solute Concentration");
    solutesVar.putAtt("units", "kg/m^3");
    m_outNetCDFVariables["solute_concentration"] = solutesVar;


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

            int heatSolverMode = -1;

            if (it != m_solverTypeFlags.end())
              heatSolverMode = it->second;

            switch (heatSolverMode)
            {
              case 1:
                m_hydHeadSolver->setSolverType(ODESolver::RK4);
                break;
              case 2:
                m_hydHeadSolver->setSolverType(ODESolver::RKQS);
                break;
              case 3:
                m_hydHeadSolver->setSolverType(ODESolver::CVODE_ADAMS);
                break;
              case 4:
                m_hydHeadSolver->setSolverType(ODESolver::CVODE_BDF);
                break;
              case 5:
                m_hydHeadSolver->setSolverType(ODESolver::EULER);
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
              m_hydHeadSolver->setAbsoluteTolerance(abs_tol);

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
              m_hydHeadSolver->setRelativeTolerance(rel_tol);

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
            std::string code = options[1].toUpper().toStdString();
            auto it = m_solverTypeFlags.find(code);

            int heatSolverMode = -1;

            if (it != m_solverTypeFlags.end())
              heatSolverMode = it->second;

            switch (heatSolverMode)
            {
              case 1:
                m_heatSolver->setSolverType(ODESolver::RK4);
                break;
              case 2:
                m_heatSolver->setSolverType(ODESolver::RKQS);
                break;
              case 3:
                m_heatSolver->setSolverType(ODESolver::CVODE_ADAMS);
                break;
              case 4:
                m_heatSolver->setSolverType(ODESolver::CVODE_BDF);
                break;
              case 5:
                m_heatSolver->setSolverType(ODESolver::EULER);
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
            errorMessage = "Temeprature solver type error";
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
            double abs_tol = options[1].toDouble(&ok);

            if (ok)
              m_heatSolver->setAbsoluteTolerance(abs_tol);

            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Heat solver absolute tolerance error";
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
            double rel_tol = options[1].toDouble(&ok);

            if (ok)
              m_heatSolver->setRelativeTolerance(rel_tol);

            foundError = !ok;
          }
          else
          {
            foundError = true;
          }

          if (foundError)
          {
            errorMessage = "Heat solver relative tolerance error";
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
      case 17:
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
      case 18:
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
      case 19:
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
      case 20:
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
      case 21:
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
      case 22:
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
      case 23:
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
      case 24:
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
      case 25:
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
      case 26:
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
      case 27:
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
      case 28:
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
      case 29:
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
      case 30:
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

  if (columns.size() == 4)
  {
    bool foundError = false;

    if (m_addedSoluteCount < (int)m_solutes.size())
    {
      m_solutes[m_addedSoluteCount] = columns[0].toStdString();

      std::string solverType = columns[1].toStdString();
      auto it = m_solverTypeFlags.find(solverType);

      if (it != m_solverTypeFlags.end())
      {
        int solverTypeCode = it->second;

        switch (solverTypeCode)
        {
          case 1:
            m_soluteSolvers[m_addedSoluteCount]->setSolverType(ODESolver::RK4);
            break;
          case 2:
            m_soluteSolvers[m_addedSoluteCount]->setSolverType(ODESolver::RKQS);
            break;
          case 3:
            m_soluteSolvers[m_addedSoluteCount]->setSolverType(ODESolver::CVODE_ADAMS);
            break;
          case 4:
            m_soluteSolvers[m_addedSoluteCount]->setSolverType(ODESolver::CVODE_BDF);
            break;
          case 5:
            m_soluteSolvers[m_addedSoluteCount]->setSolverType(ODESolver::EULER);
            break;
          default:
            foundError = true;
            break;
        }

        if (foundError)
        {
          errorMessage = "Solute error";
          return false;
        }

        bool parsed;
        double abs_tol = columns[2].toDouble(&parsed);

        if (parsed)
        {
          m_soluteSolvers[m_addedSoluteCount]->setAbsoluteTolerance(abs_tol);
        }
        else
        {
          errorMessage = "Solute absolute tolerance error";
          return false;
        }

        double rel_tol = columns[3].toDouble(&parsed);

        if (parsed)
        {
          m_soluteSolvers[m_addedSoluteCount]->setRelativeTolerance(rel_tol);
        }
        else
        {
          errorMessage = "Solute relative tolerance error";
          return false;
        }
      }
      else
      {
        errorMessage = "Solute error";
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

  if (columns.size() > 10)
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

      bool topWidthOk ;
      double topWidth = columns[7].toDouble(&topWidthOk);

      bool depthOk ;
      double depth = columns[8].toDouble(&depthOk);

      bool bedThicknessOk ;
      double bedThickness = columns[9].toDouble(&bedThicknessOk);

      bool channelKOk ;
      double channelK = columns[10].toDouble(&channelKOk);

      if (lengthOk && hydHeadOk && topWidthOk && bottomElOk && topElOk &&
          depthOk && bedThicknessOk && channelKOk)
      {
        Element *element = addElement(id.toStdString(), ej1, ej2);
        element->length = length;
        element->channelWSE = element->z + depth;
        element->channelBedThickness = bedThickness;
        element->channelWidth = topWidth;
        element->channelBedHydCond = channelK;

        for(int j = 0 ; j < m_totalCellsPerElement; j++)
        {
          ElementCell *elementCell = element->elementCells[j];
          elementCell->hydHead.value = elementCell->prevHydHead.value = hydHead;
          elementCell->bedRockElev = bottomEl;
          elementCell->topElev = topEl;
        }

        int currentColumn = 11;

        if(m_solveHeatTransport)
        {
          if(columns.size() > 11)
          {
            bool tempOk;
            double temp = columns[11].toDouble(&tempOk);

            bool chanTempOk;
            double chanTemp = columns[12].toDouble(&chanTempOk);

            if(tempOk && chanTempOk)
            {
              element->channelTemperature = chanTemp;
              currentColumn = 13;

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

        if(columns.size() ==  currentColumn - 1 + m_solutes.size() * 2)
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
                element->elementCells[j]->hydConX = values[j];
              }
            }
            break;
          case 5:
            {
              for(int j = 0; j < m_numLeftCellsPerElement + m_numRightCellsPerElement; j++)
              {
                element->elementCells[j]->hydConY = values[j];
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

bool GWModel::readInputFileElementCellBCHydHeadTag(const QString &line, QString &errorMessage)
{

}

bool GWModel::readInputFileElementCellBCHydHeadDerivTag(const QString &line, QString &errorMessage)
{

}

bool GWModel::readInputFileElementCellBCFluxTag(const QString &line, QString &errorMessage)
{

}

bool GWModel::readInputFileElementCellHeadDepFluxTag(const QString &line, QString &errorMessage)
{

}

bool GWModel::readInputFileElementCellTempInitTag(const QString &line, QString &errorMessage)
{

}

bool GWModel::readInputFileElementCellBCTempTag(const QString &line, QString &errorMessage)
{

}

bool GWModel::readInputFileElementCellSoluteInitTag(const QString &line, QString &errorMessage)
{

}

bool GWModel::readInputFileElementCellBCSoluteTag(const QString &line, QString &errorMessage)
{

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
    float *temperature = new float[m_elements.size() * m_totalCellsPerElement];
    float *solutes = new float[m_elements.size() * m_solutes.size() * m_totalCellsPerElement];


#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < (int)m_elements.size(); i++)
    {
      Element *element = m_elements[i];

      for(int j = 0; j < m_totalCellsPerElement; j++)
      {
        ElementCell *elementCell = element->elementCells[j];

        hydHead[j + i * m_totalCellsPerElement] = elementCell->hydHead.value;
        temperature[j + i * m_totalCellsPerElement] = elementCell->temperature.value;

        for (size_t k = 0; k < m_solutes.size(); k++)
        {
          solutes[j + i * m_totalCellsPerElement + k * m_totalCellsPerElement * m_elements.size()] = elementCell->soluteConcs[k].value;
        }
      }
    }

//    size_t size = sizeof(float) * m_totalCellsPerElement;
    m_outNetCDFVariables["hydraulic_head"].putVar(std::vector<size_t>({currentTime, 0, 0}), std::vector<size_t>({1, m_elements.size(), (size_t)m_totalCellsPerElement}), hydHead);

    m_outNetCDFVariables["temperature"].putVar(std::vector<size_t>({currentTime, 0, 0}), std::vector<size_t>({1, m_elements.size(), (size_t)m_totalCellsPerElement}), temperature);

    if(m_solutes.size())
    {
      m_outNetCDFVariables["solute_concentration"].putVar(std::vector<size_t>({currentTime, 0, 0, 0}), std::vector<size_t>({1, m_solutes.size(), m_elements.size(), (size_t)m_totalCellsPerElement}), solutes);
    }

    delete[] hydHead;
    delete[] temperature;
    delete[] solutes;

    if(m_flushToDisk)
    {
      m_outputNetCDF->sync();
    }
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
                                                             {"[BC_HYD_HEAD]", 8},
                                                             {"[BC_HYD_HEAD_DERIV]", 9},
                                                             {"[BC_FLUX]", 10},
                                                             {"[BC_HEAD_DEPENDENT_FLUX]", 11},
                                                             {"[TEMP_INIT]", 12},
                                                             {"[BC_TEMP]", 13},
                                                             {"[SOLUTE_INIT]", 14},
                                                             {"[BC_SOLUTE]", 15},
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
                                                           {"HYD_HEAD_SOLVER", 9},
                                                           {"HYD_HEAD_SOLVER_ABS_TOL", 10},
                                                           {"HYD_HEAD_SOLVER_REL_TOL", 11},
                                                           {"SOLVE_TEMP", 12},
                                                           {"TEMP_SOLVER", 13},
                                                           {"TEMP_SOLVER_ABS_TOL", 14},
                                                           {"TEMP_SOLVER_REL_TOL", 15},
                                                           {"WATER_DENSITY", 16},
                                                           {"WATER_SPECIFIC_HEAT_CAPACITY", 17},
                                                           {"DEFAULT_SED_DENSITY", 18},
                                                           {"DEFAULT_SED_SPECIFIC_HEAT_CAPACITY", 19},
                                                           {"DEFAULT_POROSITY", 20},
                                                           {"DEFAULT_HYD_COND_X", 21},
                                                           {"DEFAULT_HYD_COND_Y", 22},
                                                           {"DEFAULT_SPECIFIC_STORAGE", 23},
                                                           {"NUM_LEFT_CELLS", 24},
                                                           {"NUM_RIGHT_CELLS", 25},
                                                           {"NUM_SOLUTES", 26},
                                                           {"VERBOSE", 27},
                                                           {"FLUSH_TO_DISK_FREQ", 28},
                                                           {"PRINT_FREQ", 29},
                                                           {"DEFAULT_CELL_WIDTH", 30},
                                                         });

const unordered_map<string, int> GWModel::m_solverTypeFlags({{"RK4", 1},
                                                             {"RKQS", 2},
                                                             {"ADAMS", 3},
                                                             {"BDF", 4},
                                                             {"EULER", 5}
                                                            });

const QRegExp GWModel::m_dateTimeDelim("(\\,|\\t|\\\n|\\/|\\s+|\\:)");

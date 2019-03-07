/*!
*  \file    GWComponent.cpp
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
#include "gwcomponent.h"
#include "core/valuedefinition.h"
#include "spatial/linestring.h"
#include "core/dimension.h"
#include "core/abstractoutput.h"
#include "core/idbasedargument.h"
#include "gwmodel.h"
#include "progresschecker.h"
#include "temporal/timedata.h"
#include "element.h"
#include "elementjunction.h"
#include "spatial/point.h"
#include "elementinput.h"
#include "elementoutput.h"
#include "core/unit.h"
#include "core/unitdimensions.h"

using namespace HydroCouple;
using namespace HydroCouple::Temporal;

GWComponent::GWComponent(const QString &id, GWComponentInfo *modelComponentInfo)
  : AbstractTimeModelComponent(id, modelComponentInfo),
    m_GWComponentInfo(modelComponentInfo),
    m_modelInstance(nullptr),
    m_parent(nullptr),
    m_inputFilesArgument(nullptr),
    m_timeDimension(nullptr),
    m_geometryDimension(nullptr)
{
  m_timeDimension = new Dimension("TimeDimension",this);
  m_geometryDimension = new Dimension("ElementGeometryDimension", this);


  m_heatFluxUnit = new Unit(this);
  m_heatFluxUnit->setCaption("Heat Source (W or J/s)");
  m_heatFluxUnit->setConversionFactorToSI(1.0);
  m_heatFluxUnit->setOffsetToSI(0.0);
  m_heatFluxUnit->dimensionsInternal()->setPower(HydroCouple::Mass, 1.0);
  m_heatFluxUnit->dimensionsInternal()->setPower(HydroCouple::Length, 2.0);
  m_heatFluxUnit->dimensionsInternal()->setPower(HydroCouple::Time, -3.0);

  m_temperatureUnit = new Unit(this);
  m_temperatureUnit->setCaption("Temperature (°C)");
  m_temperatureUnit->setConversionFactorToSI(1.0);
  m_temperatureUnit->setOffsetToSI(273.15);
  m_temperatureUnit->dimensionsInternal()->setPower(HydroCouple::Temperature, 1.0);

  m_soluteFluxUnit = new Unit(this);
  m_soluteFluxUnit->setCaption("Mass Flux (kg/s)");
  m_soluteFluxUnit->dimensionsInternal()->setPower(HydroCouple::Mass, 1.0);
  m_soluteFluxUnit->dimensionsInternal()->setPower(HydroCouple::Time, -1.0);

  m_soluteUnit = new Unit(this);
  m_soluteUnit->setCaption("Concentration (kg/m^3)");
  m_soluteUnit->dimensionsInternal()->setPower(HydroCouple::Mass, 1.0);
  m_soluteUnit->dimensionsInternal()->setPower(HydroCouple::Length, -3.0);

  m_soluteConcQuantity = new Quantity(QVariant::Double, m_soluteUnit, this);
  m_soluteConcFluxQuantity = new Quantity(QVariant::Double, m_soluteFluxUnit, this);

  createArguments();
}

GWComponent::~GWComponent()
{

}

QList<QString> GWComponent::validate()
{
  if(isInitialized())
  {
    setStatus(IModelComponent::Validating,"Validating...");

    //check connections

    setStatus(IModelComponent::Valid,"");
  }
  else
  {
    //throw has not been initialized yet.
  }

  return QList<QString>();
}

void GWComponent::prepare()
{
  if(!isPrepared() && isInitialized() && m_modelInstance)
  {
    for(auto output :  outputsInternal())
    {
      for(auto adaptedOutput : output->adaptedOutputs())
      {
        adaptedOutput->initialize();
      }
    }

    updateOutputValues(QList<HydroCouple::IOutput*>());

    setStatus(IModelComponent::Updated ,"Finished preparing model");
    setPrepared(true);
  }
  else
  {
    setPrepared(false);
    setStatus(IModelComponent::Failed ,"Error occured when preparing model");
  }
}

void GWComponent::update(const QList<HydroCouple::IOutput *> &requiredOutputs)
{
  if(status() == IModelComponent::Updated)
  {
    setStatus(IModelComponent::Updating);

    double minConsumerTime = std::max(m_modelInstance->currentDateTime(), getMinimumConsumerTime());

    while (m_modelInstance->currentDateTime() <= minConsumerTime &&
           m_modelInstance->currentDateTime() < m_modelInstance->endDateTime())
    {
      m_modelInstance->update();

      if(progressChecker()->performStep(m_modelInstance->currentDateTime()))
      {
        setStatus(IModelComponent::Updated , "Simulation performed time-step | DateTime: " + QString::number(m_modelInstance->currentDateTime(), 'f') , progressChecker()->progress());
      }
    }

    updateOutputValues(requiredOutputs);

    currentDateTimeInternal()->setJulianDay(m_modelInstance->currentDateTime());

    if(m_modelInstance->currentDateTime() >=  m_modelInstance->endDateTime())
    {
      setStatus(IModelComponent::Done , "Simulation finished successfully", 100);
    }
    else
    {
      if(progressChecker()->performStep(m_modelInstance->currentDateTime()))
      {
        setStatus(IModelComponent::Updated , "Simulation performed time-step | DateTime: " + QString::number(m_modelInstance->currentDateTime(), 'f') , progressChecker()->progress());
      }
      else
      {
        setStatus(IModelComponent::Updated);
      }
    }
  }
}

void GWComponent::finish()
{
  if(isPrepared())
  {
    setStatus(IModelComponent::Finishing , "GWComponent with id " + id() + " is being disposed" , 100);

    std::list<std::string> errors;
    m_modelInstance->finalize(errors);
    initializeFailureCleanUp();

    m_channelSoluteInputs.clear();
    m_channelOutflowSoluteOutputs.clear();

    setPrepared(false);
    setInitialized(false);

    setStatus(IModelComponent::Finished , "GWComponent with id " + id() + " has been disposed" , 100);
    setStatus(IModelComponent::Created , "GWComponent with id " + id() + " ran successfully and has been re-created" , 100);
  }
}

GWModel *GWComponent::modelInstance() const
{
  return m_modelInstance;
}

ICloneableModelComponent *GWComponent::parent() const
{
  return m_parent;
}

ICloneableModelComponent *GWComponent::clone()
{
  return m_parent;
}

QList<ICloneableModelComponent*> GWComponent::clones() const
{
  return m_clones;
}

void GWComponent::initializeFailureCleanUp()
{
  if(m_modelInstance)
  {
    delete m_modelInstance;
    m_modelInstance = nullptr;
  }
}

void GWComponent::createArguments()
{
  createInputFileArguments();
}

void GWComponent::createInputFileArguments()
{
  QStringList fidentifiers;
  fidentifiers.append("Input File");
  fidentifiers.append("Output NetCDF File");

  Quantity *fquantity = Quantity::unitLessValues("InputFilesQuantity", QVariant::String, this);
  fquantity->setDefaultValue("");
  fquantity->setMissingValue("");

  Dimension *dimension = new Dimension("IdDimension","Dimension for identifiers",this);

  m_inputFilesArgument = new IdBasedArgumentString("FileIO", fidentifiers, dimension, fquantity, this);
  m_inputFilesArgument->setCaption("Model Input/Output Files");
  m_inputFilesArgument->addFileFilter("Input File (*.inp)");
  m_inputFilesArgument->addFileFilter("NetCDF File (*.nc)");
  m_inputFilesArgument->setMatchIdentifiersWhenReading(true);

  addArgument(m_inputFilesArgument);
}

bool GWComponent::initializeArguments(QString &message)
{
  bool initialized = true;

  initialized = initializeInputFilesArguments(message);

  if(initialized)
  {
    createGeometries();
  }

  return initialized;
}

bool GWComponent::initializeInputFilesArguments(QString &message)
{

  QString inputFilePath = (*m_inputFilesArgument)["Input File"];
  QFileInfo inputFile = getAbsoluteFilePath(inputFilePath);

  if(inputFile.exists())
  {
    initializeFailureCleanUp();

    m_modelInstance = new GWModel(this);
    m_modelInstance->setInputFile(inputFile);

    QString netCDFOutput = QString((*m_inputFilesArgument)["Output NetCDF File"]);
    if(!netCDFOutput.isEmpty() && !netCDFOutput.isNull())
      m_modelInstance->setOutputNetCDFFile(QFileInfo(netCDFOutput));

    //    QString csvOutput = QString((*m_inputFilesArgument)["Output CSV File"]);
    //    if(!csvOutput.isEmpty() && !csvOutput.isNull())
    //      m_modelInstance->setOutputCSVFile(QFileInfo(csvOutput));

    std::list<std::string> errors;
    bool initialized = m_modelInstance->initialize(errors);

    for (std::string errorMsg : errors)
    {
      message += "/n" + QString::fromStdString(errorMsg);
    }

    if(initialized)
    {
      timeHorizonInternal()->setJulianDay(m_modelInstance->startDateTime());
      timeHorizonInternal()->setDuration(m_modelInstance->endDateTime() - m_modelInstance->startDateTime());
      currentDateTimeInternal()->setJulianDay(m_modelInstance->startDateTime());
      progressChecker()->reset(m_modelInstance->startDateTime(), m_modelInstance->endDateTime());
    }

    return initialized;
  }
  else
  {
    message = "Input file does not exist: " + inputFile.absoluteFilePath();
    return false;
  }

  return true;
}

void GWComponent::createGeometries()
{

  m_elementGeometries.clear();
  m_elementJunctionGeometries.clear();

  for(int i = 0; i < m_modelInstance->numElements() ; i++)
  {
    Element *element = m_modelInstance->getElement(i);
    ElementJunction *from = element->upstreamJunction;
    ElementJunction *to   = element->downstreamJunction;

    HCLineString *lineString = new HCLineString(QString::fromStdString(element->id));
    lineString->setMarker(i);
    HCPoint *p1 = new HCPoint(from->x , from->y, QString::fromStdString(from->id), lineString);
    HCPoint *p2 = new HCPoint(to->x , to->y, QString::fromStdString(to->id), lineString);
    lineString->addPoint(p1);
    lineString->addPoint(p2);

    m_elementJunctionGeometries.push_back(QSharedPointer<HCPoint>(new HCPoint(from->x , from->y, from->z, QString::fromStdString(from->id), nullptr)));
    m_elementJunctionGeometries.push_back(QSharedPointer<HCPoint>(new HCPoint(to->x , to->y, to->z, QString::fromStdString(to->id), nullptr)));

    m_elementGeometries.push_back(QSharedPointer<HCLineString>(lineString));
  }
}

void GWComponent::createInputs()
{
  createChannelWSEInput();

  createChannelWidthInput();

  createChannelTemperatureInput();

  for(int i = 0; i < m_modelInstance->numSolutes() ; i++)
  {
    createChannelSoluteInput(i);
  }

  createWaterAgeFluxInput();
}

void GWComponent::createChannelWSEInput()
{
  Quantity *wseQuantity = Quantity::lengthInMeters(this);

  m_channelWSEInput = new ElementInput("ChannelWSEInput",
                                       m_timeDimension,
                                       m_geometryDimension,
                                       wseQuantity,
                                       ElementInput::ChannelWSE,
                                       this);

  m_channelWSEInput->setCaption("Channel Water Surface Elevation (m)");

  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  m_channelWSEInput->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime()- 1.0/1000000.0, m_channelWSEInput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), m_channelWSEInput);

  m_channelWSEInput->addTime(dt1);
  m_channelWSEInput->addTime(dt2);

  addInput(m_channelWSEInput);
}

void GWComponent::createChannelWidthInput()
{
  Quantity *widthQuantity = Quantity::lengthInMeters(this);

  m_channelWidthInput = new ElementInput("ChannelWidthInput",
                                         m_timeDimension,
                                         m_geometryDimension,
                                         widthQuantity,
                                         ElementInput::ChannelWidth,
                                         this);

  m_channelWidthInput->setCaption("Channel Width (m)");

  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  m_channelWidthInput->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime()- 1.0/1000000.0, m_channelWidthInput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), m_channelWidthInput);

  m_channelWidthInput->addTime(dt1);
  m_channelWidthInput->addTime(dt2);

  addInput(m_channelWidthInput);
}

void GWComponent::createChannelTemperatureInput()
{
  Quantity *temperatureQuantity = new Quantity(QVariant::Double, m_temperatureUnit, this);

  m_channelTemperatureInput = new ElementInput("ChannelTemperatureInput",
                                               m_timeDimension,
                                               m_geometryDimension,
                                               temperatureQuantity,
                                               ElementInput::ChannelTemperature,
                                               this);

  m_channelTemperatureInput->setCaption("Channel Temperature (m)");

  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  m_channelTemperatureInput->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime()- 1.0/1000000.0, m_channelTemperatureInput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), m_channelTemperatureInput);

  m_channelTemperatureInput->addTime(dt1);
  m_channelTemperatureInput->addTime(dt2);

  addInput(m_channelTemperatureInput);
}

void GWComponent::createChannelSoluteInput(int soluteIndex)
{
  QString soluteName = QString::fromStdString(m_modelInstance->solute(soluteIndex));

  ElementInput *channelSoluteInput = new ElementInput("Channel" + soluteName + "Input",
                                                      m_timeDimension,
                                                      m_geometryDimension,
                                                      m_soluteConcQuantity,
                                                      ElementInput::ChannelSolute,
                                                      this);

  channelSoluteInput->setCaption("Channel " + soluteName + " Concentration (kg/m^3)");
  channelSoluteInput->setSoluteIndex(soluteIndex);

  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  channelSoluteInput->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime()- 1.0/1000000.0, channelSoluteInput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), channelSoluteInput);

  channelSoluteInput->addTime(dt1);
  channelSoluteInput->addTime(dt2);

  m_channelSoluteInputs.push_back(channelSoluteInput);
  addInput(channelSoluteInput);
}

void GWComponent::createWaterAgeFluxInput()
{
  if(m_modelInstance->simulateWaterAge())
  {
    Quantity *timeQuantity = Quantity::unitLessValues("Time/Time", QVariant::Double, this);

    ElementInput *soluteFluxInput  = new ElementInput("WaterAgeFluxInput",
                                                      m_timeDimension,
                                                      m_geometryDimension,
                                                      timeQuantity,
                                                      ElementInput::ChannelSolute,
                                                      this);

    soluteFluxInput->setCaption("Element Water Age (days)");
    soluteFluxInput->setSoluteIndex(m_modelInstance->numSolutes());

    QList<QSharedPointer<HCGeometry>> geometries;

    for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
    {
      geometries.append(lineString);
    }

    soluteFluxInput->addGeometries(geometries);

    SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime() - 1.0/1000000.0, soluteFluxInput);
    SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), soluteFluxInput);

    soluteFluxInput->addTime(dt1);
    soluteFluxInput->addTime(dt2);

    addInput(soluteFluxInput);
  }
}

void GWComponent::createOutputs()
{
  createChannelOutflowOutput();

  createChannelOutflowHeatFluxOutput();

  for(int i = 0; i < m_modelInstance->numSolutes() ; i++)
  {
    createChannelOutflowSoluteFluxOutput(i);
  }

  createWaterAgeFluxOutput();
}

void GWComponent::createChannelOutflowOutput()
{
  Quantity *flowQuantity = Quantity::flowInCMS(this);

  m_channelOutflowOutput = new ElementOutput("ChannelOutflow",
                                             m_timeDimension,
                                             m_geometryDimension,
                                             flowQuantity,
                                             ElementOutput::ChannelOutflow,
                                             this);

  m_channelOutflowOutput->setCaption("Channel Outflow (m^3/s)");

  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  m_channelOutflowOutput->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime()- 1.0/1000000.0, m_channelOutflowOutput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), m_channelOutflowOutput);

  m_channelOutflowOutput->addTime(dt1);
  m_channelOutflowOutput->addTime(dt2);

  addOutput(m_channelOutflowOutput);
}

void GWComponent::createChannelOutflowHeatFluxOutput()
{
  Quantity *heatQuantity = new Quantity(QVariant::Double, m_heatFluxUnit, this);

  m_channelOutflowHeatOutput = new ElementOutput("ChannelOutflowHeat",
                                                 m_timeDimension,
                                                 m_geometryDimension,
                                                 heatQuantity,
                                                 ElementOutput::ChannelOuflowHeatFlux,
                                                 this);

  m_channelOutflowHeatOutput->setCaption("Channel Outflow Heat (J/s)");

  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  m_channelOutflowHeatOutput->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime()- 1.0/1000000.0, m_channelOutflowHeatOutput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), m_channelOutflowHeatOutput);

  m_channelOutflowHeatOutput->addTime(dt1);
  m_channelOutflowHeatOutput->addTime(dt2);

  addOutput(m_channelOutflowHeatOutput);
}

void GWComponent::createChannelOutflowSoluteFluxOutput(int soluteIndex)
{

  QString soluteName = QString::fromStdString(m_modelInstance->solute(soluteIndex));

  ElementOutput *channelOutflowSoluteOutput = new ElementOutput("ChannelOutflow" + soluteName,
                                                                m_timeDimension,
                                                                m_geometryDimension,
                                                                m_soluteConcFluxQuantity,
                                                                ElementOutput::ChannelOutflowSoluteFlux,
                                                                this);

  channelOutflowSoluteOutput->setCaption("Channel Outflow " + soluteName + " Flux (kg/s)");
  channelOutflowSoluteOutput->setSoluteIndex(soluteIndex);

  QList<QSharedPointer<HCGeometry>> geometries;

  for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
  {
    geometries.append(lineString);
  }

  channelOutflowSoluteOutput->addGeometries(geometries);

  SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime()- 1.0/1000000.0, channelOutflowSoluteOutput);
  SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), channelOutflowSoluteOutput);

  channelOutflowSoluteOutput->addTime(dt1);
  channelOutflowSoluteOutput->addTime(dt2);

  m_channelOutflowSoluteOutputs.push_back(channelOutflowSoluteOutput);
  addOutput(channelOutflowSoluteOutput);
}

void GWComponent::createWaterAgeFluxOutput()
{
  if(m_modelInstance->simulateWaterAge())
  {
    Quantity *timeQuantity = Quantity::timeInDays(this);
    ElementOutput *soluteOutput  = new ElementOutput("WaterAgeOutput",
                                                     m_timeDimension,
                                                     m_geometryDimension,
                                                     timeQuantity ,
                                                     ElementOutput::ChannelOutflowSoluteFlux,
                                                     this);
    soluteOutput->setCaption("Element Water Age Flux (days/s)");
    soluteOutput->setSoluteIndex(m_modelInstance->numSolutes());

    QList<QSharedPointer<HCGeometry>> geometries;

    for(const QSharedPointer<HCGeometry> &lineString : m_elementGeometries)
    {
      geometries.append(lineString);
    }

    soluteOutput->addGeometries(geometries);

    SDKTemporal::DateTime *dt1 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime() - 1.0/1000000.0, soluteOutput);
    SDKTemporal::DateTime *dt2 = new SDKTemporal::DateTime(m_modelInstance->currentDateTime(), soluteOutput);

    soluteOutput->addTime(dt1);
    soluteOutput->addTime(dt2);

    addOutput(soluteOutput);
  }
}

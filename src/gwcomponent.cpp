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

    setPrepared(false);
    setInitialized(false);

    setStatus(IModelComponent::Finished , "GWComponent with id " + id() + " has been disposed" , 100);
    setStatus(IModelComponent::Created , "GWComponent with id " + id() + " ran successfully and has been re-created" , 100);
  }
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

}

void GWComponent::createOutputs()
{

}

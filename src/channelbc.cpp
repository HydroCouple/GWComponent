#include "stdafx.h"
#include "channelbc.h"
#include "gwmodel.h"
#include "temporal/timeseries.h"
#include "element.h"
#include "core/datacursor.h"

ChannelBC::ChannelBC(Variable variable,
                     Element *startElement,
                     Element *endElement,
                     GWModel *model):
  QObject(model),
  m_variable(variable),
  m_startElement(startElement),
  m_endElement(endElement),
  m_timeSeries(nullptr),
  m_applyBCFunction(nullptr),
  m_parentModel(model)
{
  m_dataCursor = new DataCursor();
}

ChannelBC::~ChannelBC()
{
  delete m_dataCursor;
}

void ChannelBC::findAssociatedGeometries()
{
  m_profile.clear();
  m_parentModel->findProfile(m_startElement, m_endElement, m_profile);
}

void ChannelBC::prepare()
{
  m_match = false;

  switch (m_variable)
  {
    case WIDTH:
      {
        m_applyBCFunction = &ChannelBC::applyWidthBC;
      }
      break;
    case WSE:
      {
        m_applyBCFunction = &ChannelBC::applyWSEBC;
      }
      break;
    case TEMPERATURE:
      {
        m_applyBCFunction = &ChannelBC::applyTemperatureBC;
      }
      break;
    case SOLUTE:
      {
        m_applyBCFunction = &ChannelBC::applySoluteBC;
      }
      break;
  }

  if(m_timeSeries->numColumns() == (int)m_profile.size())
  {
    m_match = true;
  }
}

void ChannelBC::applyBoundaryConditions(double dateTime)
{
  (this->*m_applyBCFunction)(dateTime);
}

void ChannelBC::clear()
{
  m_profile.clear();
}

ChannelBC::Variable ChannelBC::variable() const
{
  return m_variable;
}

void ChannelBC::setVariable(Variable variable)
{
  m_variable = variable;
}

int ChannelBC::soluteIndex() const
{
  return m_soluteIndex;
}

void ChannelBC::setSoluteIndex(int soluteIndex)
{
  m_soluteIndex = soluteIndex;
}

Element *ChannelBC::startElement() const
{
  return m_startElement;
}

void ChannelBC::setStartElement(Element *element)
{
  m_startElement = element;
}

Element *ChannelBC::endElement() const
{
  return m_endElement;
}

void ChannelBC::setEndElement(Element *element)
{
  m_endElement = element;
}

std::vector<Element*> ChannelBC::profiles() const
{
  return m_profile;
}

QSharedPointer<TimeSeries> ChannelBC::timeSeries() const
{
  return m_timeSeries;
}

void ChannelBC::setTimeSeries(const QSharedPointer<TimeSeries> &timeseries)
{
  m_timeSeries = timeseries;
  m_dataCursor->setMin(0);
  m_dataCursor->setMax(timeseries->numRows() - 1);
}

void ChannelBC::applyWidthBC(double dateTime)
{
  double value = 0;

  if(m_match)
  {
    for(size_t j = 0; j < m_profile.size(); j++)
    {
      if(m_timeSeries->interpolate(dateTime, j, m_dataCursor, value))
      {
        m_profile[j]->channelWidth = value;
      }
    }
  }
  else
  {
    if(m_timeSeries->interpolate(dateTime, 0, m_dataCursor, value))
    {
      for(size_t j = 0; j < m_profile.size(); j++)
      {
        m_profile[j]->channelWidth = value;
      }
    }
  }
}

void ChannelBC::applyWSEBC(double dateTime)
{
  double value = 0;

  if(m_match)
  {
    for(size_t j = 0; j < m_profile.size(); j++)
    {
      if(m_timeSeries->interpolate(dateTime, j, m_dataCursor, value))
      {
        m_profile[j]->channelWSE = value;
      }
    }
  }
  else
  {
    if(m_timeSeries->interpolate(dateTime, 0, m_dataCursor, value))
    {
      for(size_t j = 0; j < m_profile.size(); j++)
      {
        m_profile[j]->channelWSE = value;
      }
    }
  }
}

void ChannelBC::applyTemperatureBC(double dateTime)
{
  double value = 0;

  if(m_match)
  {
    for(size_t j = 0; j < m_profile.size(); j++)
    {
      if(m_timeSeries->interpolate(dateTime, j, m_dataCursor, value))
      {
        m_profile[j]->channelTemperature = value;
      }
    }
  }
  else
  {
    if(m_timeSeries->interpolate(dateTime, 0, m_dataCursor, value))
    {
      for(size_t j = 0; j < m_profile.size(); j++)
      {
        m_profile[j]->channelTemperature = value;
      }
    }
  }
}

void ChannelBC::applySoluteBC(double dateTime)
{
  double value = 0;

  if(m_match)
  {
    for(size_t j = 0; j < m_profile.size(); j++)
    {
      if(m_timeSeries->interpolate(dateTime, j, m_dataCursor, value))
      {
        m_profile[j]->channelSoluteConcs[m_soluteIndex] = value;
      }
    }
  }
  else
  {
    if(m_timeSeries->interpolate(dateTime, 0, m_dataCursor, value))
    {
      for(size_t j = 0; j < m_profile.size(); j++)
      {
        m_profile[j]->channelSoluteConcs[m_soluteIndex] = value;
      }
    }
  }
}

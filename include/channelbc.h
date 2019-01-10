#ifndef CHANNELBC_H
#define CHANNELBC_H


#include <QObject>
#include <QSharedPointer>

#include "iboundarycondition.h"

struct Element;
class TimeSeries;
class GWModel;
class DataCursor;
class ChannelBC;

typedef void (ChannelBC::*ApplyBC)(double dateTime);

class GWCOMPONENT_EXPORT ChannelBC : public QObject,
    public virtual IBoundaryCondition
{

    Q_OBJECT

  public:

    enum Variable
    {
      WIDTH = 1,
      WSE = 2,
      TEMPERATURE = 3,
      SOLUTE = 4
    };

  public:

    ChannelBC(Variable variable,
              Element *startElement,
              Element *endElement,
              GWModel *model);

    virtual ~ChannelBC();

    void  findAssociatedGeometries() override final;

    void prepare() override final;

    void applyBoundaryConditions(double dateTime) override final;

    void clear() override final;

    Variable variable() const;

    void setVariable(Variable variable);

    int soluteIndex() const;

    void setSoluteIndex(int soluteIndex);

    Element *startElement() const;

    void setStartElement(Element *element);

    Element *endElement() const;

    void setEndElement(Element *element);

    std::vector<Element*> profiles() const;

    QSharedPointer<TimeSeries> timeSeries() const;

    void setTimeSeries(const QSharedPointer<TimeSeries> &timeseries);

  private:

    void applyWidthBC(double dateTime);

    void applyWSEBC(double dateTime);

    void applyTemperatureBC(double dateTime);

    void applySoluteBC(double dateTime);


  private:

    Variable m_variable;
    Element *m_startElement, *m_endElement;
    int m_soluteIndex;
    bool m_match;
    std::vector<Element*> m_profile;
    QSharedPointer<TimeSeries> m_timeSeries;
    ApplyBC m_applyBCFunction;
    DataCursor *m_dataCursor;
    GWModel *m_parentModel;
};

#endif // CHANNELBC_H

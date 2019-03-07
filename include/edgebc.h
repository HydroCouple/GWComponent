#ifndef EDGEBC_H
#define EDGEBC_H

#include <QObject>
#include <vector>
#include <QSharedPointer>

#include "iboundarycondition.h"

struct Element;
struct ElementCell;
class TimeSeries;
class GWModel;
class EdgeBC;
class TimeSeries;
class DataCursor;

typedef void (EdgeBC::*ApplyBoundaryVarBC)(double dateTime);

class GWCOMPONENT_EXPORT EdgeBC : public QObject,
    public virtual IBoundaryCondition
{

    Q_OBJECT

  public:
    enum Variable
    {
      HYD_HEAD = 1,
      GRAD_HYD_HEAD = 2,
      TEMPERATURE = 3,
      GRAD_TEMPERATURE = 4,
      SOLUTE = 5,
      GRAD_SOLUTE = 6
    };

  public:
    enum Edge
    {
      RIGHT = 0,
      DOWN = 1,
      LEFT = 2,
      UP = 3
    };

  public:

    EdgeBC(Variable variable,
           Edge edge,
           Element *startElement,
           Element *endElement,
           int startElementCell,
           int endElementCell,
           GWModel *model);

    virtual ~EdgeBC();

    void  findAssociatedGeometries() override final;

    void prepare() override final;

    void applyBoundaryConditions(double dateTime) override final;

    void clear() override final;

    Variable variable() const;

    void setVariable(Variable variable);

    Edge edge() const;

    void setEdge(Edge edge);

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

    void applyHydHeadBC(double dateTime);

    void applyGradHydHeadBC(double dateTime);

    void applyTempBC(double dateTime);

    void applyGradTempBC(double dateTime);

    void applySoluteBC(double dateTime);

    void applyGradSoluteBC(double dateTime);

  private:

    Variable m_variable;
    Edge m_edge, m_nEdge;
    Element *m_startElement, *m_endElement;
    int m_cellLength;
    int m_startElementCell, m_endElementCell;
    int m_soluteIndex;
    bool m_match;
    std::vector<Element*> m_profile;
    std::vector<ElementCell*> m_neighbourCells;
    QSharedPointer<TimeSeries> m_timeSeries;
    ApplyBoundaryVarBC m_applyBCFunction;
    DataCursor *m_dataCursor;
    GWModel *m_parentModel;
};

#endif // EDGEBC_H

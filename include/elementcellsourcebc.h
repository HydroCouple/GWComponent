#ifndef ELEMENTCELLSOURCE_H
#define ELEMENTCELLSOURCE_H

#include <QObject>
#include <QSharedPointer>

#include "iboundarycondition.h"

struct Element;
struct ElementCell;
class TimeSeries;
class GWModel;
class DataCursor;
class ElementCellSourceBC;

typedef void (ElementCellSourceBC::*ApplySourceBC)(double dateTime);
typedef double (ElementCellSourceBC::*GetGeometryMultiplier)(ElementCell *elementCell);

class GWCOMPONENT_EXPORT ElementCellSourceBC : public QObject,
    public virtual IBoundaryCondition
{

    Q_OBJECT

  public:

    enum Variable
    {
      FLOW = 1,
      HEAT = 2,
      SOLUTE = 3
    };

  public:

    enum GeometryMultiplier
    {
      None = 0,
      Area = 1,
      Volume = 2,
      Length = 3,
      Width = 4,
      SaturatedThickness = 5
    };

  public:

    ElementCellSourceBC(Variable variable,
                        GeometryMultiplier geometryMultiplier,
                        Element *startElement,
                        Element *endElement,
                        int startElementCell,
                        int endElementCell,
                        GWModel *model);

    virtual ~ElementCellSourceBC();

    void  findAssociatedGeometries() override final;

    void prepare() override final;

    void applyBoundaryConditions(double dateTime) override final;

    void clear() override final;

    Variable variable() const;

    void setVariable(Variable variable);

    GeometryMultiplier geometryMultiplier() const;

    void setGeometryMultiplier(GeometryMultiplier multiplier);

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

    void applyFlowBC(double dateTime);

    void applyHeatBC(double dateTime);

    void applySoluteBC(double dateTime);

    double getZero(ElementCell *elementCell);

    double getArea(ElementCell *elementCell);

    double getVolume(ElementCell *elementCell);

    double getLength(ElementCell *elementCell);

    double getWidth(ElementCell *elementCell);

    double getSaturatedThickness(ElementCell *elementCell);

  private:

    Variable m_variable;
    GeometryMultiplier m_geometryMultiplier;
    Element *m_startElement, *m_endElement;
    int m_cellLength;
    int m_startElementCell, m_endElementCell;
    int m_soluteIndex;
    bool m_match;
    std::vector<Element*> m_profile;
    QSharedPointer<TimeSeries> m_timeSeries;
    GetGeometryMultiplier m_getGeometryMultiplier;
    ApplySourceBC m_applyBCFunction;
    DataCursor *m_dataCursor;
    GWModel *m_parentModel;
};

#endif // ELEMENTCELLSOURCE_H

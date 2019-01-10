#ifndef ELEMENTOUTPUT_H
#define ELEMENTOUTPUT_H


#include "gwcomponent_global.h"
#include "spatiotemporal/timegeometryoutput.h"

class GWComponent;

class GWCOMPONENT_EXPORT ElementOutput: public TimeGeometryOutputDouble
{

    Q_OBJECT

  public:

    enum VariableType
    {
      ChannelOutflow,
      ChannelOuflowHeatFlux,
      ChannelOutflowSoluteFlux,
    };


  public:

    ElementOutput(const QString &id,
                  Dimension *timeDimension,
                  Dimension *geometryDimension,
                  ValueDefinition *valueDefinition,
                  VariableType variableType,
                  GWComponent *modelComponent);


    virtual ~ElementOutput();

    void updateValues(HydroCouple::IInput *querySpecifier) override;

    void updateValues() override;

    VariableType variableType() const;

    void setVariableType(VariableType variableType);

    int soluteIndex() const;

    void setSoluteIndex(int soluteIndex);

  private:

    GWComponent *m_component;
    VariableType m_variableType;
    int m_soluteIndex;
};



#endif // ELEMENTOUTPUT_H

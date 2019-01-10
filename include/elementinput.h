#ifndef ELEMENTINPUT_H
#define ELEMENTINPUT_H

#include "gwcomponent_global.h"
#include "spatiotemporal/timegeometrymultiinput.h"

#include <unordered_map>

class GWComponent;

class GWCOMPONENT_EXPORT ElementInput : public  TimeGeometryMultiInputDouble
{
  public:

    enum VariableType
    {
      ChannelWidth,
      ChannelWSE,
      ChannelTemperature,
      ChannelSolute,
    };

    ElementInput(const QString &id,
                 Dimension *timeDimension,
                 Dimension *geometryDimension,
                 ValueDefinition *valueDefinition,
                 VariableType variableType,
                 GWComponent *modelComponent);

    virtual ~ElementInput();

    bool addProvider(HydroCouple::IOutput *provider) override;

    bool removeProvider(HydroCouple::IOutput *provider) override;

    bool canConsume(HydroCouple::IOutput *provider, QString &message) const override;

    void retrieveValuesFromProvider() override;

    void applyData() override;

    VariableType variableType() const;

    void setVariableType(VariableType variableType);

    int soluteIndex() const;

    void setSoluteIndex(int soluteIndex);

  private:

    GWComponent *m_component;
    VariableType m_variableType;
    std::unordered_map<HydroCouple::IOutput*, std::unordered_map<int,int>> m_geometryMapping;
    int m_soluteIndex;
};

#endif // ELEMENTINPUT_H

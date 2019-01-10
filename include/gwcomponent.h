/*!
*  \file    GWComponent.h
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


#ifndef GWCOMPONENT_H
#define GWCOMPONENT_H

#include "gwcomponentinfo.h"
#include "temporal/abstracttimemodelcomponent.h"

class Dimension;
class GWModel;
class Unit;
class HCGeometry;
class ElementInput;
class ElementOutput;
class Quantity;

class GWCOMPONENT_EXPORT GWComponent: public AbstractTimeModelComponent,
    public virtual HydroCouple::ICloneableModelComponent
{
    Q_OBJECT

    Q_INTERFACES(HydroCouple::ICloneableModelComponent)

  public:

    /*!
     * \brief GWComponent constructor
     * \param id Unique identifier for this component instance.
     * \param modelComponentInfo the parent ModelComponentInfo that generated this component instance.
     */
    GWComponent(const QString &id, GWComponentInfo* modelComponentInfo = nullptr);

    /*!
     * \brief ~GWComponent destructor
     */
    virtual ~GWComponent();

    /*!
     * \brief validate validates this component model instance
     * \return Returns a list of error messages.
     */
    QList<QString> validate() override;

    /*!
     * \brief prepare Prepares the model component instance.
     */
    void prepare() override;

    /*!
     * \brief update
     * \param requiredOutputs
     */
    void update(const QList<HydroCouple::IOutput*> &requiredOutputs = QList<HydroCouple::IOutput*>()) override;

    /*!
     * \brief finish
     */
    void finish() override;

    /*!
     * \brief modelInstance
     * \return
     */
    GWModel *modelInstance() const;

    /*!
     * \brief parent
     * \return
     */
    HydroCouple::ICloneableModelComponent* parent() const override;

    /*!
     * \brief clone
     * \return
     */
    HydroCouple::ICloneableModelComponent* clone() override;

    /*!
     * \brief clones
     * \return
     */
    QList<HydroCouple::ICloneableModelComponent*> clones() const override;

  protected:

    bool removeClone(GWComponent *component);

    /*!
     * \brief intializeFailureCleanUp
     */
    void initializeFailureCleanUp() override;

  private:

    /*!
     * \brief createArguments
     */
    void createArguments() override;

    /*!
     * \brief createInputFileArguments
     */
    void createInputFileArguments();

    /*!
     * \brief initializeArguments
     * \param message
     * \return
     */
    bool initializeArguments(QString &message) override;

    /*!
     * \brief initializeInputFilesArguments
     * \param message
     * \return
     */
    bool initializeInputFilesArguments(QString &message);

    /*!
     * \brief createGeometriesMap
     */
    void createGeometries();

    /*!
     * \brief createInputs
     */
    void createInputs() override;

    /*!
     * \brief createChannelWSEInput
     */
    void createChannelWSEInput();

    /*!
     * \brief createChannelWidthInput
     */
    void createChannelWidthInput();

    /*!
     * \brief createChannelTemperatureInput
     */
    void createChannelTemperatureInput();

    /*!
     * \brief createChannelSoluteInput
     * \param soluteIndex
     */
    void createChannelSoluteInput(int soluteIndex);

    /*!
     * \brief createOutputs
     */
    void createOutputs() override;

    /*!
     * \brief createChannelOutflowOutput
     */
    void createChannelOutflowOutput();

    /*!
     * \brief createChannelOutflowHeatFluxOutput
     */
    void createChannelOutflowHeatFluxOutput();

    /*!
     * \brief createChannelOutflowSoluteFluxOutput
     * \param soluteIndex
     */
    void createChannelOutflowSoluteFluxOutput(int soluteIndex);

  private:

    GWComponentInfo *m_GWComponentInfo;
    GWModel *m_modelInstance;

    GWComponent *m_parent;
    QList<HydroCouple::ICloneableModelComponent*> m_clones;

    IdBasedArgumentString *m_inputFilesArgument;

    Dimension *m_timeDimension,
    *m_geometryDimension;

    Unit *m_radiationFluxUnit,
         *m_heatFluxUnit,
         *m_temperatureUnit,
         *m_soluteUnit,
         *m_soluteFluxUnit;

    Quantity *m_soluteConcQuantity, *m_soluteConcFluxQuantity;

    ElementInput *m_channelWSEInput, *m_channelWidthInput, *m_channelTemperatureInput;
    std::vector<ElementInput*> m_channelSoluteInputs;

    ElementOutput *m_channelOutflowOutput, *m_channelOutflowHeatOutput;
    std::vector<ElementOutput*> m_channelOutflowSoluteOutputs;

    std::vector<QSharedPointer<HCGeometry>> m_elementGeometries;
    std::vector<QSharedPointer<HCGeometry>> m_elementJunctionGeometries;

};

#endif // GWCOMPONENT_H

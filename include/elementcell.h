#ifndef ELEMENTCELL_H
#define ELEMENTCELL_H

#include "variable.h"
#include "gwcomponent_global.h"

#include <string>
#include <vector>

struct Element;
struct ElementCell;
struct ElementJunction;
class GWModel;

/*!
 *\brief Function pointer to calculate temperature advection to eliminate costly if else function calls
 */
typedef double (ElementCell::*ComputeTempAdv)(double dt, double T[]);

/*!
 *\brief Function pointer to calculate solute advection to eliminate costly if else function calls
 */
typedef double (ElementCell::*ComputeSoluteAdv)(double dt, double S[], int soluteIndex);

/*!
 *
 */
typedef double (ElementCell::*ComputeEdgeDeriv)(int edgeIndex, double dt, double H[]);

/*!
 *
 */
typedef double (ElementCell::*ComputeEdgeVariableDeriv)(int edgeIndex, double dt, double H[], int soluteIndex);

struct ElementCell
{

    friend struct Element;

    /*!
     * \brief ElementCell
     * \param cindex
     * \param cparent
     */
    ElementCell(int cindex,
                Element *cparent);

    ~ElementCell();

    /*!
     * \brief index
     */
    int index ;

    /*!
     * \brief elementCellIndex
     */
    int elementCellIndex;

    /*!
     * \brief isBedCell
     */
    bool isBedCell = false;

    /*!
     * \brief topBedCell
     */
    ElementCell *topBedCell = nullptr;

    /*!
     * \brief bedCells
     */
    ElementCell **bedCells = nullptr;

    /*!
     * \brief zIndex
     */
    int zIndex = 0;

    /*!
     * \brief hydHead
     */
    Variable hydHead;

    /*!
     * \brief prevHydHead
     */
    Variable prevHydHead;

    /*!
     * \brief edgeHydHead
     */
    Variable *edgeHydHead;

    /*!
     * \brief gradHydHeadX
     */
    Variable *edgeGradHydHead;

    /*!
     * \brief externalInflow
     */
    double externalInflow;

    /*!
     * \brief channelInflow
     */
    double channelInflow;

    /*!
     * \brief channelInflowFlux
     */
    double channelInflowFlux;

    /*!
     * \brief totalMassBalance
     */
    double totalMassBalance;

    /*!
     * \brief temperature (°C)
     */
    Variable temperature;

    /*!
     * \brief prevTemperature (°C)
     */
    Variable prevTemperature;

    /*!
     * \brief edgeTemperatures
     */
    Variable *edgeTemperatures;

    /*!
     * \brief gradEdgeTemperatures
     */
    Variable *edgeGradTemperatures;

    /*!
     * \brief edgeHeatFluxes
     */
    double *edgeHeatFluxes;

    /*!
     * \brief externalHeatFluxes
     */
    double externalHeatFluxes;

    /*!
     * \brief channelHeatFlux
     */
    double channelHeatFlux;

    /*!
     * \brief channelHeatRate
     */
    double channelHeatRate;

    /*!
     * \brief totalHeatBalance
     */
    double totalHeatBalance;

    /*!
     * \brief totalExternalHeatBalance
     */
    double totalExternalHeatBalance;

    /*!
     * \brief soluteConcs (kg/m^3)
     */
    Variable *soluteConcs;

    /*!
     * \brief prevSoluteConcs (kg/m^3)
     */
    Variable *prevSoluteConcs;

    /*!
     * \brief edgeSoluteConcs
     */
    Variable **edgeSoluteConcs;

    /*!
     * \brief gradEdgeSoluteConcs
     */
    Variable **edgeGradSoluteConcs;

    /*!
     * \brief edgeSoluteConcFluxes
     */
    double **edgeSoluteConcFluxes;

    /*!
     * \brief externalSoluteFluxes
     */
    double *externalSoluteFluxes;

    /*!
     * \brief channelSoluteFlux
     */
    double *channelSoluteFlux;

    /*!
     * \brief channelSoluteRate
     */
    double *channelSoluteRate;

    /*!
     * \brief totalSoluteMassBalance
     */
    double *totalSoluteMassBalance;

    /*!
     * \brief hydConX
     */
    double *hydCon;

    /*!
     * \brief sedDensity
     */
    double sedDensity;

    /*!
     * \brief setCp
     */
    double sedCp;

    /*!
     * \brief porosity
     */
    double porosity;

    /*!
     * \brief topElev
     */
    double topElev;

    /*!
     * \brief bottomElev
     */
    double bottomElev;

    /*!
     * \brief bedRockElev
     */
    double bedRockElev;

    /*!
     * \brief bedRockElevs
     */
    double *bedRockElevs;

    /*!
     * \brief width
     */
    double width;

    /*!
     * \brief wettedWidth
     */
    double wettedWidth;

    /*!
     * \brief specificYield
     */
    double specificYield = 0.2;

    /*!
     * \brief specificStorage
     */
    double specificStorage = 6.2E-3;

    /*!
     * \brief dispersivityX (m)
     */
    double *dispersivity;

    /*!
     * \brief sedThermalConductivity
     */
    double sedThermalConductivity;

    /*!
     * \brief edgeFlows
     */
    double *edgeFlows;


    double **edgeFlowsComp;

    /*!
     * \brief depth
     */
    double depth;

    /*!
     * \brief satDepth
     */
    double satDepth;

    /*!
     * \brief edgeHydCons
     */
    double *edgeHydCons;

    /*!
     * \brief effectiveKe
     */
    double *edgeEffectiveKe;

    /*!
     * \brief edgePorosity
     */
    double *edgePorosity;

    /*!
     * \brief edgeHeatDispersionCoeff
     */
    double *edgeMechDispersionCoeff;

    /*!
     * \brief edgeHeatPecletNumbers
     */
    double *edgeHeatPecletNumbers;

    /*!
     * \brief edgeSolutePecletNumbers
     */
    double *edgeSolutePecletNumbers;

    /*!
     * \brief centerY
     */
    double centerY;

    /*!
     * \brief volume
     */
    double volume;

    /*!
     * \brief prevVolume
     */
    double prev_volume;

    /*!
     * \brief dvolume_dt
     */
    double dvolume_dt;


    double dsolute_dt;

    /*!
     * \brief start
     */
    bool start;

    /*!
     * \brief parentElement
     */
    Element *parentElement;

    /*!
     * \brief neighbour
     * \param edge
     * \param top
     * \return
     */
    ElementCell *neighbour(int edge, bool top = false);

    /*!
     * \brief flowLengthP
     * \param edge
     * \return
     */
    double flowLengthP(int edge);

    /*!
     * \brief flowLengthN
     * \param edge
     * \param neighbour
     * \return
     */
    double flowLengthN(int edge, ElementCell *n = nullptr);

    /*!
     * \brief flowArea
     * \param edge
     * \return
     */
    double flowArea(int edge, int layer);

    /*!
     * \brief layer
     * \param level
     * \return
     */
    ElementCell *layer(int level);

    /*!
     * \brief initializeBedCells
     */
    void initializeBedCells();

    /*!
     * \brief deleteBedCells
     */
    void deleteBedCells();

    /*!
     * \brief isTopCell
     * \return
     */
    bool isTopCell();

    /*!
     * \brief initialize
     */
    void initialize();

    /*!
     * \brief initializeSolutes
     */
    void initializeSolutes();

    /*!
     * \brief computeDHydHeadDt
     * \param dt
     * \param H
     * \return
     */
    double computeDHydHeadDt(double dt, double H[]);

    /*!
     * \brief computeDTDt
     * \param dt
     * \param T
     * \return
     */
    double computeDTDt(double dt, double T[]);

    /*!
     * \brief computeDTDtDispersion
     * \param dt
     * \param T
     * \return
     */
    double computeDTDtDispersion(double dt, double T[]);

    /*!
     * \brief computeDTDtUpwind
     * \param dt
     * \param T
     * \return
     */
    double computeDTDtUpwind(double dt, double T[]);

    /*!
     * \brief computeDTDtCentral
     * \param dt
     * \param T
     * \return
     */
    double computeDTDtCentral(double dt, double T[]);

    /*!
     * \brief computeDTDtHybrid
     * \param dt
     * \param T
     * \return
     */
    double computeDTDtHybrid(double dt, double T[]);

    /*!
     * \brief computeDTDtUltimate
     * \param dt
     * \param T
     * \return
     */
    double computeDTDtTVD(double dt, double T[]);

    /*!
     * \brief computeDSoluteDt
     * \param dt
     * \param S
     * \return
     */
    double computeDSoluteDt(double dt, double S[], int soluteIndex);

    /*!
     * \brief computeDSoluteDtDispersion
     * \param dt
     * \param S
     * \param soluteIndex
     * \return
     */
    double computeDSoluteDtDispersion(double dt, double S[], int soluteIndex);

    /*!
     * \brief computeDSoluteDtUpwind
     * \param dt
     * \param S
     * \param soluteIndex
     * \return
     */
    double computeDSoluteDtUpwind(double dt, double S[], int soluteIndex);

    /*!
     * \brief ElementCell::computeDSoluteDtCentral
     * \param dt
     * \param S
     * \param soluteIndex
     * \return
     */
    double computeDSoluteDtCentral(double dt, double S[], int soluteIndex);

    /*!
     * \brief computeDSoluteDtHybrid
     * \param dt
     * \param S
     * \param soluteIndex
     * \return
     */
    double computeDSoluteDtHybrid(double dt, double S[], int soluteIndex);

    /*!
     * \brief computeDSoluteDtTVD
     * \param dt
     * \param S
     * \param soluteIndex
     * \return
     */
    double computeDSoluteDtTVD(double dt, double S[], int soluteIndex);

    /*!
     * \brief computeCourantNumber
     * \return
     */
    double computeCourantFactor() const;

    /*!
     * \brief computeDiffusionFactor
     * \return
     */
    double computeDiffusionFactor() ;

    /*!
     * \brief computeDerivedHydraulics
     */
    void calculatePreComputedHydHeadVariables();

    /*!
     * \brief calculatePreComputedTempVariables
     */
    void calculatePreComputedTempVariables();

    /*!
     * \brief calculatePreComputedSoluteVariables
     */
    void calculatePreComputedSoluteVariables();

    /*!
     * \brief computePecletNumbers
     */
    void computeHeatPecletNumbers();

    /*!
     * \brief computeSolutePecletNumbers
     */
    void computeSolutePecletNumbers();

    /*!
     * \brief computeMassBalance
     * \param timeStep
     */
    void computeMassBalance(double timeStep);

    /*!
     * \brief computeHeatBalance
     */
    void computeHeatBalance(double timeStep);

    /*!
     * \brief computeSoluteBalance
     * \param soluteIndex
     */
    void computeSoluteBalance(double timeStep, int soluteIndex);

    /*!
     * \brief computeEdgeDepths
     */
    void computeDepth();

    /*!
     * \brief computeEdgeDepth
     * \param edge
     * \param layer
     */
    double computeEdgeDepth(int edge, int layerIndex);

    /*!
     * \brief computeVolumeDerivative
     */
    void computeVolumeDerivative();

  private:

    /*!
     * \brief deleteSoluteVariables
     */
    void deleteSoluteVariables();

    /*!
     * \brief computeEdgeHydCons
     */
    void computeEdgeHydCons();

    /*!
     * \brief computeEdgeDispersionCoefficients
     */
    void computeEdgeDispersionCoefficients();

    /*!
     * \brief computeChannelMassFlux
     * \param dt
     * \param T
     */
    void computeChannelMassFlux();

    /*!
     * \brief computeChannelHeatFlux
     * \param dt
     * \param T
     */
    void computeChannelHeatFlux();

    /*!
         * \brief computeChannelSoluteFlux
         * \param dt
         * \param S
         * \param soluteIndex
         */
    void computeChannelSoluteFlux(int soluteIndex);

    /*!
     * \brief computeHydHeadEdgeHeadBC
     * \param edgeIndex
     * \param dt
     * \param H
     * \return
     */
    double computeEdgeHydHeadBC(int edgeIndex, double dt, double H[]);

    /*!
     * \brief computeHydHeadEdgeGradHeadBC
     * \param edgeIndex
     * \param dt
     * \param H
     * \return
     */
    double computeEdgeGradHydHeadBC(int edgeIndex, double dt, double H[]);

    /*!
     * \brief computeHydHeadNeighborBC
     * \param edgeIndex
     * \param dt
     * \param H
     * \return
     */
    double computeNeighborHydHeadBC(int edgeIndex, double dt, double H[]);

    /*!
     * \brief computeNeighborLayersHydHeadBC
     * \param edgeIndex
     * \param dt
     * \param H
     * \return
     */
    double computeNeighborLayersHydHeadBC(int edgeIndex, double dt, double H[]);

    /*!
     * \brief computeHydHeadZeroBC
     * \param edgeIndex
     * \param dt
     * \param H
     * \return
     */
    double computeZeroHydHeadBC(int edgeIndex, double dt, double H[]);

    /*!
     * \brief computeDTDtEdgeBC
     * \param edgeIndex
     * \param dt
     * \param T
     * \return
     */
    double computeEdgeTempBC(int edgeIndex, double dt, double T[]);

    /*!
     * \brief computeEdgeGradTempBC
     * \param edgeIndex
     * \param dt
     * \param T
     * \return
     */
    double computeEdgeGradTempBC(int edgeIndex, double dt, double T[]);

    /*!
     * \brief computeDTDtNeighborBC
     * \param edgeIndex
     * \param dt
     * \param T
     * \return
     */
    double computeNeighborTempBC(int edgeIndex, double dt, double T[]);

    /*!
     * \brief computeNeighborLayersTempBC
     * \param edgeIndex
     * \param dt
     * \param T
     * \return
     */
    double computeNeighborLayersTempBC(int edgeIndex, double dt, double T[]);

    /*!
     * \brief computeDTDtZeroBC
     * \param edgeIndex
     * \param dt
     * \param T
     * \return
     */
    double computeZeroTempBC(int edgeIndex, double dt, double T[]);

    /*!
     * \brief computeDSoluteDtEdgeBC
     * \param edgeIndex
     * \param dt
     * \param T
     * \param soluteIndex
     * \return
     */
    double computeEdgeSoluteBC(int edgeIndex, double dt, double S[], int soluteIndex);

    /*!
     * \brief computeEdgeGradSoluteBC
     * \param edgeIndex
     * \param dt
     * \param S
     * \param soluteIndex
     * \return
     */
    double computeEdgeGradSoluteBC(int edgeIndex, double dt, double S[], int soluteIndex);

    /*!
     * \brief computeNeighborSoluteBC
     * \param edgeIndex
     * \param dt
     * \param T
     * \param soluteIndex
     * \return
     */
    double computeNeighborSoluteBC(int edgeIndex, double dt, double S[], int soluteIndex);

    /*!
     * \brief computeNeighborLayersSoluteBC
     * \param edgeIndex
     * \param dt
     * \param S
     * \param soluteIndex
     * \return
     */
    double computeNeighborLayersSoluteBC(int edgeIndex, double dt, double S[], int soluteIndex);

    /*!
     * \brief computeZeroSoluteBC
     * \param edgeIndex
     * \param dt
     * \param T
     * \return
     */
    double computeZeroSoluteBC(int edgeIndex, double dt, double S[], int soluteIndex);

    /*!
     * \brief computeTempAdv Pointer to function to compute temperature advection.
     */
    ComputeTempAdv computeTempAdv;

    /*!
     * \brief computeSoluteAdv Pointer to function to compute solute advection.
     */
    ComputeSoluteAdv computeSoluteAdv;

    /*!
     * \brief computeEdgeHeadDeriv
     */
    ComputeEdgeDeriv *computeEdgeHeadDerivs;

    /*!
     * \brief computeEdgeHeadDerivs
     */
    ComputeEdgeDeriv *computeEdgeTempDerivs;


    ComputeEdgeVariableDeriv **computeEdgeSoluteDerivs;

  private:

    double rhom_Cm;

    const static int nEdgeIndex[];

    const static double ef1[];

    const static double ef2[];

    double *retardationFactor;

};


#endif // ELEMENTCELL_H


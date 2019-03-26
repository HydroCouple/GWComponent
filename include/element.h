#ifndef ELEMENT_H
#define ELEMENT_H

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


typedef double (ElementCell::*ComputeEdgeDeriv)(int edgeIndex, double dt, double H[]);

typedef double (ElementCell::*ComputeEdgeVariableDeriv)(int edgeIndex, double dt, double H[], int soluteIndex);


struct ElementCell
{

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
     * \brief externalHeatFluxes
     */
    double externalHeatFluxes;

    /*!
     * \brief channelHeatFlux
     */
    double channelHeatFlux;

    /*!
     * \brief totalHeatBalance
     */
    double totalHeatBalance;

    /*!
     * \brief totalExternalHeatBalance
     */
    double totalExternalHeatBalance;

    /*!
     * \brief numSolutes
     */
    int numSolutes = 0;

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
     * \brief externalSoluteFluxes
     */
    double *externalSoluteFluxes;

    /*!
     * \brief channelSoluteFlux
     */
    double *channelSoluteFlux;

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
     * \brief bedRockElev
     */
    double bedRockElev;

    /*!
     * \brief edgeBedRockElev
     */
    double *edgeBedRockElevs;

    /*!
     * \brief width
     */
    double width;

    /*!
     * \brief wettedWidth
     */
    double wettedWidth;

    /*!
     * \brief specificStorage
     */
    double specificStorage;

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

    /*!
     * \brief edgeBedRockElevs
     */
    //    double *edgeBedRockElevs;

    /*!
     * \brief edgeDepths
     */
    double *edgeDepths;

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

    /*!
     * \brief start
     */
    bool start;

    /*!
     * \brief parentElement
     */
    Element *parentElement;

    /*!
     * \brief neighbors
     */
    ElementCell **neighbors;


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
    double computeDiffusionFactor() const;

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
    void computeEdgeDepths();

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
     * \brief computeEdgeDepths
     * \param H
     */
    void computeEdgeDepths(double H[]);

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

    //    /*!
    //     * \brief computeChannelMassFlux
    //     * \param dt
    //     * \param T
    //     */
    //    void computeChannelMassFlux(double dt, double H[]);

    //    /*!
    //     * \brief computeChannelHeatFlux
    //     * \param dt
    //     * \param T
    //     */
    //    void computeChannelHeatFlux(double dt, double T[]);

    //    /*!
    //     * \brief computeChannelSoluteFlux
    //     * \param dt
    //     * \param S
    //     * \param soluteIndex
    //     */
    //    void computeChannelSoluteFlux(double dt, double S[], int soluteIndex);

    /*!
     * \brief computeHydHeadEdgeHeadBC
     * \param edgeIndex
     * \param dt
     * \param H
     * \return
     */
    double computeEdgeHydHeadHeadBC(int edgeIndex, double dt, double H[]);

    /*!
     * \brief computeHydHeadEdgeGradHeadBC
     * \param edgeIndex
     * \param dt
     * \param H
     * \return
     */
    double computeEdgeGradHydHeadHeadBC(int edgeIndex, double dt, double H[]);

    /*!
     * \brief computeHydHeadNeighborBC
     * \param edgeIndex
     * \param dt
     * \param H
     * \return
     */
    double computeNeighborHydHeadBC(int edgeIndex, double dt, double H[]);

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

    double *flowLengthP;

    double *flowLengthN;

    double *flowWidth;

    const static int nEdgeIndex[];

    const static double ef1[];

    const static double ef2[];

    double *retardationFactor;

};

/*!
 * \brief This struct represents the channel control volume
 */
struct GWCOMPONENT_EXPORT Element
{
    /*!
    * \brief Element - Creates an instance of the control volume element used to represent a computational
    * element in a reach.
    * \param numSolutes - Number of solutes that are going to be transported in addition to temperature.
    * \param from - The upstream junction of this element.
    * \param to - The downstream junction of this element.
    * \param project
    */
   Element(const std::string &id, ElementJunction *upstream, ElementJunction *downstream,  GWModel *model);

   /*!
    * \brief ~Element - Destructor for this class.
    */
   ~Element();

   /*!
    * \brief index unique identifier for element
    */
   int index;

   /*!
    * \brief id
    */
   std::string id;

   /*!
    * \brief x
    */
   double x;

   /*!
    * \brief y
    */
   double y;

   /*!
    * \brief z
    */
   double z;

   /*!
    * \brief length
    */
   double length;

   /*!
    * \brief upstreamElement
    */
   Element *upstreamElement;

   /*!
    * \brief downstreamElement
    */
   Element *downstreamElement;

   /*!
    * \brief fromJunction
    */
   ElementJunction *upstreamJunction;

   /*!
    * \brief toJunction
    */
   ElementJunction *downstreamJunction;

   /*!
    * \brief channelWidth
    */
   double channelWidth;

   /*!
    * \brief channelWSE
    */
   double channelWSE;

   /*!
    * \brief channelTemperature
    */
   double channelTemperature;

   /*!
    * \brief channelConductance
    */
   double channelBedHydCond;

   /*!
    * \brief hydConZ
    */
   double hydConZ;

   /*!
    * \brief channelBedThickness
    */
   double channelBedThickness;

   /*!
    * \brief channelSoluteConcs
    */
   double *channelSoluteConcs;

   /*!
    * \brief channelInflow
    */
   double channelInflow;

   /*!
    * \brief channelInflowFlux
    */
   double channelInflowFlux;

   /*!
    * \brief channelHeat
    */
   double channelHeatRate;

   /*!
    * \brief channelHeatFlux
    */
   double channelHeatFlux;

   /*!
    * \brief channelSoluteRate
    */
   double *channelSoluteRate;

   /*!
    * \brief channelSoluteFlux
    */
   double *channelSoluteFlux;

   /*!
    * \brief elementCells
    */
   std::vector<ElementCell*> elementCells;

   /*!
    * \brief model
    */
   GWModel *model;

   /*!
    * \brief initialize
    */
   void initialize();

   /*!
    * \brief initializeElementCells
    */
   void initializeElementCells();


  private:

   /*!
    * \brief initializeSolutes
    */
   void initializeSolutes();

   /*!
    * \brief setUpstreamElement
    */
   void setUpstreamElement();

   /*!
    * \brief setDownStreamElement
    */
   void setDownStreamElement();

};

#endif // ELEMENT_H

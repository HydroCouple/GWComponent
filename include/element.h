#ifndef ELEMENT_H
#define ELEMENT_H

#include "variable.h"
#include "gwcomponent_global.h"

#include <string>
#include <vector>

struct Element;
struct ElementJunction;
class GWModel;

struct ElementCell
{


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
    Variable *gradHydHead;

    /*!
     * \brief externalInflow
     */
    double externalInflow;

    /*!
     * \brief totalMassBalance
     */
    double totalMassBalance;

    /*!
     * \brief totalExternalInflow
     */
    double totalExternalInflow;

    /*!
     * \brief temperature (°C)
     */
    Variable temperature;

    /*!
     * \brief prevTemperature (°C)
     */
    Variable prevTemperature;

    /*!
     * \brief externalHeatFluxes
     */
    double externalHeatFluxes;

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
     * \brief externalSoluteFluxes
     */
    double *externalSoluteFluxes;

    /*!
     * \brief totalSoluteMassBalance
     */
    double *totalSoluteMassBalance;

    /*!
     * \brief hydConX
     */
    double hydConX;

    /*!
     * \brief hydConY
     */
    double hydConY;

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
     * \brief width
     */
    double width;

    /*!
     * \brief specificStorage
     */
    double specificStorage;

    /*!
     * \brief dispersionX
     */
    double dispersionX;

    /*!
     * \brief dispersionY
     */
    double dispersionY;

    /*!
     * \brief edgeFlows
     */
    double *edgeFlows;

    /*!
     * \brief satDepth
     */
    double satDepth;

    /*!
     * \brief edgeDepths
     */
    double *edgeDepths;

    /*!
     * \brief depth
     */
    double depth;

    /*!
     * \brief edgeBottomElevs
     */
    double *edgeBottomElevs;

    /*!
     * \brief edgeHydCons
     */
    double *edgeHydCons;

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
    double prevVolume;

    /*!
     * \brief dVolumedt
     */
    double dVolumedt;

    /*!
     * \brief start
     */
    bool start;

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
     * \brief computeDSoluteDt
     * \param dt
     * \param S
     * \return
     */
    double computeDSoluteDt(double dt, double S[], int soluteIndex);

    /*!
     * \brief computeCourantNumber
     * \return
     */
    double computeCourantFactor() const;

    /*!
     * \brief computeDerivedHydraulics
     */
    void computeDerivedHydraulics();

    /*!
     * \brief computeDVolumeDt
     */
    void computeDVolumeDt();

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

  private:

    void computeEdgeHydCons();

    void computeEdgeDepths();

    void computeEdgeDepths(double H[]);

  private:

    Element *parentElement;

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
    * \brief channelBedThickness
    */
   double channelBedThickness;


   /*!
    * \brief channelSoluteConcs
    */
   double *channelSoluteConcs;

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

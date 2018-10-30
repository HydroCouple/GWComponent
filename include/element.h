#ifndef ELEMENT_H
#define ELEMENT_H

#include "variable.h"
#include "gwcomponent_global.h"

struct Element;
struct ElementJunction;
class GWModel;

struct ElementCell
{
    enum Orientation
    {
      Left,
      Right
    };

    ElementCell(int index,
                double widthFactor,
                Orientation orientation,
                Element *parent);

    ~ElementCell();

    int index ;

    Orientation orientation;

    /*!
     * \brief hydHead
     */
    Variable hydHead;

    /*!
     * \brief prevHydHead
     */
    Variable prevHydHead;

    /*!
     * \brief gradHydHeadX
     */
    Variable gradHydHeadX;

    /*!
     * \brief gradHydHeadY
     */
    Variable gradHydHeadY;

    /*!
     * \brief temperature (°C)
     */
    Variable temperature;

    /*!
     * \brief prevTemperature (°C)
     */
    Variable prevTemperature;

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
     * \brief dispersionX
     */
    double dispersionX;

    /*!
     * \brief dispersionY
     */
    double dispersionY;

    /*!
     * \brief volume
     */
    double volume;

    /*!
     * \brief prevVolume
     */
    double prevVolume;

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
    double computeDSoluteDt(double dt, double S[]);

    /*!
     * \brief computeCourantNumber
     * \return
     */
    double computeCourantFactor() const;

    /*!
     * \brief commputeDispersionFactor
     * \return
     */
    double computeDispersionFactor() const;

    /*!
     * \brief computeDerivedHydraulics
     */
    void computeDerivedHydraulics();

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

   int numLeftCells;

   double leftWidth;

   double *leftCellWidthFactors;

   int numRightCells;

   double rightWidth;

   double *rightCellWidthFactors;

   double length;

   ElementCell *elementCells;

   void initialize();

   void initializeElementCells();

   void initializeSolutes();

}

#endif // ELEMENT_H

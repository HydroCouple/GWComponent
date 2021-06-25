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
 * \brief This struct represents the channel control volume
 */
struct GWCOMPONENT_EXPORT Element
{
    friend struct ElementCell;

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

   /*!
    * \brief computeChannelMassFlux
    */
   void computeChannelMassFlux();

   /*!
    * \brief computeChannelHeatFlux
    */
   void computeChannelHeatFlux();

   /*!
    * \brief computeChannelSoluteFlux
    * \param soluteIndex
    */
   void computeChannelSoluteFlux(int soluteIndex);


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

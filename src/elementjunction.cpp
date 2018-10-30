/*!
*  \file    ElementJunction.cpp
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


#include "stdafx.h"
#include "elementjunction.h"
#include "element.h"

#include <math.h>

ElementJunction::ElementJunction(const std::string &id, double x, double y, double z, HTSModel *model)
  :id(id), x(x), y(y), z(z),
    model(model)
{

}

ElementJunction::~ElementJunction()
{

  while (outgoingElements.size())
  {
    Element *element = *outgoingElements.begin();
    delete element;
  }

  while (incomingElements.size())
  {
    Element *element = *incomingElements.begin();
    delete element;
  }
}

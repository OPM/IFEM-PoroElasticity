// $Id$
//==============================================================================
//!
//! \file SIMPoroElasticity.C
//!
//! \date May 29 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Simulation driver for poroelasticity problems.
//!
//==============================================================================

#include "SIMPoroElasticity.h"
#include "PoroSolutions.h"
#include "PoroMaterial.h"

#include "AnaSol.h"
#include "Functions.h"


template<> bool SIMPoroEl2D::parseDimSpecific (const TiXmlElement* elem)
{
  std::string type;
  utl::getAttribute(elem,"type",type,true);
  if (type == "terzhagi")
  {
    double height = 1.0, load = 1.0;
    utl::getAttribute(elem,"height",height);
    utl::getAttribute(elem,"load",load);
    IFEM::cout <<"\tAnalytical solution: Terzhagi, height = "<< height
               <<" load = "<< load << std::endl;
    Elasticity* elp = this->getIntegrand();
    PoroMaterial* pmat = static_cast<PoroMaterial*>(elp->getMaterial());
    double gacc = elp->getGravity().length();
    SIM2D::mySol = new AnaSol(new TerzhagiPressure(pmat,gacc,height,load));
  }
  else if (type == "terzhagi-stationary")
  {
    double load = 1.0;
    utl::getAttribute(elem,"load",load);
    IFEM::cout <<"\tAnalytical solution: Terzhagi (stationary),"
               <<" load = "<< load << std::endl;
    Elasticity* elp = this->getIntegrand();
    PoroMaterial* pmat = static_cast<PoroMaterial*>(elp->getMaterial());
    SIM2D::mySol = new AnaSol(new ConstFunc(0.0), nullptr,
                              new StationaryTerzhagiDisplacement(pmat,load));
  }
  else
    return false;

  return true;
}

// $Id$
//==============================================================================
//!
//! \file PoroMaterialTexture.C
//!
//! \date Mar 8 2019
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for poroelastic material models with textures.
//!
//==============================================================================

#include "PoroMaterialTexture.h"
#include "Utilities.h"
#include "Functions.h"
#include "Tensor.h"
#include "Vec3Oper.h"
#include "MatVec.h"
#include "IFEM.h"
#include "tinyxml.h"


void PoroMaterialTexture::parse (const TiXmlElement* elem)
{
  props.parse(elem);
  this->PoroMaterial::parse(elem);
}


void PoroMaterialTexture::printLog () const
{
  IFEM::cout <<"\tConstitutive Properties:";
  props.printLog();
  IFEM::cout << std::endl;
}


double PoroMaterialTexture::getPorosity (const Vec3& X) const
{
  double result;
  if (!props.getProperty("poro",X,result))
    return this->PoroMaterial::getPorosity(X);

  return result;
}


Vec3 PoroMaterialTexture::getPermeability (const Vec3& X) const
{
  Vec3 result;
  if (!props.getProperty("permx",X,result[0]))
    return this->PoroMaterial::getPermeability(X);

  if (!props.getProperty("permy",X,result[1]))
    result[1] =  result[0];
  if (!props.getProperty("permz",X,result[2]))
    result[2] = result[0];

  return result;
}


double PoroMaterialTexture::getBulkFluid(const Vec3& X) const
{
  double result;
  if (!props.getProperty("bulkfluid",X,result))
    return this->PoroMaterial::getBulkFluid(X);

  return result;
}


double PoroMaterialTexture::getBulkSolid(const Vec3& X) const
{
  double result;
  if (!props.getProperty("bulksolid",X,result))
    return this->PoroMaterial::getBulkSolid(X);

  return result;
}


double PoroMaterialTexture::getBulkMedium(const Vec3& X) const
{
  double result;
  if (!props.getProperty("bulkmedium",X,result))
    return this->PoroMaterial::getBulkMedium(X);

  return result;
}


double PoroMaterialTexture::getSolidDensity (const Vec3& X) const
{
  double result;
  if (!props.getProperty("rhos",X,result))
    return this->PoroMaterial::getFluidDensity(X);

  return result;
}


double PoroMaterialTexture::getFluidDensity (const Vec3& X) const
{
  double result;
  if (!props.getProperty("rhof",X,result))
    return this->PoroMaterial::getFluidDensity(X);

  return result;
}


double PoroMaterialTexture::getViscosity (const Vec3& X) const
{
  double result;
  if (!props.getProperty("mu",X,result))
    return this->PoroMaterial::getViscosity(X);

  return result;
}


double PoroMaterialTexture::getStiffness (const Vec3& X) const
{
  double result;
  if (!props.getProperty("stiffness",X,result))
    return this->PoroMaterial::getStiffness(X);

  return result;
}


double PoroMaterialTexture::getPoisson (const Vec3& X) const
{
  double result;
  if (!props.getProperty("poisson",X,result))
    return this->PoroMaterial::getPoisson(X);

  return result;
}

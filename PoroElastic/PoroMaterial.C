// $Id$
//==============================================================================
//!
//! \file PoroMaterial.C
//!
//! \date Apr 29 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for poroelastic material models.
//!
//==============================================================================

#include "PoroMaterial.h"
#include "Utilities.h"
#include "Functions.h"
#include "Tensor.h"
#include "Vec3Oper.h"
#include "MatVec.h"
#include "IFEM.h"
#include "tinyxml.h"


template<>
RealFunc* PoroMaterial::FuncConstPair<RealFunc>::parse(const char* val,
                                                       const std::string& type)
{
  return utl::parseRealFunc(val, type);
}


template<>
ScalarFunc* PoroMaterial::FuncConstPair<ScalarFunc>::parse(const char* val,
                                                           const std::string& type)
{
  return utl::parseTimeFunc(val, type);
}


template<>
VecFunc* PoroMaterial::FuncConstPair<VecFunc>::parse(const char* val,
                                                     const std::string& type)
{
  return utl::parseVecFunc(val, type);
}


/*!
  \brief Template function to parse a property value from and XML-element.
*/

template<class T>
static bool propertyParse(PoroMaterial::FuncConstPair<T>& data,
                          const TiXmlElement* elem,
                          const char* attr, const char* tag)
{
  std::string constant;
  if (utl::getAttribute(elem,attr,constant)) {
    std::stringstream str;
    str << constant;
    str >> data.constant;
    return true;
  }

  const TiXmlElement* child = elem->FirstChildElement(tag);
  if (!child) return false;
  const TiXmlNode* aval = child->FirstChild();
  if (!aval) return false;

  IFEM::cout <<" ";
  std::string type;
  utl::getAttribute(child,"type",type,true);
  data.function = data.parse(aval->Value(),type);

  return data.function != nullptr;
}


void PoroMaterial::parse(const TiXmlElement* elem)
{
  propertyParse(Emod, elem, "E", "stiffness");
  propertyParse(nu, elem, "nu", "poisson");
  propertyParse(rhof, elem, "rhof", "fluiddensity");
  propertyParse(rhos, elem, "rhos", "soliddensity");

  propertyParse(thermalexpansion, elem, "alpha", "thermalexpansion");
  propertyParse(heatcapacity, elem, "cp", "heatcapacity");
  propertyParse(conductivity, elem, "kappa", "conductivity");

  propertyParse(porosity, elem, "poro", "porosity");
  propertyParse(permeability, elem, "perm", "permeability");
  propertyParse(bulkw, elem, "Kw", "waterbulk");
  propertyParse(bulks, elem, "Ks", "solidbulk");
  propertyParse(bulkm, elem, "Ko", "mediumbulk");
}


void PoroMaterial::printLog() const
{
  IFEM::cout <<"\tConstitutive Properties: "
             <<"\n\t\tYoung's Modulus, E = "<< Emod.constant
             <<"\n\t\tPoisson's Ratio, nu = "<< nu.constant;
  IFEM::cout <<"\n\tDensities: "
             <<"\n\t\tDensity of Fluid, rhof = "<< rhof.constant
             <<"\n\t\tDensity of Solid, rhos = "<< rhos.constant;
  IFEM::cout <<"\n\tBulk Moduli: "
             <<"\n\t\tBulk Modulus of Water, Kw = "<< bulkw.constant
             <<"\n\t\tBulk Modulus of Solid, Ks = "<< bulks.constant
             <<"\n\t\tBulk Modulus of Medium, Ko = "<< bulkm.constant;
  IFEM::cout <<"\n\tPorosity, n = "<< porosity.constant << std::endl;
}


double PoroMaterial::getThermalExpansion (double T) const
{
  return thermalexpansion.evaluate(T);
}


double PoroMaterial::getHeatCapacity (double T) const
{
  return heatcapacity.evaluate(T);
}


double PoroMaterial::getThermalConductivity(double T) const
{
  return conductivity.evaluate(T);
}


double PoroMaterial::getPorosity(const Vec3& X) const
{
  return porosity.evaluate(X);
}


Vec3 PoroMaterial::getPermeability(const Vec3& X) const
{
  return permeability.evaluate(X);
}


double PoroMaterial::getFluidDensity(const Vec3& X) const
{
  return rhof.evaluate(X);
}


double PoroMaterial::getSolidDensity(const Vec3& X) const
{
  return rhos.evaluate(X);
}


double PoroMaterial::getMassDensity(const Vec3& X) const
{
  double poro = porosity.evaluate(X);
  return rhos.evaluate(X)*(1.0-poro) + rhof.evaluate(X)*poro;
}


double PoroMaterial::getBulkWater(const Vec3& X) const
{
  return bulkw.evaluate(X);
}


double PoroMaterial::getBulkSolid(const Vec3& X) const
{
  return bulks.evaluate(X);
}


double PoroMaterial::getBulkMedium(const Vec3& X) const
{
  return bulkm.evaluate(X);
}


double PoroMaterial::getStiffness(const Vec3& X) const
{
  return Emod.evaluate(X);
}


double PoroMaterial::getPoisson(const Vec3& X) const
{
  return nu.evaluate(X);
}


bool PoroMaterial::evaluate (Matrix& Cmat, SymmTensor& sigma, double& U,
                             const FiniteElement&, const Vec3& X,
                             const Tensor&, const SymmTensor& eps, char iop,
                             const TimeDomain*, const Tensor*) const
{
  double E = Emod.evaluate(X);
  double v = nu.evaluate(X);

  const size_t nsd = eps.dim();
  const size_t nstrc = nsd*(nsd+1)/2;
  Cmat.resize(nstrc,nstrc,true);

  if (nsd == 1)
  {
    // Constitutive matrix for 1D
    Cmat(1,1) = E;
    if (iop > 0)
    {
      sigma = eps; sigma *= E;
      if (iop == 3)
        U = 0.5*sigma(1,1)*eps(1,1);
    }
    return true;
  }
  else if (v < 0.0 || v >= 0.5)
  {
    std::cerr <<" *** PoroMaterial::evaluate: Poisson's ratio "<< v
              <<" out of range [0,0.5>."<< std::endl;
    return false;
  }

  double C3 = E/(2.0+2.0*v);
  double C2 = C3*v/(0.5-v);
  double C1 = C2 + 2.0*C3;

  // Constitutive matrix for 2D plane-strain and 3D
  for (size_t i = 1; i <= nstrc; i++)
    for (size_t j = 1; j <= nstrc; j++)
      if (i == j)
        Cmat(i,j) = i <= nsd ? C1 : C3;
      else if (i <= nsd && j <= nsd)
        Cmat(i,j) = C2;

  if (iop > 0) // Calculate the stress tensor, sigma = C*eps
    if (!Cmat.multiply(eps,sigma))
      return false;

  if (iop == 3) // Calculate strain energy density, // U = 0.5*sigma:eps
    U = 0.5*sigma.innerProd(eps);

  return true;
}


bool PoroMaterial::evaluate (double& lambda, double& mu,
                             const FiniteElement& fe, const Vec3& X) const
{
  double E = Emod.evaluate(X);
  double v = nu.evaluate(X);

  if (v < 0.0 || v >= 0.5)
  {
    std::cerr <<" *** PoroMaterial::evaluate: Poisson's ratio "<< v
              <<" out of range [0,0.5>."<< std::endl;
    return false;
  }

  // Evaluate the Lame parameters
  mu = 0.5*E/(1.0+v);
  lambda = mu*v/(0.5-v);

  return true;
}

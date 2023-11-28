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


/*!
  \brief Parses an XML element. Specialization for RealFunc.
*/

template<>
RealFunc* PoroMaterial::FuncConstPair<RealFunc>::parse (const char* val,
                                                        const std::string& type)
{
  return utl::parseRealFunc(val, type);
}


/*!
  \brief Parses an XML element. Specialization for ScalarFunc.
*/

template<>
ScalarFunc* PoroMaterial::FuncConstPair<ScalarFunc>::parse (const char* val,
                                                            const std::string& type)
{
  return utl::parseTimeFunc(val, type);
}


/*!
  \brief Parses an XML element. Specialization for VecFunc.
*/

template<>
VecFunc* PoroMaterial::FuncConstPair<VecFunc>::parse (const char* val,
                                                      const std::string& type)
{
  return utl::parseVecFunc(val, type);
}


void PoroMaterial::parse (const tinyxml2::XMLElement* elem)
{
  Emod.propertyParse(elem, "E", "stiffness");
  nu.propertyParse(elem, "nu", "poisson");

  rhof.propertyParse(elem, "rhof", "fluiddensity");
  rhos.propertyParse(elem, "rhos", "soliddensity");
  viscosity.propertyParse(elem, "mu", "viscosity");

  porosity.propertyParse(elem, "poro", "porosity");
  permeability.propertyParse(elem, "perm", "permeability");
  bulkf.propertyParse(elem, "Kf", "fluidbulk");
  bulks.propertyParse(elem, "Ks", "solidbulk");
  bulkm.propertyParse(elem, "Ko", "mediumbulk");

  utl::getAttribute(elem, "alpha", alpha);
  utl::getAttribute(elem, "Minv", Minv);
}


void PoroMaterial::printLog () const
{
  IFEM::cout <<"\tConstitutive Properties:"
             <<"\n\t\tYoung's Modulus, E = "<< Emod.constant
             <<"\n\t\tPoisson's Ratio, nu = "<< nu.constant;
  IFEM::cout <<"\n\tDensities:"
             <<"\n\t\tDensity of Fluid, rhof = "<< rhof.constant
             <<"\n\t\tDensity of Solid, rhos = "<< rhos.constant;
  IFEM::cout <<"\n\tBulk Moduli:";
  if (alpha >= 0.0 && Minv >= 0.0)
    IFEM::cout <<"\n\t\tBiot's coefficient, alpha = "<< alpha
               <<"\n\t\tBiot's inverse modulus, M^-1 = "<< Minv;
  else
    IFEM::cout <<"\n\t\tBulk Modulus of Fluid, Kf = "<< bulkf.constant
               <<"\n\t\tBulk Modulus of Solid, Ks = "<< bulks.constant
               <<"\n\t\tBulk Modulus of Medium, Ko = "<< bulkm.constant;
  IFEM::cout <<"\n\tPorosity, n = "<< porosity.constant
             <<"\n\tPermeability, K = "<< permeability.constant << std::endl;
}


double PoroMaterial::getPorosity (const Vec3& X) const
{
  return porosity.evaluate(X);
}


Vec3 PoroMaterial::getPermeability (const Vec3& X) const
{
  return permeability.evaluate(X);
}


double PoroMaterial::getFluidDensity (const Vec3& X) const
{
  return rhof.evaluate(X);
}


double PoroMaterial::getSolidDensity (const Vec3& X) const
{
  return rhos.evaluate(X);
}


double PoroMaterial::getMassDensity (const Vec3& X) const
{
  double poro = porosity.evaluate(X);
  return rhos.evaluate(X)*(1.0-poro) + rhof.evaluate(X)*poro;
}


double PoroMaterial::getViscosity (const Vec3& X) const
{
  return viscosity.evaluate(X);
}


double PoroMaterial::getBulkFluid (const Vec3& X) const
{
  return bulkf.evaluate(X);
}


double PoroMaterial::getBulkSolid (const Vec3& X) const
{
  return bulks.evaluate(X);
}


double PoroMaterial::getBulkMedium (const Vec3& X) const
{
  return bulkm.evaluate(X);
}


double PoroMaterial::getBiotCoeff (const Vec3& X) const
{
  if (alpha >= 0.0 && alpha <= 1.0)
    return alpha;

  return 1.0 - bulkm.evaluate(X)/bulks.evaluate(X);
}


double PoroMaterial::getBiotModulus (const Vec3& X, double al, double po) const
{
  if (Minv >= 0.0)
    return Minv;

  return (al-po)/bulks.evaluate(X) + po/bulkf.evaluate(X);
}


double PoroMaterial::getStiffness (const Vec3& X) const
{
  return Emod.evaluate(X);
}


double PoroMaterial::getPoisson (const Vec3& X) const
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
                             const FiniteElement&, const Vec3& X) const
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

// $Id$
//==============================================================================
//!
//! \file PoroMaterial.h
//!
//! \date Apr 29 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for poroelastic material models.
//!
//==============================================================================

#ifndef _PORO_MATERIAL_H
#define _PORO_MATERIAL_H

#include "MaterialBase.h"
#include "Function.h"

class TiXmlElement;


/*!
  \brief Class representing a material model for a poroelastic problem.
*/

class PoroMaterial : public Material
{
public:
  //! \brief Helper template for wrapping a constant/function pair.
  template<class Function> struct FuncConstPair
  {
    Function* function;                 //!< Function definition
    typename Function::Output constant; //!< Constant

    //! \brief Default constructor.
    FuncConstPair() { function = nullptr; constant = 0.0; }

    //! \brief Parses an XML element. Specialized per function type.
    Function* parse(const char*, const std::string&) { return nullptr; }

    //! \brief Evaluates the function at the given point \b X.
    typename Function::Output evaluate(const typename Function::Input& X) const
    {
      return function ? (*function)(X) : constant;
    }
  };

  //! \brief The constructor initializes Biot's parameters to -1 (undefined).
  PoroMaterial() { alpha = Minv = -1.0; }
  //! \brief Empty destructor.
  virtual ~PoroMaterial() {}

  //! \brief Parses material parameters from an XML element.
  virtual void parse(const TiXmlElement* elem);

  //! \brief Prints out material parameters to the log stream.
  virtual void printLog() const;

  //! \brief Evaluates the mass density of the fluid at current point.
  double getFluidDensity(const Vec3&) const;
  //! \brief Evaluates the mass density of the solid at current point.
  double getSolidDensity(const Vec3&) const;
  //! \brief Evaluates the mass density at current point.
  virtual double getMassDensity(const Vec3&) const;
  //! \brief Returns the dynamic viscosity at current point.
  double getViscosity(const Vec3& X) const;
  //! \brief Returns porosity at the current point.
  double getPorosity(const Vec3& X) const;
  //! \brief Returns permeability at the current point.
  Vec3 getPermeability(const Vec3& X) const;
  //! \brief Returns bulk modulus of the fluid at the current point.
  double getBulkFluid(const Vec3& X) const;
  //! \brief Returns bulk modulus of the solid at the current point.
  double getBulkSolid(const Vec3& X) const;
  //! \brief Returns bulk modulus of the medium at the current point.
  double getBulkMedium(const Vec3& X) const;
  //! \brief Returns Biot's coefficient at the current point.
  double getBiotCoeff(const Vec3& X) const;
  //! \brief Returns the inverse Biot's modulus at the current point.
  double getBiotModulus(const Vec3& X, double al, double po) const;
  //! \brief Returns stiffness at the current point.
  virtual double getStiffness(const Vec3& X) const;
  //! \brief Returns Poisson's ratio at the current point.
  virtual double getPoisson(const Vec3& X) const;

  //! \brief Evaluates the constitutive relation at an integration point.
  //! \param[out] Cmat Constitutive matrix at current point
  //! \param[out] sigma Stress tensor at current point
  //! \param[out] U Strain energy density at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] eps Strain tensor at current point
  //! \param[in] iop Calculation option;
  //!   0 : Calculate the constitutive matrix only,
  //!   1 : Calculate Cauchy stresses and the constitutive matrix,
  //!   3 : Calculate the strain energy density only.
  virtual bool evaluate(Matrix& Cmat, SymmTensor& sigma, double& U,
                        const FiniteElement&, const Vec3& X,
                        const Tensor&, const SymmTensor& eps,
                        char iop = 1, const TimeDomain* = nullptr,
                        const Tensor* = nullptr) const;

  //! \brief Evaluates the Lame-parameters at an integration point.
  //! \param[out] lambda Lame's first parameter
  //! \param[out] mu Lame's second parameter (shear modulus)
  //! \param[in] X Cartesian coordinates of current point
  virtual bool evaluate(double& lambda, double& mu,
                        const FiniteElement&, const Vec3& X) const;

protected:
  FuncConstPair<RealFunc> Emod;         //!< Young's modulus
  FuncConstPair<RealFunc> nu;           //!< Poisson's ratio
  FuncConstPair<RealFunc> rhof;         //!< Fluid density
  FuncConstPair<RealFunc> rhos;         //!< Solid density
  FuncConstPair<RealFunc> viscosity;    //!< Dynamic viscosity

  FuncConstPair<RealFunc> porosity;     //!< Porosity
  FuncConstPair<VecFunc>  permeability; //!< Permeability
  FuncConstPair<RealFunc> bulkf;        //!< Bulk modulus of fluid
  FuncConstPair<RealFunc> bulks;        //!< Bulk modulus of solid
  FuncConstPair<RealFunc> bulkm;        //!< Bulk modulus of medium

  double alpha; //!< Biot's coefficient
  double Minv;  //!< Inverse Biot's modulus (compressibility modulus)
};

#endif

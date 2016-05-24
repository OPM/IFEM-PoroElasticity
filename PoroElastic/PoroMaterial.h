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
#include "MatVec.h"
#include "Vec3.h"

class TiXmlElement;


/*!
  \brief Class representing a material model for a poroelastic problem.
*/

class PoroMaterial : public Material
{
public:
  // \brief Helper template for wrapping a constant/function pair.
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

  //! \brief Empty constructor.
  PoroMaterial() {}
  //! \brief Empty destructor.
  virtual ~PoroMaterial() {}

  //! \brief Parses material parementers from an XML element.
  virtual void parse(const TiXmlElement* elem);

  //! \brief Prints out material parameters to the log stream.
  virtual void printLog() const;

  //! \brief Evaluates the mass density of the fluid at current point.
  double getFluidDensity(const Vec3&) const;
  //! \brief Evaluates the mass density of the solid at current point.
  double getSolidDensity(const Vec3&) const;
  //! \brief Evaluates the mass density at current point.
  virtual double getMassDensity(const Vec3&) const;
  //! \brief Evaluates the effective heat capacity at the current point
  virtual double getHeatCapacity(double T) const;
  //! \brief Evaluates the heat capacity at the current point
  double getFluidHeatCapacity(double T) const;
  //! \brief Evaluates the heat capacity at the current point
  double getSolidHeatCapacity(double T) const;
  //! \brief Evaluates the effective thermal conductivity at the current point
  virtual double getThermalConductivity(double T) const;
  //! \brief Evaluates the thermal conductivity of the fluid at the current point
  double getFluidThermalConductivity(double T) const;
  //! \brief Evaluates the thermal conductivity of the solid at the current point
  double getSolidThermalConductivity(double T) const;
  //! \brief Evaluates the thermal expansion of the solid at the current point
  double getSolidThermalExpansion(double T) const;
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
  //! \param[in] fe Finite element quantities at current point
  //! \param[in] X Cartesian coordinates of current point
  virtual bool evaluate(double& lambda, double& mu,
                        const FiniteElement& fe, const Vec3& X) const;

  //! \brief Get scaling factor.
  //! \param X Point to evaluate parameters in
  //! \param time Time stepping parameters
  //! \param gacc Gravity acceleration
  double getScaling(const Vec3& X, const TimeDomain& time,
                    double gacc) const;

protected:
  FuncConstPair<RealFunc> Emod; //!< Young's modulus
  FuncConstPair<RealFunc> nu;   //!< Poisson's ratio
  FuncConstPair<RealFunc> rhof; //!< Fluid density
  FuncConstPair<RealFunc> rhos; //!< Solid density

  FuncConstPair<ScalarFunc> fheatcapacity; //!< Specific heat capacity for fluid
  FuncConstPair<ScalarFunc> sheatcapacity; //!< Specific heat capacity for solid
  FuncConstPair<ScalarFunc> fconductivity; //!< Thermal conductivity
  FuncConstPair<ScalarFunc> sconductivity; //!< Thermal conductivity
  FuncConstPair<ScalarFunc> sexpansion;    //!< Thermal expansion

  FuncConstPair<RealFunc> porosity;     //!< Porosity
  FuncConstPair<VecFunc>  permeability; //!< Permeability
  FuncConstPair<RealFunc> bulkf;        //!< Bulk modulus of fluid 
  FuncConstPair<RealFunc> bulks;        //!< Bulk modulus of solid
  FuncConstPair<RealFunc> bulkm;        //!< Bulk modulus of medium
};

#endif

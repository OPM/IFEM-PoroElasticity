// $Id$
//==============================================================================
//!
//! \file PoroMaterialTexture.h
//!
//! \date Mar 7 2019
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for poroelastic material models with textures.
//!
//==============================================================================

#ifndef _PORO_MATERIAL_TEXTURE_H
#define _PORO_MATERIAL_TEXTURE_H

#include "PoroMaterial.h"
#include "TextureProperties.h"

class TiXmlElement;


/*!
  \brief Class representing a material model for a poroelastic problem.
*/

class PoroMaterialTexture : public PoroMaterial
{
public:
  //! \brief The constructor initializes Biot's parameters to -1 (undefined).
  PoroMaterialTexture() = default;
  //! \brief Empty destructor.
  virtual ~PoroMaterialTexture() {}

  //! \brief Parses material parameters from an XML element.
  void parse(const TiXmlElement* elem) override;

  //! \brief Prints out material parameters to the log stream.
  void printLog() const override;

  //! \brief Evaluates the mass density of the fluid at current point.
  double getFluidDensity(const Vec3&) const override;
  //! \brief Evaluates the mass density of the solid at current point.
  double getSolidDensity(const Vec3&) const override;
  //! \brief Returns the dynamic viscosity at current point.
  double getViscosity(const Vec3& X) const override;
  //! \brief Returns porosity at the current point.
  double getPorosity(const Vec3& X) const override;
  //! \brief Returns permeability at the current point.
  Vec3 getPermeability(const Vec3& X) const override;
  //! \brief Returns bulk modulus of the fluid at the current point.
  double getBulkFluid(const Vec3& X) const override;
  //! \brief Returns bulk modulus of the solid at the current point.
  double getBulkSolid(const Vec3& X) const override;
  //! \brief Returns bulk modulus of the medium at the current point.
  double getBulkMedium(const Vec3& X) const override;
  //! \brief Returns stiffness at the current point.
  virtual double getStiffness(const Vec3& X) const override;
  //! \brief Returns Poisson's ratio at the current point.
  virtual double getPoisson(const Vec3& X) const override;

protected:
  TextureProperties props; //!< Map of textured properties
};

#endif

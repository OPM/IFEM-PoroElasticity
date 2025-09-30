//==============================================================================
//!
//! \file TestSIMPoroMaterial.C
//!
//! \date May 13 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for class for poro-elastic material models.
//!
//==============================================================================

#include "PoroMaterial.h"

#include "Catch2Support.h"
#include <tinyxml2.h>


TEST_CASE("TestPoroMaterial.Parse")
{
  tinyxml2::XMLDocument doc;
  doc.LoadFile("Plaxis1DVerif.xinp");
  REQUIRE(doc.RootElement());

  const tinyxml2::XMLElement* elem = doc.RootElement()->FirstChildElement("poroelasticity");
  REQUIRE(elem != nullptr);

  const tinyxml2::XMLElement* iso = elem->FirstChildElement("isotropic");

  PoroMaterial mat;
  mat.parse(iso);

  Vec3 X;
  REQUIRE_THAT(mat.getFluidDensity(X), WithinRel(1000.0));
  REQUIRE_THAT(mat.getSolidDensity(X), WithinRel(2700.0));
  REQUIRE_THAT(mat.getViscosity(X), WithinRel(9810.0));
  REQUIRE_THAT(mat.getPorosity(X), WithinRel(0.5));
  REQUIRE_THAT(mat.getStiffness(X), WithinRel(1000000.0));
  REQUIRE_THAT(mat.getPoisson(X), WithinAbs(0.0, 1e-13));
  REQUIRE_THAT(mat.getBulkFluid(X), WithinRel(1e99));
  REQUIRE_THAT(mat.getBulkSolid(X), WithinRel(1e99));
  REQUIRE_THAT(mat.getBulkMedium(X), WithinAbs(0.0, 1e-13));
  Vec3 perm = mat.getPermeability(X);
  REQUIRE_THAT(perm[0], WithinRel(0.0000000115741));
  REQUIRE_THAT(perm[1], WithinRel(0.0000000115741));
  REQUIRE_THAT(perm[2], WithinAbs(0.0, 1e-13));
}

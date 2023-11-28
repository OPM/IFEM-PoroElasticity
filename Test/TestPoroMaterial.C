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

#include "gtest/gtest.h"
#include "tinyxml2.h"

TEST(TestPoroMaterial, Parse)
{
  tinyxml2::XMLDocument doc;
  doc.LoadFile("Plaxis1DVerif.xinp");
  if (!doc.RootElement())
    ASSERT_TRUE(false);

  const tinyxml2::XMLElement* elem = doc.RootElement()->FirstChildElement("poroelasticity");
  ASSERT_TRUE(elem != nullptr);

  const tinyxml2::XMLElement* iso = elem->FirstChildElement("isotropic");

  PoroMaterial mat;
  mat.parse(iso);

  Vec3 X;
  ASSERT_FLOAT_EQ(mat.getFluidDensity(X), 1000.0);
  ASSERT_FLOAT_EQ(mat.getSolidDensity(X), 2700.0);
  ASSERT_FLOAT_EQ(mat.getViscosity(X), 9810.0);
  ASSERT_FLOAT_EQ(mat.getPorosity(X), 0.5);
  ASSERT_FLOAT_EQ(mat.getStiffness(X), 1000000.0);
  ASSERT_FLOAT_EQ(mat.getPoisson(X), 0.0);
  ASSERT_FLOAT_EQ(mat.getBulkFluid(X), 1e99);
  ASSERT_FLOAT_EQ(mat.getBulkSolid(X), 1e99);
  ASSERT_FLOAT_EQ(mat.getBulkMedium(X), 0.0);
  Vec3 perm = mat.getPermeability(X);
  ASSERT_FLOAT_EQ(perm[0], 0.0000000115741);
  ASSERT_FLOAT_EQ(perm[1], 0.0000000115741);
  ASSERT_FLOAT_EQ(perm[2], 0.0);

}

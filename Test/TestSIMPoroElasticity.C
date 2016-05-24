//==============================================================================
//!
//! \file TestSIMPoroElasticity.C
//!
//! \date May 13 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for simulation driver for poroelasticity problems.
//!
//==============================================================================

#include "SIMPoroElasticity.h"
#include "SIM2D.h"

#include "gtest/gtest.h"


TEST(TestSIMPoroElasticity, Parse)
{
  SIMPoroElasticity<SIM2D> sim;
  EXPECT_TRUE(sim.read("Plaxis1DVerif.xinp"));

  const PoroElasticity* poro = static_cast<const PoroElasticity*>(sim.getProblem());

  Vec3 grav = poro->getGravity();
  ASSERT_FLOAT_EQ(grav.x,0.0);
  ASSERT_FLOAT_EQ(grav.y,0.0);
  ASSERT_FLOAT_EQ(grav.z,0.0);
  TimeDomain t;
  t.dt = 1.0;
  ASSERT_FLOAT_EQ(poro->getScaling(Vec3(), t), 920642222.26448686158677615343);
}

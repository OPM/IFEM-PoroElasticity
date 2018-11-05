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
#include "TimeStep.h"

#include "gtest/gtest.h"


TEST(TestSIMPoroElasticity, Parse)
{
  SIMPoroElasticity<SIM2D> sim;
  ASSERT_TRUE(sim.read("Plaxis1DVerif.xinp"));
  ASSERT_TRUE(sim.init(TimeStep()));

  const PoroElasticity* poro = static_cast<const PoroElasticity*>(sim.getProblem());

  Vec3 grav = poro->getGravity();
  ASSERT_FLOAT_EQ(grav.x,0.0);
  ASSERT_FLOAT_EQ(grav.y,9.81);
  ASSERT_FLOAT_EQ(grav.z,0.0);
}

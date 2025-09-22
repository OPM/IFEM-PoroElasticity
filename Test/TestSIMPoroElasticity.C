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
#include "PoroElasticity.h"
#include "SIM2D.h"
#include "TimeStep.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;


TEST_CASE("TestSIMPoroElasticity.Parse")
{
  SIMPoroElasticity<SIM2D> sim;
  REQUIRE(sim.read("Plaxis1DVerif.xinp"));
  REQUIRE(sim.init(TimeStep()));

  const PoroElasticity* poro = static_cast<const PoroElasticity*>(sim.getProblem());

  Vec3 grav = poro->getGravity();
  REQUIRE_THAT(grav.x, WithinAbs(0.0, 1e-14));
  REQUIRE_THAT(grav.y, WithinRel(9.81));
  REQUIRE_THAT(grav.z, WithinAbs(0.0, 1e-14));
}

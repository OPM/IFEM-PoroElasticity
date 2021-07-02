// $Id$
//==============================================================================
//!
//! \file SIMStatPoroElasticity.h
//!
//! \date Sep 10 2018
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Static simulation driver for poroelasticity problems.
//!
//==============================================================================

#ifndef _SIM_STAT_PORO_ELASTICITY_H_
#define _SIM_STAT_PORO_ELASTICITY_H_

#include "SIMPoroElasticity.h"
#include "PoroElasticity.h"

/*!
  \brief Driver class for quasi-static poroelasticity problems.
*/

template<class Dim, bool staticFlow> class SIMStatPoroElasticity : public SIMPoroElasticity<Dim>
{
public:
  //! \brief Default constructor.
  SIMStatPoroElasticity() {}

  //! \brief Constructor for mixed problems.
  explicit SIMStatPoroElasticity(const std::vector<unsigned char>& flds)
    : SIMPoroElasticity<Dim>(flds) {}

  //! \brief Empty destructor.
  virtual ~SIMStatPoroElasticity() {}

  //! \brief Computes the solution for the current load step.
  virtual bool solveStep(TimeStep& tp)
  {
    this->printStep(tp.step,tp.time);

    Vector empty;
    if (!this->updateDirichlet(tp.time.t,&empty))
      return false;

    double oldtol = utl::zero_print_tol;
    utl::zero_print_tol = 1.0e-8;
    bool ok = this->SIMPoroElasticity<Dim>::solveStep(tp);
    utl::zero_print_tol = oldtol;

    return ok;
  }

protected:
  //! \brief Returns the actual integrand.
  virtual Elasticity* getIntegrand()
  {
    if (!Dim::myProblem)
      Dim::myProblem = new PoroElasticity(Dim::dimension, Dim::nf.size() > 1, staticFlow);
    return static_cast<Elasticity*>(Dim::myProblem);
  }
};

#endif

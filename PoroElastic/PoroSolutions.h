// $Id$
//==============================================================================
//!
//! \file PoroSolutions.h
//!
//! \date Jun 22 2015
//!
//! \author Arne Morten Kvarving
//!
//! \brief Analytic solutions for PoroElasticity problems.
//!
//==============================================================================

#ifndef _PORO_SOLUTIONS_H
#define _PORO_SOLUTIONS_H

#include "Function.h"


class PoroElasticity;


class TerzhagiPressure : public RealFunc
{
public:
  //! \brief Empty constructor.
  TerzhagiPressure(const PoroElasticity* mat, double h, double l)
    : m_mat(mat), height(h), load(l) {}
  //! \brief Empty destructor.
  virtual ~TerzhagiPressure() {}

protected:
  //! \brief Evaluates the analytic pressure field at the point \a X.
  virtual double evaluate(const Vec3& X) const;

  const PoroElasticity* m_mat; //!< Pointer to integrand with material parameters

private:
  double height, load;
};

#endif
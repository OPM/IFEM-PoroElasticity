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

class PoroMaterial;


/*!
  \brief Analytic pressure solution for the Terzhagi problem.
*/

class TerzhagiPressure : public RealFunc
{
public:
  //! \brief The constructor initializes the data members.
  TerzhagiPressure(PoroMaterial* m, double g, double h, double l)
    : myMat(m), gacc(g), height(h), load(l) {}
  //! \brief Empty destructor.
  virtual ~TerzhagiPressure() {}

protected:
  //! \brief Evaluates the analytic pressure field at the point \a X.
  virtual double evaluate(const Vec3& X) const;

private:
  PoroMaterial* myMat;  //!< Pointer to material parameters object
  double        gacc;   //!< Gravitation constant
  double        height; //!< Domain height
  double        load;   //!< Load factor
};


/*!
  \brief Analytic displacement solution for the Terzhagi problem.
*/

class StationaryTerzhagiDisplacement : public VecFunc
{
public:
  //! \brief The constructor initializes the data members.
  StationaryTerzhagiDisplacement(PoroMaterial* m, double l)
    : myMat(m), load(l) {}
  //! \brief Empty destructor.
  virtual ~StationaryTerzhagiDisplacement() {}

protected:
  //! \brief Evaluates the analytic displacement field at the point \a X.
  virtual Vec3 evaluate(const Vec3& X) const;

private:
  PoroMaterial* myMat; //!< Pointer to material parameters object
  double        load;  //!< Load factor
};

#endif

// $Id$
//==============================================================================
//!
//! \file PoroSolutions.C
//!
//! \date Jun 22 2015
//!
//! \author Arne Morten Kvarving
//!
//! \brief Analytic solutions for PoroElasticity problems.
//!
//==============================================================================

#include "PoroSolutions.h"
#include "PoroMaterial.h"
#include "Vec3.h"


double TerzhagiPressure::evaluate (const Vec3& X) const
{
  double k = myMat->getPermeability(X).x;
  double rhof = myMat->getFluidDensity(X);
  double E = myMat->getStiffness(X);
  double nu = myMat->getPoisson(X);
  double n = myMat->getPorosity(X);
  double alpha = myMat->getBiotCoeff(X);
  double Minv = myMat->getBiotModulus(X,alpha,n);

  double K  = E / (3.0 - 6.0*nu);
  double G  = E / (2.0 + 2.0*nu);
  double cv = k / (rhof * gacc * (alpha*alpha / (K + G*4.0/3.0) + Minv));

  const Vec4& Xt = static_cast<const Vec4&>(X);
  double stime = cv * Xt.t / (height * height);

  int i, j;
  double p_ret = 0.0;
  double p_add = 1.0;
  for (i = j = 1; fabs(p_add) > 1.0e-15; i++, j += 2) {
    double trig_arg = j * M_PI * X.y * 0.5 / height;
    double exp_c = exp(-j*j * M_PI * M_PI * 0.25 * stime);
    double p_fac = 4.0 * load / (j * M_PI);
    p_add = cos(trig_arg) * exp_c * p_fac;
    p_ret += (i%2) ? p_add : -p_add;
  }

  return p_ret;
}


Vec3 StationaryTerzhagiDisplacement::evaluate (const Vec3& X) const
{
  double E = myMat->getStiffness(X);
  double nu = myMat->getPoisson(X);
  double c = -load/E * (1.0+nu) * (1.0-2.0*nu) / (1.0-nu);

  return Vec3(0.0, c*X.y, 0.0);
}

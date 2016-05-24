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
#include "PoroElasticity.h"
#include "PoroMaterial.h"
#include "Vec3.h"


double TerzhagiPressure::evaluate (const Vec3& X) const
{
  double y = X.y;

  const PoroMaterial* mat = m_mat->getMaterial();

  const Vec4& Xt = static_cast<const Vec4&>(X);

  double k = mat->getPermeability(X)[0];
  double rhof = mat->getFluidDensity(X);
  double E = mat->getStiffness(X);
  double nu = mat->getPoisson(X);
  double Kf = mat->getBulkFluid(X);
  double Ks = mat->getBulkSolid(X);
  double Ko = mat->getBulkMedium(X);
  double n = mat->getPorosity(X);

  double alpha = 1.0 - (Ko/Ks);
  double Minv = ((alpha - n)/Ks) + (n/Kf);
  double K = E / 3 / (1 - 2*nu);
  double G = E / 2 / (1 + nu);
  double mv = 1 / (K + 4*G/3);
  double cv = k / rhof / m_mat->getGAcc() / (alpha * alpha * mv + Minv);

  double stime = cv * Xt.t / height / height;

  int i = 1;
  double ret = 0.0;
  while (true) {
    double trig_arg = (2*i - 1) * M_PI * y / 2 / height;
    double exp_c = exp(-(2*i-1) * (2*i-1) * M_PI * M_PI * stime / 4);
    double sign = i % 2 == 0 ? -1 : 1;
    double p_fac = 4 / (2.0*i - 1) / M_PI * load;

    double p_add = cos(trig_arg) * exp_c * p_fac * sign;
    ret += p_add;

    if (fabs(p_add) < 1e-15)
      break;
    i++;
  }

  return ret;
}


Vec3 StationaryTerzhagiDisplacement::evaluate (const Vec3& X) const
{
  double y = X.y;

  const PoroMaterial* mat = m_mat->getMaterial();
  double E = mat->getStiffness(X);
  double nu = mat->getPoisson(X);

  double c = -load/E * (1+nu) * (1-2*nu) / (1-nu);

  return Vec3(0, c*y, 0);
}

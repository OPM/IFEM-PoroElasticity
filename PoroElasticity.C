// $Id$
//==============================================================================
//!
//! \file PoroElasticity.C
//!
//! \date April 16 2015
//!
//! \author Yared Bekele
//!
//! \brief Integrand implementations for time-dependent poroelasticity problems.
//!
//==============================================================================

#include "PoroElasticity.h"
#include "PoroMaterial.h"
#include "FiniteElement.h"
#include "TimeDomain.h"
#include "Utilities.h"
#include "Tensor.h"
#include "Vec3Oper.h"

typedef std::vector<int> IntVec;  //!< General integer vector


PoroElasticity::PoroElasticity (unsigned short int n, int order) : Elasticity(n)
{
  primsol.resize(1+order); // Current and previous timestep solutions required
  sc = 0.0;
  gacc = 9.81; //kmo: Danger! hard-coded physical property. Why not derive this one from gravity.length() instead ???
}


LocalIntegral* PoroElasticity::getLocalIntegral (const std::vector<size_t>& nen,
                                                 size_t, bool neumann) const
{
  ElmMats* result;
  if (m_mode == SIM::DYNAMIC)
    result = new NewmarkMats<MixedElmMats>(nsd * nen[0], nen[1], neumann, 0.25, 0.5);
  else
    result = new MixedElmMats(nsd * nen[0], nen[1], neumann);
  return result;
}


LocalIntegral* PoroElasticity::getLocalIntegral (size_t nen,
                                                 size_t, bool neumann) const
{
  ElmMats* result;
  if (m_mode == SIM::DYNAMIC)
    result = new NewmarkMats<NonMixedElmMats>(nsd * nen, nen, neumann, 0.25, 0.5);
  else
    result = new NonMixedElmMats(nsd * nen, nen, neumann);
  return result;
}


bool PoroElasticity::initElement (const std::vector<int>& MNPC,
                                  const std::vector<size_t>& elem_sizes,
                                  const std::vector<size_t>& basis_sizes,
                                  LocalIntegral& elmInt)
{
  if (primsol.front().empty()) return true;

  // Extract the element level solution vectors
  elmInt.vec.resize(NSOL);
  std::vector<int>::const_iterator fstart = MNPC.begin() + elem_sizes[0];
  int ierr = utl::gather(IntVec(MNPC.begin(),fstart),nsd,primsol.front(),elmInt.vec[U])
           + utl::gather(IntVec(fstart,MNPC.end()),0,1,primsol.front(),elmInt.vec[P],nsd*basis_sizes[0],basis_sizes[0]);

  if (ierr == 0) return true;

  std::cerr << " *** PoroElasticity::initElement: Detected " << ierr/3
            << " node numbers out of range." << std::endl;

  return false;
}


bool PoroElasticity::initElementBou (const std::vector<int>& MNPC,
                                     const std::vector<size_t>& elem_sizes,
                                     const std::vector<size_t>& basis_sizes,
                                     LocalIntegral& elmInt)
{
  return this->initElement(MNPC,elem_sizes,basis_sizes,elmInt);
}


bool PoroElasticity::initElement (const std::vector<int>& MNPC,
                                  LocalIntegral& elmInt)
{
  if (primsol.empty() || primsol.front().empty())
    return true;

  // Extract the element level solution vectors
  elmInt.vec.resize(NSOL);
  int ierr = 0;
  Matrix temp(nsd+1, MNPC.size());
  ierr += utl::gather(MNPC, nsd+1, primsol.front(), temp);
  Matrix temp2(nsd, MNPC.size());
  for (size_t k = 1; k <= nsd; ++k)
    temp2.fillRow(k, temp.getRow(k).ptr());
  elmInt.vec[U] = temp2;
  elmInt.vec[P] = temp.getRow(nsd+1);

  if (ierr == 0) return true;

  std::cerr << " *** PoroElasticity::initElement: Detected " << ierr/3
            << " node numbers out of range." << std::endl;

  return false;
}


bool PoroElasticity::initElementBou (const std::vector<int>& MNPC,
                                     LocalIntegral& elmInt)
{
  return this->initElement(MNPC,elmInt);
}


bool PoroElasticity::evalIntMx (LocalIntegral& elmInt,
                                const MxFiniteElement& fe,
                                const TimeDomain& time, const Vec3& X) const
{
  return evalInt(elmInt, fe, time, X);
}


bool PoroElasticity::evalStiffnessMatrix(Matrix& mx, const Matrix &B, const Matrix &C, double detJxW) const
{
  Matrix CB;
  CB.multiply(C, B, false, false);
  CB *= detJxW;
  mx.multiply(B, CB, true, false, true);

  return true;
}


bool PoroElasticity::evalMassMatrix(Matrix &mx, const Vector &basis, double rho, double detJxW) const
{
  Matrix temp(basis.size(), basis.size());
  temp.outer_product(basis, basis);
  temp *= rho * detJxW;

  for (size_t i = 0; i < basis.size(); i++)
    for (size_t j = 0; j < basis.size(); j++)
      for (size_t k = 1; k <= nsd; k++)
        mx(i*nsd+k, j*nsd+k) = temp(i+1,j+1);

  return true;
}


bool PoroElasticity::evalCouplingMatrix(Matrix &mx, const Matrix &B, const Vector &basis,
                                        double scl, double alpha, const Vector &m, double detJxW) const
{
  const size_t nstrc = nsd * (nsd + 1) / 2;
  Matrix K(basis.size(), nstrc);

  for (size_t i = 1; i <= basis.size(); i++)
    for (size_t j = 1; j <= nstrc; j++)
      K(i,j) += scl * m(j) * alpha * basis(i) * detJxW;

  mx.multiply(B, K, true, true, true);

  return true;
}


bool PoroElasticity::evalCompressibilityMatrix(Matrix &mx, const Vector& basis,
                                               double scl, double Minv, double detJxW) const
{
  Matrix temp(mx.rows(), mx.cols());
  temp.outer_product(basis, basis);
  temp *= scl * scl * Minv * detJxW;
  mx += temp;

  return true;
}


bool PoroElasticity::evalPermeabilityMatrix(Matrix &mx, const Matrix &grad, double scl,
                                            const Vec3 &permeability, double acc_dens, double detJxW) const
{
  for (size_t i = 1; i <= grad.rows(); i++)
    for (size_t j = 1; j <= grad.rows(); j++)
      for (size_t k = 1; k <= nsd; k++)
        mx(i,j) += scl * scl * permeability[k-1] / acc_dens * grad(i,k) * grad(j,k) * detJxW;

  return true;
}


bool PoroElasticity::evalInt (LocalIntegral& elmInt,
                              const FiniteElement& fe,
                              const TimeDomain& time, const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);
  const PoroMaterial* pmat = dynamic_cast<const PoroMaterial*>(material);
  if (!pmat) {
    std::cerr << __FUNCTION__ << ": No material data." << std::endl;
    return false;
  }

  Matrix Bmat, Cmat;
  if (!this->formBmatrix(Bmat,fe.dNdX))
    return false;

  SymmTensor eps(nsd), sigma(nsd); double U = 0.0;
  if (!material->evaluate(Cmat,sigma,U,fe,X,eps,eps,0))
    return false;

  Vec3 permeability = pmat->getPermeability(X);

  double scl(sc);
  if (scl == 0.0)
    scl = sqrt(pmat->getStiffness(X) * pmat->getFluidDensity(X) * gacc / permeability[0] / time.dt);

  // Biot's coefficient
  double Ko = pmat->getBulkMedium(X);
  double Ks = pmat->getBulkSolid(X);
  double Kw = pmat->getBulkWater(X);
  double poro = pmat->getPorosity(X);

  double alpha = 1.0 - (Ko/Ks);
  // Inverse of the compressibility modulus
  double Minv = ((alpha - poro)/Ks) + (poro/Kw);
  double rho = pmat->getMassDensity(X);

  // Define the unit Voigt vector
  Vector m(Cmat.rows());
  for (size_t i = 1; i <= nsd; i++)
    m(i) = 1.0;

  if (!evalStiffnessMatrix(elMat.A[uu_K], Bmat, Cmat, fe.detJxW))
    return false;
  if (!evalCouplingMatrix(elMat.A[up], Bmat, fe.basis(2), scl, alpha, m, fe.detJxW))
    return false;
  if (!evalCompressibilityMatrix(elMat.A[pp_S], fe.basis(2), scl, Minv, fe.detJxW))
    return false;
  if (!evalPermeabilityMatrix(elMat.A[pp_P], fe.grad(2), scl, permeability,
                              pmat->getFluidDensity(X) * gacc, fe.detJxW))
    return false;
  if (m_mode == SIM::DYNAMIC && !evalMassMatrix(elMat.A[uu_M], fe.basis(1), rho, fe.detJxW))
    return false;

  return true;
}


bool PoroElasticity::evalBouMx (LocalIntegral& elmInt,
                                const MxFiniteElement& fe,
                                const TimeDomain& time, const Vec3& X,
                                const Vec3& normal) const
{
  return evalBou(elmInt, fe, time, X, normal);
}


bool PoroElasticity::evalBou (LocalIntegral& elmInt,
                              const FiniteElement& fe,
                              const TimeDomain& time, const Vec3& X,
                              const Vec3& normal) const
{
  if (!tracFld && !fluxFld)
  {
    std::cerr << " *** PoroElasticity::evalBouMx: No fluxes/tractions." << std::endl;
    return false;
  }

  const PoroMaterial* pmat = dynamic_cast<const PoroMaterial*>(material);
  if (!pmat) {
    std::cerr << __FUNCTION__ << ": No material data." << std::endl;
    return false;
  }

  // Evaluate the surface traction
  Vec4 Xt = static_cast<const Vec4&>(X);
  Xt.t = time.t;
  Vec3 trac = this->getTraction(Xt, normal);
  Vec3 permeability = pmat->getPermeability(X);

  double scl(sc);
  if (scl == 0.0)
    scl = sqrt(pmat->getStiffness(X) * pmat->getFluidDensity(X) * gacc / permeability[0] / time.dt);

  // Integrate the force vector fu
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);
  for (size_t i = 1; i <= fe.basis(1).size(); i++)
    for (unsigned short int j = 1; j <= nsd; j++)
      elMat.b[Fu](nsd*(i-1)+j) += trac[j-1] * fe.basis(1)(i) * fe.detJxW;

  return true;
}


bool PoroElasticity::finalizeElement (LocalIntegral& elmInt, const TimeDomain& time, size_t)
{
  Mats& elMat = static_cast<Mats&>(elmInt);

  // Construct the C-matrix (needed by both)
  elMat.add_pu(up, sys_C);
  elMat.add_pp(pp_S, sys_C);

  if (m_mode != SIM::DYNAMIC) {
    // Construct the system matrix
    elMat.add_uu(uu_K, sys);
    elMat.add_up(up, sys, -1.0);
    elMat.add_pu(up, sys);
    elMat.add_pp(pp_S, sys);
    elMat.add_pp(pp_P, sys, time.dt);

    // Contribution to RHS from previous timestep
    elMat.form_vector(elmInt.vec[U], elmInt.vec[P], Fprev);
    elMat.b[Fprev] = elMat.A[sys_C] * elMat.b[Fprev];

    // Contribution to RHS from current timestep
    elMat.b[Fp] *= time.dt;
    elMat.form_vector(elMat.b[Fu], elMat.b[Fp], Fsys);

    elMat.b[Fsys] += elMat.b[Fprev];
  } else {
    // Construct the M-matrix
    elMat.add_uu(uu_M, sys_M, -1.0);

    // Construct the C-matrix
    elMat.add_pu(up, sys_C);
    elMat.add_pp(pp_S, sys_C);

    // Construct the K-matrix
    elMat.add_uu(uu_K, sys_K);
    elMat.add_up(up, sys_K, -1.0);
    elMat.add_pp(pp_P, sys_K);

    // In case of dynamic mode, we add the zero-order contribution
    // already here, due to (possibly) nonlinear terms
    elMat.form_vector(elMat.b[Fu], elMat.b[Fp], Fsys);
    elMat.form_vector(elmInt.vec[U], elmInt.vec[P], Fprev);
    elMat.b[Fsys] -= elMat.A[sys_K] * elMat.b[Fprev];
  }

  return true;
}


bool PoroElasticity::finalizeElementBou(LocalIntegral& elmInt, const FiniteElement&, const TimeDomain& time)
{
  Mats& elMat = static_cast<Mats&>(elmInt);

  if (m_mode != SIM::DYNAMIC)
    elMat.b[Fp] *= time.dt;
  elMat.form_vector(elMat.b[Fu], elMat.b[Fp], Fsys);

  return true;
}


size_t PoroElasticity::getNoFields (int fld) const
{
  if (fld < 2)
    return nsd+1;
  return nsd * (nsd + 1);
}


std::string PoroElasticity::getField1Name (size_t i, const char* prefix) const
{
  if (i == 11)
    return "Displacements";
  if (i == 12)
    return "Pressure";

  if (i >= nsd)
    i = 3;

  static const char* s[5] = {"u_x", "u_y", "u_z", "p^w"};

  if (!prefix)
    return s[i];
  return prefix + std::string(" ") + s[i];
}


std::string PoroElasticity::getField2Name (size_t i, const char* prefix) const
{
  static const char* s[][6] = {{"x", "y", "xy"},
                               {"x", "y", "z", "yz", "xz", "xy"}};
  size_t ncomps = nsd * (nsd + 1) / 2;

  std::string name = (i < ncomps ? "eps" : "sig") + std::string("_") + s[nsd-2][i % ncomps];
  if (!prefix)
    return name;

  return prefix + std::string(" ") + name;
}


bool PoroElasticity::evalSol(Vector& s, const MxFiniteElement& fe,
                             const Vec3& X, const std::vector<int>& MNPC,
                             const std::vector<size_t>& elem_sizes) const
{
  Vector eV;
  std::vector<int>::const_iterator fstart = MNPC.begin() + elem_sizes[0];
  utl::gather(IntVec(MNPC.begin(),fstart),nsd,primsol.front(),eV);

  return this->evalSolCommon(s,fe,X,Vector(&eV[0],nsd*fe.N.size()));
}


bool PoroElasticity::evalSol(Vector& s, const FiniteElement& fe,
                             const Vec3& X,
                             const std::vector<int>& MNPC) const
{
  Vector eV;
  utl::gather(MNPC,nsd+1,primsol.front(),eV);

  Vector disp(nsd * fe.N.size());
  for (size_t i = 0; i < nsd; i++)
    for (size_t bfun = 0; bfun < fe.N.size(); bfun++)
      disp[nsd*bfun+i] = eV[(nsd+1)*bfun+i];

  return this->evalSolCommon(s,fe,X,disp);
}


bool PoroElasticity::evalSolCommon (Vector& s,
                                    const FiniteElement& fe, const Vec3& X,
                                    const Vector& disp) const
{
  if (!material)
  {
    std::cerr << __FUNCTION__ <<": No material data."<< std::endl;
    return false;
  }

  Matrix Bmat;
  if (!this->formBmatrix(Bmat,fe.dNdX))
    return false;

  SymmTensor eps(nsd), sigma(nsd);
  if (!Bmat.multiply(disp,eps))
    return false;

  Matrix Cmat; double U = 0.0;
  if (!material->evaluate(Cmat,sigma,U,fe,X,eps,eps))
    return false;

  s = eps;
  const RealArray& sig = sigma;
  s.insert(s.end(),sig.begin(),sig.end());

  return true;
}

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


//! \brief Enum for element level solution vectors
enum SolutionVectors
{
  U = 0,                        // Displacement
  P = 1,                        // Pore pressure
  NSOL = 2
};


//! \brief Enum for element level right-hand-side vectors
enum ResidualVectors
{
  // System vectors, "global" size
  Fsys = 0,                     // Final RHS vector
  Fcur = 1,                     // RHS contribution from this timestep
  Fprev = 2,                    // RHS contribution from previous timestep

  // Sub-vectors, sized according to the bases in question
  Fu = 3,                       // Traction and body forces
  Fp = 4,                       // Flux and body forces

  Fres = 5,
  NVEC = 6
};


//! \brief Enum for element level left-hand-side matrices
enum TangentMatrices
{
  // System matrices, "global" size (note: order is fixed according to Newmark API)
  sys = 0,                      // Final newton matrix
  sys_mass = 1,                 // Mass (acceleration term in Newmark)
  sys_gstiff = 2,               // Geometric stiffness (velocity term in Newmark)
  sys_mstiff = 3,               // Material stiffness (zero-order term in Newmark)

  // Sub-matrices, sized according to the bases in question
  uu_K = 4,                     // Stiffness matrix
  uu_M = 5,                     // Mass matrix
  up = 6,                       // Coupling matrix
  pp_S = 7,                     // Compressibility matrix
  pp_P = 8,                     // Permeability matrix

  NMAT = 9
};


PoroElasticity::Mats::Mats(size_t ndof_displ, size_t ndof_press, bool neumann)
{
  resize(NMAT, NVEC);
  this->ndof_displ = ndof_displ;
  this->ndof_press = ndof_press;

  size_t ndof_tot = ndof_displ + ndof_press;

  rhsOnly = neumann;
  withLHS = !neumann;
  b[Fsys].resize(ndof_tot);
  b[Fprev].resize(ndof_tot);
  b[Fu].resize(ndof_displ);
  b[Fp].resize(ndof_press);
  b[Fres].resize(ndof_tot);

  if (!neumann)
  {
    A[sys].resize(ndof_tot, ndof_tot);
    A[sys_mass].resize(ndof_tot, ndof_tot);
    A[sys_gstiff].resize(ndof_tot, ndof_tot);
    A[sys_mstiff].resize(ndof_tot, ndof_tot);
    A[uu_K].resize(ndof_displ, ndof_displ);
    A[uu_M].resize(ndof_displ, ndof_displ);
    A[up].resize(ndof_displ, ndof_press);
    A[pp_S].resize(ndof_press, ndof_press);
    A[pp_P].resize(ndof_press, ndof_press);
  }
}

const Matrix& PoroElasticity::Mats::getNewtonMatrix () const
{
  return A[sys];
}


const Vector& PoroElasticity::Mats::getRHSVector () const
{
  return b[Fsys];
}


void PoroElasticity::MixedElmMats::add_uu(size_t source, size_t target, double scale)
{
  A[target].addBlock(A[source], scale, 1, 1);
}


void PoroElasticity::MixedElmMats::add_up(size_t source, size_t target, double scale)
{
  A[target].addBlock(A[source], scale, 1, 1 + ndof_displ);
}


void PoroElasticity::MixedElmMats::add_pu(size_t source, size_t target, double scale)
{
  A[target].addBlock(A[source], scale, 1 + ndof_displ, 1, true);
}


void PoroElasticity::MixedElmMats::add_pp(size_t source, size_t target, double scale)
{
  A[target].addBlock(A[source], scale, 1 + ndof_displ, 1 + ndof_displ);
}


void PoroElasticity::MixedElmMats::form_vector(const Vector &u, const Vector &p, size_t target)
{
  b[target] = u;
  b[target].insert(b[target].end(), p.begin(), p.end());
}


void PoroElasticity::NonMixedElmMats::add_uu(size_t source, size_t target, double scale)
{
  Matrix &T = A[target];
  Matrix &S = A[source];
  size_t nsd = ndof_displ / ndof_press;
  size_t nf = nsd + 1;

  for (size_t i = 1; i <= ndof_press; ++i)
    for (size_t j = 1; j <= ndof_press; ++j)
      for (size_t l = 1; l <= nsd; ++l)
        for (size_t k = 1; k <= nsd; ++k)
          T(nf*(i-1)+l, nf*(j-1)+k) += scale * S(nsd*(i-1)+l, nsd*(j-1)+k);
}


void PoroElasticity::NonMixedElmMats::add_up(size_t source, size_t target, double scale)
{
  Matrix &T = A[target];
  Matrix &S = A[source];
  size_t nsd = ndof_displ / ndof_press;
  size_t nf = nsd + 1;

  for (size_t i = 1; i <= ndof_press; ++i)
    for (size_t j = 1; j <= ndof_press; ++j)
      for (size_t l = 1; l <= nsd; ++l)
        T(nf*(i-1)+l, j*nf) += scale * S(nsd*(i-1)+l, j);
}


void PoroElasticity::NonMixedElmMats::add_pu(size_t source, size_t target, double scale)
{
  Matrix &T = A[target];
  Matrix &S = A[source];
  size_t nsd = ndof_displ / ndof_press;
  size_t nf = nsd + 1;

  for (size_t i = 1; i <= ndof_press; ++i)
    for (size_t j = 1; j <= ndof_press; ++j)
      for (size_t l = 1; l <= nsd; ++l)
        T(j*nf, nf*(i-1)+l) += scale * S(nsd*(i-1)+l, j);
}


void PoroElasticity::NonMixedElmMats::add_pp(size_t source, size_t target, double scale)
{
  Matrix &T = A[target];
  Matrix &S = A[source];
  size_t nf = ndof_displ / ndof_press + 1;

  for (size_t i = 1; i <= ndof_press; ++i)
    for (size_t j = 1; j <= ndof_press; ++j)
      T(i*nf, j*nf) += scale * S(i, j);
}


void PoroElasticity::NonMixedElmMats::form_vector(const Vector &u, const Vector &p, size_t target)
{
  utl::interleave(u, p, b[target], ndof_displ / ndof_press, 1);
}


PoroElasticity::PoroElasticity (unsigned short int n, int order) : Elasticity(n)
{
  primsol.resize(1+order); // Current and previous timestep solutions required
  sc = 0.0;
  gacc = 9.81; //kmo: Danger! hard-coded physical property. Why not derive this one from gravity.length() instead ???
}


LocalIntegral* PoroElasticity::getLocalIntegral (const std::vector<size_t>& nen,
                                                 size_t, bool neumann) const
{
  ElmMats* result = new MixedElmMats(nsd * nen[0], nen[1], neumann);
  return result;
}


LocalIntegral* PoroElasticity::getLocalIntegral (size_t nen,
                                                 size_t, bool neumann) const
{
  ElmMats* result = new NonMixedElmMats(nsd * nen, nen, neumann);
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

  // Construct the system matrix
  elMat.add_uu(uu_K, sys);
  elMat.add_up(up, sys, -1.0);
  elMat.add_pu(up, sys);
  elMat.add_pp(pp_S, sys);
  elMat.add_pp(pp_P, sys, time.dt);

  // Construct the geometric stiffness matrix
  elMat.add_pu(up, sys_gstiff);
  elMat.add_pp(pp_S, sys_gstiff);

  // Contribution to RHS from previous timestep
  elMat.form_vector(elmInt.vec[U], elmInt.vec[P], Fprev);
  elMat.b[Fprev] = elMat.A[sys_gstiff] * elMat.b[Fprev];

  // Contribution to RHS from current timestep
  elMat.b[Fp] *= time.dt;
  elMat.form_vector(elMat.b[Fu], elMat.b[Fp], Fcur);

  elMat.b[Fsys] = elMat.b[Fcur] + elMat.b[Fprev];

  return true;
}


bool PoroElasticity::finalizeElementBou(LocalIntegral& elmInt, const FiniteElement&, const TimeDomain& time)
{
  Mats& elMat = static_cast<Mats&>(elmInt);

  // Contribution to RHS from current timestep
  elMat.b[Fp] *= time.dt;
  elMat.form_vector(elMat.b[Fu], elMat.b[Fp], Fcur);

  elMat.b[Fsys] = elMat.b[Fcur] + elMat.b[Fprev];

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

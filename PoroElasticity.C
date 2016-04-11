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
  U = 0, // displacement
  P = 1, // pore pressure
  NSOL = 2
};


//! \brief Enum for element level right-hand-side vectors
enum ResidualVectors
{
  Fu = 0,
  Fp = 1,
  Fres = 2,
  Fprev = 3,
  NVEC = 4
};


//! \brief Enum for element level left-hand-side matrices
enum TangentMatrices
{
  uu = 0,                       // Stiffness matrix
  up = 1,                       // Coupling matrix
  pp_S = 2,                     // Compressibility matrix
  pp_P = 3,                     // Permeability matrix
  pp = 4,                       // S + dt P
  Ktan = 5,
  Kprev = 6,
  NMAT = 7
};


PoroElasticity::MixedElmMats::MixedElmMats ()
{
  this->resize(NMAT,NVEC); // Number of element matrices and vectors
}


void PoroElasticity::MixedElmMats::makeNewtonMatrix_P(Matrix& N, size_t pp_idx) const
{
  size_t n = A[uu].rows();
  N.fillBlock(A[up], 1+n, 1, true);
  N.fillBlock(A[pp_idx], 1+n, 1+n);
}


void PoroElasticity::MixedElmMats::makeNewtonMatrix_U(Matrix& N) const
{
  size_t n = A[uu].rows();
  N.fillBlock(A[uu], 1, 1);
  N.addBlock(A[up], -1.0, 1, 1+n);
}


const Matrix& PoroElasticity::MixedElmMats::getNewtonMatrix () const
{
  Matrix& N = const_cast<Matrix&>(A[Ktan]);

  makeNewtonMatrix_U(N);
  makeNewtonMatrix_P(N, pp);

  return A[Ktan];
}


const Vector& PoroElasticity::MixedElmMats::getRHSVector () const
{
  Vector& F = const_cast<Vector&>(b[Fres]);
  size_t n = b[Fu].size();

  std::copy(b[Fu].begin(), b[Fu].end(), F.begin());
  std::copy(b[Fp].begin(), b[Fp].end(), F.begin()+n);

  F += b[Fprev];

  return b[Fres];
}


PoroElasticity::NonMixedElmMats::NonMixedElmMats ()
{
  this->resize(NMAT,NVEC); // Number of element matrices and vectors
}


void PoroElasticity::NonMixedElmMats::makeNewtonMatrix_P(Matrix& N, size_t pp_idx) const
{
  size_t n = A[pp].rows();
  size_t nsd = A[uu].rows()/A[pp].rows();
  size_t nf = nsd + 1;

  for (size_t i = 1; i <= n; ++i) {
    for (size_t j = 1; j <= n; ++j) {
      // PP-coupling
      N(nf*i, nf*j) = A[pp_idx](i,j);
      for (size_t l = 1; l <= nsd; ++l) {
        // PU-coupling
        N(j*nf, nf*(i-1)+l) = A[up](nsd*(i-1)+l, j);
      }
    }
  }
}


void PoroElasticity::NonMixedElmMats::makeNewtonMatrix_U(Matrix& N) const
{
  size_t n = A[pp].rows();
  size_t nsd = A[uu].rows()/A[pp].rows();
  size_t nf = nsd + 1;

  for (size_t i = 1; i <= n; ++i) {
    for (size_t j = 1; j <= n; ++j) {
      for (size_t l = 1; l <= nsd; ++l) {
        // UP-coupling
        N(nf*(i-1)+l, j*nf) = -A[up](nsd*(i-1)+l, j);
        for (size_t k = 1; k <= nsd; ++k) {
          // UU-coupling
          N(nf*(i-1)+l, (j-1)*nf+k) = A[uu](nsd*(i-1)+l, (j-1)*nsd+k);
        }
      }
    }
  }
}


const Matrix& PoroElasticity::NonMixedElmMats::getNewtonMatrix () const
{
  Matrix& N = const_cast<Matrix&>(A[Ktan]);

  makeNewtonMatrix_U(N);
  makeNewtonMatrix_P(N, pp);

  return A[Ktan];
}


const Vector& PoroElasticity::NonMixedElmMats::getRHSVector () const
{
  Vector& F = const_cast<Vector&>(b[Fres]);

  size_t nsd = b[Fu].size()/b[Fp].size();
  utl::interleave(b[Fu], b[Fp], F, nsd, 1);

  F += b[Fprev];

  return b[Fres];
}


PoroElasticity::PoroElasticity (unsigned short int n, int order) : Elasticity(n)
{
  primsol.resize(1+order); // Current and previous timestep solutions required
  sc = 0.0;
  gacc = 9.81; //kmo: Danger! hard-coded physical property. Why not derive this one from gravity.length() instead ???
}


void PoroElasticity::initLocalIntegral(ElmMats *result, size_t ndof_displ,
                                       size_t ndof_press, bool neumann) const
{
  size_t ndof_tot = ndof_displ + ndof_press;

  result->rhsOnly = neumann;
  result->withLHS = !neumann;
  result->b[Fres].resize(ndof_tot);
  result->b[Fprev].resize(ndof_tot);
  result->b[Fu].resize(ndof_displ);
  result->b[Fp].resize( ndof_press);

  if (!neumann)
  {
    result->A[uu].resize(ndof_displ, ndof_displ);
    result->A[up].resize(ndof_displ, ndof_press);
    result->A[pp_S].resize(ndof_press, ndof_press);
    result->A[pp_P].resize(ndof_press, ndof_press);
    result->A[pp].resize(ndof_press, ndof_press);
    result->A[Ktan].resize(ndof_tot, ndof_tot);
    result->A[Kprev].resize(ndof_tot, ndof_tot);
  }
}


LocalIntegral* PoroElasticity::getLocalIntegral (const std::vector<size_t>& nen,
                                                 size_t, bool neumann) const
{
  ElmMats* result = new MixedElmMats();
  initLocalIntegral(result, nsd * nen[0], nen[1], neumann);
  return result;
}


LocalIntegral* PoroElasticity::getLocalIntegral (size_t nen,
                                                 size_t, bool neumann) const
{
  ElmMats* result = new NonMixedElmMats();
  initLocalIntegral(result, nsd * nen, nen, neumann);
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

  if (!evalStiffnessMatrix(elMat.A[uu], Bmat, Cmat, fe.detJxW))
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
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  elMat.A[pp].add(elMat.A[pp_S]);
  elMat.A[pp].add(elMat.A[pp_P], time.dt);

  Vector prevSol;

  MixedElmMats* mMat = dynamic_cast<MixedElmMats*>(&elmInt);
  if (mMat) {
    mMat->makeNewtonMatrix_P(elMat.A[Kprev], pp_S);
    prevSol = elmInt.vec[U];
    prevSol.insert(prevSol.end(), elmInt.vec[P].begin(), elmInt.vec[P].end());
  } else {
    NonMixedElmMats* mMat = dynamic_cast<NonMixedElmMats*>(&elmInt);
    mMat->makeNewtonMatrix_P(elMat.A[Kprev], pp_S);
    utl::interleave(elmInt.vec[U], elmInt.vec[P], prevSol, nsd, 1);
  }

  elMat.b[Fprev] = elMat.A[Kprev] * prevSol;

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

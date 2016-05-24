// $Id$
//==============================================================================
//!
//! \file ThermoPoroElasticity.C
//!
//! \date
//!
//! \author Yared Bekele
//!
//! \brief Integrand implementations for non-isothermal PoroElasticity problems
//!
//==============================================================================

#include "ThermoPoroElasticity.h"
#include "ASMbase.h"
#include "ASMmxBase.h"
#include "FiniteElement.h"
#include "TimeDomain.h"
#include "Utilities.h"
#include "Tensor.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "Vec3Oper.h"
#include "VTF.h"
#include "StabilizationUtils.h"

//! \brief Enum for element level solution vectors
enum SolutionVectors
{
  Uo = 0,         // Previous displacement
  Po = 1,         // Previous pore pressure
  To = 2,         // Previous temperature
  Uc = 3,         // Current displacement
  Pc = 4,         // Current pore pressure
  Tc = 5,         // Current temperature
  NSOL = 6
};


//! \brief Enum for element level right-hand-side vectors
enum ResidualVectors
{
  Fu = 0,
  Fp = 1,
  FT = 2,
  Fprev = 3,
  Fnow = 4,
  Fres = 5,
  Fc = 6,
  NVEC = 7
};


//! \brief Enum for element level left-hand-side matrices
enum TangentMatrices
{
  uu = 0,
  up = 1,
  uT = 2,
  pp = 3,
  TT = 4,
  Kprev = 5,
  Know = 6,
  Ktan = 7,
  NMAT = 8
};


ThermoPoroElasticity::MixedElmMats::MixedElmMats()
{
  this->resize(NMAT,NVEC);
}


const Matrix& ThermoPoroElasticity::MixedElmMats::getNewtonMatrix() const
{
  Matrix& N = const_cast<Matrix&>(A[Ktan]);

  size_t i,j;
  size_t ru = A[uu].rows();
  size_t rp = A[pp].rows();

  for (i = 1; i <= ru; i++)
  {
    for (j = 1; j <= ru; j++)
    {
      N(i,j) = A[uu](i,j);
    }
    for (j = 1; j <= rp; j++)
    {
      size_t k = ru+2*j-1;
      N(i,k) = A[up](i,j);
      N(k,i) = A[up](i,j);
      size_t l = ru+2*j;
      N(i,l) = A[uT](i,j);
      //N(l,i) = A[uT](i,j);
    }
  }

  for (i = 1; i <= rp; i++)
  {
    for (j = 1; j <= rp; j++)
    {
      size_t ki = ru+2*i-1;
      size_t kj = ru+2*j-1;
      N(ki,kj) = A[pp](i,j);
      size_t li = ru+2*i;
      size_t lj = ru+2*j;
      N(li,lj) = A[TT](i,j);
    }
  }

  return A[Ktan];
}


const Vector& ThermoPoroElasticity::MixedElmMats::getRHSVector() const
{
  Vector& F = const_cast<Vector&>(b[Fres]);

  size_t ru = b[Fu].size();
  size_t rp = b[Fp].size();

  for (size_t i = 1; i <= ru; i++)
    F(i) = b[Fu](i);

  for (size_t i = 1; i <= rp; i++)
  {
    F(ru+2*i-1) = b[Fp](i);
    F(ru+2*i  ) = b[FT](i);
  }

  F += b[Fprev];
  F -= b[Fnow];
  // Robin boundary contribution to the residual
  F -= A[Know]*b[Fc];

  return b[Fres];
}


ThermoPoroElasticity::ThermoPoroElasticity(unsigned short int n, int order, bool stab) :
  Elasticity(n), gacc(9.81), mat(nullptr), SUPG(stab)
{
  primsol.resize(1+order);
  tracFld = nullptr;
  fluxFld = nullptr;
  eS = 1;
}


Vec3 ThermoPoroElasticity::getTraction(const Vec3& X, const Vec3& n) const
{
  if (fluxFld)
    return (*fluxFld)(X);
  else if (tracFld)
    return (*tracFld)(X,n);
  else
    return Vec3();
}


LocalIntegral* ThermoPoroElasticity::getLocalIntegral(const std::vector<size_t>& nen,
                                                      size_t, bool neumann) const
{
  const size_t nedof1 = nsd*nen[0];
  const size_t nedof = nedof1 + 2*nen[1];

  ElmMats* result = new MixedElmMats();

  result->rhsOnly = neumann;
  result->withLHS = !neumann;
  result->b[Fu].resize(nedof1);
  result->b[Fp].resize(nen[1]);
  result->b[FT].resize(nen[1]);
  result->b[Fprev].resize(nedof);
  result->b[Fnow].resize(nedof);
  result->b[Fres].resize(nedof);
  result->b[Fc].resize(nedof);

  /*if(!neumann)
    {*/
  result->A[uu].resize(nedof1,nedof1);
  result->A[up].resize(nedof1,nen[1]);
  result->A[uT].resize(nedof1,nen[1]);
  result->A[pp].resize(nen[1],nen[1]);
  result->A[TT].resize(nen[1],nen[1]);
  result->A[Kprev].resize(nedof,nedof);
  result->A[Know].resize(nedof,nedof);
  result->A[Ktan].resize(nedof,nedof);
  //}

  return result;
}


bool ThermoPoroElasticity::initElement(const std::vector<int>& MNPC,
                                       const std::vector<size_t>& elem_sizes,
                                       const std::vector<size_t>& basis_sizes,
                                       LocalIntegral& elmInt)
{
  if (primsol.front().empty()) return true;

  // Extract the element level solution vectors
  elmInt.vec.resize(NSOL);
  std::vector<int>::const_iterator fstart = MNPC.begin() + elem_sizes[0];
  int ierr = utl::gather(IntVec(MNPC.begin(), fstart),nsd,primsol[0],elmInt.vec[Uc])
    + utl::gather(IntVec(fstart,MNPC.end()),0,2,primsol[0],elmInt.vec[Pc],nsd*basis_sizes[0],basis_sizes[0])|
    + utl::gather(IntVec(fstart,MNPC.end()),1,2,primsol[0],elmInt.vec[Tc],nsd*basis_sizes[0],basis_sizes[0])
    + utl::gather(IntVec(MNPC.begin(),fstart),nsd,primsol[1],elmInt.vec[Uo])
    + utl::gather(IntVec(fstart,MNPC.end()),0,2,primsol[1],elmInt.vec[Po],nsd*basis_sizes[0],basis_sizes[0]);
  + utl::gather(IntVec(fstart,MNPC.end()),1,2,primsol[1],elmInt.vec[To],nsd*basis_sizes[0],basis_sizes[0]);

  if (ierr == 0) return true;

  std::cerr << " *** ThermoPoroElasticity::initElement: Detected " << ierr/3
    << " node numbers out of range." << std::endl;

  return false;
}


bool ThermoPoroElasticity::initElementBou(const std::vector<int>& MNPC,
                                          const std::vector<size_t>& elem_sizes,
                                          const std::vector<size_t>& basis_sizes,
                                          LocalIntegral& elmInt)
{
  return this->IntegrandBase::initElementBou(MNPC,elem_sizes,basis_sizes,elmInt);
}


bool ThermoPoroElasticity::evalIntMx(LocalIntegral& elmInt,
                                     const MxFiniteElement& fe,
                                     const TimeDomain& time, const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  if (!mat)
  {
    std::cerr << __FUNCTION__ << ":No material data." << std::endl;
    return false;
  }

  size_t i,j,k;

  Matrix Bmat, Cmat, CB;
  if(!this->formBmatrix(Bmat,fe.grad(1)))
    return false;

  SymmTensor eps(nsd), sigma(nsd);
  double U = 0.0;
  if (!mat->evaluate(Cmat,sigma,U,fe,X,eps,eps,0))
    return false;

  Vec3 hydCond = mat->getPermeability(X);
  double rhof = mat->getFluidDensity(X);
  double Ko = mat->getBulkMedium(X);
  double Ks = mat->getBulkSolid(X);
  double Kf = mat->getBulkFluid(X);
  double poro = mat->getPorosity(X);
  // Biot's coefficient
  double alpha = 1.0 - (Ko/Ks);
  // Inverse of the compressibility modulus
  double Minv = ((alpha - poro)/Ks) + (poro/Kf);
  // m vector for 2D
  Vec3 m(1,1,0);

  // Integration of Cuu
  CB.multiply(Cmat,Bmat,false,false);
  CB *= -1.0 * fe.detJxW;
  elMat.A[uu].multiply(Bmat,CB,true,false,true);

  // Integration of Cup
  Matrix Cuptmp;
  const size_t nstrc = nsd*(nsd+1)/2;
  Cuptmp.resize(nstrc,fe.basis(2).size());
  for (i = 1; i <= nstrc; i++)
    for (j = 1; j <= fe.basis(2).size(); j++)
      Cuptmp(i,j) += scl1*m[i-1]*alpha*fe.basis(2)(j)*fe.detJxW;

  elMat.A[up].multiply(Bmat,Cuptmp,true,false,true);

  double T = elMat.vec[Tc].dot(fe.basis(2));
  double rhoc = mat->getHeatCapacity(T);
  double lambda = mat->getThermalConductivity(T);

  // Integration of CuT
  Matrix BC, CuTtmp;
  BC.multiply(Bmat,Cmat,true,false);
  CuTtmp.resize(nstrc,fe.basis(2).size());
  double alpha_s = mat->getSolidThermalExpansion(T);
  for (i = 1; i <= nstrc; i++)
    for (j = 1; j <= fe.basis(2).size(); j++)
      CuTtmp(i,j) += scl2*m[i-1]*(alpha_s/3.0)*fe.basis(2)(j)*fe.detJxW;

  elMat.A[uT].multiply(BC,CuTtmp,false,false,true);

  // Integration of Cpp
  Matrix Cpp;
  Cpp.resize(fe.basis(2).size(),fe.basis(2).size());
  for (i = 1; i <= fe.basis(2).size(); i++)
    for (j = 1; j <= fe.basis(2).size(); j++)
      Cpp(i,j) += scl1*scl1*fe.basis(2)(i)*Minv*fe.basis(2)(j)*fe.detJxW;

  // Integration of Kpp
  Matrix Kpp;
  Kpp.resize(fe.basis(2).size(),fe.basis(2).size());
  for (i = 1; i <= fe.basis(2).size(); i++)
    for (j = 1; j <= fe.basis(2).size(); j++)
      for (k = 1; k <= nsd; k++)
        Kpp(i,j) += scl1*scl1*fe.grad(2)(i,k)*(hydCond[k-1]/(rhof*gacc))*fe.grad(2)(j,k)*fe.detJxW;

  elMat.A[pp] += Cpp;
  elMat.A[pp].add(Kpp,time.dt);

  // Integration of CTT;
  Matrix CTT;
  CTT.resize(fe.basis(2).size(),fe.basis(2).size());
  for (i = 1; i <= fe.basis(2).size(); i++)
    for (j = 1; j <= fe.basis(2).size(); j++)
      CTT(i,j) += scl2*scl2*fe.basis(2)(i)*rhoc*fe.basis(2)(j)*fe.detJxW;

  // Integration of KTT
  double rhofcf = rhof * mat->getFluidHeatCapacity(T);
  Vector gradP;
  fe.grad(2).multiply(elMat.vec[Pc],gradP,true);
  Matrix KTT;
  KTT.resize(fe.basis(2).size(),fe.basis(2).size());

  Vec3 vel;
  for (k = 1; k <= nsd; k++)
    vel[k-1] = -1.0*gradP(k)*hydCond[k-1]/(rhof*gacc);

  for (i = 1; i <= fe.basis(2).size(); i++) {
    for (j = 1; j <= fe.basis(2).size(); j++) {
      double laplace = 0.0, convection = 0.0;
      for (k = 1; k <= nsd; k++) {
        laplace += fe.grad(2)(i,k)*fe.grad(2)(j,k);
        convection += fe.basis(2)(i)*rhofcf*vel[k-1]*fe.grad(2)(j,k);
      }
      KTT(i,j) += scl2*scl2*(laplace*lambda + convection)*fe.detJxW;
    }
  }

  elMat.A[TT] += CTT;
  elMat.A[TT].add(KTT,time.dt);

  size_t ru = elMat.A[uu].rows();
  size_t rp = elMat.A[pp].rows();

  for (i = 1; i <= rp; i++) {
    for (j = 1; j <= rp; j++) {
      size_t ki = ru+2*i-1;
      size_t kj = ru+2*j-1;
      elMat.A[Kprev](ki,kj) += Cpp(i,j);
      size_t li = ru+2*i;
      size_t lj = ru+2*j;
      elMat.A[Kprev](li,lj) += CTT(i,j);
    }
  }

  if (SUPG)
  {
    Matrix KTTs, CTTs;
    KTTs.resize(fe.basis(2).size(),fe.basis(2).size());
    CTTs.resize(fe.basis(2).size(),fe.basis(2).size());
    double h = StabilizationUtils::getElementSize(fe.XC,nsd);
    double alpha_e = h*fabs(nsd*vel[0])/(2*lambda);
    double taue = h/(2*fabs(nsd*vel[0])) * (1/tanh(alpha_e) - 1/alpha_e);
    std::cout << "h = " << h << ", alpha_e = " << alpha_e << ", taue = " << taue << std::endl;
    for (i = 1; i <= fe.basis(2).size(); i++) {
      for (j = 1; j <= fe.basis(2).size(); j++) {
        double a = 0.0, b = 0.0, c = 0.0;
        for (k = 1; k <= nsd; k++) {
          a += fe.grad(2)(i,k)*vel[k-1]*vel[k-1]*fe.grad(2)(j,k);
          b += fe.grad(2)(i,k)*vel[k-1]*fe.hess(2)(i,k,k);
          c += fe.grad(2)(i,k)*vel[k-1]*fe.basis(2)(j);
        }
        KTTs(i,j) += taue * rhoc * (rhoc * a + lambda * b) * fe.detJxW;
        CTTs(i,j) += taue * rhoc * c * fe.detJxW;
      }
    }
    elMat.A[TT] += CTTs;
    elMat.A[TT].add(KTTs,time.dt);
    for (i = 1; i <= rp; i++) {
      for (i = 1; i <= rp; i++) {
        size_t li = ru+2*i;
        size_t lj = ru+2*j;
        elMat.A[Kprev](li,lj) += CTTs(i,j);
      }
    }
  }

  return true;
}


bool ThermoPoroElasticity::evalBouMx(LocalIntegral& elmInt,
                                     const MxFiniteElement& fe,
                                     const TimeDomain& time,
                                     const Vec3& X, const Vec3& normal) const
{
  if (!tracFld && !fluxFld)
  {
    std::cerr << " *** ThermoPoroElasticity::evalBouMx: No fluxes/tractions." << std::endl;
    return false;
  }
  else if (!eS)
  {
    std::cerr << " *** ThermoPoroElasticity::evalBouMx: No load vector." << std::endl;
    return false;
  }

  // Evaluate the surface traction
  Vec4 Xt = static_cast<const Vec4&>(X);
  Xt.t = time.t;
  Vec3 tr2 = this->getTraction(Xt,normal);
  Xt.t -= time.dt;
  Vec3 tr1 = this->getTraction(Xt,normal);
  Vec3 dtr;
  dtr = tr2 - tr1;

  // Integration of vector Fu
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);
  for (size_t i = 1; i <= fe.basis(1).size(); i++)
    for (unsigned short int j = 1; j <= nsd; j++)
      elMat.b[Fu](nsd*(i-1)+j) += -1.0*dtr[j-1]*fe.basis(1)(i)*fe.detJxW;

  // Integration of vector Fp
  for (size_t i = 1; i <= fe.basis(2).size(); i++)
    for (size_t k = 1; k <= nsd; k++)
      elMat.b[Fp](i) += 0.0*scl1; // TO DO

  // Integration of vector FT
  for (size_t i = 1; i <= fe.basis(2).size(); i++)
    for (size_t k = 1; k <= nsd; k++)
      elMat.b[FT](i) += 0.0*scl2; // TO DO

  return true;
}


bool ThermoPoroElasticity::finalizeElement(LocalIntegral& elmInt,
                                           const TimeDomain&, size_t)
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  size_t i,j;

  size_t ru = elMat.A[uu].rows();
  size_t rp = elMat.A[pp].rows();

  for (i = 1; i <= ru; i++)
  {
    for (j = 1; j <= ru; j++) {
      elMat.A[Kprev](i,j) = elMat.A[uu](i,j);
      elMat.A[Know](i,j)  = elMat.A[uu](i,j);
    }
    for (j = 1; j <= rp; j++)
    {
      size_t k = ru+2*j-1;
      elMat.A[Kprev](i,k) = elMat.A[up](i,j);
      elMat.A[Know](i,k)  = elMat.A[up](i,j);
      elMat.A[Kprev](k,i) = elMat.A[up](i,j);
      elMat.A[Know](k,i)  = elMat.A[up](i,j);
      size_t l = ru+2*j;
      elMat.A[Kprev](i,l) = elMat.A[uT](i,j);
      elMat.A[Know](i,l)  = elMat.A[uT](i,j);
      //elMat.A[Kprev](l,i) = elMat.A[uT](i,j);
      //elMat.A[Know](l,i)  = elMat.A[uT](i,j);
    }
  }

  for (i = 1; i <= rp; i++) {
    for (j = 1; j <= rp; j++) {
      size_t ki = ru+2*i-1;
      size_t kj = ru+2*j-1;
      elMat.A[Know](ki,kj)  = elMat.A[pp](i,j);
      size_t li = ru+2*i;
      size_t lj = ru+2*j;
      elMat.A[Know](li,lj)  = elMat.A[TT](i,j);
    }
  }

  Vector PoTo, PcTc, prevSol, currSol;
  PoTo.resize(2*elMat.vec[Po].size());
  PcTc.resize(2*elMat.vec[Pc].size());
  for (i = 1; i <= elMat.vec[Po].size(); i++) {
    PoTo(2*i-1) = elMat.vec[Po](i);
    PoTo(2*i  ) = elMat.vec[To](i);
    PcTc(2*i-1) = elMat.vec[Pc](i);
    PcTc(2*i  ) = elMat.vec[Tc](i);
  }

  prevSol = elMat.vec[Uo];
  prevSol.insert(prevSol.end(),PoTo.begin(),PoTo.end());

  currSol = elMat.vec[Uc];
  currSol.insert(currSol.end(),PcTc.begin(),PcTc.end());

  elMat.b[Fprev] = elMat.A[Kprev]*prevSol;
  elMat.b[Fnow]  = elMat.A[Know]*currSol;

  return true;
}

bool ThermoPoroElasticity::evalSol(Vector& s, const MxFiniteElement& fe,
                                   const Vec3& X, const std::vector<int>& MNPC1,
                                   const std::vector<size_t>& elem_sizes) const
{
  return true;
}

size_t ThermoPoroElasticity::getNoFields(int fld) const
{
  if (fld < 2)
    return nsd+2;
  else
    return nsd+5;
}


std::string ThermoPoroElasticity::getField1Name(size_t i, const char* prefix) const
{
  static const char* s[4] = { "u_x", "u_y", "p^w", "T" };
  if (i == 11)
    return "first_basis";
  if (i == 12)
    return "second_basis";

  if(!prefix) return s[i];

  return prefix + std::string(" ") + s[i];
}


std::string ThermoPoroElasticity::getField2Name (size_t i, const char* prefix) const
{
  static const char* s2[] = {"eps_x","eps_y","eps_xy","sig_x",
    "sig_y","sig_z","sig_xy"};
  if (!prefix) return s2[i];

  return prefix + std::string(" ") + s2[i];
}


ThermoPoroElasticity::WeakDirichlet::WeakDirichlet(unsigned short int n) :
  nsd(n), flux(nullptr)
{
  primsol.resize(2);
}


bool ThermoPoroElasticity::WeakDirichlet::initElementBou(const std::vector<int>& MNPC,
                                                         const std::vector<size_t>& elem_sizes,
                                                         const std::vector<size_t>& basis_sizes,
                                                         LocalIntegral& elmInt)
{
  if (primsol.front().empty()) return true;

  // Extract the element level solution vectors
  elmInt.vec.resize(NSOL);
  std::vector<int>::const_iterator fstart = MNPC.begin() + elem_sizes[0];
  int ierr = utl::gather(IntVec(MNPC.begin(), fstart),nsd,primsol[0],elmInt.vec[Uc])
    + utl::gather(IntVec(fstart,MNPC.end()),0,2,primsol[0],elmInt.vec[Pc],nsd*basis_sizes[0],basis_sizes[0])|
    + utl::gather(IntVec(fstart,MNPC.end()),1,2,primsol[0],elmInt.vec[Tc],nsd*basis_sizes[0],basis_sizes[0])
    + utl::gather(IntVec(MNPC.begin(),fstart),nsd,primsol[1],elmInt.vec[Uo])
    + utl::gather(IntVec(fstart,MNPC.end()),0,2,primsol[1],elmInt.vec[Po],nsd*basis_sizes[0],basis_sizes[0]);
  + utl::gather(IntVec(fstart,MNPC.end()),1,2,primsol[1],elmInt.vec[To],nsd*basis_sizes[0],basis_sizes[0]);

  if (ierr == 0) return true;

  std::cerr << " *** ThermoPoroElasticity::initElement: Detected " << ierr/3
    << " node numbers out of range." << std::endl;

  return false;
}


LocalIntegral* ThermoPoroElasticity::WeakDirichlet::getLocalIntegral(const std::vector<size_t>& nen,
                                                                     size_t,
                                                                     bool neumann) const
{
  const size_t nedof1 = nsd*nen[0];           //!< Number of DOFs on basis 1
  const size_t nedof  = nedof1 + 2*nen[1];    //!< Total number of DOFs

  ElmMats* result = new MixedElmMats();

  //result->rhsOnly = false;
  result->withLHS = true;
  result->b[Fu].resize(nedof1);
  result->b[Fp].resize(nen[1]);
  result->b[FT].resize(nen[1]);
  result->b[Fprev].resize(nedof);
  result->b[Fnow].resize(nedof);
  result->b[Fres].resize(nedof);
  result->b[Fc].resize(nedof);

  result->A[uu].resize(nedof1,nedof1);
  result->A[up].resize(nedof1,nen[1]);
  result->A[uT].resize(nedof1,nen[1]);
  result->A[pp].resize(nen[1],nen[1]);
  result->A[TT].resize(nen[1],nen[1]);
  result->A[Kprev].resize(nedof,nedof);
  result->A[Know].resize(nedof,nedof);
  result->A[Ktan].resize(nedof,nedof);

  return result;
}


bool ThermoPoroElasticity::WeakDirichlet::evalBouMx(LocalIntegral& elmInt,
                                                    const MxFiniteElement& fe,
                                                    const TimeDomain& time,
                                                    const Vec3& X,
                                                    const Vec3& normal) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  // Evaluate the Neumann heat flux on the boundary
  double qT = 0.0;
  if (flux)
    qT = (*flux)(X);

  size_t ru = elMat.A[uu].rows();

  for (size_t i = 1; i <= fe.basis(2).size(); i++)
  {
    for (size_t j = 1; j <= fe.basis(2).size(); j++)
    {
      size_t li = ru+2*i;
      size_t lj = ru+2*j;
      elMat.A[Know](li,lj)  += time.dt*fe.basis(2)(i)*lambdae*fe.basis(2)(j)*fe.detJxW;
      elMat.A[TT](i,j) += time.dt*fe.basis(2)(i)*lambdae*fe.basis(2)(j)*fe.detJxW;
    }
    elMat.b[FT](i) += time.dt*(qT*fe.basis(2)(i) + lambdae*Te*fe.basis(2)(i))*fe.detJxW;
  }

  Vector PoTo, PcTc, prevSol, currSol;
  PoTo.resize(2*elMat.vec[Po].size());
  PcTc.resize(2*elMat.vec[Pc].size());
  for (size_t i = 1; i <= elMat.vec[Po].size(); i++) {
    PoTo(2*i-1) = elMat.vec[Po](i);
    PoTo(2*i  ) = elMat.vec[To](i);
    PcTc(2*i-1) = elMat.vec[Pc](i);
    PcTc(2*i  ) = elMat.vec[Tc](i);
  }

  elMat.b[Fc] = elMat.vec[Uc];
  elMat.b[Fc].insert(elMat.b[Fc].end(),PcTc.begin(),PcTc.end());

  return true;
}

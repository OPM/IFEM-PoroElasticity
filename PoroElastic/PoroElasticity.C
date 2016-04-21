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
#include "IFEM.h"
#include "tinyxml.h"


PoroElasticity::PoroElasticity (unsigned short int n, int order) : Elasticity(n)
{
  eS = Fu + 1;               // Elasticity::evalBou stores traction here
  primsol.resize(1+order);   // Current and previous timestep solutions required
  sc = 0.0;
  gacc = 9.81;               // kmo: Danger! hard-coded physical property. Why
                             // not derive this one from gravity.length()
                             // instead ???
}


void PoroElasticity::setMode (SIM::SolutionMode mode)
{
  Elasticity::setMode(mode);
  eS = Fu + 1;
}


bool PoroElasticity::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"scaling"))
    return this->Elasticity::parse(elem);
  else if (utl::getAttribute(elem,"value",sc))
    IFEM::cout <<"\tScaling: sc = "<< sc << std::endl;

  return true;
}


Material* PoroElasticity::parseMatProp (const TiXmlElement* elem, bool)
{
  IFEM::cout <<" Poroelastic material, see \"Problem definition:\" below."
             << std::endl;
  material = new PoroMaterial();
  material->parse(elem);
  return material;
}


void PoroElasticity::printLog () const
{
  IFEM::cout <<"PoroElasticity: scaling = "<< sc << std::endl;
  this->Elasticity::printLog();
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
  if (primsol.empty() || primsol.front().empty())
    return true;

  // Split the nodal correspondance array into one for each solution field
  std::vector<int>::const_iterator pStart = MNPC.begin() + elem_sizes.front();
  std::vector<int> MNPCu(MNPC.begin(),pStart);
  std::vector<int> MNPCp(pStart,MNPC.end());

  if (elmInt.vec.size() < NSOL)
    elmInt.vec.resize(NSOL);

  // Extract the element level solution vectors
  int ierr = utl::gather(MNPCu, nsd, primsol.front(), elmInt.vec[Vu]);
  ierr += utl::gather(MNPCp, 0, 1, primsol.front(), elmInt.vec[Vp],
                      nsd*basis_sizes.front(), basis_sizes.front());

  if (m_mode == SIM::DYNAMIC) {
    // Extract the element level velocity vectors
    ierr += utl::gather(MNPCu, nsd, primsol[1], elmInt.vec[Vuvel]);
    ierr += utl::gather(MNPCp, 0, 1, primsol[1], elmInt.vec[Vpvel],
                        nsd*basis_sizes.front(), basis_sizes.front());

    // Extract the element level acceleration vector for displacement
    ierr += utl::gather(MNPCu, nsd, primsol[2], elmInt.vec[Vuacc]);
  }

  if (ierr == 0) return true;

  std::cerr <<" *** PoroElasticity::initElement: Detected "<< ierr/2
            <<" node numbers out of range."<< std::endl;

  return false;
}


bool PoroElasticity::initElement (const std::vector<int>& MNPC,
                                  LocalIntegral& elmInt)
{
  if (primsol.empty() || primsol.front().empty())
    return true;

  if (elmInt.vec.size() < NSOL)
    elmInt.vec.resize(NSOL);

  Matrix temp(nsd+1, MNPC.size());

  // Extract the element level solution vectors
  int ierr = utl::gather(MNPC, nsd+1, primsol.front(), temp);
  elmInt.vec[Vp] = temp.getRow(nsd+1);
  elmInt.vec[Vu] = temp.expandCols(-1);

  if (m_mode == SIM::DYNAMIC) {
    temp.resize(nsd+1, MNPC.size());
    ierr += utl::gather(MNPC, nsd+1, primsol[1], temp);
    elmInt.vec[Vpvel] = temp.getRow(nsd+1);
    elmInt.vec[Vuvel] = temp.expandCols(-1);

    temp.resize(nsd+1, MNPC.size());
    ierr += utl::gather(MNPC, nsd+1, primsol[2], temp);
    elmInt.vec[Vuacc] = temp.expandCols(-1);
  }

  if (ierr == 0)
    return true;

  std::cerr <<" *** PoroElasticity::initElement: Detected "<< ierr
            <<" node numbers out of range."<< std::endl;

  return false;
}


bool PoroElasticity::evalIntMx (LocalIntegral& elmInt,
                                const MxFiniteElement& fe,
                                const TimeDomain& time, const Vec3& X) const
{
  return this->evalInt(elmInt, fe, time, X);
}


bool PoroElasticity::evalElasticityMatrices (ElmMats& elMat, const Matrix& B,
                                             const FiniteElement& fe,
                                             const Vec3& X) const
{
  Matrix C;
  SymmTensor eps(nsd), sigma(nsd); double U = 0.0;
  if (!material->evaluate(C,sigma,U,fe,X,eps,eps,0))
    return false;

  Matrix CB;
  CB.multiply(C,B).multiply(fe.detJxW); // CB = dSdE*B*|J|*w
  elMat.A[uu_K].multiply(B,CB,true,false,true); // EK += B^T * CB

  return true;
}


bool PoroElasticity::evalMassMatrix (Matrix& mx, const Vector& N, double scl) const
{
  Matrix temp(N.size(), N.size());
  temp.outer_product(N, N);
  temp *= scl;

  for (size_t i = 0; i < N.size(); i++)
    for (size_t j = 0; j < N.size(); j++)
      for (size_t k = 1; k <= nsd; k++)
        mx(i*nsd+k, j*nsd+k) = temp(i+1, j+1);

  return true;
}


bool PoroElasticity::evalCouplingMatrix (Matrix& mx, const Matrix& B,
                                         const Vector& N, double scl) const
{
  Matrix K(N.size(), nsd * (nsd + 1) / 2);
  for (size_t i = 1; i <= N.size(); i++)
    for (size_t j = 1; j <= nsd; j++)
      K(i,j) = scl * N(i);

  mx.multiply(B, K, true, true, true);

  return true;
}


bool PoroElasticity::evalCompressibilityMatrix (Matrix& mx, const Vector& N,
                                                double scl) const
{
  Matrix temp(mx.rows(), mx.cols());
  if (!temp.outer_product(N,N))
    return false;

  mx.add(temp,scl);

  return true;
}


bool PoroElasticity::evalPermeabilityMatrix (Matrix& mx, const Matrix& dNdX,
                                             const Vec3& permeability,
                                             double scl) const
{
  for (size_t i = 1; i <= dNdX.rows(); i++)
    for (size_t j = 1; j <= dNdX.rows(); j++)
      for (size_t k = 1; k <= nsd; k++)
        mx(i,j) += scl * permeability[k-1]  * dNdX(i,k) * dNdX(j,k);

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

  Matrix Bmat;
  if (!this->formBmatrix(Bmat,fe.dNdX))
    return false;

  Vec3 permeability = pmat->getPermeability(X);

  double rhot = pmat->getMassDensity(X);
  double rhog = pmat->getFluidDensity(X) * gacc;
  double scl = sc;
  if (scl == 0.0)
    scl = sqrt(pmat->getStiffness(X)*rhog/(permeability.x*time.dt));

  // Biot's coefficient
  double Ko = pmat->getBulkMedium(X);
  double Ks = pmat->getBulkSolid(X);
  double Kw = pmat->getBulkWater(X);
  double poro = pmat->getPorosity(X);

  double alpha = 1.0 - (Ko/Ks);
  // Inverse of the compressibility modulus
  double Minv = ((alpha - poro)/Ks) + (poro/Kw);

  if (m_mode == SIM::DYNAMIC && !this->evalMassMatrix(elMat.A[uu_M], fe.basis(1),
                                                      rhot*fe.detJxW))
    return false;

  if (!this->evalElasticityMatrices(elMat, Bmat, fe, X))
    return false;

  if (!this->evalCouplingMatrix(elMat.A[up], Bmat, fe.basis(2),
                                scl*alpha*fe.detJxW))
    return false;

  if (!this->evalCompressibilityMatrix(elMat.A[pp_S], fe.basis(2), 
                                       scl*scl*Minv*fe.detJxW))
    return false;

  return this->evalPermeabilityMatrix(elMat.A[pp_P], fe.grad(2),
                                      permeability, scl*scl/rhog*fe.detJxW);
}


bool PoroElasticity::evalBouMx (LocalIntegral& elmInt,
                                const MxFiniteElement& fe,
                                const TimeDomain&, const Vec3& X,
                                const Vec3& normal) const
{
  // Using the inherited standard method from Elasticity
  return this->evalBou(elmInt, fe, X, normal);
}


bool PoroElasticity::finalizeElement (LocalIntegral& elmInt, const TimeDomain& time, size_t)
{
  Mats& elMat = static_cast<Mats&>(elmInt);
  elMat.setStepSize(time.dt);

  return true;
}


bool PoroElasticity::finalizeElementBou (LocalIntegral& elmInt, const FiniteElement&,
                                         const TimeDomain& time)
{
  Mats& elMat = static_cast<Mats&>(elmInt);
  elMat.setStepSize(time.dt);

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


bool PoroElasticity::evalSol (Vector& s, const MxFiniteElement& fe,
                              const Vec3& X, const std::vector<int>& MNPC,
                              const std::vector<size_t>& elem_sizes) const
{
  Vector eV;
  std::vector<int>::const_iterator fstart = MNPC.begin() + elem_sizes.front();
  utl::gather(std::vector<int>(MNPC.begin(),fstart),nsd,primsol.front(),eV);

  return this->evalSolCommon(s,fe,X,Vector(&eV[0],nsd*fe.N.size()));
}


bool PoroElasticity::evalSol (Vector& s, const FiniteElement& fe,
                              const Vec3& X, const std::vector<int>& MNPC) const
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

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
#include "ElmNorm.h"
#include "TimeDomain.h"
#include "Utilities.h"
#include "Functions.h"
#include "AnaSol.h"
#include "Tensor.h"
#include "Vec3Oper.h"
#include "IFEM.h"
#include "tinyxml.h"


PoroElasticity::PoroElasticity (unsigned short int n, bool mix, bool staticFlow) :
  Elasticity(n), staticFlow(staticFlow)
{
  gravity[n-1] = 9.81; // Default gravity acceleration
  npv = mix ? 0 : nsd+1; // Number of primary unknowns per node (non-mixed only)

  nondim = calculateEnergy = useDynCoupling = residual = false;

  volumeFlux = fluxFld = nullptr;
}


PoroElasticity::~PoroElasticity ()
{
  delete volumeFlux;
}


bool PoroElasticity::parse (const TiXmlElement* elem)
{
  if (!strncasecmp(elem->Value(),"calcEn",6))
    calculateEnergy = true;
  else if (!strncasecmp(elem->Value(),"useDynC",7))
    useDynCoupling = true;
  else if (!strncasecmp(elem->Value(),"volumeflux",10))
  {
    std::string type;
    utl::getAttribute(elem,"type",type);
    volumeFlux = utl::parseRealFunc(elem->GetText(),type);
  }
  else if (!strcasecmp(elem->Value(),"nondimensionalize"))
  {
    nondim = true;
    IFEM::cout << "\tNondimensionalized mode enabled" << std::endl;

    // The default values of characteristics are all 1.0.
    // Later, we will compute some of them automatically if they haven't been set in the input file.
    // Therefore, set them to zero here (acts like a sentintel value).
    cars.E = cars.alpha = cars.perm = 0.0;

    if (utl::getAttribute(elem, "E", cars.E))
      IFEM::cout << "\t\tCharacteristic E = " << cars.E << std::endl;
    if (utl::getAttribute(elem, "alpha", cars.alpha))
      IFEM::cout << "\t\tCharacteristic alpha = " << cars.alpha << std::endl;
    if (utl::getAttribute(elem, "perm", cars.perm))
      IFEM::cout << "\t\tCharacteristic permeability = " << cars.perm << std::endl;
  }
  else
    return this->Elasticity::parse(elem);

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


void PoroElasticity::getNodalDofTypes (std::vector<char>& dType) const
{
  dType.resize(npv,'D');
  if (npv > nsd)
    dType.back() = 'P';
}


void PoroElasticity::printLog () const
{
  IFEM::cout <<"PoroElasticity: useDynCoupling = "
             << std::boolalpha << useDynCoupling << std::endl;
  this->Elasticity::printLog();
}


void PoroElasticity::setMode (SIM::SolutionMode mode)
{
  this->Elasticity::setMode(mode);
  eS = Fu + 1; // Elasticity::evalBou stores traction force vectors here
  iS = mode == SIM::DYNAMIC ? 1 : 0; // Flag calculation of internal forces
}


bool PoroElasticity::init (const TimeDomain& time)
{
  if (!nondim)
    return true;

  IFEM::cout << "PoroElasticity:" << std::endl;

  const PoroMaterial* pmat = dynamic_cast<const PoroMaterial*>(material);
  if (pmat)
  {
    Vec3 X;

    if (cars.E == 0.0) {
      cars.E = pmat->getStiffness(X);
      IFEM::cout << "\tComputed characteristic E = " << cars.E << std::endl;
    }

    if (cars.alpha == 0.0) {
      cars.alpha = pmat->getBiotCoeff(X);
      IFEM::cout << "\tComputed characteristic alpha = " << cars.alpha << std::endl;
    }

    if (cars.perm == 0.0) {
      Vec3 perm = pmat->getPermeability(X);
      for (int d = 1; d <= nsd; d++)
        cars.perm += perm(d);
      cars.perm /= pmat->getViscosity(X) * nsd;
      IFEM::cout << "\tComputed characteristic permeability = " << cars.perm << std::endl;
    }
  }

  cars.normalize();

  IFEM::cout << "\tComputed characteristic pressure = " << cars.p << " [Pa]" << std::endl
             << "\tComputed characteristic time = " << cars.t << " [s]" << std::endl;

  return true;
}


LocalIntegral* PoroElasticity::getLocalIntegral (const std::vector<size_t>& nen,
                                                 size_t, bool neumann) const
{
  if (m_mode == SIM::DYNAMIC)
    return new NewmarkMats<MixedElmMats>(nen[0], nen[1], neumann,
                                         intPrm[2], intPrm[3],
                                         intPrm[0], intPrm[1], useDynCoupling, nsd);
  else
    return new MixedElmMats(nen[0], nen[1], neumann, staticFlow ? 0 : 1, residual, nsd);
}


LocalIntegral* PoroElasticity::getLocalIntegral (size_t nen,
                                                 size_t, bool neumann) const
{
  if (m_mode == SIM::DYNAMIC)
    return new NewmarkMats<StdElmMats>(nen, nen, neumann,
                                       intPrm[2], intPrm[3],
                                       intPrm[0], intPrm[1], useDynCoupling, nsd);
  else
    return new StdElmMats(nen, nen, neumann, staticFlow ? 0 : 1, residual, nsd);
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

  size_t nsol = m_mode == SIM::DYNAMIC ? NSOL : 2;
  if (elmInt.vec.size() < nsol)
    elmInt.vec.resize(nsol);

  // Extract the element level solution vectors
  int ierr = utl::gather(MNPCu, nsd, primsol.front(), elmInt.vec[Vu]);
  ierr += utl::gather(MNPCp, 0, 1, primsol.front(), elmInt.vec[Vp],
                      nsd*basis_sizes.front(), basis_sizes.front());

  if (m_mode == SIM::DYNAMIC)
  {
    // Extract the element level velocity vectors
    ierr += utl::gather(MNPCu, nsd, primsol[1], elmInt.vec[Vuvel]);
    ierr += utl::gather(MNPCp, 0, 1, primsol[1], elmInt.vec[Vpvel],
                        nsd*basis_sizes.front(), basis_sizes.front());
    // Extract the element level acceleration vector for displacement
    ierr += utl::gather(MNPCu, nsd, primsol[2], elmInt.vec[Vuacc]);
  }

  if (ierr == 0)
    return true;

  std::cerr <<" *** PoroElasticity::initElement: Detected "<< ierr
            <<" node numbers out of range."<< std::endl;
  return false;
}


bool PoroElasticity::initElement (const std::vector<int>& MNPC,
                                  LocalIntegral& elmInt)
{
  if (primsol.empty() || primsol.front().empty())
    return true;

  size_t nsol = m_mode == SIM::DYNAMIC ? NSOL : 2;
  if (elmInt.vec.size() < nsol)
    elmInt.vec.resize(nsol);

  Matrix temp(npv, MNPC.size());

  // Extract the element level solution vectors
  int ierr = utl::gather(MNPC, npv, primsol.front(), temp);
  elmInt.vec[Vp] = temp.getRow(npv);
  elmInt.vec[Vu] = temp.expandRows(-1);

  size_t ip = primsol.size();
  if (m_mode == SIM::DYNAMIC && ip > 2)
  {
    // Extract the element level acceleration vector
    temp.resize(npv, MNPC.size());
    ierr += utl::gather(MNPC, npv, primsol[--ip], temp);
    elmInt.vec[Vuacc] = temp.expandRows(-1);
    // Extract the element level velocity vectors
    temp.resize(npv, MNPC.size());
    ierr += utl::gather(MNPC, npv, primsol[--ip], temp);
    elmInt.vec[Vpvel] = temp.getRow(npv);
    elmInt.vec[Vuvel] = temp.expandRows(-1);
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
  SymmTensor eps(nsd), sigma(nsd);
  if (!B.multiply(elMat.vec[Vu],eps))
    return false;

  Matrix C; double U = 0.0;
  if (!material->evaluate(C,sigma,U,fe,X,eps,eps,iS))
    return false;

  // Integrate the (material) stiffness
  Matrix CB;
  if (CB.multiply(C,B).multiply(fe.detJxW).empty()) // CB = dSdE*B*|J|*w
    return false;

  if (elMat.A[uu_K].multiply(B,CB,true,false,true).empty()) // K_uu += B^T * CB
    return false;

  // Integrate the internal forces, if deformed configuration
  if (iS == 0 || eps.isZero(1.0e-16))
    return true;

  sigma *= fe.detJxW;
  return B.multiply(sigma,elMat.b[Fu],true,-1); // F_u -= B^T*sigma
}


bool PoroElasticity::evalCouplingMatrix (Matrix& mx, const Matrix& B,
                                         const Vector& N, double scale) const
{
  Matrix K(N.size(), nsd * (nsd + 1) / 2);
  for (size_t i = 1; i <= N.size(); i++)
    for (size_t j = 1; j <= nsd; j++)
      K(i,j) = scale * N(i);

  return !mx.multiply(B,K,true,true,true).empty();
}


bool PoroElasticity::evalCompressibilityMatrix (Matrix& mx, const Vector& N,
                                                double scale) const
{
  return mx.outer_product(N,N,true,scale);
}


void PoroElasticity::evalPermeabilityMatrix (Matrix& mx, const Matrix& dNdX,
                                             const SymmTensor& K,
                                             double scale) const
{
  for (size_t i = 1; i <= dNdX.rows(); i++)
    for (size_t j = 1; j <= dNdX.rows(); j++)
      for (size_t k = 1; k <= nsd; k++)
        for (size_t l = 1; l <= nsd; l++)
          mx(i,j) += scale * dNdX(i,k) * K(k,l) * dNdX(j,l);
}


void PoroElasticity::evalDynCouplingMatrix (Matrix& mx, const Vector& Nu,
                                            const Matrix& dNpdx,
                                            const SymmTensor& K,
                                            double scale) const
{
  for (size_t i = 1; i <= dNpdx.rows(); i++)
    for (size_t j = 1; j <= Nu.size(); j++)
      for (size_t k = 1; k <= nsd; k++)
        for (size_t l = 1; l <= nsd; l++)
          mx(nsd*(j-1) + k, i) += dNpdx(i,k) * K(k,l) * Nu(j) * scale;
}


bool PoroElasticity::formPermeabilityTensor (SymmTensor& K,
                                             const Vectors&,
                                             const FiniteElement&,
                                             const Vec3& X) const
{
  const PoroMaterial* pmat = dynamic_cast<const PoroMaterial*>(material);
  if (!pmat) return false;

  Vec3 permeability = pmat->getPermeability(X);

  K.zero();
  for (size_t i = 1; i <= K.dim(); i++)
    K(i,i) = permeability(i);

  return true;
}


bool PoroElasticity::evalInt (LocalIntegral& elmInt,
                              const FiniteElement& fe,
                              const TimeDomain& time, const Vec3& X) const
{
  const PoroMaterial* pmat = dynamic_cast<const PoroMaterial*>(material);
  if (!pmat)
  {
    std::cerr <<" *** PoroElasticity::evalInt: No material data."<< std::endl;
    return false;
  }

  Matrix Bmat;
  if (!this->formBmatrix(Bmat,fe.dNdX))
    return false;

  SymmTensor Kperm(nsd); // Evaluate the permeability tensor
  if (!this->formPermeabilityTensor(Kperm,elmInt.vec,fe,X))
    return false;

  // Evaluate other material parameters
  double visc  = pmat->getViscosity(X);
  double poro  = pmat->getPorosity(X);
  double alpha = pmat->getBiotCoeff(X);
  double Minv  = pmat->getBiotModulus(X,alpha,poro);

  // Integrate the element matrices, depending on solution mode
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  if (!elMat.A[uu_M].empty())
    this->formMassMatrix(elMat.A[uu_M], fe.basis(1), X, fe.detJxW);

  if (!elMat.A[up_D].empty())
    this->evalDynCouplingMatrix(elMat.A[up_D], fe.basis(1), fe.grad(2),
                                Kperm, 1.0 / visc * fe.detJxW);

  this->evalPermeabilityMatrix(elMat.A[pp_P], fe.grad(2),
                               Kperm, 1.0 / visc * fe.detJxW);

  if (!this->evalElasticityMatrices(elMat, Bmat, fe, X))
    return false;

  if (!this->evalCouplingMatrix(elMat.A[up_Q], Bmat, fe.basis(2),
                                alpha * fe.detJxW))
    return false;

  // In the fully static formulation, we don't need the S-matrix
  if (m_mode == SIM::DYNAMIC || !staticFlow)
    if (!this->evalCompressibilityMatrix(elMat.A[pp_S], fe.basis(2),
                                         Minv * fe.detJxW))
        return false;

  if (volumeFlux)
    elMat.b[Fp].add(fe.basis(2),(*volumeFlux)(X)*fe.detJxW);

  return true;
}


bool PoroElasticity::evalBou (LocalIntegral& elmInt,
                              const FiniteElement& fe,
                              const TimeDomain&, const Vec3& X,
                              const Vec3& normal) const
{
  if (fluxFld)
  {
    double hbar = -(*fluxFld)(X);
    static_cast<ElmMats&>(elmInt).b[Fp].add(fe.basis(2),hbar*fe.detJxW);
    return true;
  }

  // Using the inherited standard method from Elasticity for tractions
  return this->Elasticity::evalBou(elmInt,fe,X,normal);
}


bool PoroElasticity::evalBouMx (LocalIntegral& elmInt,
                                const MxFiniteElement& fe,
                                const TimeDomain& time, const Vec3& X,
                                const Vec3& normal) const
{
  return this->evalBou(elmInt,fe,time,X,normal);
}


bool PoroElasticity::finalizeElement (LocalIntegral& elmInt,
                                      const TimeDomain& time, size_t)
{
  Mats& mats = static_cast<Mats&>(elmInt);

  mats.setStepSize(time.dt);

  // All scales are 1.0 if nondimensional mode is not set, so this is safe
  mats.nondimensionalize(cars);

  return true;
}


bool PoroElasticity::finalizeElementBou (LocalIntegral& elmInt,
                                         const FiniteElement&,
                                         const TimeDomain& time)
{
  Mats& mats = static_cast<Mats&>(elmInt);

  mats.setStepSize(time.dt);

  // All scales are 1.0 if nondimensional mode is not set, so this is safe
  mats.nondimensionalize(cars);

  return true;
}


size_t PoroElasticity::getNoFields (int fld) const
{
  if (fld < 2)
    return nsd+1;

  // Strain, stress and fluid content
  return nsd * (nsd + 1) + 1;
}


std::string PoroElasticity::getField1Name (size_t i, const char* prefix) const
{
  if (i == 11 && npv != nsd+1)
    return "Displacement";
  else if (i == 11 && npv == nsd+1)
    return nsd == 2 ? "Displacement_x&&Displacement_y&&Pressure" :
                      "Displacement_x&&Displacement_y&&Displacement_z&&Pressure";
  if (i == 12)
    return "Pressure";
  else if (i >= nsd)
    i = 3;

  static const char* s[4] = {"u_x", "u_y", "u_z", "p^w"};

  if (!prefix)
    return s[i];

  return prefix + std::string(" ") + s[i];
}


std::string PoroElasticity::getField2Name (size_t i, const char* prefix) const
{
  static const char* s[][6] = {{"x", "y", "xy", "", "", ""},
                               {"x", "y", "z", "yz", "xz", "xy"}};

  size_t ncmp = nsd * (nsd + 1) / 2;
  std::string name;
  if (i < 2 * ncmp)
    name = std::string(i < ncmp ? "eps_" : "sig_") + s[nsd-2][i%ncmp];
  else
    name = "FluidContent";

  if (!prefix)
    return name;
  return prefix + std::string(" ") + name;
}


bool PoroElasticity::evalSol (Vector& s, const MxFiniteElement& fe,
                              const Vec3& X, const std::vector<int>& MNPC,
                              const std::vector<size_t>& elem_sizes,
                              const std::vector<size_t>&) const
{
  auto split = MNPC.begin() + elem_sizes.front();
  std::vector<int> MNPC1(MNPC.begin(), split);
  std::vector<int> MNPC2(split, split + elem_sizes[1]);

  Vector displacement;
  utl::gather(MNPC1, nsd, primsol.front(), displacement);

  Vector pressure;
  utl::gather(MNPC2, 1, primsol.front(), pressure);

  return this->evalSol(s,fe,X,displacement,pressure.dot(fe.basis(2)));
}


bool PoroElasticity::evalSol (Vector& s, const FiniteElement& fe,
                              const Vec3& X, const std::vector<int>& MNPC) const
{
  Vector eV;
  utl::gather(MNPC,npv,primsol.front(),eV);

  Vector displacement(nsd * fe.N.size());
  for (size_t i = 0; i < nsd; i++)
    for (size_t a = 0; a < fe.N.size(); a++)
      displacement[nsd*a+i] = eV[npv*a+i];

  double pressure = eV.dot(fe.N, nsd, nsd+1);
  return this->evalSol(s,fe,X,displacement,pressure);
}


bool PoroElasticity::evalSol (Vector& s, const FiniteElement& fe,
                              const Vec3& X, const Vector& disp, double pressure) const
{
  if (!material)
  {
    std::cerr <<" *** PoroElasticity::evalSol: No material data."<< std::endl;
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

  const PoroMaterial* pmat = dynamic_cast<const PoroMaterial*>(material);
  if (!pmat) {
    std::cerr <<" *** PoroElasticity::evalSol: Wrong material data."<< std::endl;
    return false;
  }

  double alpha = pmat->getBiotCoeff(X);
  double porosity = pmat->getPorosity(X);
  double Minv = pmat->getBiotModulus(X, alpha, porosity);

  // Fluid content is given implicitly, see (10) in Bekele2017omi
  double fluidContent = eps.trace() * pmat->getBiotCoeff(X) + pressure * Minv;
  s.insert(s.end(), fluidContent);

  return true;
}


NormBase* PoroElasticity::getNormIntegrand (AnaSol* sol) const
{
  if (!calculateEnergy)
    return nullptr;

  if (sol)
    return new PoroNorm(*const_cast<PoroElasticity*>(this),
                        sol->getVectorSol(), sol->getVectorSecSol(),
                        sol->getScalarSol(), sol->getScalarSecSol());

  return new PoroNorm(*const_cast<PoroElasticity*>(this));
}


PoroNorm::PoroNorm (PoroElasticity& poroel,
                    VecFunc* disp, TensorFunc* /*d_disp*/,
                    RealFunc* press, VecFunc* /*d_press*/) : NormBase(poroel)
{
  nrcmp = 2;
  displacement = disp;
//  d_displacement = d_disp;
  pressure = press;
//  d_pressure = d_press;
}


size_t PoroNorm::getNoFields (int group) const
{
  // FIXME: Projected solutions
  if (group == 0)
    return 1;

  size_t ret = 2;
  if (displacement) ret++;
  if (pressure) ret++;
  return ret;
}


std::string PoroNorm::getName (size_t i, size_t j, const char* prefix) const
{
  if (i == 0 || j == 0 || j > 4)
    return this->NormBase::getName(i, j, prefix);

  static const char* s[4] = {
    "|u^h|_l2",
    "|p^h|_l2",
    "|u-u^h|_l2",
    "|p-p^h|_l2",
  };

  if (!displacement && j > 2)
    j++;

  if (!prefix)
    return s[j-1];

  return prefix + std::string(" ") + s[j-1];
}


bool PoroNorm::evalIntMx (LocalIntegral& elmInt, const MxFiniteElement& fe,
                          const TimeDomain& time, const Vec3& X) const
{
  return this->evalInt(elmInt, fe, time, X);
}


bool PoroNorm::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                        const TimeDomain& time, const Vec3& X) const
{
  ElmNorm& norms = static_cast<ElmNorm&>(elmInt);
  size_t nsd = fe.grad(1).cols();

  // Numerical displacement and pressure
  Vec3   disp_h;
  double press_h;
  if (fe.getNoBasis() == 1)
  {
    for (size_t i = 0; i < nsd; i++)
      disp_h[i] = norms.vec[PoroElasticity::Vu].dot(fe.N,i,nsd+1);
    press_h = norms.vec[PoroElasticity::Vp].dot(fe.N,nsd,nsd+1);
  }
  else
  {
    for (size_t i = 0; i < nsd; i++)
      disp_h[i] = norms.vec[PoroElasticity::Vu].dot(fe.basis(1),i,nsd);
    press_h = norms.vec[PoroElasticity::Vp].dot(fe.basis(2));
  }

  // Norm index counter
  size_t np = 0;

  // Displacement L2-norm
  norms[np++] += disp_h * disp_h * fe.detJxW;

  // Pressure L2-norm
  norms[np++] += press_h * press_h * fe.detJxW;

  // Displacement error L2-norm
  if (displacement)
  {
    Vec3 disp_a = (*displacement)(Vec4(X,time.t));
    norms[np++] += (disp_a - disp_h) * (disp_a - disp_h) * fe.detJxW;
  }

  // Pressure error L2-norm
  if (pressure)
  {
    double press_a = (*pressure)(Vec4(X,time.t));
    norms[np] += (press_a - press_h) * (press_a - press_h) * fe.detJxW;
  }

  return true;
}

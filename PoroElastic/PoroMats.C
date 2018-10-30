// $Id$
//==============================================================================
//!
//! \file PoroMats.C
//!
//! \date April 20 2015
//!
//! \author Eivind Fonn
//!
//! \brief Element matrix implementation for PoroElasticity.
//!
//==============================================================================

#include "PoroElasticity.h"
#include "Utilities.h"

#include <cassert>


PoroElasticity::Mats::Mats (size_t ndof_displ, size_t ndof_press,
                            bool neumann, char dynamic, int nbasis, int nsd)
  : BlockElmMats(2, nbasis), dynamic(dynamic)
{
  this->resize(NMAT, NVEC);

  rhsOnly = neumann;
  withLHS = !neumann;
  this->redim(1, ndof_displ, nsd, 1);
  this->redim(2, ndof_press, 1, nbasis);
  this->redimOffDiag(3,0);
  this->redimOffDiag(4,0);
  this->redimNewtonMat();

  if (withLHS) {
    if (dynamic > 1)
      A[uu_M].resize(A[uu_K].rows(), A[uu_K].rows());
    if (dynamic == 3)
      A[up_D].resize(A[up_Q].rows(), A[up_Q].cols());
    A[pp_P].resize(A[pp_S].rows(), A[pp_S].cols());
  }

  h = lambda = 0.0;
}


// NOTE: This method assumes full- or half-static formulation.
// It should be overridden in subclasses for dynamic integration.
const Matrix& PoroElasticity::Mats::getNewtonMatrix () const
{
  assert(dynamic < 2);

  if (dynamic == 1)
    const_cast<Matrix&>(A[pu_Q]).addBlock(A[up_Q], 1.0, 1, 1, true);

  // NOTE: We mutate up_Q on the assumption that getRHSVector is called
  // before getNewtonMatrix, as the former assumes that up_Q has the sign
  // it expects. The const-ness of these two methods is highly dubious.
  const_cast<Matrix&>(A[up_Q]) *= -1.0;

#if INT_DEBUG > 2
  std::cout <<"\nPoroElasticity::Mats::getNewtonMatrix:"
            <<"\nElement stiffness matrix, K_uu"<< A[uu_K]
            <<"\nElement coupling matrix, Q_up"<< A[up_Q]
            <<"\nElement compressibility matrix, S_pp"<< A[pp_S]
            <<"\nElement permeability matrix, P_pp"<< A[pp_P];
#endif

  const_cast<Matrix&>(A[pp_S]).add(A[pp_P], (dynamic == 1) ? h : 1.0);

  const Matrix& result = this->BlockElmMats::getNewtonMatrix();
#if INT_DEBUG > 2
  std::cout <<"\nElement coefficient matrix" << result;
#endif
  return result;
}


// NOTE: This method assumes full- or half-static formulation.
// It should be overridden in subclasses for dynamic integration.
const Vector& PoroElasticity::Mats::getRHSVector () const
{
  assert(dynamic < 2);

#if INT_DEBUG > 2
  std::cout <<"\nPoroElasticity::Mats::getRHSVector:"<< std::endl;
  if (vec.size() > Vu) std::cout <<"Element displacement, Vu"<< vec[Vu];
  if (vec.size() > Vp) std::cout <<"Element pressure, Vp"<< vec[Vp];
  std::cout <<"S_ext-S_int"<< b[Fu] <<"S_p"<< b[Fp];
#endif

  // In case of half-static formulation, add some terms to the pressure
  // RHS arising from the previous timestep.
  // Fp' = Fp * h + Q_pu * u + S_pp * p
  if (dynamic == 1) {
    Vector& fp = const_cast<Vector&>(b[Fp]);
    fp *= h;
    if (A.size() > up_Q && vec.size() > Vu)
      A[up_Q].multiply(vec[Vu], fp, true, true);
    if (A.size() > pp_S && vec.size() > Vp)
      A[pp_S].multiply(vec[Vp], fp, false, true);
  }

  const Vector& result = this->BlockElmMats::getRHSVector();
#if INT_DEBUG > 2
  std::cout <<"\nElement right-hand-side vector"<< result;
#endif
  return result;
}

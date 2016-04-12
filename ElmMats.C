// $Id$
//==============================================================================
//!
//! \file ElmMats.C
//!
//! \date April 12 2016
//!
//! \author Eivind Fonn
//!
//! \brief Element matrix implementation for PoroElasticity
//!
//==============================================================================

#include "PoroElasticity.h"
#include "Utilities.h"


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

  if (!neumann)
  {
    A[sys].resize(ndof_tot, ndof_tot);
    A[sys_M].resize(ndof_tot, ndof_tot);
    A[sys_C].resize(ndof_tot, ndof_tot);
    A[sys_K].resize(ndof_tot, ndof_tot);
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

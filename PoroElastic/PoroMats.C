// $Id$
//==============================================================================
//!
//! \file PoroMats.C
//!
//! \date April 20 2015
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
  this->resize(NMAT, NVEC);
  this->ndof_displ = ndof_displ;
  this->ndof_press = ndof_press;

  size_t ndof_tot = ndof_displ + ndof_press;

  rhsOnly = neumann;
  withLHS = !neumann;
  b[Fsys].resize(ndof_tot);
  b[Fu].resize(ndof_displ);
  b[Fp].resize(ndof_press);

  if (!neumann) {
    A[sys].resize(ndof_tot, ndof_tot);
    A[uu_K].resize(ndof_displ, ndof_displ);
    A[uu_M].resize(ndof_displ, ndof_displ);
    A[up].resize(ndof_displ, ndof_press);
    A[pp_S].resize(ndof_press, ndof_press);
    A[pp_P].resize(ndof_press, ndof_press);
  }
}


const Matrix& PoroElasticity::Mats::getNewtonMatrix() const
{
  Matrix& res = const_cast<Matrix&>(A[sys]); res.fill(0.0);
  this->add_uu(A[uu_K], res);
  this->add_up(A[up], res, -1.0);
  this->add_pu(A[up], res);
  this->add_pp(A[pp_S], res);
  this->add_pp(A[pp_P], res, h);
  return A[sys];
}


const Vector& PoroElasticity::Mats::getRHSVector() const
{
  Vector temp(b[Fp].size());
  temp.add(b[Fp], h);
  if (vec.size() > 1) {
    A[up].multiply(vec[Vu], temp, true, true);
    A[pp_S].multiply(vec[Vp], temp, false, true);
  }

  Vector& res = const_cast<Vector&>(b[Fsys]); res.fill(0.0);
  this->form_vector(b[Fu], temp, res);

  return b[Fsys];
}


void PoroElasticity::MixedElmMats::add_uu(const Matrix& source, Matrix& target, double scale) const
{
  target.addBlock(source, scale, 1, 1);
}


void PoroElasticity::MixedElmMats::add_up(const Matrix& source, Matrix& target, double scale) const
{
  target.addBlock(source, scale, 1, 1 + ndof_displ);
}


void PoroElasticity::MixedElmMats::add_pu(const Matrix& source, Matrix& target, double scale) const
{
  target.addBlock(source, scale, 1 + ndof_displ, 1, true);
}


void PoroElasticity::MixedElmMats::add_pp(const Matrix& source, Matrix& target, double scale) const
{
  target.addBlock(source, scale, 1 + ndof_displ, 1 + ndof_displ);
}


void PoroElasticity::MixedElmMats::form_vector(const Vector &u, const Vector &p, Vector& target) const
{
  target = u;
  target.insert(target.end(), p.begin(), p.end());
}


void PoroElasticity::NonMixedElmMats::add_uu(const Matrix& source, Matrix& target, double scale) const
{
  size_t nsd = ndof_displ / ndof_press;
  size_t nf = nsd + 1;

  for (size_t i = 1; i <= ndof_press; ++i)
    for (size_t j = 1; j <= ndof_press; ++j)
      for (size_t l = 1; l <= nsd; ++l)
        for (size_t k = 1; k <= nsd; ++k)
          target(nf*(i-1)+l, nf*(j-1)+k) += scale * source(nsd*(i-1)+l, nsd*(j-1)+k);
}


void PoroElasticity::NonMixedElmMats::add_up(const Matrix& source, Matrix& target, double scale) const
{
  size_t nsd = ndof_displ / ndof_press;
  size_t nf = nsd + 1;

  for (size_t i = 1; i <= ndof_press; ++i)
    for (size_t j = 1; j <= ndof_press; ++j)
      for (size_t l = 1; l <= nsd; ++l)
          target(nf*(i-1)+l, j*nf) += scale * source(nsd*(i-1)+l, j);
}


void PoroElasticity::NonMixedElmMats::add_pu(const Matrix& source, Matrix& target, double scale) const
{
  size_t nsd = ndof_displ / ndof_press;
  size_t nf = nsd + 1;

  for (size_t i = 1; i <= ndof_press; ++i)
    for (size_t j = 1; j <= ndof_press; ++j)
      for (size_t l = 1; l <= nsd; ++l)
        target(j*nf, nf*(i-1)+l) += scale * source(nsd*(i-1)+l, j);
}


void PoroElasticity::NonMixedElmMats::add_pp(const Matrix& source, Matrix& target, double scale) const
{
  size_t nsd = ndof_displ / ndof_press;
  size_t nf = nsd + 1;

  for (size_t i = 1; i <= ndof_press; ++i)
    for (size_t j = 1; j <= ndof_press; ++j)
      target(i*nf, j*nf) += scale * source(i, j);
}


void PoroElasticity::NonMixedElmMats::form_vector(const Vector& u, const Vector& p, Vector& target) const
{
  utl::interleave(u, p, target, ndof_displ / ndof_press, 1);
}

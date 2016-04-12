// $Id$
//==============================================================================
//!
//! \file PoroElasticity.h
//!
//! \date April 16 2015
//!
//! \author Yared Bekele
//!
//! \brief Integrand implementations for time-dependent PoroElasticity problems.
//!
//==============================================================================

#ifndef _PORO_ELASTICITY_H_
#define _PORO_ELASTICITY_H_

#include "Elasticity.h"
#include "ElmMats.h"


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
  Fprev = 1,                    // RHS contribution from previous timestep

  // Sub-vectors, sized according to the bases in question
  Fu = 2,                       // Traction and body forces
  Fp = 3,                       // Flux and body forces

  // Reserved for Newmark integration
  Fv = 4,                       // Velocity
  Fa = 5,                       // Acceleration

  NVEC = 6
};


//! \brief Enum for element level left-hand-side matrices
enum TangentMatrices
{
  // System matrices, "global" size (note: order is fixed according to Newmark API)
  sys = 0,                      // Final newton matrix
  sys_M = 1,                    // Acceleration term in Newmark
  sys_C = 2,                    // Velocity term in Newmark
  sys_K = 3,                    // Zero-order term in Newmark

  // Sub-matrices, sized according to the bases in question
  uu_K = 4,                     // Stiffness matrix
  uu_M = 5,                     // Mass matrix
  up = 6,                       // Coupling matrix
  pp_S = 7,                     // Compressibility matrix
  pp_P = 8,                     // Permeability matrix

  NMAT = 9
};


/*!
  \brief Class representing the integrand of the PoroElasticity problem.
*/

class PoroElasticity : public Elasticity
{
  /*!
   * \brief Superclass for PoroElasticity element matrices
   */
  class Mats : public ElmMats
  {
  public:
    //! \brief Default constructor
    //! \param[in] ndof_displ Number of dofs in displacement
    //! \param[in] ndof_press Number of dofs in pressure
    //! \param[in] neumann Whether or not we are assembling Neumann BCs
    Mats(size_t ndof_displ, size_t ndof_press, bool neumann);
    //! \brief Empty destructor
    virtual ~Mats() {}
    //! \brief Returns the element level Newton matrix
    virtual const Matrix& getNewtonMatrix() const;
    //! \brief Returns the element level RHS vector
    virtual const Vector& getRHSVector() const;
    //! \brief Adds in a UU-matrix to a system matrix
    virtual void add_uu(size_t source, size_t target, double scale = 1.0) = 0;
    //! \brief Adds in a UP-matrix to a system matrix
    virtual void add_up(size_t source, size_t target, double scale = 1.0) = 0;
    //! \brief Adds in the transpose of a UP-matrix to a system matrix
    virtual void add_pu(size_t source, size_t target, double scale = 1.0) = 0;
    //! \brief Adds in a PP-matrix to a system matrix
    virtual void add_pp(size_t source, size_t target, double scale = 1.0) = 0;
    //! \brief Forms a system vector out of two sub-vectors
    virtual void form_vector(const Vector &u, const Vector &p, size_t target) = 0;
  protected:
    size_t ndof_displ, ndof_press, nsd;
  };

  /*!
   * \brief Class representing an element matrix for the mixed PoroElasticity problem
   */
  class MixedElmMats : public Mats
  {
  public:
    //! \brief Default constructor
    MixedElmMats(size_t ndof_displ, size_t ndof_press, bool neumann)
      : Mats(ndof_displ, ndof_press, neumann) {}
    //! \brief Empty destructor
    virtual ~MixedElmMats() {}
    //! \brief Adds in a UU-matrix to a system matrix
    virtual void add_uu(size_t source, size_t target, double scale = 1.0);
    //! \brief Adds in a UP-matrix to a system matrix
    virtual void add_up(size_t source, size_t target, double scale = 1.0);
    //! \brief Adds in the transpose of a UP-matrix to a system matrix
    virtual void add_pu(size_t source, size_t target, double scale = 1.0);
    //! \brief Adds in a PP-matrix to a system matrix
    virtual void add_pp(size_t source, size_t target, double scale = 1.0);
    //! \brief Forms a system vector out of two sub-vectors
    virtual void form_vector(const Vector &u, const Vector &p, size_t target);
  };

  /*!
   * \brief Class representing an element matrix for the non-mixed PoroElasticity problem
   */
  class NonMixedElmMats : public Mats
  {
  public:
    //! \brief Default constructor
    NonMixedElmMats(size_t ndof_displ, size_t ndof_press, bool neumann)
      : Mats(ndof_displ, ndof_press, neumann) {}
    //! \brief Empty destructor
    virtual ~NonMixedElmMats() {}
    //! \brief Adds in a UU-matrix to a system matrix
    virtual void add_uu(size_t source, size_t target, double scale = 1.0);
    //! \brief Adds in a UP-matrix to a system matrix
    virtual void add_up(size_t source, size_t target, double scale = 1.0);
    //! \brief Adds in the transpose of a UP-matrix to a system matrix
    virtual void add_pu(size_t source, size_t target, double scale = 1.0);
    //! \brief Adds in a PP-matrix to a system matrix
    virtual void add_pp(size_t source, size_t target, double scale = 1.0);
    //! \brief Forms a system vector out of two sub-vectors
    virtual void form_vector(const Vector &u, const Vector &p, size_t target);
  };

  /*!
   * \brief Newmark element matrices for PoroElasticity
   */
  template<class M>
  class NewmarkMats : public M
  {
  public:
    //! \brief Default constructor
    NewmarkMats(size_t ndof_displ, size_t ndof_press, bool neumann, double b, double c)
      : M(ndof_displ, ndof_press, neumann), beta(b), gamma(c) {}
    //! \brief Empty destructor
    virtual ~NewmarkMats() {}
    //! \brief Updates the time step size and the \a isPredictor flag.
    //! \param[in] dt New time step size
    void setStepSize(double dt, int) { h = dt; }
    //! \brief Returns the element level Newton matrix
    virtual const Matrix& getNewtonMatrix() const
    {
      Matrix& N = const_cast<Matrix&>(M::A.front());
      N = M::A[sys_M];
      N.add(M::A[sys_C], gamma*h);
      N.add(M::A[sys_K], beta*h*h);
      return M::A.front();
    }
    //! \brief Returns the element level RHS vector
    virtual const Vector& getRHSVector() const
    {
      Vector &F = const_cast<Vector&>(M::b.front());
      F.add(M::A[sys_M] * M::b[Fa], -1.0);
      F.add(M::A[sys_C] * M::b[Fv], -1.0);
      return M::b.front();
    }
  protected:
    double beta, gamma, h;
  };

public:
  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] order Order of the time-integration scheme
  PoroElasticity(unsigned short int n = 3, int order = 1);
  //! \brief Empty destructor.
  virtual ~PoroElasticity() {}

  //! \brief Defines the scaling factor.
  void setScaling(double scaling) { sc = scaling; }
  //! \brief Obtain current scaling factor.
  double getScaling() const { return sc; }

  //! \brief Returns the current gravity vector.
  const Vec3 getGravity() const { return gravity; }

  //! \brief Computes the stiffness matrix for a quadrature point.
  bool evalStiffnessMatrix(Matrix& mx, const Matrix &B, const Matrix &C, double detJxW) const;

  //! \brief Computes the n-dimensional mass matrix for a quadrature point.
  bool evalMassMatrix(Matrix &mx, const Vector &basis, double rho, double detJxW) const;

  //! \brief Computes the coupling matrix for a quadrature point.
  bool evalCouplingMatrix(Matrix &mx, const Matrix &B, const Vector &basis,
                          double scl, double alpha, const Vector &m, double detJxW) const;

  //! \brief Computes the compressibility matrix for a quadrature point.
  bool evalCompressibilityMatrix(Matrix &mx, const Vector &basis,
                                 double scl, double Minv, double detJxW) const;

  //! \brief Computes the permeability matrix for a quadrature point.
  bool evalPermeabilityMatrix(Matrix &mx, const Matrix &grad, double scl,
                              const Vec3 &permeability, double acc_dens, double detJxW) const;

  using IntegrandBase::getLocalIntegral;
  //! \brief Returns a local integral container for the given element
  //! \param[in] nen1 Number of nodes on element for basis 1
  //! \param[in] nen2 Number of nodes on element for basis 2
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  virtual LocalIntegral* getLocalIntegral(const std::vector<size_t>& nen,
                                          size_t, bool neumann) const;
  //! \brief Returns a local integral contribution object for the given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] iEl Global element number (1-based)
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  virtual LocalIntegral* getLocalIntegral(size_t nen, size_t iEl,
                                          bool neumann = false) const;

  //! \brief Initializes current element for numerical integration
  //! \param[in] MNPC1 Nodal point correspondence for basis 1
  //! \param[in] MNPC2 Nodal point correspondence for basis 2
  //! \param[in] n1 Number of nodes in basis 1 on this patch
  //! \param elmInt The local integral object for current element
  virtual bool initElement(const std::vector<int>& MNPC,
                           const std::vector<size_t>& elem_sizes,
                           const std::vector<size_t>& basis_sizes,
                           LocalIntegral& elmInt);

  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param elmInt Local integral for element
  //!
  //! \details This method is invoked once before starting the numerical
  //! integration loop over the Gaussian quadrature points over an element.
  //! It is supposed to perform all the necessary internal initializations
  //! needed before the numerical integration is started for current element.
  virtual bool initElement(const std::vector<int>& MNPC, LocalIntegral& elmInt);

  //! \brief Initializes current element for numerical boundary integration (mixed)
  //! \param[in] MNPC1 Nodal point correspondence for basis 1
  //! \param[in] MNPC2 Nodal point correspondence for basis 2
  //! \param elmInt The local integral object for current element
  virtual bool initElementBou(const std::vector<int>& MNPC,
                              const std::vector<size_t>& elem_sizes,
                              const std::vector<size_t>& basis_sizes,
                              LocalIntegral& elmInt);

  //! \brief Initializes current element for boundary integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param elmInt Local integral for element
  virtual bool initElementBou(const std::vector<int>& MNPC,
                              LocalIntegral& elmInt);
  //! \brief Evaluates the integrand at an interior point
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalIntMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
                         const TimeDomain& time, const Vec3& X) const;

  //! \brief Evaluates the integrand at an interior point
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const TimeDomain& time, const Vec3& X) const;

  //! \brief Evaluates the integrand at a boundary point
  //! \param elmInt The local interal object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBouMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
                         const TimeDomain& time, const Vec3& X,
                         const Vec3& normal) const;

  //! \brief Evaluates the integrand at a boundary point
  //! \param elmInt The local interal object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                       const TimeDomain& time, const Vec3& X,
                       const Vec3& normal) const;

  using IntegrandBase::evalSol;
  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s The solution field values at current point
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  virtual bool evalSol(Vector& s, const FiniteElement& fe, const Vec3& X,
		       const std::vector<int>& MNPC) const;
  //! \brief Evaluates the secondary solution at a result point (mixed problem).
  //! \param[out] s The solution field values at current point
  //! \param[in] fe Mixed finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the bases
  //! \param[in] elem_sizes Size of each basis on the element
  virtual bool evalSol(Vector& s, const MxFiniteElement& fe, const Vec3& X,
                       const std::vector<int>& MNPC,
                       const std::vector<size_t>& elem_sizes) const;

  //! \brief Finalizes the element quantities after the numerical integration
  //! \details This method is invoked once for each element, after the numerical
  //! integration loop over interior points is finished and before the resulting
  //! element quantities are assembled into their system level equivalents
  virtual bool finalizeElement(LocalIntegral&, const TimeDomain&, size_t);

  //! \brief Finalizes the element quantities after boundary integration
  //! \details This method is invoked once for each element, after the numerical
  //! integration loop over boundary points is finished and before the resulting
  //! element quantities are assembled into their system level equivalents
  virtual bool finalizeElementBou(LocalIntegral&, const FiniteElement&, const TimeDomain&);

  //! \brief Returns whether a mixed formulation is used
  virtual bool mixedFormulation() const { return true; }

  //! \brief Returns the number of primary/secondary solution field components
  //! \param[in] fld Which field set to consider (1=primary,2=secondary)
  virtual size_t getNoFields(int fld = 1) const;

  //! \brief Returns the name of a primary solution field component
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual std::string getField1Name(size_t i, const char* prefix = 0) const;

  //! \brief Returns the name of a secondary solution field component
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual std::string getField2Name(size_t i, const char* prefix = 0) const;

private:
  //! \brief Evaluates the secondary solution at a result point
  //! (shared code between mixed and non-mixed)
  //! \param[out] s The solution field values at current point
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] disp The displacement coefficients
  bool evalSolCommon(Vector& s,
                     const FiniteElement& fe, const Vec3& X,
                     const Vector& disp) const;

private:
  double sc;   //!< Scaling factor
  double gacc; //!< Gravitational acceleration
};

#endif

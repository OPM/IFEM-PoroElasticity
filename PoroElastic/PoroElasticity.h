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
#include "PoroMaterial.h"


//! \brief Enum for element level solution vectors
enum SolutionVectors
{
  Vu = 0,                        // Displacement
  Vp = 1,                        // Pore pressure

  // Newmark integration stuff
  Vuvel = 2,                    // Displacement velocity
  Vpvel = 3,                    // Pressure velocity
  Vuacc = 4,                    // Displacement acceleration

  NSOL = 5
};


//! \brief Enum for element level right-hand-side vectors
enum ResidualVectors
{
  Fsys = 0,                     // Final RHS vector

  // Sub-vectors, sized according to the bases in question
  Fu = 1,                       // Traction and body forces
  Fp = 2,                       // Flux and body forces

  NVEC = 3
};


//! \brief Enum for element level left-hand-side matrices
enum TangentMatrices
{
  sys = 0,                      // Final Newton matrix

  // Sub-matrices, sized according to the bases in question
  uu_K = 1,                     // Stiffness matrix
  uu_M = 2,                     // Mass matrix
  up_Q = 3,                     // Coupling matrix
  up_D = 4,                     // Dynamic coupling matrix
  pp_S = 5,                     // Compressibility matrix
  pp_P = 6,                     // Permeability matrix

  NMAT = 7
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
    //! \brief Updates the time step size
    //! \param[in] dt New time step size
    void setStepSize(double dt) { h = dt; }
    //! \brief Returns the element level Newton matrix
    virtual const Matrix& getNewtonMatrix() const;
    //! \brief Returns the element level RHS vector
    virtual const Vector& getRHSVector() const;
    //! \brief Adds in a UU-matrix to a system matrix
    virtual void add_uu(const Matrix& source, Matrix& target, double scale = 1.0) const = 0;
    //! \brief Adds in a UP-matrix to a system matrix
    virtual void add_up(const Matrix& source, Matrix& target, double scale = 1.0) const = 0;
    //! \brief Adds in the transpose of a UP-matrix to a system matrix
    virtual void add_pu(const Matrix& source, Matrix& target, double scale = 1.0) const = 0;
    //! \brief Adds in a PP-matrix to a system matrix
    virtual void add_pp(const Matrix& source, Matrix& target, double scale = 1.0) const = 0;
    //! \brief Forms a system vector out of two sub-vectors
    virtual void form_vector(const Vector &u, const Vector &p, Vector& target) const = 0;
  protected:
    size_t ndof_displ, ndof_press, nsd;
    double h;
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
    virtual void add_uu(const Matrix& source, Matrix& target, double scale = 1.0) const;
    //! \brief Adds in a UP-matrix to a system matrix
    virtual void add_up(const Matrix& source, Matrix& target, double scale = 1.0) const;
    //! \brief Adds in the transpose of a UP-matrix to a system matrix
    virtual void add_pu(const Matrix& source, Matrix& target, double scale = 1.0) const;
    //! \brief Adds in a PP-matrix to a system matrix
    virtual void add_pp(const Matrix& source, Matrix& target, double scale = 1.0) const;
    //! \brief Forms a system vector out of two sub-vectors
    virtual void form_vector(const Vector &u, const Vector &p, Vector& target) const;
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
    virtual void add_uu(const Matrix& source, Matrix& target, double scale = 1.0) const;
    //! \brief Adds in a UP-matrix to a system matrix
    virtual void add_up(const Matrix& source, Matrix& target, double scale = 1.0) const;
    //! \brief Adds in the transpose of a UP-matrix to a system matrix
    virtual void add_pu(const Matrix& source, Matrix& target, double scale = 1.0) const;
    //! \brief Adds in a PP-matrix to a system matrix
    virtual void add_pp(const Matrix& source, Matrix& target, double scale = 1.0) const;
    //! \brief Forms a system vector out of two sub-vectors
    virtual void form_vector(const Vector &u, const Vector &p, Vector& target) const;
  };

  /*!
   * \brief Newmark element matrices for PoroElasticity
   */
  template<class P> class NewmarkMats : public P
  {
  public:
    //! \brief Default constructor
    NewmarkMats(size_t ndof_d, size_t ndof_p, bool neumann,
                double b, double c, double md = 0, double sd = 0)
      : P(ndof_d, ndof_p, neumann)
      , beta(b), gamma(c), m_damp(md), s_damp(sd) {}
    //! \brief Empty destructor
    virtual ~NewmarkMats() {}
    //! \brief Returns the element level Newton matrix
    virtual const Matrix& getNewtonMatrix() const
    {
      double gh = gamma * P::h;
      double bh2 = beta * P::h * P::h;

      Matrix& res = const_cast<Matrix&>(P::A.front()); res.fill(0.0);
      this->add_uu(P::A[uu_M], res, 1.0);
      this->add_uu(P::A[uu_K], res, bh2);
      this->add_up(P::A[up_Q], res, -bh2);
      this->add_pu(P::A[up_Q], res, gh);
      this->add_pp(P::A[pp_S], res, gh);
      this->add_pp(P::A[pp_P], res, bh2);

      if (m_damp > 0.0)
        this->add_uu(P::A[uu_M], res, gh * m_damp);
      if (s_damp > 0.0)
        this->add_uu(P::A[uu_K], res, gh * s_damp);

      return P::A.front();
    }
    //! \brief Returns the element level RHS vector
    virtual const Vector& getRHSVector() const
    {
      Vector tu(P::b[Fu]), tp(P::b[Fp]);
      if (P::A.size() > up_Q && P::vec.size() > Vp)
        P::A[up_Q].multiply(P::vec[Vp],    tu, false, 1);
      if (P::A.size() > pp_P && P::vec.size() > Vp)
        P::A[pp_P].multiply(P::vec[Vp],    tp, false,-1);
      if (P::A.size() > uu_M && P::vec.size() > Vuacc)
        P::A[uu_M].multiply(P::vec[Vuacc], tu, false,-1);
      if (P::A.size() > pp_S && P::vec.size() > Vpvel)
        P::A[pp_S].multiply(P::vec[Vpvel], tp, false,-1);
      if (P::A.size() > up_Q && P::vec.size() > Vuvel)
        P::A[up_Q].multiply(P::vec[Vuvel], tp, true, -1);

      if (P::A.size() >> uu_M && P::vec.size() > Vuvel && m_damp > 0.0) {
        Vector auvel(P::vec[Vuvel]);
        auvel *= m_damp;
        P::A[uu_M].multiply(auvel, tu, false, -1);
      }

      if (P::A.size() > uu_K && P::vec.size() > Vuvel && s_damp > 0.0) {
        Vector buvel(P::vec[Vuvel]);
        buvel *= s_damp;
        P::A[uu_K].multiply(buvel, tu, false, -1);
      }

      this->form_vector(tu, tp, const_cast<Vector&>(P::b.front()));
      return P::b.front();
    }
  protected:
    double beta, gamma, m_damp, s_damp;
  };


public:
  //! \brief Default constructor.
  //! \param[in] n Number of spatial dimensions
  PoroElasticity(unsigned short int n = 3);
  //! \brief Empty destructor.
  virtual ~PoroElasticity() {}

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  virtual void setMode(SIM::SolutionMode mode);

  using Elasticity::parseMatProp;
  //! \brief Parses material properties from an XML-element.
  virtual Material* parseMatProp(const TiXmlElement* elem, bool);
  //! \brief Parses a data section from an XML-element.
  virtual bool parse(const TiXmlElement* elem);

  //! \brief Prints out the problem definition to the log stream.
  virtual void printLog() const;

  //! \brief Returns current scaling factor (for unit testing).
  double getScaling(const Vec3& X, const TimeDomain& time) const;

  //! \brief Returns the current gravity vector (for unit testing).
  const Vec3 getGravity() const { return gravity; }

  //! \brief Returns the acceleration due to gravity.
  double getGAcc() const { return gacc; }

  //! \brief Returns a pointer to the material.
  const PoroMaterial* getMaterial() const { return static_cast<PoroMaterial*>(material); }

  //! \brief Returns a local integral contribution object for the given element.
  //! \param[in] nen Number of nodes on element for each basis
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  virtual LocalIntegral* getLocalIntegral(const std::vector<size_t>& nen,
                                          size_t, bool neumann) const;
  //! \brief Returns a local integral contribution object for the given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  virtual LocalIntegral* getLocalIntegral(size_t nen,
                                          size_t, bool neumann) const;

  using Elasticity::initElement;
  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Nodal point correspondence for each basis
  //! \param[in] elem_sizes Size of each basis on the element
  //! \param[in] basis_sizes Size of each basis on the patch level
  //! \param elmInt The local integral object for current element
  virtual bool initElement(const std::vector<int>& MNPC,
                           const std::vector<size_t>& elem_sizes,
                           const std::vector<size_t>& basis_sizes,
                           LocalIntegral& elmInt);

  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param elmInt Local integral for element
  virtual bool initElement(const std::vector<int>& MNPC, LocalIntegral& elmInt);

  //! \brief Initializes current element for boundary integration.
  //! \details Does nothing, solution vectors not needed.
  virtual bool initElementBou(const std::vector<int>&,
                              const std::vector<size_t>&,
                              const std::vector<size_t>&,
                              LocalIntegral&) { return true; }

  //! \brief Initializes current element for boundary integration.
  //! \details Does nothing, solution vectors not needed.
  virtual bool initElementBou(const std::vector<int>&,
                              LocalIntegral&) { return true; }

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
 
  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local interal object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBouMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
                         const TimeDomain&, const Vec3& X,
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
                       const std::vector<size_t>& elem_sizes,
                       const std::vector<size_t>& basis_sizes) const;

  //! \brief Finalizes the element quantities after the numerical integration.
  virtual bool finalizeElement(LocalIntegral&, const TimeDomain&, size_t);

  //! \brief Finalizes the element quantities after boundary integration.
  virtual bool finalizeElementBou(LocalIntegral&, const FiniteElement&,
                                  const TimeDomain&);

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  virtual NormBase* getNormIntegrand(AnaSol*) const;

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

  //! \brief Computes the coupling matrix for a quadrature point.
  bool evalCouplingMatrix(Matrix& mx, const Matrix& B, const Vector& N,
                          double scl) const;

  //! \brief Computes the compressibility matrix for a quadrature point.
  bool evalCompressibilityMatrix(Matrix& mx, const Vector& N, double scl) const;

  //! \brief Computes the permeability matrix for a quadrature point.
  bool evalPermeabilityMatrix(Matrix& mx, const Matrix& dNdX,
                              const Vec3& permeability, double scl) const;

  //! \brief Computes the dynamic coupling matrix for a quadrature point.
  bool evalDynamicCouplingMatrix(Matrix& mx, const Vector& Nu, const Matrix& dNpdx,
                                 const Vec3& permeability, double scl) const;

protected:
  //! \brief Computes the elasticity matrices for a quadrature point.
  //! \param elmInt The element matrix object to receive the contributions
  //! \param[in] B Strain-displacement matrix of current integration point
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalElasticityMatrices(ElmMats& elMat, const Matrix& B,
                                      const FiniteElement& fe,
                                      const Vec3& X) const;

private:
  double sc;   //!< Scaling factor
  double gacc; //!< Gravitational acceleration
};


/*!
 * \brief Class representing the integrand of the poroelasticity norms.
 */
class PoroNorm : public NormBase
{
public:
  //! \brief Initialize members.
  //! \param[in] poroel The poroelasticity problem to evaluate norms for
  //! \param[in] disp Displacement solution (optional)
  //! \param[in] d_disp Derivative of displacement (optional)
  //! \param[in] press Pressure solution (optional)
  //! \param[in] d_press Derivative of pressure (optional)
  PoroNorm(PoroElasticity& poroel,
           VecFunc* disp = 0, TensorFunc* d_disp = 0,
           RealFunc* press = 0, VecFunc* d_press = 0);

  //! Empty destructor
  virtual ~PoroNorm() {}

  //! \brief Returns the number of norm groups or the size of a specified group.
  //! \param[in] group The norm group to return the size of
  //! (if zero, return the number of groups)
  virtual size_t getNoFields(int group = 0) const;

  //! \brief Returns the name of a norm quantity.
  //! \param[in] i The norm group (one-based index)
  //! \param[in] j The norm number (one-based index)
  //! \param[in] prefix Common prefix for all norm names
  virtual std::string getName(size_t i, size_t j, const char* prefix = 0) const;

  //! \brief Evaluates the integrand at an interior point (mixed).
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalIntMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
                         const TimeDomain& time, const Vec3& X) const;

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const TimeDomain& time, const Vec3& X) const;

  //! \brief Initializes the current element for numerical integration (mixed).
  virtual bool initElement(const std::vector<int>& MNPC,
                           const std::vector<size_t>& elem_sizes,
                           const std::vector<size_t>& basis_sizes,
                           LocalIntegral& elmInt);

  //! \brief Initializes the current element for numerical integration.
  virtual bool initElement(const std::vector<int>& MNPC, LocalIntegral& elmInt);

private:
  VecFunc* displacement;
  TensorFunc* d_displacement;
  RealFunc* pressure;
  VecFunc* d_pressure;
};

#endif

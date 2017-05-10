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


/*!
  \brief Class representing the integrand of the PoroElasticity problem.
*/

class PoroElasticity : public Elasticity
{
public:
  //! \brief Enum for element-level solution vectors.
  enum SolutionVectors
  {
    Vu = 0,    // Displacement
    Vp = 1,    // Pore pressure

    // Newmark integration stuff
    Vuvel = 2, // Velocity
    Vpvel = 3, // Pore pressure rate
    Vuacc = 4, // Acceleration

    NSOL = 5
  };

  //! \brief Enum for element-level right-hand-side vectors.
  enum ResidualVectors
  {
    Fsys = 0,  // Final RHS vector

    // Sub-vectors, sized according to the bases in question
    Fu = 1,    // Internal elastic forces, external traction and body forces
    Fp = 2,    // Pressure flux and body forces

    NVEC = 3
  };

  //! \brief Enum for element-level left-hand-side matrices.
  enum TangentMatrices
  {
    Nsys = 0,  // Final Newton matrix

    // Sub-matrices, sized according to the bases in question
    uu_K = 1,  // Stiffness matrix
    uu_M = 2,  // Mass matrix
    up_Q = 3,  // Coupling matrix
    up_D = 4,  // Dynamic coupling matrix
    pp_S = 5,  // Compressibility matrix
    pp_P = 6,  // Permeability matrix

    NMAT = 7
  };

protected:
  //! \brief Superclass for the PoroElasticity element matrices.
  class Mats : public ElmMats
  {
  public:
    //! \brief Default constructor.
    //! \param[in] ndof_displ Number of dofs in displacement
    //! \param[in] ndof_press Number of dofs in pressure
    //! \param[in] neumann Whether or not we are assembling Neumann BCs
    //! \param[in] dynamic Option for dynamic analysis
    //! (0: static analysis, 1: allocate M, 2: allocate M and D)
    Mats(size_t ndof_displ, size_t ndof_press, bool neumann, char dynamic);
    //! \brief Empty destructor.
    virtual ~Mats() {}
    //! \brief Updates the time step size.
    void setStepSize(double dt) { h = dt; }
    //! \brief Updates the perpendicular crack stretch.
    void setCrackStretch(double cs) { lambda = cs; }
    //! \brief Returns the element level Newton matrix.
    virtual const Matrix& getNewtonMatrix() const;
    //! \brief Returns the element level right-hand-side vector.
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
    double h;      //!< Current time step size
    double lambda; //!< Perpendicular crack stretch at current location
  };

private:
  //! \brief Class representing the element matrices for mixed formulation.
  class MixedElmMats : public Mats
  {
  public:
    //! \brief The constructor forwards to the parent class constructor.
    MixedElmMats(size_t ndof_d, size_t ndof_p, bool neumann, char dyn = 0)
      : Mats(ndof_d, ndof_p, neumann, dyn) {}
    //! \brief Empty destructor.
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

  //! \brief Class representing the element matrices for standard formulation.
  class NonMixedElmMats : public Mats
  {
  public:
    //! \brief The constructor forwards to the parent class constructor.
    NonMixedElmMats(size_t ndof_d, size_t ndof_p, bool neumann, char dyn = 0)
      : Mats(ndof_d, ndof_p, neumann, dyn) {}
    //! \brief Empty destructor.
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

  //! \brief Class representing the element matrices for dynamic problems.
  template<class P> class NewmarkMats : public P
  {
  public:
    //! \brief The constructor initializes the time integration parameters.
    NewmarkMats(size_t ndof_d, size_t ndof_p, bool neumann,
                double b, double c, double m = 0.0, double s = 0.0,
                bool useDynCpl = false)
      : P(ndof_d, ndof_p, neumann, useDynCpl ? 2 : 1),
        beta(fabs(b)), gamma(c), m_damp(m), s_damp(s), slvDisp(b < 0.0) {}
    //! \brief Empty destructor.
    virtual ~NewmarkMats() {}
    //! \brief Returns the element-level Newton matrix.
    virtual const Matrix& getNewtonMatrix() const
    {
      Matrix& res = const_cast<Matrix&>(P::A.front()); res.fill(0.0);
      double betah2 = beta * P::h * P::h;
      double gammah = gamma * P::h;
      this->add_uu(P::A[uu_M], res, 1.0    + gammah * m_damp);
      this->add_uu(P::A[uu_K], res, betah2 + gammah * s_damp);
      this->add_up(P::A[up_Q], res,-betah2);
      this->add_pu(P::A[up_Q], res, gammah);
      if (!P::A[up_D].empty())
        this->add_pu(P::A[up_D], res, -1.0);
      this->add_pp(P::A[pp_S], res, gammah);
      this->add_pp(P::A[pp_P], res, betah2);
      if (slvDisp) res.multiply(1.0/betah2);
#if INT_DEBUG > 2
      std::cout <<"\nPoroElasticity::NewmarkMats::getNewtonMatrix:"
                <<"\nElement stiffness matrix, K_uu"<< P::A[uu_K]
                <<"\nElement mass matrix, M_uu"<< P::A[uu_M]
                <<"\nElement coupling matrix, Q_up"<< P::A[up_Q]
                <<"\nElement dynamic coupling matrix, D_up"<< P::A[up_D]
                <<"\nElement compressibility matrix, S_pp"<< P::A[pp_S]
                <<"\nElement permeability matrix, P_pp"<< P::A[pp_P]
                <<"\nElement coefficient matrix"<< P::A[Nsys];
#endif
      return P::A.front();
    }
    //! \brief Returns the element-level right-hand-side vector.
    virtual const Vector& getRHSVector() const
    {
      Vector tu(P::b[Fu]), tp(P::b[Fp]);
      if (P::A.size() > up_Q && P::vec.size() > Vp)
        P::A[up_Q].multiply(P::vec[Vp],    tu, false, 1);
      if (P::A.size() > pp_P && P::vec.size() > Vp)
        P::A[pp_P].multiply(P::vec[Vp],    tp, false,-1);
      if (P::A.size() > uu_M && P::vec.size() > Vuacc)
        P::A[uu_M].multiply(P::vec[Vuacc], tu, false,-1);
      if (P::A.size() > up_D && P::vec.size() > Vuacc && !P::A[up_D].empty())
        P::A[up_D].multiply(P::vec[Vuacc], tp, true,  1);
      if (P::A.size() > pp_S && P::vec.size() > Vpvel)
        P::A[pp_S].multiply(P::vec[Vpvel], tp, false,-1);
      if (P::A.size() > up_Q && P::vec.size() > Vuvel)
        P::A[up_Q].multiply(P::vec[Vuvel], tp, true, -1);
      if (P::A.size() > uu_M && P::vec.size() > Vuvel && m_damp > 0.0)
        P::A[uu_M].multiply(P::vec[Vuvel]*m_damp, tu, false, -1);
      if (P::A.size() > uu_K && P::vec.size() > Vuvel && s_damp > 0.0)
        P::A[uu_K].multiply(P::vec[Vuvel]*s_damp, tu, false, -1);

      this->form_vector(tu, tp, const_cast<Vector&>(P::b.front()));
#if INT_DEBUG > 2
      std::cout <<"\nPoroElasticity::NewmarkMats::getRHSVector:";
      if (P::vec.size() > Vu)
        std::cout <<"Element displacement, Vu"<< P::vec[Vu];
      if (P::vec.size() > Vp)
        std::cout <<"Element pressure, Vp"<< P::vec[Vp];
      if (P::vec.size() > Vuacc)
        std::cout <<"Element acceleration, Vuacc"<< P::vec[Vuacc];
      if (P::vec.size() > Vuvel)
        std::cout <<"Element velocity, Vuvel"<< P::vec[Vuvel];
      if (P::vec.size() > Vpvel)
        std::cout <<"Element pressure rate, Vpvel"<< P::vec[Vpvel];
      std::cout <<"\nElement right-hand-side vector"<< P::b[Fsys];
#endif
      return P::b.front();
    }
  protected:
    double beta;   //!< Time integration parameter
    double gamma;  //!< Time integration parameter
    double m_damp; //!< Mass-proportional damping
    double s_damp; //!< Stiffness-proportional damping
    bool  slvDisp; //!< If \e true, solve for incremental displacements
  };

public:
  //! \brief Default constructor.
  //! \param[in] n Number of spatial dimensions
  PoroElasticity(unsigned short int n = 3);
  //! \brief Empty destructor.
  virtual ~PoroElasticity() {}

  using Elasticity::parseMatProp;
  //! \brief Parses material properties from an XML-element.
  virtual Material* parseMatProp(const TiXmlElement* elem, bool);
  //! \brief Parses a data section from an XML-element.
  virtual bool parse(const TiXmlElement* elem);

  //! \brief Prints out the problem definition to the log stream.
  virtual void printLog() const;

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  virtual void setMode(SIM::SolutionMode mode);

  //! \brief Returns the scaling factor at given location.
  double getScaling(const Vec3& X, double dt = 0.0) const;

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

  using Elasticity::evalIntMx;
  //! \brief Evaluates the integrand at an interior point
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalIntMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
                         const TimeDomain& time, const Vec3& X) const;

  using Elasticity::evalInt;
  //! \brief Evaluates the integrand at an interior point
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const TimeDomain& time, const Vec3& X) const;

  using Elasticity::evalBouMx;
  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local interal object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBouMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
                         const TimeDomain&, const Vec3& X,
                         const Vec3& normal) const;

  using Elasticity::evalSol;
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
                       const std::vector<size_t>&) const;

  using Elasticity::finalizeElement;
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
  void evalPermeabilityMatrix(Matrix& mx, const Matrix& dNdX,
                              const SymmTensor& K, double scl) const;

  //! \brief Computes the dynamic coupling matrix for a quadrature point.
  void evalDynCouplingMatrix(Matrix& mx, const Vector& Nu, const Matrix& dNpdx,
                             const SymmTensor& K, double scl) const;

protected:
  //! \brief Computes the elasticity matrices for a quadrature point.
  //! \param elmInt The element matrix object to receive the contributions
  //! \param[in] B Strain-displacement matrix of current integration point
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalElasticityMatrices(ElmMats& elMat, const Matrix& B,
                                      const FiniteElement& fe,
                                      const Vec3& X) const;

  //! \brief Evaluates the permeability tensor at a quadrature point.
  //! \param[out] K The permeability tensor
  //! \param[in] eV Element solution vectors
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool formPermeabilityTensor(SymmTensor& Kperm,
                                      const Vectors& eV,
                                      const FiniteElement& fe,
                                      const Vec3& X) const;

private:
  double sc; //!< Scaling factor

  bool calculateEnergy; //!< If \e true, perform energy norm calculation
  bool useDynCoupling;  //!< If \e true, include the dynamic coupling matrix

  RealFunc *volumeFlux; //!< Applied volumetric flux

  friend class PoroNorm;
};


/*!
  \brief Class representing the integrand of the poroelasticity norms.
*/

class PoroNorm : public NormBase
{
public:
  //! \brief The only constructor initializes its data members.
  //! \param[in] poroel The poroelasticity problem to evaluate norms for
  //! \param[in] disp Displacement solution (optional)
  //! \param[in] d_disp Derivative of displacement (optional)
  //! \param[in] press Pressure solution (optional)
  //! \param[in] d_press Derivative of pressure (optional)
  PoroNorm(PoroElasticity& poroel,
           VecFunc* disp = nullptr, TensorFunc* d_disp = nullptr,
           RealFunc* press = nullptr, VecFunc* d_press = nullptr);
  //! \brief Empty destructor.
  virtual ~PoroNorm() {}

  //! \brief Returns the number of norm groups or the size of a specified group.
  //! \param[in] group The norm group to return the size of
  //! (if zero, return the number of groups)
  virtual size_t getNoFields(int group = 0) const;

  //! \brief Returns the name of a norm quantity.
  //! \param[in] i The norm group (one-based index)
  //! \param[in] j The norm number (one-based index)
  //! \param[in] prefix Common prefix for all norm names
  virtual std::string getName(size_t i, size_t j, const char* prefix) const;

  using NormBase::evalIntMx;
  //! \brief Evaluates the integrand at an interior point (mixed).
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalIntMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
                         const TimeDomain& time, const Vec3& X) const;

  using NormBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const TimeDomain& time, const Vec3& X) const;

private:
  VecFunc*    displacement;   //!< Analytical displacement field
//  TensorFunc* d_displacement; //!< Analytical displacement gradient
  RealFunc*   pressure;       //!< Analytical pressure field
//  VecFunc*    d_pressure;     //!< Analytical pressure gradient
};

#endif

// $Id$
//==============================================================================
//!
//! \file ThermoPoroElasticity.h
//!
//! \date
//!
//! \author Yared Bekele
//!
//! \brief Integrand implementations for non-isothermal PoroElasticity problems
//!
//==============================================================================

#ifndef _THERMO_POROELASTICITY_H_
#define _THERMO_POROELASTICITY_H_

#include "BDF.h"
#include "Vec3.h"
#include "ElmMats.h"
#include "Elasticity.h"
#include "IntegrandBase.h"
#include "PoroMaterial.h"

/*!
 * \brief Class representing the integrand of a fully
          coupled Thermo-PoroElasticity problem.
*/

class ThermoPoroElasticity : public Elasticity
{
  /*!
   * \brief Class representing an element matrix for a
   *        coupled Thermo-PoroElasticity problem
   */
  class MixedElmMats : public ElmMats
  {
    public:
      //! \brief Default constructor
      MixedElmMats();
      //! \brief Empty destructor
      virtual ~MixedElmMats() {}
      //! \brief Returns the element level Newton matrix
      virtual const Matrix& getNewtonMatrix() const;
      //! \brief Returns the element level RHS vector
      virtual const Vector& getRHSVector() const;
  };

public:
  /*!
   * \brief Class for integrating Robin boundary conditions
  */
  class WeakDirichlet : public IntegrandBase
  {
  public:
    //! \brief Default constructor.
    //! \param[in] n Number of spatial dimensions
    WeakDirichlet(unsigned short int n);

    //! \brief Empty destructor.
    virtual ~WeakDirichlet() {}

    //! \brief Defines the flux function
    void setFlux(RealFunc* f) { flux = f; }
    //! \brief Defines the temperature of the environment.
    void setEnvtTemperature(double envT) { Te = envT; }
    //! \brief Defines the thermal conductivity of the environment.
    void setEnvtConductivity(double envCond) { lambdae = envCond; }

    //! \brief Returns that this integrand has no interior contributions.
    virtual bool hasInteriorTerms() const { return false; }

    //! \brief Returns a local integral container for the given element.
    //! \param[in] nen1 Number of nodes on element for basis 1
    //! \param[in] nen2 Number of nodes on element for basis 2
    //! \param[in] neumann Whether or not we are assembling Neumann BCs
    virtual LocalIntegral* getLocalIntegral(const std::vector<size_t>& nen,
                                            size_t, bool neumann) const;

    //! \brief Initializes current element for numerical boundary integration (mixed)
    //! \param[in] MNPC1 Nodal point correspondence for basis 1
    //! \param[in] MNPC2 Nodal point correspondence for basis 2
    //! \param elmInt The local integral object for current element
    virtual bool initElementBou(const std::vector<int>& MNPC,
                                const std::vector<size_t>& basis_sizes,
                                const std::vector<size_t>& elem_sizes,
                                LocalIntegral& elmInt);

    //! \brief Evaluates the integrands at a boundary point.
    //! \param elmInt The local integral object to receive the contributions
    //! \param[in] fe Finite element data of current integration point
    //! \param[in] time Parameters for nonlineat and time-dependent simulations
    //! \param[in] X Cartesian coordinates of current integration point
    //! \param[in] normal Boundary normal vector at integration point
    virtual bool evalBouMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
                           const TimeDomain& time, const Vec3& X,
                           const Vec3& normal) const;

    //! \brief Finalizes the element quantities after the numerical integration.
    //! \details This method is invoked once for each element, afte the numerical
    //! integration loop over interior points is finished and before the resulting
    //! element quantities are assembled into their system level equivalents
    //virtual bool finalizeElement(LocalIntegral&, const TimeDomain&, size_t);

  protected:
    unsigned short int nsd;       //!< Number of space dimensions
    RealFunc* flux;               //!< Flux function
    double Te;                    //!< Temperature of environment
    double lambdae;               //!< Thermal conductivity of environment
  };

  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] order The order of the time integration
  ThermoPoroElasticity(unsigned short int n = 2, int order = 1, bool stab = false);

  //! \brief The destructor frees the dynamically allocated data objects.
  virtual ~ThermoPoroElasticity() {}

  //! \brief Defines the traction field to use in Neumann boundary conditions
  void setTraction(TractionFunc* tf) { tracFld = tf; }

  //! \brief Defines the trction field to use in Neumann boundary conditions
  void setTraction(VecFunc* tf) { fluxFld = tf; }

  //! \brief Defines the scaling values between U and P and U and T
  void setScalingValues(double scl_u_p, double scl_u_T)
    { scl1 = scl_u_p; scl2 = scl_u_T; }

  //! \brief Evaluates the boundary traction field (if any) at specified point
  //! \param[in] X Cartesian coordinate of the current integration point
  //! \param[in] n Outward-directed unit normal vector at current point
  virtual Vec3 getTraction(const Vec3& X, const Vec3& n) const;

  //! \brief Defines the material properties.
  virtual void setMaterial(PoroMaterial* material) { mat = material; }

  virtual int getIntegrandType() const
  { return SUPG ? ELEMENT_CORNERS | SECOND_DERIVATIVES : STANDARD; }

  //! \brief Returns a local integral container for the given element
  //! \param[in] nen1 Number of nodes on element for basis 1
  //! \param[in] nen2 Number of nodes on element for basis 2
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  virtual LocalIntegral* getLocalIntegral(const std::vector<size_t>& nen,
                                          size_t, bool neumann) const;

  //! \brief Initializes current element for numerical integration
  //! \param[in] MNPC1 Nodal point correspondence for basis 1
  //! \param[in] MNPC2 Nodal point correspondence for basis 2
  //! \param[in] n1 Number of nodes in basis 1 on this patch
  //! \param elmInt The local integral object for current element
  virtual bool initElement(const std::vector<int>& MNPC,
                           const std::vector<size_t>& basis_sizes,
                           const std::vector<size_t>& elem_sizes,
                           LocalIntegral& elmInt);

  //! \brief Initializes current element for numerical boundary integration (mixed)
  //! \param[in] MNPC1 Nodal point correspondence for basis 1
  //! \param[in] MNPC2 Nodal point correspondence for basis 2
  //! \param elmInt The local integral object for current element
  virtual bool initElementBou(const std::vector<int>& MNPC,
                           const std::vector<size_t>& basis_sizes,
                           const std::vector<size_t>& elem_sizes,
                           LocalIntegral& elmInt);

  //! \brief Evaluates the integrand at an interior point
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalIntMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
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

  //! \brief Evaluates the secondary solution at a result point (mixed problem).
  //! \param[out] s The solution field values at current point
  //! \param[in] fe Mixed finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC1 Nodal point correspondance for the basis 1
  //! \param[in] MNPC2 Nodal point correspondance for the basis 2
  virtual bool evalSol(Vector& s, const MxFiniteElement& fe, const Vec3& X,
                       const std::vector<int>& MNPC1,
                       const std::vector<size_t>& elem_sizes) const;

  //! \brief Finalizes the element quantities after the numerical integration
  //! \details This method is invoked once for each element, after the numerical
  //! integration loop over interior points is finished and before the resulting
  //! element quantities are assembled into their system level equivalents
  virtual bool finalizeElement(LocalIntegral&, const TimeDomain&, size_t);

  //! \brief Returns whether a mixed formulation is used
  //virtual bool mixedFormulation() const { return true; }

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

protected:
  unsigned short int eS;     //!< Index to element load vector
  double gacc;               //!< Gravitational acceleration
  double scl1;               //!< Scaling value for displacement and pressure
  double scl2;               //!< Scaling value for displacement and temperature
  TractionFunc* tracFld;     //!< Pointer to implicit boundary traction field
  VecFunc* fluxFld;          //!< Pointer to explicit boundary traction field
  TimeIntegration::BDF bdf;  //!< BDF time discretization parameters
  PoroMaterial* mat;         //!< Material data
  bool SUPG; //!< SUPG for advection. Currently broken
};

#endif  // _THERMO_POROELASTICITY_H_

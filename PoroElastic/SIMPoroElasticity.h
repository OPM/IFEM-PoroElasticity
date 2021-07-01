// $Id$
//==============================================================================
//!
//! \file SIMPoroElasticity.h
//!
//! \date April 16 2015
//!
//! \author Yared Bekele
//!
//! \brief Simulation driver for poroelasticity problems.
//!
//==============================================================================

#ifndef _SIM_PORO_ELASTICITY_H_
#define _SIM_PORO_ELASTICITY_H_

#include "SIMElasticityWrap.h"


/*!
  \brief Driver class for poroelastic simulators.
*/

template<class Dim> class SIMPoroElasticity : public SIMElasticityWrap<Dim>
{
public:
  //! \brief The default constructor sets the solution dimension for each basis.
  SIMPoroElasticity();

  //! \brief Empty destructor.
  virtual ~SIMPoroElasticity() {}

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  virtual std::string getName() const { return "PoroElasticity"; }

  //! \brief Computes the solution for the current time step.
  virtual bool solveStep(TimeStep& tp);

  //! \brief Evaluates some iteration norms for convergence assessment.
  //! \param[in] u Global primary solution vector
  //! \param[in] r Global residual vector associated with the solution vector
  //! \param[out] eNorm Energy norm of solution increment
  //! \param[out] rNorm Residual norm of solution increment
  //! \param[out] dNorm Displacement norm of solution increment
  virtual void iterationNorms(const Vector& u, const Vector& r,
                              double& eNorm, double& rNorm, double& dNorm) const;

  //! \brief Prints a summary of the calculated solution to std::cout.
  virtual void printSolutionSummary(const Vector& solution, int = 0,
                                    const char* = nullptr, std::streamsize = 0);

  //! \brief Computes energy norms on the converged solution.
  bool postSolve(TimeStep& tp);

  //! \brief Initializes for integration of Neumann terms for a given property.
  virtual bool initNeumann(size_t propInd);

protected:
  //! \brief Returns the actual integrand.
  virtual Elasticity* getIntegrand();

  using SIMElasticityWrap<Dim>::parse;
  //! \brief Parses a data section from an XML element
  virtual bool parse(const TiXmlElement* elem);

  using SIMElasticityWrap<Dim>::parseAnaSol;
  //! \brief Parses the analytical solution from an XML element.
  virtual bool parseAnaSol(const TiXmlElement* elem);

  //! \brief Parses dimension-specific analytical solution from an XML element.
  bool parseDimSpecific(const TiXmlElement*) { return false; }

private:
  double scaleD; //!< Displacement DOF scaling in convergence checks
  double scaleP; //!< Pressure DOF scaling in convergence checks
};

#endif

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
  explicit SIMPoroElasticity(const std::vector<unsigned char>& = {});

  //! \brief Empty destructor.
  virtual ~SIMPoroElasticity() {}

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  std::string getName() const override { return "PoroElasticity"; }

  //! \brief Computes the solution for the current time step.
  bool solveStep(TimeStep& tp) override;

  //! \brief Evaluates some iteration norms for convergence assessment.
  //! \param[in] u Global primary solution vector
  //! \param[in] r Global residual vector associated with the solution vector
  //! \param[out] eNorm Energy norm of solution increment
  //! \param[out] rNorm Residual norm of solution increment
  //! \param[out] dNorm Displacement norm of solution increment
  void iterationNorms(const Vector& u, const Vector& r,
                      double& eNorm, double& rNorm, double& dNorm) const override;

  //! \brief Prints a summary of the calculated solution to std::cout.
  void printSolutionSummary(const Vector& solution, int = 0,
                            const char* = nullptr, std::streamsize = 0) override;

  //! \brief Computes energy norms on the converged solution.
  bool postSolve(TimeStep& tp);

  //! \brief Initializes for integration of Neumann terms for a given property.
  bool initNeumann(size_t propInd) override;

protected:
  //! \brief Returns the actual integrand.
  Elasticity* getIntegrand() override;

  using SIMElasticityWrap<Dim>::parse;
  //! \brief Parses a data section from an XML element
  bool parse(const TiXmlElement* elem) override;

  using SIMElasticityWrap<Dim>::parseAnaSol;
  //! \brief Parses the analytical solution from an XML element.
  bool parseAnaSol(const TiXmlElement* elem) override;

  //! \brief Parses dimension-specific analytical solution from an XML element.
  bool parseDimSpecific(const TiXmlElement*) { return false; }

private:
  double scaleD; //!< Displacement DOF scaling in convergence checks
  double scaleP; //!< Pressure DOF scaling in convergence checks
};

#endif

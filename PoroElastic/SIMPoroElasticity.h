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
#include "SIM2D.h"
#include "SIM3D.h"
#include "PoroElasticity.h"
#include "ASMmxBase.h"


/*!
  \brief Driver class for poroelastic simulators.
*/

template<class Dim> class SIMPoroElasticity : public SIMElasticityWrap<Dim>
{
public:
  //! \brief The default constructor sets the solution dimension for each basis.
  SIMPoroElasticity()
  {
    if (ASMmxBase::Type > ASMmxBase::NONE)
      Dim::nf = { Dim::dimension, 1 }; // mixed formulation
    else
      Dim::nf = { Dim::dimension+1 }; // standard formulation

    Dim::myHeading = "Poroelasticity solver";
    SIMElasticity<Dim>::myContext = "poroelasticity";
  }

  //! \brief Empty destructor.
  virtual ~SIMPoroElasticity() {}

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  virtual std::string getName() const { return "PoroElasticity"; }

  //! \brief Computes the solution for the current time step.
  virtual bool solveStep(TimeStep& tp)
  {
    this->setMode(SIM::STATIC);
    this->setQuadratureRule(Dim::opt.nGauss[0],true);
    if (!this->assembleSystem(tp.time,SIMsolution::solution))
      return false;

    if (!this->solveSystem(SIMsolution::solution.front(),Dim::msgLevel-1))
      return false;

    this->printSolutionSummary(SIMsolution::solution.front());

    return this->postSolve(tp);
  }

  //! \brief Prints a summary of the calculated solution to std::cout.
  virtual void printSolutionSummary(const Vector& solution, int = 0,
                                    const char* = nullptr, std::streamsize = 0)
  {
    const size_t nsd = this->getNoSpaceDim();
    size_t iMax[nsd+1];
    double dMax[nsd+1];
    double dNorm = this->solutionNorms(solution,dMax,iMax,this->getNoFields(1));
    double pNorm = 0.0;
    if (this->getNoFields(2) > 0)
      pNorm = this->solutionNorms(solution,dMax+nsd,iMax+nsd,
                                  this->getNoFields(2),'P');

    IFEM::cout <<"  Primary solution summary: L2-norm            : "
               << utl::trunc(dNorm);
    if (pNorm != 0.0)
      IFEM::cout <<"\n                   Pressure L2-norm            : "<< pNorm;

    char D = 'X';
    for (size_t d = 0; d < nsd; d++, D++)
      if (utl::trunc(dMax[d]) != 0.0)
        IFEM::cout <<"\n                            Max "<< char('X'+d)
                   <<"-displacement : "<< dMax[d] <<" node "<< iMax[d];
    if (utl::trunc(dMax[nsd]) != 0.0)
      IFEM::cout <<"\n                            Max pressure       : "
                 << dMax[nsd] <<" node "<< iMax[nsd] <<"\n";
  }

  //! \brief Computes energy norms on the converged solution.
  bool postSolve(TimeStep& tp)
  {
    NormBase* norm = this->getNormIntegrand();
    if (!norm) return true;

    Vectors gNorms;
    this->setMode(SIM::RECOVERY);
    this->setQuadratureRule(Dim::opt.nGauss[1]);
    bool ok = this->solutionNorms(tp.time,SIMsolution::solution,gNorms);
    if (ok && !gNorms.empty())
      for (size_t i = 1; i <= gNorms.front().size(); i++)
        if (utl::trunc(gNorms.front()(i)) != 0.0)
          IFEM::cout << utl::adjustRight(33,norm->getName(1,i))
                     << gNorms.front()(i) << std::endl;

    delete norm;
    return ok;
  }

protected:
  //! \brief Returns the actual integrand.
  virtual Elasticity* getIntegrand()
  {
    if (!Dim::myProblem)
      Dim::myProblem = new PoroElasticity(Dim::dimension, Dim::nf.size() > 1);
    return static_cast<Elasticity*>(Dim::myProblem);
  }

  using SIMElasticityWrap<Dim>::parseDimSpecific;
  //! \brief Parses a dimension-specific data section from an XML element.
  virtual bool parseDimSpecific(const TiXmlElement*);
};


typedef SIMPoroElasticity<SIM2D> SIMPoroEl2D; //!< 2D specific driver
typedef SIMPoroElasticity<SIM3D> SIMPoroEl3D; //!< 3D specific driver

//! \brief Template specialization - 2D specific input parsing.
template<> bool SIMPoroEl2D::parseDimSpecific(const TiXmlElement* elem);
//! \brief Template specialization - 3D specific input parsing.
template<> bool SIMPoroEl3D::parseDimSpecific(const TiXmlElement* elem);

#endif

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
#include "SAM.h"
#include "PoroElasticity.h"
#include "ASMmxBase.h"
#include "tinyxml.h"


/*!
  \brief Driver class for poroelastic simulators.
*/

template<class Dim> class SIMPoroElasticity : public SIMElasticityWrap<Dim>
{
public:
  //! \brief The default constructor sets the solution dimension for each basis.
  SIMPoroElasticity(const std::vector<unsigned char>& = {})
  {
    if (ASMmxBase::Type > ASMmxBase::NONE)
      Dim::nf = { Dim::dimension, 1 }; // mixed formulation
    else
      Dim::nf = { Dim::dimension+1 }; // standard formulation

    Dim::myHeading = "Poroelasticity solver";
    SIMElasticity<Dim>::myContext = "poroelasticity";
    scaleD = scaleP = 0.0;
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

  //! \brief Evaluates some iteration norms for convergence assessment.
  //! \param[in] u Global primary solution vector
  //! \param[in] r Global residual vector associated with the solution vector
  //! \param[out] eNorm Energy norm of solution increment
  //! \param[out] rNorm Residual norm of solution increment
  //! \param[out] dNorm Displacement norm of solution increment
  virtual void iterationNorms(const Vector& u, const Vector& r,
                              double& eNorm, double& rNorm, double& dNorm) const
  {
    if (scaleD > 0.0)
    {
      eNorm = Dim::mySam->dot(r,u,'D')*scaleD;
      rNorm = Dim::mySam->norm2(r,'D')*scaleD;
      dNorm = Dim::mySam->norm2(u,'D')*scaleD;
    }
    else if (scaleP > 0.0)
      eNorm = rNorm = dNorm = 0.0;
    else
      this->Dim::iterationNorms(u,r,eNorm,rNorm,dNorm);

    if (scaleP > 0.0)
    {
      eNorm += Dim::mySam->dot(r,u,'P')*scaleP;
      rNorm += Dim::mySam->norm2(r,'P')*scaleP;
      dNorm += Dim::mySam->norm2(u,'P')*scaleP;
    }

    if (scaleD+scaleP > 0.0)
    {
      eNorm /= (scaleD+scaleP);
      rNorm /= (scaleD+scaleP);
      dNorm /= (scaleD+scaleP);
    }
  }

  //! \brief Prints a summary of the calculated solution to std::cout.
  virtual void printSolutionSummary(const Vector& solution, int = 0,
                                    const char* = nullptr, std::streamsize = 0)
  {
    const size_t nsd = this->getNoSpaceDim();
    size_t iMax[nsd+1];
    double dMax[nsd+1];
    double dNorm = this->solutionNorms(solution,dMax,iMax,nsd);
    double pNorm = this->solutionNorms(solution,dMax+nsd,iMax+nsd,1,'P');

    IFEM::cout <<"  Primary solution summary: L2-norm            : "
               << utl::trunc(dNorm)
               <<"\n                   Pressure L2-norm            : "<< pNorm;

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
    bool ok = this->solutionNorms(tp.time,this->getSolutions(),gNorms);
    if (ok && !gNorms.empty())
      for (size_t i = 1; i <= gNorms.front().size(); i++)
        if (utl::trunc(gNorms.front()(i)) != 0.0)
          IFEM::cout << utl::adjustRight(33,norm->getName(1,i))
                     << gNorms.front()(i) << std::endl;

    delete norm;
    return ok;
  }

  //! \brief Initializes for integration of Neumann terms for a given property.
  virtual bool initNeumann(size_t propInd)
  {
    PoroElasticity* prob = dynamic_cast<PoroElasticity*>(Dim::myProblem);
    if (!prob) return false;

    // If we find a real function, then it's a flux boundary term
    typename Dim::SclFuncMap::const_iterator rit = Dim::myScalars.find(propInd);
    if (rit != Dim::myScalars.end())
    {
      prob->setBoundaryFlux(rit->second);
      return true;
    }

    // If not, let parent class try to find a traction or vector term
    prob->setBoundaryFlux(nullptr);
    return this->SIMElasticityWrap<Dim>::initNeumann(propInd);
  }

protected:
  //! \brief Returns the actual integrand.
  virtual Elasticity* getIntegrand()
  {
    if (!Dim::myProblem)
      Dim::myProblem = new PoroElasticity(Dim::dimension, Dim::nf.size() > 1);
    return static_cast<Elasticity*>(Dim::myProblem);
  }

  using SIMElasticityWrap<Dim>::parse;
  //! \brief Parses a data section from an XML element
  virtual bool parse(const TiXmlElement* elem)
  {
    if (!strcasecmp(elem->Value(),"poroelasticity"))
    {
      const TiXmlElement* child = elem->FirstChildElement();
      for (; child; child = child->NextSiblingElement())
        if (!strcasecmp(child->Value(),"normscaling"))
        {
          utl::getAttribute(child,"displacement",scaleD);
          utl::getAttribute(child,"pressure",scaleP);
          IFEM::cout <<"\tDisplacement scaling: "<< scaleD
                     <<"\n\tPressure scaling: "<< scaleP << std::endl;
        }
    }
    return this->SIMElasticityWrap<Dim>::parse(elem);
  }

  using SIMElasticityWrap<Dim>::parseAnaSol;
  //! \brief Parses the analytical solution from an XML element.
  virtual bool parseAnaSol(const TiXmlElement* elem)
  {
    if (Dim::mySol || this->parseDimSpecific(elem))
      return true;

    IFEM::cout <<"\tAnalytical solution: Expression"<< std::endl;
    Dim::mySol = new AnaSol(elem);
    return true;
  }

  //! \brief Parses dimension-specific analytical solution from an XML element.
  bool parseDimSpecific(const TiXmlElement*) { return false; }

private:
  double scaleD; //!< Displacement DOF scaling in convergence checks
  double scaleP; //!< Pressure DOF scaling in convergence checks
};


typedef SIMPoroElasticity<SIM2D> SIMPoroEl2D; //!< 2D specific driver

//! \brief Template specialization - 2D specific input parsing.
template<> bool SIMPoroEl2D::parseDimSpecific(const TiXmlElement* elem);

#endif

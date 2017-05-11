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
#ifdef HAS_CEREAL
#include <cereal/cereal.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#endif


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

  //! \brief Initializes the solution vectors.
  virtual bool init(const TimeStep&)
  {
    bool ok = this->setMode(SIM::STATIC);

    solution.resize(this->getNoSolutions());
    for (size_t i = 0; i < solution.size(); i++)
      solution[i].resize(this->getNoDOFs(),true);

    this->setQuadratureRule(Dim::opt.nGauss[0],true);
    return ok;
  }

  //! \brief Advances the time step one step forward.
  virtual bool advanceStep(TimeStep& tp)
  {
    // Update vectors between time steps
    const int nNusols = solution.size();
    for (int n = nNusols-1; n > 0; n--)
      solution[n] = solution[n-1];

    return this->SIMElasticity<Dim>::advanceStep(tp);
  }

  //! \brief Computes the solution for the current time step.
  virtual bool solveStep(TimeStep& tp)
  {
    if (Dim::msgLevel >= 0)
      IFEM::cout <<"\n  step = "<< tp.step
                 <<"  time = "<< tp.time.t << std::endl;

    if (!this->assembleSystem(tp.time,solution))
      return false;

    if (!this->solveSystem(solution.front()))
      return false;

    this->printSolutionSummary(solution.front());

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
    this->setQuadratureRule(Dim::opt.nGauss[1]);
    bool ok = this->solutionNorms(tp.time,solution,gNorms);
    if (ok && !gNorms.empty())
      for (size_t i = 1; i < gNorms.front().size(); i++)
        if (utl::trunc(gNorms.front()(i)) != 0.0)
          IFEM::cout << utl::adjustRight(36,norm->getName(1,i))
                     << gNorms.front()(i) << std::endl;

    delete norm;
    return ok;
  }

  //! \brief Serialize internal state for restarting purposes.
  //! \param data Container for serialized data
  bool serialize(DataExporter::SerializeData& data)
  {
#ifdef HAS_CEREAL
    std::ostringstream str;
    cereal::BinaryOutputArchive ar(str);
    for (size_t i = 0; i < solution.size(); ++i)
      ar(solution[i]);
    data.insert(std::make_pair(this->getName(), str.str()));
    return true;
#else
    return false;
#endif
  }

  //! \brief Set internal state from a serialized state.
  //! \param[in] data Container for serialized data
  bool deSerialize(const DataExporter::SerializeData& data)
  {
#ifdef HAS_CEREAL
    std::stringstream str;
    auto it = data.find(this->getName());
    if (it != data.end()) {
      str << it->second;
      cereal::BinaryInputArchive ar(str);
      for (size_t i = 0; i < solution.size(); ++i)
        ar(solution[i]);
      return true;
    }
#endif
    return false;
  }

protected:
  //! \brief Returns the actual integrand.
  virtual Elasticity* getIntegrand()
  {
    if (!Dim::myProblem)
      Dim::myProblem = new PoroElasticity(Dim::dimension);
    return static_cast<Elasticity*>(Dim::myProblem);
  }

  using SIMElasticityWrap<Dim>::parseDimSpecific;
  //! \brief Parses a dimension-specific data section from an XML element.
  virtual bool parseDimSpecific(const TiXmlElement*);

  //! \brief Returns a const reference to current solution vector.
  virtual const Vector& getSolution(int idx = 0) const { return solution[idx]; }

private:
  Vectors solution; //!< Solution vectors
};


typedef SIMPoroElasticity<SIM2D> SIMPoroEl2D; //!< 2D specific driver
typedef SIMPoroElasticity<SIM3D> SIMPoroEl3D; //!< 3D specific driver

//! \brief Template specialization - 2D specific input parsing.
template<> bool SIMPoroEl2D::parseDimSpecific(const TiXmlElement* elem);
//! \brief Template specialization - 3D specific input parsing.
template<> bool SIMPoroEl3D::parseDimSpecific(const TiXmlElement* elem);

#endif

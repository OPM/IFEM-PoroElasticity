// $Id$
//==============================================================================
//!
//! \file SIMDynPoroElasticity.h
//!
//! \date Apr 23 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Dynamic simulation driver for poroelasticity problems.
//!
//==============================================================================

#ifndef _SIM_DYN_PORO_ELASTICITY_H_
#define _SIM_DYN_PORO_ELASTICITY_H_

#include "SIMPoroElasticity.h"


/*!
  \brief Driver class for dynamic poroelasticity problems.
*/

template<class Dim, class DynSIM>
class SIMDynPoroElasticity : public SIMPoroElasticity<Dim>
{
public:
  //! \brief Default constructor.
  SIMDynPoroElasticity() : dSim(*this) {}

  //! \brief Constructor for mixed problems.
  explicit SIMDynPoroElasticity(const std::vector<unsigned char>& flds)
    : SIMPoroElasticity<Dim>(flds), dSim(*this) {}

  //! \brief Empty destructor.
  virtual ~SIMDynPoroElasticity() {}

  //! \brief Prints out problem-specific data to the log stream.
  virtual void printProblem() const
  {
    static short int ncall = 0;
    if (++ncall == 1) // Avoiding infinite recursive calls
      dSim.printProblem();
    else
      this->SIMPoroElasticity<Dim>::printProblem();
    --ncall;
  }

  //! \brief Initializes the problem.
  virtual bool init(const TimeStep& tp, bool withRF)
  {
    if (!this->initSystem(Dim::opt.solver,1,1,0,withRF))
      return false;

    dSim.initPrm();
    dSim.initSol(3);

    bool ok = this->setMode(SIM::INIT) && this->getIntegrand()->init(tp.time);
    this->setQuadratureRule(Dim::opt.nGauss[0],true);
    return ok;
  }

  //! \brief Advances the time step one step forward.
  virtual bool advanceStep(TimeStep& tp) { return dSim.advanceStep(tp,false); }

  //! \brief Computes the solution for the current time step.
  virtual bool solveStep(TimeStep& tp)
  {
    if (dSim.solveStep(tp) != SIM::CONVERGED)
      return false;

    if (this->getNoFields(2) > 0)
    {
      // Calculate and print the pressure norms in case of mixed problem.
      // The solution norms involving displacement variables
      // are printed by the dynamic solution driver.
      size_t iMax = 0;
      double pMax = 0.0;
      double pNorm = this->solutionNorms(dSim.getSolution(),&pMax,&iMax,
                                         this->getNoFields(2),'P');
      IFEM::cout <<"  Pressure L2-norm                : "<< pNorm
                 <<"\n               Max pressure       : "<< pMax
                 <<" node "<< iMax << std::endl;
    }

    return this->postSolve(tp);
  }

  //! \brief Solves the linearized system of current iteration.
  SIM::ConvStatus solveIteration(TimeStep& p) { return dSim.solveIteration(p); }

  //! \brief Returns the maximum number of iterations.
  int getMaxit() const { return dSim.getMaxit(); }

  //! \brief Returns a const reference to current solution vector.
  virtual const Vector& getSolution(int i) const { return dSim.getSolution(i); }

  //! \brief Returns a const reference to the solution vectors.
  virtual const Vectors& getSolutions() const { return dSim.getSolutions(); }

protected:
  //! \brief Returns a reference to the solution vectors (for assignment).
  virtual Vectors& theSolutions() { return dSim.theSolutions(); }

  using SIMPoroElasticity<Dim>::parse;
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const tinyxml2::XMLElement* elem)
  {
    bool result = true;
    static short int ncall = 0;
    if (++ncall == 1) // Avoiding infinite recursive calls
      result = dSim.parse(elem);
    else
      result = this->SIMPoroElasticity<Dim>::parse(elem);
    --ncall;
    return result;
  }

private:
  DynSIM dSim; //!< Dynamic solution driver
};

#endif

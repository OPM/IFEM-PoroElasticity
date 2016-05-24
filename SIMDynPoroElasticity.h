// $Id$
//==============================================================================
//!
//! \file SIMDynPoroElasticity.h
//!
//! \date Apr 23 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Dynamic simulation driver for poroelasticity problems with.
//!
//==============================================================================

#ifndef _SIM_DYN_PORO_ELASTICITY_H_
#define _SIM_DYN_PORO_ELASTICITY_H_

#include "NewmarkSIM.h"
#include "SIMPoroElasticity.h"


/*!
  \brief Driver class for dynamic elasticity problems with fracture.
*/

template<class Dim, class DynSIM=NewmarkSIM>
class SIMDynPoroElasticity : public SIMPoroElasticity<Dim>
{
public:
  //! \brief Default constructor.
  SIMDynPoroElasticity() : dSim(*this) {}

  //! \brief Constructor for mixed problems.
  SIMDynPoroElasticity(const std::vector<unsigned char>& flds)
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
  virtual bool init(const TimeStep&)
  {
    dSim.initPrm();
    dSim.initSol(3);

    bool ok = this->setMode(SIM::INIT);
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

    if (SIMPoroElasticity<Dim>::doPrintNorms) {
      Matrix eNorm;
      Vectors gNorm;
      this->solutionNorms(tp.time, Vectors(1, dSim.getSolution()), gNorm, &eNorm);
      this->printNorms(gNorm);
    }

    this->printSolutionSummary(dSim.getSolution());
    return true;
  }

  //! \brief Solves the linearized system of current iteration.
  //! \param[in] tp Time stepping parameters
  SIM::ConvStatus solveIteration(TimeStep& tp)
  {
    return dSim.solveIteration(tp);
  }

  //! \brief Returns the maximum number of iterations.
  int getMaxit() const { return dSim.getMaxit(); }

  //! \brief Returns a const reference to current solution vector.
  virtual const Vector& getSolution(int i) const { return dSim.getSolution(i); }

protected:
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem)
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

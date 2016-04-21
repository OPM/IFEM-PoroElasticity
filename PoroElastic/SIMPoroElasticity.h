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

#include "SIMElasticity.h"
#include "PoroElasticity.h"
#include "DataExporter.h"


/*!
  \brief Driver class for poroelastic simulators.
*/

template<class Dim> class SIMPoroElasticity : public SIMElasticity<Dim>
{
public:
  //! \brief Default constructor.
  SIMPoroElasticity()
  {
    Dim::msgLevel = 1;
    Dim::myHeading = "Elasticity solver";
  }

  //! \brief This constructor sets the number of solution fields for each basis.
  SIMPoroElasticity(const std::vector<unsigned char>& fields)
  {
    Dim::nf = fields;
    Dim::msgLevel = 1;
    Dim::myHeading = "Poroelasticity solver";
    SIMElasticity<Dim>::myContext = "poroelasticity";
  }

  //! \brief Empty destructor.
  virtual ~SIMPoroElasticity() {}

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  virtual std::string getName() const { return "PoroElasticity"; }

  //! \brief Registers solution fields for data output.
  void registerFields(DataExporter& exporter)
  {
    int flag = DataExporter::PRIMARY;
    if (!Dim::opt.pSolOnly)
      flag |= DataExporter::SECONDARY;
    exporter.registerField("u","solution",DataExporter::SIM,flag);
    exporter.setFieldValue("u",this,&this->getSolution());
  }

  //! \brief Initializes the solution vectors.
  virtual bool init(const TimeStep&)
  {
    solution.resize(this->getNoSolutions());
    for (size_t i = 0; i < solution.size(); i++)
      solution[i].resize(this->getNoDOFs(),true);

    this->registerField("init1", solution[0]);
    this->registerField("init2", solution[1]);

    bool ok = this->setMode(SIM::STATIC);
    this->setQuadratureRule(Dim::opt.nGauss[0],true);
    return ok;
  }

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  //! \param geoBlk Running geometry block counter
  bool saveModel(char* fileName, int& geoBlk, int&)
  {
    return Dim::opt.format < 0 ? true : this->writeGlvG(geoBlk,fileName);
  }

  //! \brief Saves the converged results of a given time step to VTF file.
  //! \param[in] tp Time stepping parameters
  //! \param nBlock Running result block counter
  virtual bool saveStep(const TimeStep& tp, int& nBlock)
  {
    if (tp.step%Dim::opt.saveInc > 0 || Dim::opt.format < 0)
      return true;

    int iDump = 1 + tp.step/Dim::opt.saveInc;
    if (!this->writeGlvS(solution.front(),iDump,nBlock,tp.time.t,"vector",89))
      return false;

    return this->writeGlvStep(iDump,tp.time.t);
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

    if (!this->assembleSystem(tp.time, solution))
      return false;

    if (!this->solveSystem(solution.front()))
      return false;

    this->printSolutionSummary(solution.front());
    return true;
  }

  //! \brief Prints a summary of the calculated solution to std::cout.
  virtual void printSolutionSummary(const Vector& solution, int = 0,
                                    const char* = nullptr, std::streamsize = 0)
  {
    if (SIMElasticity<Dim>::myContext != "poroelasticity")
    {
      this->Dim::printSolutionSummary(solution);
      return;
    }

    const size_t nsd = this->getNoSpaceDim();
    size_t iMax[nsd+1];
    double dMax[nsd+1];
    double dNorm = this->solutionNorms(solution,dMax,iMax,this->getNoFields(1));
    double pNorm = 0.0;
    if (this->getNoFields(2) > 0)
      pNorm = this->solutionNorms(solution,dMax+nsd,iMax+nsd,
                                  this->getNoFields(2),'P');

    std::stringstream str;
    if (this->adm.getProcId() == 0)
    {
      str <<"  Primary solution summary: L2-norm            : "
          << utl::trunc(dNorm);
      if (pNorm != 0.0)
        str <<"\n                   Pressure L2-norm            : "<< pNorm;

      char D = 'X';
      for (size_t d = 0; d < nsd; d++, D++)
        if (utl::trunc(dMax[d]) != 0.0)
          str <<"\n                            Max "<< char('X'+d)
              <<"-displacement : "<< dMax[d] <<" node "<< iMax[d];
      if (utl::trunc(dMax[nsd]) != 0.0)
        str <<"\n                            Max pressure       : "
            << dMax[nsd] <<" node "<< iMax[nsd] <<"\n";
    }

    utl::printSyncronized(std::cout,str,this->adm.getProcId());
  }

protected:
  //! \brief Returns the actual integrand.
  virtual Elasticity* getIntegrand()
  {
    if (!Dim::myProblem)
      Dim::myProblem = new PoroElasticity(Dim::dimension);
    return static_cast<Elasticity*>(Dim::myProblem);
  }

  //! \brief Returns a const reference to current solution vector.
  virtual const Vector& getSolution(int idx = 0) const { return solution[idx]; }

private:
  Vectors solution; //!< Solution vectors
};

#endif

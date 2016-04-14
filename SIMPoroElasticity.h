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
#include "PoroMaterial.h"
#include "DataExporter.h"
#include "InitialConditionHandler.h"
#include "SIMSolver.h"
#include "Vec3Oper.h"
#include "Profiler.h"


/*!
  \brief Driver class for poro-elastic simulators.
*/

template<class Dim> class SIMPoroElasticity : public SIMElasticity<Dim>
{
public:
  //! \brief Dummy declaration, no setup properties needed
  typedef bool SetupProps;

  //! \brief The constructor allocates the problem integrand.
  SIMPoroElasticity(const std::vector<unsigned char>& fields)
  {
    Dim::myProblem = new PoroElasticity(Dim::dimension);
    Dim::nf = fields;
  }

  //! \brief Empty destructor.
  virtual ~SIMPoroElasticity() {}

  using SIMElasticity<Dim>::parse;
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem)
  {
    if (strcasecmp(elem->Value(),"poroelasticity"))
      return this->SIMElasticity<Dim>::parse(elem);

    std::vector<Material*>& mDat = SIMElasticity<Dim>::mVec;

    const TiXmlElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement())
      if (!strcasecmp(child->Value(),"scaling")) {
        double sc = 1.0e14;
        utl::getAttribute(child,"value",sc);
        static_cast<PoroElasticity*>(Dim::myProblem)->setScaling(sc);
        IFEM::cout <<"\tScaling: sc = "<< sc << std::endl;
      }
      else if (!strcasecmp(child->Value(),"isotropic")) {
        int code = this->parseMaterialSet(child,mDat.size());
        std::cout <<"\tMaterial code "<< code <<":" << std::endl;
        mDat.push_back(new PoroMaterial);
        mDat.back()->parse(child);
      }
      else
        this->SIMElasticity<Dim>::parse(child);

    if (!mDat.empty())
      static_cast<Elasticity*>(Dim::myProblem)->setMaterial(mDat.front());

    return true;
  }

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  virtual std::string getName() const { return "PoroElasticity"; }

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  //! \param[out] geoBlk Running geometry block counter
  //! \param[out] nBlock Running result block counter
  bool saveModel(char* fileName, int& geoBlk, int& nBlock)
  {
    if (Dim::opt.format < 0)
      return true;

    nBlock = 0;
    return this->writeGlvG(geoBlk,fileName);
  }

  //! \brief Saves the converged results of a given time step to VTF file.
  //! \param[in] tp Time stepping parameters
  //! \param nBlock Running result block counter
  bool saveStep(const TimeStep& tp, int& nBlock)
  {
    if (tp.step%Dim::opt.saveInc > 0 || Dim::opt.format < 0)
      return true;

    int iDump = 1 + tp.step/Dim::opt.saveInc;
    if (!this->writeGlvS(solution.front(),iDump,nBlock,tp.time.t,"vector",89))
      return false;

    return this->writeGlvStep(iDump,tp.time.t);
  }

  //! \brief Initializes the simulator time stepping loop.
  bool init()
  {
    solution.resize(this->getNoSolutions());
    for (size_t i=0;i<solution.size();++i)
      solution[i].resize(this->getNoDOFs(),true);

    this->registerField("init1", solution[0]);
    this->registerField("init2", solution[1]);

    SIM::setInitialConditions(*this);
    return true;
  }

  //! \brief Advances the time step one step forward.
  bool advanceStep(TimeStep&)
  {
    // Update vectors between time steps
    const int nNusols = solution.size();
    for (int n = nNusols-1; n > 0; n--)
      solution[n] = solution[n-1];

    return true;
  }

  //! \brief Computes the solution for the current time step.
  bool solveStep(TimeStep& tp)
  {
    if (Dim::msgLevel >= 0)
      IFEM::cout <<"\n  step = "<< tp.step <<"  time = "<< tp.time.t << std::endl;

    if (!this->assembleSystem(tp.time, solution))
      return false;

    if (!this->solveSystem(solution.front(), Dim::msgLevel-1,"displacement+pressure"))
      return false;

    return true;
  }

  //! \brief Register fields for data export
  void registerFields(DataExporter& exporter)
  {
    exporter.registerField("u,p","primary",DataExporter::SIM,DataExporter::PRIMARY);
    exporter.setFieldValue("u,p",this,&solution.front());
  }

  //! \brief Prints a summary of the calculated solution to std::cout.
  //! \param[in] solution The solution vector
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  virtual void printSolutionSummary(const Vector& solution, int, const char*,
                                    std::streamsize outPrec)
  {
    const size_t nsd = this->getNoSpaceDim();
    size_t iMax[nsd+1];
    double dMax[nsd+1];
    double dNorm = this->solutionNorms(solution,dMax,iMax,this->getNoFields(1));
    if (this->getNoFields(2) > 0) {
      double dpNorm = this->solutionNorms(solution,dMax+nsd,
                                          iMax+nsd,this->getNoFields(2), 'P');
      dNorm = sqrt(pow(dNorm,2.0)+pow(dpNorm,2.0));
    }

    std::stringstream str;
    if (this->adm.getProcId() == 0)
    {
      if (outPrec > 0) str.precision(outPrec);

      str <<"  Primary solution summary: L2-norm            : "<< utl::trunc(dNorm);

      char D = 'X';
      for (size_t d = 0; d < nsd; d++, D++)
        if (utl::trunc(dMax[d]) != 0.0)
          str <<"\n                            Max "<< char('X'+d)
              <<"-displacement : "<< dMax[d] <<" node "<< iMax[d];
        str <<"\n                            Max pressure       : "
            << dMax[nsd] <<" node "<< iMax[nsd] << std::endl;
    }

    utl::printSyncronized(std::cout,str,this->adm.getProcId());
  }

private:
  Vectors solution; //!< Solution vectors
};


//! \brief Partial specialization for configurator
template<class Dim>
struct SolverConfigurator< SIMPoroElasticity<Dim> > {
  int setup(SIMPoroElasticity<Dim>& pe, const bool&, char* infile)
  {
    utl::profiler->start("Model input");

    // Read input file
    if(!pe.read(infile))
      return 1;

    utl::profiler->stop("Model input");

    // Configure finite element library
    if(!pe.preprocess())
      return 2;

    // Setup integration
    pe.setQuadratureRule(pe.opt.nGauss[0],true);
    pe.initSystem(pe.opt.solver);
    pe.setMode(SIM::STATIC);
    pe.init();

    return 0;
  }
};

#endif

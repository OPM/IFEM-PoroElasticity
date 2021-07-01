// $Id$
//==============================================================================
//!
//! \file main.C
//!
//! \date April 16 2015
//!
//! \author Yared Bekele
//!
//! \brief Main program for the isogeometric poroelasticity solver.
//!
//==============================================================================

#include "IFEM.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "SIMStatPoroElasticity.h"
#include "SIMDynPoroElasticity.h"
#include "SIMSolver.h"
#include "ASMmxBase.h"
#include "GenAlphaSIM.h"
#include "Profiler.h"


/*!
  \brief Enumerates the different integration options.
*/
enum Integrator {
  FULL_STATIC,                  //!< Quasistatic in both displacements and flow (aka fully static)
  HALF_STATIC,                  //!< Quasistatic in displacements, dynamic in flow (aka half static)
  NEWMARK,                      //!< Dynamic, Newmark timestepping
  GENALPHA                      //!< Dynamic, Generalized alpha timestepping
};


/*!
  \brief Dynamic simulation driver.
*/

template<class T> class SIMDriver : public SIMSolver<T>
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  explicit SIMDriver(T& s) : SIMSolver<T>(s) {}
  //! \brief Empty destructor.
  virtual ~SIMDriver() {}
  //! \brief Overrides the stop time that was read from the input file.
  void setStopTime(double t) { SIMSolver<T>::tp.stopTime = t; }
};


/*!
  \brief Creates the poroelastic simulator and launches the simulation.
  \param[in] infile The input file to parse
  \param[in] stopTime Stop time of the simulation (if non-negative)
*/

template<class Dim, class Sim> int runSimulator (char* infile, double stopTime)
{
  utl::profiler->start("Model input");
  IFEM::cout <<"\n\n0. Parsing input file(s)."
             <<"\n========================="<< std::endl;

  // Establish the poroelastic FE model
  Sim model;
  if (!model.read(infile))
    return 1;

  model.opt.print(IFEM::cout) << std::endl;

  // Establish the simulation driver
  SIMDriver<Sim> solver(model);
  if (!solver.read(infile))
    return 1;

  utl::profiler->stop("Model input");

  if (stopTime >= 0.0)
    solver.setStopTime(stopTime);

  IFEM::cout <<"\n\n10. Preprocessing the finite element model:"
             <<"\n==========================================="<< std::endl;

  // Preprocess the model and establish data structures for the algebraic system
  if (!model.preprocess())
    return 2;

  // Initialize the linear solvers and solution vectors,
  // include assembly of reaction forces
  if (!model.init(solver.getTimePrm(),true))
    return 2;

  // Load initial conditions unless a restart
  if (model.opt.restartFile.empty())
    model.setInitialConditions();
  else if (solver.restart(model.opt.restartFile,model.opt.restartStep) < 0)
    return 2;

  // HDF5 output
  if (model.opt.dumpHDF5(infile))
    solver.handleDataOutput(model.opt.hdf5,model.getProcessAdm(),
                            model.opt.saveInc,model.opt.restartInc);

  return solver.solveProblem(infile,"100. Starting the simulation");
}


/*!
  \brief Creates the dynamic poroelastic simulator and launches the simulation.
  \param[in] infile The input file to parse
  \param[in] integrator The time integrator to use
  \param[in] T1 Stop time of the simulation (if non-negative)
*/

template<class Dim> int runSimulator (char* infile, Integrator integrator, double T1)
{
  switch (integrator) {
  case FULL_STATIC: return runSimulator<Dim, SIMStatPoroElasticity<Dim, true>>(infile, T1);
  case HALF_STATIC: return runSimulator<Dim, SIMStatPoroElasticity<Dim, false>>(infile, T1);
  case NEWMARK: return runSimulator<Dim, SIMDynPoroElasticity<Dim, NewmarkSIM>>(infile, T1);
  case GENALPHA: return runSimulator<Dim, SIMDynPoroElasticity<Dim, GenAlphaSIM>>(infile, T1);
  default: return 1;
  }
}


/*!
  \brief Main program for the isogeometric poroelasticity solver.
*/

int main (int argc, char** argv)
{
  Profiler prof(argv[0]);
  utl::profiler->start("Initialization");

  char* infile = nullptr;
  Integrator integrator = FULL_STATIC;
  bool twoD = false;
  double stopTime = -1.0;
  ASMmxBase::Type = ASMmxBase::NONE;

  IFEM::Init(argc,argv,"Poroelasticity solver");

  for (int i = 1; i < argc; i++)
    if (SIMoptions::ignoreOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!strcmp(argv[i],"-2D"))
      twoD = Elastic::planeStrain = true;
    else if (!strcmp(argv[i],"-mixed-full"))
      ASMmxBase::Type = ASMmxBase::FULL_CONT_RAISE_BASIS1;
    else if (!strncmp(argv[i],"-mixed",6))
      ASMmxBase::Type = ASMmxBase::REDUCED_CONT_RAISE_BASIS1;
    else if (!strcmp(argv[i],"-dyn2"))
      integrator = GENALPHA;
    else if (!strncmp(argv[i],"-dyn",4))
      integrator = NEWMARK;
    else if (!strcmp(argv[i],"-halfstatic"))
      integrator = HALF_STATIC;
    else if (!strcmp(argv[i],"-fullstatic"))
      integrator = FULL_STATIC;
    else if (!strncmp(argv[i],"-stopT",6) && i < argc-1)
      stopTime = atof(argv[++i]);
    else if (!infile)
      infile = argv[i];
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (!infile)
  {
    std::cout <<"Usage: "<< argv[0]
              <<" <inputfile> [-dense|-spr|-superlu[<nt>]|-samg|-petsc]\n"
              <<"       [-lag|-spec|-LR] [-2D] [-nGauss <n>] [-mixed]\n"
              <<"       [-dyn[1|2]|-halfstatic|-fullstatic]\n"
              <<"       [-vtf <format> [-nviz <nviz>] [-nu <nu> [-nv <nv>]\n"
              <<"       [-nw <nw>]] [-hdf5] [-stopTime <t>]\n";
    return 0;
  }

  IFEM::cout <<"\nInput file: "<< infile;
  IFEM::getOptions().print(IFEM::cout);
  if (stopTime >= 0.0)
    IFEM::cout <<"\nSimulation stop time: "<< stopTime;
  IFEM::cout << std::endl;

  if (twoD)
    return runSimulator<SIM2D>(infile,integrator,stopTime);
  else
    return runSimulator<SIM3D>(infile,integrator,stopTime);
}

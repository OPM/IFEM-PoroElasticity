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
#include "SIMDynPoroElasticity.h"
#include "SIMSolver.h"
#include "GenAlphaSIM.h"
#include "InitialConditionHandler.h"
#include "AppCommon.h"
#include "Profiler.h"
#include "SIMThermoPoroElasticity.h"

/*!
  \brief Creates the thermo-poroelastic simulator and launches the simulation.
  \param[in] infile The input file to parse
*/
template<class Dim> int runSimulatorThermalCoupled (char* infile, bool supg)
{
  SIMThermoPoroElasticity<Dim> model(supg);
  SIMSolver<SIMThermoPoroElasticity<Dim>> solver(model);

  // Read input file
  if(!model.read(infile) || !solver.read(infile))
    return 1;

  // Configure finite element library
  if(!model.preprocess())
    return 2;

  // Setup integration
  model.setQuadratureRule(model.opt.nGauss[0],true);
  model.initSystem(model.opt.solver,1,1,false);
  model.init(solver.getTimePrm());
  model.setInitialConditions();
  model.setAssociatedRHS(0,0);
  model.setMode(SIM::DYNAMIC);

  // HDF5 output
  DataExporter* exporter = NULL;
  if (model.opt.dumpHDF5(infile))
  {
    exporter = new DataExporter(true);
    exporter->registerWriter(new HDF5Writer(model.opt.hdf5,model.getProcessAdm()));
    exporter->registerWriter(new XMLWriter(model.opt.hdf5,model.getProcessAdm()));
  }

  int res = solver.solveProblem(infile,exporter);

  delete exporter;
  return res;
}


/*!
  \brief Creates the poroelastic simulator and launches the simulation.
  \param[in] infile The input file to parse
*/

template<class Dim, class Sim> int runSimulatorIsoThermal (char* infile)
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
  SIMSolver<Sim> solver(model);
  if (!solver.read(infile))
    return 1;

  utl::profiler->stop("Model input");
  IFEM::cout <<"\n\n10. Preprocessing the finite element model:"
             <<"\n==========================================="<< std::endl;

  // Preprocess the model and establish data structures for the algebraic system
  if (!model.preprocess())
    return 2;

  // Initialize the linear solvers
  if (!model.initSystem(model.opt.solver))
    return 2;

  // Initialize the solution fields
  model.init(TimeStep());
  SIM::setInitialConditions(model);

  // HDF5 output
  DataExporter* exporter = nullptr;
  if (model.opt.dumpHDF5(infile))
    exporter = SIM::handleDataOutput(model,solver,model.opt.hdf5,false,1,1);

  int res = solver.solveProblem(infile,exporter,"100. Starting the simulation");

  delete exporter;
  return res;
}


/*!
  \brief Creates the dynamic poroelastic simulator and launches the simulation.
  \param[in] infile The input file to parse
  \param[in] integrator The time integrator to use (0=linear quasi-static,
             1=linear Newmark, 2=Generalized alpha)
  \param[in] thermal If true include temperature effects
*/

template<class Dim> int runSimulator (char* infile,
                                      char integrator, bool thermal)
{
  if (integrator == 2)
    return runSimulatorIsoThermal<Dim, SIMDynPoroElasticity<Dim,GenAlphaSIM> >(infile);
  else if (integrator > 0)
    return runSimulatorIsoThermal<Dim, SIMDynPoroElasticity<Dim> >(infile);
  else {
    if (thermal)
      return runSimulatorThermalCoupled<Dim>(infile,false);
    else
      return runSimulatorIsoThermal<Dim, SIMPoroElasticity<Dim>>(infile);
  }
}


/*!
  \brief Main program for NURBS-based poroelasticity solver.
*/

int main (int argc, char ** argv)
{
  Profiler prof(argv[0]);
  utl::profiler->start("Initialization");

  int i;
  char* infile = 0;
  char integrator = 0;
  bool twoD = false;
  ASMmxBase::Type = ASMmxBase::NONE;
  bool thermal = false;

  IFEM::Init(argc,argv);

  for (i = 1; i < argc; i++)
    if (SIMoptions::ignoreOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!strcmp(argv[i],"-2D"))
      twoD = SIMElasticity<SIM2D>::planeStrain = true;
    else if (!strcmp(argv[i],"-mixed"))
      ASMmxBase::Type = ASMmxBase::FULL_CONT_RAISE_BASIS1;
    else if (!strcmp(argv[i],"-thermal"))
      thermal = true;
    else if (!strcmp(argv[i],"-dyn2"))
      integrator = 2;
    else if (!strncmp(argv[i],"-dyn",4))
      integrator = 1;
    else if (!infile)
      infile = argv[i];
    else
      std::cerr <<"*** Unknown option ignored: "<< argv[i] << std::endl;

  if (!infile)
  {
    std::cout <<"Usage: "<< argv[0]
              <<" <inputfile> [-dense|-spr|-superlu[<nt>]|-samg|-petsc]\n      "
              <<" [-lag|-spec|-LR] [-2D] [-nGauss <n>] [-mixed] [-dyn[1|2]]\n"
              <<"       [-vtf <format> [-nviz <nviz>] [-nu <nu> [-nv <nv>]"
              <<" [-nw <nw>]] [-hdf5]\n";
    return 0;
  }

  IFEM::cout <<"\n >>> IFEM Poroelasticity Solver <<<"
             <<"\n =================================="
             <<"\n Executing command:\n";
  for (i = 0; i < argc; i++) IFEM::cout <<" "<< argv[i];
  IFEM::cout <<"\n\n Input file: "<< infile;
  IFEM::getOptions().print(IFEM::cout);
  IFEM::cout << std::endl;

  if (twoD)
    return runSimulator<SIM2D>(infile,integrator,thermal);
  else
    return runSimulator<SIM3D>(infile,integrator,thermal);
}

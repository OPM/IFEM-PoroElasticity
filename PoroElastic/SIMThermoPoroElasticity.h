// $Id$
//==============================================================================
//!
//! \file SIMThermoPoroElasticity.h
//!
//! \date
//!
//! \author Yared Bekele
//!
//! \brief Simulation driver for non-isothermal PoroElasticity problems
//!
//==============================================================================

#ifndef _SIM_THERMO_POROELASTICITY_H_
#define _SIM_THERMO_POROELASTICITY_H_

#include "ThermoPoroElasticity.h"
#include "SIMSolver.h"
#include "TimeStep.h"
#include "PoroMaterial.h"
#include "NonLinSIM.h"
#include "Utilities.h"
#include "InitialConditionHandler.h"
#include "DataExporter.h"
#include "SAM.h"
#include "ASMs2D.h"
#include "SIMNodalConstraint.h"

/*!
  \brief Driver class for non-isothermal poroelastic simulators.
*/

template<class Dim> class SIMThermoPoroElasticity : public SIMNodalConstraint<Dim>
{
public:
  //! \brief The constructor initializes the references to the integrand.
  SIMThermoPoroElasticity(bool supg) :
    SIMNodalConstraint<Dim>({Dim::dimension,2}), tpe(Dim::dimension, 1, supg),
    tpewd(Dim::dimension), nSim(*this,NonLinSIM::L2)
  {
    Dim::myProblem = &tpe;
  }

  //! \brief Destructor.
  virtual ~SIMThermoPoroElasticity()
  {
    Dim::myProblem = NULL;
    Dim::myInts.clear();
  }

  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem)
  {
    if (strcasecmp(elem->Value(),"thermoporoelasticity"))
      return SIMNodalConstraint<Dim>::parse(elem);

    const TiXmlElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement())
    {
      if (!strcasecmp(child->Value(),"materialdata")) {
        int code = this->parseMaterialSet(child,mVec.size());
        std::cout <<"\tMaterial code "<< code <<":" << std::endl;
        mVec.push_back(std::unique_ptr<PoroMaterial>(new PoroMaterial));
        mVec.back()->parse(child);
        mVec.back()->printLog();
      }
      else if (!strcasecmp(child->Value(),"envtproperties"))
      {
        double Te = 0.0, lambdae = 0.0;
        utl::getAttribute(child,"Te",Te);
        utl::getAttribute(child,"lambdae",lambdae);
        tpewd.setEnvtTemperature(Te);
        tpewd.setEnvtConductivity(lambdae);
        std::cout << "\tEnvironment thermal properties:"
                  << "\n\t\tTe = " << Te << "\tlambdae = " << lambdae << std::endl;
      }
      else if (!strcasecmp(child->Value(),"scaling"))
      {
        double scl1 = 1.0, scl2 = 1.0;
        utl::getAttribute(child,"scl1",scl1);
        utl::getAttribute(child,"scl2",scl2);
        tpe.setScalingValues(scl1,scl2);
        std::cout << "\tScaling values:"
                  << "\n\t\tScale U-P = " << scl1 << "\tScale U-T = " << scl2 << std::endl;
      }
      else if (!strcasecmp(child->Value(),"nonlinearsolver"))
        nSim.parse(child);
      else
        this->SIMNodalConstraint<Dim>::parse(child);
    }

    if (!mVec.empty())
      tpe.setMaterial(mVec.front().get());

    return true;
  }

  //! \brief Performs some pre-processing tasks on the FE model.
  //! \details This method is reimplemented to couple the weak Dirichlet
  //! integrand to the generic Neumann property codes.
  virtual void preprocessA()
  {
    Dim::myInts.insert(std::make_pair(0,Dim::myProblem));

    PropertyVec::iterator p;
    for (p = Dim::myProps.begin(); p != Dim::myProps.end(); p++)
      if (p->pcode == Property::NEUMANN_GENERIC)
      {
        if (Dim::myInts.find(p->pindx) == Dim::myInts.end())
          Dim::myInts.insert(std::make_pair(p->pindx,&tpewd));
      }
  }

  //! \brief Initializes for integration of Neumann terms for a given property.
  //! \param[in] propInd Physical property index
  virtual bool initNeumann(size_t propInd)
  {
    typename Dim::SclFuncMap::const_iterator sit = Dim::myScalars.find(propInd);
    typename Dim::VecFuncMap::const_iterator vit = Dim::myVectors.find(propInd);
    typename Dim::TracFuncMap::const_iterator tit = Dim::myTracs.find(propInd);

    if (sit != Dim::myScalars.end())
      tpewd.setFlux(sit->second);
    else if (vit != Dim::myVectors.end())
      tpe.setTraction(vit->second);
    else if (tit != Dim::myTracs.end())
      tpe.setTraction(tit->second);
    else
      return false;

    return true;
  }

  //! \brief Evaluates some iteration norms for convergence assessment.
  //! \param[in] u Global primary solution vector
  //! \param[in] r Global residual vector associated with the solution vector
  //! \param[out] eNorm Energy norm of solution increment
  //! \param[out] rNorm Residual norm of solution increment
  //! \param[out] dNorm Displacement norm of solution increment
  void iterationNorms(const Vector& u, const Vector& r,
                      double& eNorm, double& rNorm, double& dNorm) const
  {
    eNorm = this->mySam->dot(r,u,'A');
    rNorm = this->mySam->norm2(r,'D') + this->mySam->norm2(r,'P');
    dNorm = this->mySam->norm2(u,'D') + this->mySam->norm2(u,'P');
  }

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  virtual std::string getName() const
  {
    return "ThermoPoroElasticity";
  }

  //! \brief Obtain const reference to solution vector.
  //! \param[in] i Solution vector to get reference to.
  //! \return Const reference to requested vector.
  const Vector& getSolution(int i)
  {
    return nSim.getSolution(i);
  }

  //! \brief Initializes the simulator time stepping loop
  bool init(const TimeStep& tp)
  {
    // Initialize solution vectors
    //nSim.init(2,RealArray(this->getNoDOFs(),0.0));
    nSim.initSol();
    size_t n, nSols = this->getNoSolutions();
    std::string str = "vector1";
    for (n = 0; n < nSols; n++, str[6]++)
      this->registerField(str,nSim.getSolution(n));

    return true;
  }

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  //! \param[out] geoBlk Running geometry block counter
  //! \param[out] nBlock Running result block counter
  bool saveModel(char* fileName, int& geoBlk, int& nBlock)
  {
    if (Dim::opt.format < 0)
      return true;

    geoBlk = nBlock = 0;
    return this->writeGlvG(geoBlk,fileName);
  }

  //! \brief Saves the converged results to VTF file of a given time step.
  //! \param[in] tp Time step identifier
  //! \param[in] nBlock Running VTF block counter
  bool saveStep(const TimeStep& tp, int& nBlock)
  {
     if (tp.step > 0) {
       double old = utl::zero_print_tol;
       utl::zero_print_tol = 1e-16;
       this->savePoints(nSim.getSolution(),tp.time.t,tp.step);
       utl::zero_print_tol = old;
     }

    if (tp.step%Dim::opt.saveInc > 0 || Dim::opt.format < 0)
      return true;

    int iDump = 1 + tp.step/Dim::opt.saveInc;

    // Write solution fields
    bool result = this->writeGlvS(nSim.getSolution(), iDump, nBlock,
                                  tp.time.t, "vector", 79);

    return result && this->writeGlvStep(iDump, tp.time.t);
  }

  //! \brief Advances the time step one step forward.
  bool advanceStep(TimeStep& tp)
  {
    return nSim.advanceStep(tp,false);
  }

  //! \brief Computes the solution for the current time step.
  bool solveStep(TimeStep& tp)
  {
    if (nSim.solveStep(tp,SIM::DYNAMIC) != SIM::CONVERGED)
      return false;

    return this->postSolve(tp);
  }

  //! \brief Post processing of solution if needed
  bool postSolve(const TimeStep&, bool = false)
  {
    return true;
  }

  //! \brief Register fields for data export
  void registerFields(DataExporter& exporter)
  {
    exporter.registerField("u-p-T","primary",DataExporter::SIM,DataExporter::PRIMARY);
    exporter.setFieldValue("u-p-T",this,&nSim.getSolution());
  }

  //! \brief Sets initial conditions.
  void setInitialConditions()
  {
    SIM::setInitialConditions(*this);
  }

  //! \brief Initializes material properties for integration of interior terms.
  //! \param[in] propInd Physical property index
  virtual bool initMaterial(size_t propInd)
  {
    if (propInd >= mVec.size())
      propInd = mVec.size()-1;

    tpe.setMaterial(mVec[propInd].get());
    return true;
  }

  //! \brief Prints a summary of the calculated solution.
  //! \param[in] solution The solution vector
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  virtual void printSolutionSummary(const Vector& solution, int, const char*,
                                    std::streamsize outPrec)
  {
    const size_t nsd = this->getNoSpaceDim();
    size_t iMax[nsd+2];
    double dMax[nsd+2];
    double dNorm;
    double dNormu = this->solutionNorms(solution, dMax, iMax, nsd, 'D');
    double dNormpt = this->solutionNorms(solution, dMax+nsd, iMax+nsd, 2, 'P');
    dNorm = sqrt(pow(dNormu,2.0)+pow(dNormpt,2.0));

    if (outPrec > 0) IFEM::cout.precision(outPrec);

    IFEM::cout <<"  Primary solution summary: L2-norm            : "<< utl::trunc(dNorm);

    char D = 'X';
    for (size_t d = 0; d < nsd; d++, D++)
      if (utl::trunc(dMax[d]) != 0.0)
        IFEM::cout <<"\n                            Max "<< char('X'+d)
                   <<"-displacement : "<< dMax[d] <<" node "<< iMax[d];
    if (utl::trunc(dMax[nsd]) != 0.0)
      IFEM::cout <<"\n                            Max pressure       : "
                 << dMax[nsd] <<" node "<< iMax[nsd] << std::endl;
    if (utl::trunc(dMax[nsd+1]) != 0.0)
      IFEM::cout <<"                            Max temperature    : "
                 << dMax[nsd+1] <<" node "<< iMax[nsd+1] << std::endl;
  }

private:
  ThermoPoroElasticity tpe; //!< ThermoPoroElasticity integrand
  ThermoPoroElasticity::WeakDirichlet tpewd; //!< ThermoPoroElasticity Robin integrand
  std::vector<std::unique_ptr<PoroMaterial>> mVec;  //!< Material data
  NonLinSIM nSim; //!< Nonlinear simulator
};

#endif  // _SIM_THERMO_POROELASTICITY_H_

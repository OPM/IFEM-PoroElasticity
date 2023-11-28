// $Id$
//==============================================================================
//!
//! \file SIMPoroElasticity.C
//!
//! \date May 29 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Simulation driver for poroelasticity problems.
//!
//==============================================================================

#include "SIMPoroElasticity.h"

#include "PoroElasticity.h"
#include "PoroMaterial.h"
#include "PoroSolutions.h"

#include "AnaSol.h"
#include "ASMmxBase.h"
#include "Functions.h"
#include "SAM.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "TimeStep.h"


template<class Dim>
SIMPoroElasticity<Dim>::SIMPoroElasticity (const std::vector<unsigned char>&)
{
  if (ASMmxBase::Type > ASMmxBase::NONE)
    Dim::nf = { Dim::dimension, 1 }; // mixed formulation
  else
    Dim::nf = { Dim::dimension+1 }; // standard formulation

  Dim::myHeading = "Poroelasticity solver";
  SIMElasticity<Dim>::myContext = "poroelasticity";
  scaleD = scaleP = 0.0;
}


template<class Dim>
bool SIMPoroElasticity<Dim>::solveStep (TimeStep& tp)
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


template<class Dim>
void SIMPoroElasticity<Dim>::iterationNorms (const Vector& u, const Vector& r,
                                             double& eNorm, double& rNorm,
                                             double& dNorm) const
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


template<class Dim>
void SIMPoroElasticity<Dim>::printSolutionSummary(const Vector& solution, int,
                                                  const char*, std::streamsize)
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


template<class Dim>
bool SIMPoroElasticity<Dim>::postSolve (TimeStep& tp)
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


template<class Dim>
bool SIMPoroElasticity<Dim>::initNeumann (size_t propInd)
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


template<class Dim>
Elasticity* SIMPoroElasticity<Dim>::getIntegrand ()
{
  if (!Dim::myProblem)
    Dim::myProblem = new PoroElasticity(Dim::dimension, Dim::nf.size() > 1);
  return static_cast<Elasticity*>(Dim::myProblem);
}


template<class Dim>
bool SIMPoroElasticity<Dim>::parse (const tinyxml2::XMLElement* elem)
{
  if (!strcasecmp(elem->Value(),"poroelasticity"))
  {
    const tinyxml2::XMLElement* child = elem->FirstChildElement();
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


template<class Dim>
bool SIMPoroElasticity<Dim>::parseAnaSol (const tinyxml2::XMLElement* elem)
{
  if (Dim::mySol || this->parseDimSpecific(elem))
    return true;

  IFEM::cout <<"\tAnalytical solution: Expression"<< std::endl;
  Dim::mySol = new AnaSol(elem);
  return true;
}


//! \brief Template specialization - 2D specific input parsing.
template<>
bool SIMPoroElasticity<SIM2D>::parseDimSpecific (const tinyxml2::XMLElement* elem)
{
  std::string type;
  utl::getAttribute(elem,"type",type,true);
  if (type == "terzhagi")
  {
    double height = 1.0, load = 1.0;
    utl::getAttribute(elem,"height",height);
    utl::getAttribute(elem,"load",load);
    IFEM::cout <<"\tAnalytical solution: Terzhagi, height = "<< height
               <<" load = "<< load << std::endl;
    Elasticity* elp = this->getIntegrand();
    PoroMaterial* pmat = static_cast<PoroMaterial*>(elp->getMaterial());
    double gacc = elp->getGravity().length();
    SIM2D::mySol = new AnaSol(new TerzhagiPressure(pmat,gacc,height,load));
  }
  else if (type == "terzhagi-stationary")
  {
    double load = 1.0;
    utl::getAttribute(elem,"load",load);
    IFEM::cout <<"\tAnalytical solution: Terzhagi (stationary),"
               <<" load = "<< load << std::endl;
    Elasticity* elp = this->getIntegrand();
    PoroMaterial* pmat = static_cast<PoroMaterial*>(elp->getMaterial());
    SIM2D::mySol = new AnaSol(new ConstFunc(0.0), nullptr,
                              new StationaryTerzhagiDisplacement(pmat,load));
  }
  else
    return false;

  return true;
}


template class SIMPoroElasticity<SIM2D>;
template class SIMPoroElasticity<SIM3D>;

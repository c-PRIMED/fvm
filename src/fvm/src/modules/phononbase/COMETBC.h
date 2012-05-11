#ifndef _COMETBC_H_
#define _COMETBC_H_

#include "misc.h"
#include "FloatVarDict.h"
#include "KSConnectivity.h"

template<class T>
struct COMETBC : public FloatVarDict<T>
{
  COMETBC()
    {
      this->defineVar("specifiedTemperature",T(300.0));
      this->defineVar("specifiedReflection",T(0.0));
      this->defineVar("FullyImplicit",false);
      this->InterfaceModel="DMM";
      this->bcType="Interface";   //-- if it is not specified, default is interface
    }
  string bcType;
  string InterfaceModel;
};

template<class T>
struct COMETModelOptions : public FloatVarDict<T>
{
  COMETModelOptions()
  {
    this->defineVar("timeStep",T(0.1));
    this->defineVar("initialTemperature",T(310.0));
    this->defineVar("Tref",T(299.0));
    this->printNormalizedResiduals = true;
    this->transient = false;
    this->timeDiscretizationOrder=1;
    this->absTolerance=1e-8;
    this->relTolerance=1-4;
    this->showResidual=5;
    this->maxLevels=2;
    this->AgglomerationMethod="FaceArea";
    this->preSweeps=0;
    this->postSweeps=2;
    this->relFactor=.1;
    this->withNormal=false;
    this->NewtonTol=1e-4;
    this->BinsPerBand=10;
  }
  
  bool printNormalizedResiduals;
  bool transient;
  int timeDiscretizationOrder;
  double absTolerance;
  double relTolerance;
  int showResidual;
  int maxLevels;
  string AgglomerationMethod;
  int preSweeps;
  int postSweeps;
  double relFactor;
  bool withNormal;
  double NewtonTol;
  int BinsPerBand;

};

template<class T>
struct COMETIC
{
  typedef KSConnectivity<T> TKConnectivity;
  typedef TKConnectivity* TKCptr;
  typedef vector<TKCptr> TKClist;
  
  COMETIC(const int mID0, const int fID0, const int mID1, const int fID1, const int faceNumber)
  {
    this->InterfaceModel="DMM";
    this->MeshID0=mID0;
    this->FgID0=fID0;
    this->MeshID1=mID1;
    this->FgID1=fID1;
  }

  const TKClist& getKConnectivity(const int fgid) const
  {
    if(fgid==FgID0)
      return Connect0to1;
    else if (fgid==FgID1)
      return Connect1to0;
 
    throw CException("FaceGroup not apart of this interface!");
  };

  TKClist& getKConnectivity(const int fgid)
  {
    if(fgid==FgID0)
      return Connect0to1;
    else if (fgid==FgID1)
      return Connect1to0;
 
    throw CException("FaceGroup not apart of this interface!");
  };

  int getSelfFaceID(const int mesh)
  {
    if(mesh==MeshID0)
      return FgID0;
    else if(mesh==MeshID1)
      return FgID1;

    throw CException("Mesh not apart of this interface!");
  }
  
  string InterfaceModel;
  TKClist Connect0to1;
  TKClist Connect1to0;
  int MeshID0;
  int FgID0;
  int MeshID1;
  int FgID1;
  
};

#endif

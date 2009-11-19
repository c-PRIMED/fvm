#ifndef _QUADRATURE_H_
#define _QUADRATURE_H_

#include <math.h>
#include "Array.h"
#include "Vector.h"

class quadrature
{
 public:
  typedef Array<double> dArray;
  typedef Vector<double,3> Vectord3;
  //typedef Array<Vector<double,3>> Vectord3Array;
  dArray* absci1Ptr;
  dArray* absci2Ptr;
  dArray* absci3Ptr;
  dArray* wts1Ptr;
  dArray* wts2Ptr; 
  dArray* wts3Ptr;
  
  // quadrature(int,int)
  quadrature(int N1, int N2, int N3, double T2)
    {
      absci1Ptr=new dArray(N1);
      dArray & absci1= *absci1Ptr;
      absci2Ptr=new dArray(N2);
      dArray & absci2= *absci2Ptr;
      absci3Ptr=new dArray(N3);
      dArray & absci3= *absci3Ptr; 
      wts1Ptr=new dArray(N1);
      dArray & wts1= *wts1Ptr;
      wts2Ptr=new dArray(N2);
      dArray & wts2= *wts2Ptr; 
       wts3Ptr=new dArray(N3);
      dArray & wts3= *wts3Ptr;
     
      //Array<Vector<double,3>> cxyz;
      double cxmin,cxmax,cymin,cymax,czmin,czmax;
      double dcx,dcy,dcz;
      double clim = 5.5;
      cxmin=-clim*sqrt(0.5*T2);cxmax=clim*sqrt(0.5*T2);
      cymin=-clim*sqrt(0.5*T2);cymax=clim*sqrt(0.5*T2);
      czmin=-clim*sqrt(0.5*T2);czmax=clim*sqrt(0.5*T2);
      dcx=(cxmax-cxmin)/(N1-1.0);
      dcy=(cymax-cymin)/(N2-1.0);  
      dcz=(czmax-czmin)/(N3-1.0);
      absci1[0]=cxmin;absci2[0]=cymin; absci3[0]=czmin;
      for  (int j3=1;j3<N3;j3++){
	absci3[j3]=absci3[j3-1]+dcz;
	wts3[j3]=dcz;
      }
      for (int j2=1;j2<N2;j2++){
	absci2[j2]=absci2[j2-1]+dcy;
	wts2[j2]=dcy;
      }
      for  (int j1=1;j1<N1;j1++){
	absci1[j1]=absci1[j1-1]+dcx;	
	wts1[j1]=dcx;
      }
    }
  /*quadrature(int& N1, int& N2, double& T2) :
    _var1(N1),
    _var2(N2),
    _var3(T2)
    {
    logCtor();
    }*/
  
  //virtual ~quadrature() {}
};

#endif

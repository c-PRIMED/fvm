// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _QUADRATURE_H_
#define _QUADRATURE_H_

#include <math.h>
#include "Array.h"
#include "Array2D.h"
#include "Vector.h"
#include <stdio.h>

template<class T>
/**
 * Class quadrature for ESBGK simulations
 * two types of constructors:
 * 1) cartesian type of discrete velocity grid (int N1, int N2, int N3, double T2)
 * 2) spherical type with Gauss-Hermite quadrature in velocity magnitude and
 *    either constant or 3/8th rule in theta and phi angles.
 */

class Quadrature 
{
 public:
  /**
   * A constructor.
   * A constructor for the discrete ordinate in velocity space based on cartesian-type and spherical-type coordinates.
   */ 
  typedef Array<T> TArray;
  typedef  Array2D<T> TArray2D;
  typedef  Quadrature<T> TQuad;
   /**
   *  Cx pointer. 
   *  Pointer to discrete velocity in x-direction.
   */
  TArray* cxPtr;  
  /**
   * Cy pointer.  
   * Pointer to discrete velocity in y-direction.
   */ 
  TArray* cyPtr; 
  /**
   * Cz pointer.  
   * Pointer to discrete velocity in z-direction.
   */  
  TArray* czPtr;   
  /**
   * dcxyz pointer. 
   * Pointer to weights associated with each direction in velocity space.
   */ 
  TArray* dcxyzPtr; 
  TArray2D* malphaBGKPtr;
  TArray2D* malphaESBGKPtr;

  int getDirCount() const {return N123;}
  int getNVCount()const {return _NV;}
  int getNthetaCount()const {return _Ntheta;}
  int getNphiCount()const {return _Nphi;}
  T get_dcx()const{return _dcx;}
  T get_dcy()const{return _dcy;}
  T get_dcz() const{return _dcz;}
  TArray& get_absci1()const{return *absci1Ptr;}
  TArray& get_absci2()const{return *absci2Ptr;}
  TArray& get_absci3()const{return *absci3Ptr;}
  TArray& get_cx()const{return *cxPtr;}
  TArray& get_cy()const{return *cyPtr;}
  TArray& get_cz()const{return *czPtr;}
  TArray& get_wts1()const{return *wts1Ptr;}
  TArray& get_wts2()const{return *wts2Ptr;}
  TArray& get_wts3()const{return *wts3Ptr;}
  TArray& get_dcxyz()const{return *dcxyzPtr;}

   // bool printCx;
  /**
   * cartesian-type member taking in 5 argumetns
   * @param N1 -Number of ordinates in x-velocity. 
   * @param N2 -Number of ordinates in y-velocity.  
   * @param N3 -Number of ordinates in z-velocity. 
   * @param clim -cut-off range in velocity/sqrt(T2/2). 
   * @param T2 -Lowest non-dimensional temperature used to set the limit on discrete velocities in each direction . 
 */
Quadrature(int N1,  int N2,  int N3, double clim,  double T2)
    {
      absci1Ptr=new TArray(N1);
      TArray & absci1= *absci1Ptr;
      absci2Ptr=new TArray(N2);
      TArray & absci2= *absci2Ptr;
      absci3Ptr=new TArray(N3);
      TArray & absci3= *absci3Ptr; 
      wts1Ptr=new TArray(N1);
      TArray & wts1= *wts1Ptr;
      wts2Ptr=new TArray(N2);
      TArray & wts2= *wts2Ptr; 
      wts3Ptr=new TArray(N3);
      TArray & wts3= *wts3Ptr;
      /**
       * integer N123
       * total number of velocity directions.
       */
      N123=N1*N2*N3;
      // BGK and ESBGK Equilibrium distribution functions
      malphaBGKPtr=new TArray2D(N123,5);
      TArray2D & malphaBGK= *malphaBGKPtr;
      
      malphaESBGKPtr=new TArray2D(N123,10);
      TArray2D & malphaESBGK= *malphaESBGKPtr;

      _NV=N1;_Ntheta=N2;_Nphi=N3;


      cxPtr= new TArray(N123); cyPtr= new TArray(N123); czPtr= new TArray(N123);
      dcxyzPtr= new TArray(N123);
      TArray & cx= *cxPtr;TArray & cy= *cyPtr;TArray & cz= *czPtr;
      TArray & dcxyz= *dcxyzPtr;
      //Array<Vector<double,3>> * cxyz;// Cxyz Vector pointer.  /* Pointer to discrete velocity in xyz-directions. cxyz= (cx,cy,cz)*/
      T cxmin,cxmax,cymin,cymax,czmin,czmax;
      T dcx,dcy,dcz;
      //const T clim (5.5);
      cxmin=-clim*sqrt(0.5*T2);cxmax=clim*sqrt(0.5*T2);
      cymin=-clim*sqrt(0.5*T2);cymax=clim*sqrt(0.5*T2);
      czmin=-clim*sqrt(0.5*T2);czmax=clim*sqrt(0.5*T2);
      dcx=(cxmax-cxmin)/(N1-1.0);
      dcy=(cymax-cymin)/(N2-1.0);  
      dcz=(czmax-czmin)/(N3-1.0);

      _dcx=dcx;_dcy=dcy;_dcz=dcz;
      absci1[0]=cxmin;absci2[0]=cymin; absci3[0]=czmin;
      wts1[0]=dcx;	wts2[0]=dcy;	wts3[0]=dcz;
      for  (int j3=1;j3<N3;j3++){
	absci3[j3]=absci3[j3-1]+dcz;
	wts3[j3]=dcz;
	//printf("%4.2f", absci3[j3]);
      }
      //printf("\n");
      for (int j2=1;j2<N2;j2++){
	absci2[j2]=absci2[j2-1]+dcy;
	wts2[j2]=dcy;
	//printf("%4.2f", absci2[j2]);
      }
      //printf("\n");
      for  (int j1=1;j1<N1;j1++){
	absci1[j1]=absci1[j1-1]+dcx;	
	wts1[j1]=dcx;
	//printf("%4.2f", absci1[j1]);
      }
      //printf("\n");
     
      int j;
      j=0;
      for(int j1=0;j1<N1;j1++){
	for (int j2=0;j2<N2;j2++){
	  for (int j3=0;j3<N3;j3++){
	    cx[j]=absci1[j1];	
	    cy[j]=absci2[j2];
	    cz[j]=absci3[j3];
	    dcxyz[j]=wts1[j1]*wts2[j2]*wts3[j3];
	    
	    malphaBGK(j,0)=1.0;
	    malphaBGK(j,1)=cx[j];
	    malphaBGK(j,2)=cy[j]; 
	    malphaBGK(j,3)=cz[j];
	    malphaBGK(j,4)=pow(cx[j],2)+pow(cy[j],2)+pow(cz[j],2); 
	    
	    malphaESBGK(j,0)=1.0;
	    malphaESBGK(j,1)=cx[j];
	    malphaESBGK(j,2)=cy[j]; 
	    malphaESBGK(j,3)=cz[j];
	    malphaESBGK(j,4)=pow(cx[j],2); 
	    malphaESBGK(j,5)=pow(cy[j],2);
	    malphaESBGK(j,6)=pow(cz[j],2);
	    malphaESBGK(j,7)=cx[j]*cy[j]; 
	    malphaESBGK(j,8)=cy[j]*cz[j];
	    malphaESBGK(j,9)=cz[j]*cx[j];
	    
	    j++;
	  }
	}
      }

      /*
      FILE * pFile;
      pFile = fopen ("cxyz.txt","w");
      for(int j=0;j<N123;j++){	  
	fprintf(pFile,"%12.6f %12.6f %12.6f %12.6f \n", cx[j],cy[j],cz[j],dcxyz[j]);}
      fclose (pFile);
      */
    }

/**
   * spherical-type member taking in 6 argumetns
   * @param option_ur =0 for constant; =2,4,8,16 for Gauss-Hermite quadrature in velocity magnitude . 
   * @param Nr =number of ordinates in velocity magnitude if option_ur=0.  
   * @param option_theta =0 for constant; =1 for 3/8th rule discretization of azimuthal angle(theta). 
   * @param n_int =number of ordinates in theta if option_theta =0;
   *              =number of coarse intervals for 3/8th rule if option_theta=1(total no. of angles = 3*n_int).
   * @param option_phi  =0 for constant; =1 for 3/8th rule discretization of polar angle(phi).
   * @param nphi_int  =number of ordinates in phi if option_phi =0; =number of coarse intervals for 3/8th rule if option_phi=1 (total no. of angles = 3*n_int+1)
 */
  Quadrature(int option_ur, int Nr, int option_theta, int n_int, int option_phi, int nphi_int)  
       {
      int N1,N2,N3;
      if(option_ur ==0){
	N1=Nr;}
      else{N1=option_ur;}
      if(option_theta ==0){
	N2=n_int;}
      else{N2=3*n_int;}
      if(option_phi ==0){
	N3=nphi_int;}
      else{N3=3*nphi_int+1;}
      _NV=N1;_Ntheta=N2;_Nphi=N3;
      absci1Ptr=new TArray(N1);
      TArray & absci1= *absci1Ptr;
      absci2Ptr=new TArray(N2);
      TArray & absci2= *absci2Ptr;
      absci3Ptr=new TArray(N3);
      TArray & absci3= *absci3Ptr; 
      wts1Ptr=new TArray(N1);
      TArray & wts1= *wts1Ptr;
      wts2Ptr=new TArray(N2);
      TArray & wts2= *wts2Ptr; 
      wts3Ptr=new TArray(N3);
      TArray & wts3= *wts3Ptr;

      N123=N1*N2*N3;
      cxPtr= new TArray(N123); cyPtr= new TArray(N123); czPtr= new TArray(N123);
      dcxyzPtr= new TArray(N123);
      TArray & cx= *cxPtr;TArray & cy= *cyPtr;TArray & cz= *czPtr;
      TArray & dcxyz= *dcxyzPtr;

  // BGK and ESBGK Equilibrium distribution functions
      malphaBGKPtr=new TArray2D(N123,5);
      TArray2D & malphaBGK= *malphaBGKPtr;
      
      malphaESBGKPtr=new TArray2D(N123,10);
      TArray2D & malphaESBGK= *malphaESBGKPtr;
      _NV=N1;_Ntheta=N2;_Nphi=N3;

      switch(option_ur){
	int j1;
      case 0:   // constant difference for Ur
	{double dh1=sqrt(3.0)*3.889/(N1); //neglect ur=.0; 
	  for (j1=0;j1<N1;j1++){
	    absci1[j1]=(j1+1)*dh1;
	    wts1[j1]=dh1*pow(absci1[j1],2.0);
	  }
	}break; 
      case 2:{
	absci1[0]=0.7539869079898871;
	absci1[1]=0.1734055298879163E+1;
	wts1[0]=0.2738413467266824; 
	wts1[1]=0.1692721159996965;
	for (j1=0;j1<N1;j1++){
	  wts1[j1]=wts1[j1]*exp(pow(absci1[j1],2.0));
	}
      }break;
      case 4:            //4-gaussian for y^2*exp(-y^2)
	{ absci1[0]=0.4238628193900528;
	  absci1[1]=0.1014332104566760E+1;
	  absci1[2]=0.1742437375162050E+1;
	  absci1[3]=0.2639813333635586E+1;
	  wts1[0]=0.7649092266787873E-1;
	  wts1[1]= 0.2435439494642453;
	  wts1[2]= 0.1162953035510695;
	  wts1[3]= 0.6783287043185401E-2;
	  for (j1=0;j1<N1;j1++){
	    wts1[j1]=wts1[j1]*exp(pow(absci1[j1],2.0));
	  }
	}break; 
      case 8:   //8-gaussian for y^2*exp(-y^2)
	{ absci1[0]=0.1990000637984294;
	  absci1[1]=0.5059526450205794;
	  absci1[2]=0.9041682182040568;
	  absci1[3]=0.1372615723971598E+1;
	  absci1[4]=0.1900969572329702E+1;
	  absci1[5]=0.2490479841967435E+1;
	  absci1[6]=0.3158780677105240E+1;
	  absci1[7]=0.3966720403265353E+1;
	  
	  wts1[0]=0.9599144336400067E-2;
	  wts1[1]=0.7072944976303661E-1;
	  wts1[2]=0.157366887003943;
	  wts1[3]=0.1429322724003870;
	  wts1[4]=0.5431444004253597E-1;
	  wts1[5]=0.7835224153141577E-2;
	  wts1[6]=0.3338952597020048E-3;
	  wts1[7]=0.2149767232664775E-5;
	  for(j1=0;j1<N1;j1++){
	    wts1[j1]=wts1[j1]*exp(pow(absci1[j1],2.0));
	  }
	}break;
      case 16: //weights & absci for y^2 exp(-y^2)
	{ absci1[0]= 0.8174913389984520E-1;
	  absci1[1]=0.2154761962759740;
	  absci1[2]=0.4003530517087630;
	  absci1[3]=0.6298538771405607;
	  absci1[4]=0.8976124329697087;
	  absci1[5]=0.1198149260240153E+1;
	  absci1[6]=0.1527188184719962E+1;
	  absci1[7]=0.1881745606015598E+1;
	  absci1[8]=0.2260132964654088E+1;
	  absci1[9]=0.2661980315279350E+1;
	  absci1[10]=0.3088376381635592E+1;
	  absci1[11]=0.3542256017930265E+1;
	  absci1[12]=0.4029312272760483E+1;
	  absci1[13]=0.4560203031431090E+1;
	  absci1[14]=0.5156826768007481E+1;
	  absci1[15]=0.58781144889155572E+1;
	  
	  wts1[0]=0.7050727473210895E-3;
	  wts1[1]=0.7107111654073120E-2;
	  wts1[2]=0.2844188515941899E-1;
	  wts1[3]=0.6660235171398239E-1;
	  wts1[4]=0.1025785712747278;
	  wts1[5]=0.1077502032531791;
	  wts1[6]=0.7747156370638879E-1;
	  wts1[7]=0.3763106373385135E-1;
	  wts1[8]=0.1204873635560290E-1;
	  wts1[9]=0.2453208613776865E-2;
	  wts1[10]=0.3020309847850189E-3;
	  wts1[11]=0.2092121075871870E-4;
	  wts1[12]=0.7314637349679360E-6;
	  wts1[13]=0.1080646863902574E-7;
	  wts1[14]=0.4828081616137754E-10;
	  wts1[15]=0.2840126937112534E-13;
	  for (j1=0;j1<N1;j1++){
	    wts1[j1]=wts1[j1]*exp(pow(absci1[j1],2.0));
	  }
	}break;
      }
      switch(option_theta){
	
      case 0:   //constant difference for Theta
	{ double pi=3.14159;
	  double dh2=2*pi/(N2); //number of intervals=no. of ordinates
	  for (int j2=0;j2<N2;j2++){
	    absci2[j2]= dh2*j2+dh2/2.0; //no points on axes
	    wts2[j2]=dh2;
	    //cout<<"j2= "<<j2<<" absci2= "<<absci2[j2]<<endl;
	  }
	}break;
      case 1:  //three-eights for theta
	{ double pi=3.14159;
	  double h=pi;                     //the range of quadrature from [-h,h];
	  double dh=2*h/n_int;             //anglular increment
	  for (int i=0;i<n_int;i++){
	    absci2[3*i]=-h+dh/3.0*(3.0*i);       
	    absci2[3*i+1]=-h+dh/3.0*(3.0*i+1);   
	    absci2[3*i+2]=-h+dh/3.0*(3.0*i+2);   
	    if (i == 0){
	      wts2[i]=dh/4.0;  //here 0 -pi and pi are the same number
	    }
	    else {
	      wts2[3*i]=dh/4.0;
	    }
	    wts2[3*i+1]=dh/8.0*3.0;
	    wts2[3*i+2]=dh/8.0*3.0;
	  }
	  //absci2(N2)=h;wts2(N2)=dh/8.0;
	}break;
      }
      
      switch(option_phi){
	
      case 0: //constant difference for Phi
	{
	  double pi=3.14159;
	  //double dh3=pi/(N3-1); //number of intervals < number of ordinates (old)
	  double dh3=pi/N3;
	  for (int j3=0;j3<N3;j3++){
	    //absci3[j3]= dh3*(j3); //old
	    absci3[j3]= dh3*(j3+0.5); //no points on axes
	    wts3[j3]=dh3*sin(absci3[j3]); 
	  
	  }
	}break;
      case 1: //three-eights for phi
	{ double pi=3.14159;
	  double dh_phi=pi/nphi_int;  //anglular increment
	  for (int i=0;i<nphi_int;i++){
	    absci3[3*i]=dh_phi/3.0*(3.0*i);      
	    absci3[3*i+1]=dh_phi/3.0*(3.0*i+1);  
	    absci3[3*i+2]=dh_phi/3.0*(3.0*i+2);  
	    if (i == 0){  
	      wts3[i]=dh_phi/8.0*sin(absci3[i]);
	    }
	    else{
	      wts3[3*i]=dh_phi/4.0*sin(absci3[3*i]);
	    }      
	    wts3[3*i+1]=dh_phi/8.0*3.0*sin(absci3[3*i+1]);
	    wts3[3*i+2]=dh_phi/8.0*3.0*sin(absci3[3*i+2]);
	  }
	  absci3[N3-1]=pi;
	  wts3[N3-1]=dh_phi/8.0*sin(absci3[N3-1]);
	}break;
      }
       int j;
       j=0;
       for (int j1=0;j1<N1;j1++){ //Ur
	 for (int j2=0;j2<N2;j2++){ //theta
	   for (int j3=0;j3<N3;j3++){ //phi
	     cx[j]=absci1[j1]*cos(absci2[j2])*sin(absci3[j3]); //cx=Ur*cos(theta)*sin(phi)
	     cy[j]=absci1[j1]*sin(absci2[j2])*sin(absci3[j3]); //cy=Ur*sin(theta)*sin(phi)
	     cz[j]=absci1[j1]*cos(absci3[j3]);                 //Ur*cos(phi)
	     dcxyz[j]=wts1[j1]*wts2[j2]*wts3[j3];
	     malphaBGK(j,0)=1.0;
	     malphaBGK(j,1)=cx[j];
	     malphaBGK(j,2)=cy[j]; 
	     malphaBGK(j,3)=cz[j];
	     malphaBGK(j,4)=pow(cx[j],2)+pow(cy[j],2)+pow(cz[j],2); 
	     
	     malphaESBGK(j,0)=1.0;
	     malphaESBGK(j,1)=cx[j];
	     malphaESBGK(j,2)=cy[j]; 
	     malphaESBGK(j,3)=cz[j];
	     malphaESBGK(j,4)=pow(cx[j],2); 
	     malphaESBGK(j,5)=pow(cy[j],2);
	     malphaESBGK(j,6)=pow(cz[j],2);
	     malphaESBGK(j,7)=cx[j]*cy[j]; 
	     malphaESBGK(j,8)=cy[j]*cz[j];
	     malphaESBGK(j,9)=cz[j]*cx[j];
	    
	     
	     j++;       
	   }
	 }
       }
       // for (int j=10;j<21;j++){cout << j <<" : "<<cx[j] << " , " << dcxyz[j] <<endl; };cout<<endl;
        
      FILE * pFile;
      pFile = fopen ("cxyz.txt","w");
      for(int j=0;j<N123;j++){	  
	fprintf(pFile,"%d %12.6E %12.6E %12.6E %12.6E \n", j,cx[j],cy[j],cz[j],dcxyz[j]);}
      fclose (pFile);
       
    }

  Quadrature()
    {}

  void CopyQuad(TQuad& copyFromQuad)
  {
    N123=copyFromQuad.getDirCount();
    _NV=copyFromQuad.getNVCount();
    _Ntheta=copyFromQuad.getNthetaCount();
    _Nphi=copyFromQuad.getNphiCount();
    _dcx=copyFromQuad.get_dcx();
    _dcy=copyFromQuad.get_dcy();
    _dcz=copyFromQuad.get_dcz();

    absci1Ptr=new TArray(_NV);
    TArray& absci1= *absci1Ptr;
    absci2Ptr=new TArray(_Ntheta);
    TArray& absci2= *absci2Ptr;
    absci3Ptr=new TArray(_Nphi);
    TArray& absci3= *absci3Ptr;

    cxPtr=new TArray(N123);
    TArray& cx= *cxPtr;
    cyPtr=new TArray(N123);
    TArray& cy= *cyPtr;
    czPtr=new TArray(N123);
    TArray& cz= *czPtr;

    wts1Ptr=new TArray(_NV);
    TArray& wts1= *wts1Ptr;
    wts2Ptr=new TArray(_Ntheta);
    TArray& wts2= *wts2Ptr;
    wts3Ptr=new TArray(_Nphi);
    TArray& wts3= *wts3Ptr;

    dcxyzPtr=new TArray(N123);
    TArray& dcxyz= *dcxyzPtr;

    malphaBGKPtr=new TArray2D(N123,5);
    TArray2D& malphaBGK= *malphaBGKPtr;
    malphaESBGKPtr=new TArray2D(N123,10);
    TArray2D& malphaESBGK= *malphaESBGKPtr;

    absci1=copyFromQuad.get_absci1();
    absci2=copyFromQuad.get_absci2();
    absci3=copyFromQuad.get_absci3();

    cx=copyFromQuad.get_cx();	
    cy=copyFromQuad.get_cy();	
    cz=copyFromQuad.get_cz();

    wts1=copyFromQuad.get_wts1();
    wts2=copyFromQuad.get_wts2();
    wts3=copyFromQuad.get_wts3();

    dcxyz=copyFromQuad.get_dcxyz();

    int j;
    j=0;
    for(int j1=0;j1<_NV;j1++){
      for (int j2=0;j2<_Ntheta;j2++){
	for (int j3=0;j3<_Nphi;j3++){
	  malphaBGK(j,0)=1.0;
	  malphaBGK(j,1)=cx[j];
	  malphaBGK(j,2)=cy[j];
	  malphaBGK(j,3)=cz[j];
	  malphaBGK(j,4)=pow(cx[j],2)+pow(cy[j],2)+pow(cz[j],2);

	  malphaESBGK(j,0)=1.0;
	  malphaESBGK(j,1)=cx[j];
	  malphaESBGK(j,2)=cy[j];
	  malphaESBGK(j,3)=cz[j];
	  malphaESBGK(j,4)=pow(cx[j],2);
	  malphaESBGK(j,5)=pow(cy[j],2);
	  malphaESBGK(j,6)=pow(cz[j],2);
	  malphaESBGK(j,7)=cx[j]*cy[j];
	  malphaESBGK(j,8)=cy[j]*cz[j];
	  malphaESBGK(j,9)=cz[j]*cx[j];

	  j++;
	}
      }
    }
  }
 
  virtual ~Quadrature() {}

 private:
  TArray* absci1Ptr; 
  TArray* absci2Ptr;
  TArray* absci3Ptr;
  TArray* wts1Ptr;
  TArray* wts2Ptr; 
  TArray* wts3Ptr;
  int N123;
  int _NV;
  int _Ntheta;
  int _Nphi;
  T _dcx;
  T _dcy;
  T _dcz;
};


#endif

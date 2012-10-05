// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "MPM_Particles.h"

#include "StorageSite.h"
#include "CRConnectivity.h"



typedef MPM::VecD3 VecD3;
typedef MPM::VecD3Array VecD3Array;

MPM::MPM(string fileName):
  _particles(0),
  _coordinates(),
  _velocities(),
  _types(),
  _temperatures()
{
  char* file;
  file = &fileName[0];

  // MPM::setandwriteParticles(file);
  //get coordinate
  const shared_ptr<VecD3Array> MPM_Coordinates = MPM::readCoordinates(file);

  //get velocity
  const shared_ptr<VecD3Array> MPM_Velocities = MPM::readVelocities(file);
  
  //get type
  const shared_ptr<Array<int> > MPM_Types = MPM::readTypes(file);

  //get temperature
  const shared_ptr<Array<double> > MPM_Temperatures = MPM::readTemperatures(file);
  //store all the information in MPM class
  MPM::Init(MPM_Coordinates, MPM_Velocities, MPM_Types, MPM_Temperatures);}

MPM::MPM( ):
  _particles(0),
  _coordinates(),
  _velocities(),
  _types(),
  _temperatures()
{}


MPM::~MPM() { }


void MPM::setandwriteParticles(const char *file)
{
 
    FILE *fp;
           
    VecD3 center;
    center[0]=0.0;
    center[1]=0.0;
    center[2]=0.0;

    int count=0;   
        
#if 0
    //set up particle cartesian coordinate
    VecD3 temp;
    
    const double innerSide = 1.0; 
    const double midSide = 3.0;
    const double outerSide = 4.2;

    int nX=20, nY=20, nZ=1;
    
    double gapX=outerSide/(nX), gapY=outerSide/(nY), gapZ=0;
    const int nMPM = (nX+1)*(nY+1)*nZ;
    Array<VecD3> solidPoint(nMPM);
    Array<VecD3> solidVelocity(nMPM);
    Array<int> type(nMPM);
    Array<double> solidTemperature(nMPM);  
    for(int i=0; i<=nX; i++){
      for(int j=0; j<=nY; j++){
	for(int k=0; k<nZ; k++){
	  temp[0]=i*gapX-outerSide/2.;
	  temp[1]=j*gapY-outerSide/2.;
	  temp[2]=k*gapZ;
	  //inner square
	  if( temp[0]>=(-innerSide/2.) && temp[0]<=(innerSide/2.) && temp[1]>=(-innerSide/2.) && temp[1]<=(innerSide/2.))
	    {
	     //surface particles
	      type[count] = 0;  
	      if( temp[0]> (innerSide/2.0-gapX) || temp[0]<(-innerSide/2.0+gapX))
		type[count]=1;
	      if( temp[1]> (innerSide/2.0-gapY) || temp[1]<(-innerSide/2.0+gapY))
		type[count]=1;
	      //rotate
	      double alfa = atan(1.);
	      solidPoint[count][0]=temp[0]*cos(alfa)-temp[1]*sin(alfa);
	      solidPoint[count][1]=temp[1]*cos(alfa)+temp[0]*sin(alfa);
	      solidPoint[count][2]=temp[2];
	          
	      count+=1;
	    }
	  
	  //outer square
	  if( !(temp[0]>(-midSide/2.) && temp[0]<(midSide/2.) && temp[1]>(-midSide/2.) && temp[1]<(midSide/2.)))
	    {
	      type[count] = 0;  
	      //surface particles
	      
	      if( temp[0]< (midSide/2.0+gapX) && temp[0]>(midSide/2.0-gapX) && temp[1]<(midSide/2.0+gapY) && temp[1]>(-midSide/2.0-gapY))
		type[count]=1;
	      if( temp[0]> (-midSide/2.0-gapX) && temp[0]<(-midSide/2.0+gapX) && temp[1]<(midSide/2.0+gapY) && temp[1]>(-midSide/2.0-gapY))
		type[count]=1;
	      if( temp[1]< (midSide/2.0+gapY) && temp[1]>(midSide/2.0-gapY) && temp[0]<(midSide/2.0+gapX) && temp[0]>(-midSide/2.0-gapX))
		type[count]=1;
	      if( temp[1]> (-midSide/2.0-gapY) && temp[1]<(-midSide/2.0+gapY) && temp[0]<(midSide/2.0+gapX) && temp[0]>(-midSide/2.0-gapX))
		type[count]=1;
	      
	      double alfa = atan(1);
	      solidPoint[count][0]=temp[0]*cos(alfa)-temp[1]*sin(alfa);
	      solidPoint[count][1]=temp[1]*cos(alfa)+temp[0]*sin(alfa);
	      solidPoint[count][2]=temp[2];
	      count+=1;
	    }
	  
	}
      }
    }
	    

#endif

#if 0
    int nX=81, nY=400, nZ=1;
    double radius1=0., radius2=0.2;
    double gapR=(radius2-radius1)/(nX-1), gapAngle=2*3.1415926/nY, gapZ=1.0/nZ;
    
    const int nMPM = nX*nY*nZ;
    Array<VecD3> solidPoint(nMPM);
    Array<VecD3> solidVelocity(nMPM);
    Array<int> type(nMPM);
    Array<double> solidTemperature(nMPM);


    //polar coordinate 
     for(int i=0; i<nX; i++){
      for(int j=0; j<nY; j++){
	for(int k=0; k<nZ; k++){	 
	  solidPoint[count][0]=(radius1+i*gapR)*(cos(j*gapAngle))+center[0];
	  solidPoint[count][1]=(radius1+i*gapR)*(sin(j*gapAngle))+center[1];
	  solidPoint[count][2]=k*gapZ+center[2];
	  
	  if(i!=(nX-1)) type[count] = 0;  //internal particles
	  if(i==(nX-1)){
	    type[count] = 1;       //surface particles	   
	  }
	  count+=1;
	}
      }
    }

     /* 
    radius1=0.9, radius2=1.5;
    gapR=(radius2-radius1)/nX, gapAngle=2*3.1415926/nY, gapZ=1.0/nZ;

      //polar coordinate 
     for(int i=0; i<nX; i++){
      for(int j=0; j<nY; j++){
	for(int k=0; k<nZ; k++){	 
	  solidPoint[count][0]=(radius1+i*gapR)*(cos(j*gapAngle))+center[0];
	  solidPoint[count][1]=(radius1+i*gapR)*(sin(j*gapAngle))+center[1];
	  solidPoint[count][2]=k*gapZ+center[2];
	  
	  if(i!=0) type[count] = 0;  //internal particles
	  if(i==0){
	    type[count] = 1;       //surface particles	   
	  }
	  count+=1;
	}
      }
    }
     */
#endif
#if 1
    int nX=21,  nY=300, nZ=400;
    const double pi = atan(1.0)*4.0;
    double radius1=0., radius2=10;
    double gapR=(radius2-radius1)/(nX-1), gapAlfa=pi/nY, gapBeta=2*pi/nZ;
    const int nMPM = nX*nY*nZ;
    Array<VecD3> solidPoint(nMPM);
    Array<VecD3> solidVelocity(nMPM);
    Array<int> type(nMPM);
    Array<double> solidTemperature(nMPM);  
    //polar coordinate 
     for(int i=0; i<nX; i++){
      for(int j=0; j<nY; j++){
	for(int k=0; k<nZ; k++){	 
	  solidPoint[count][0]=(radius1+i*gapR)*sin(j*gapAlfa)*(cos(k*gapBeta))+center[0];
	  solidPoint[count][1]=(radius1+i*gapR)*sin(j*gapAlfa)*(sin(k*gapBeta))+center[1];
	  solidPoint[count][2]=(radius1+i*gapR)*cos(j*gapAlfa)+center[2];
	  
	  type[count] = 0;  //internal particles
	  if(i==nX-1){
	    type[count] = 1;       //surface particles	   
	  }
	  count+=1;
	}
      }
    }

#endif
    //set up particle velocity
#if 1
    for(int p=0; p<count; p++){
      solidVelocity[p][0]=0.0;
      solidVelocity[p][1]=0.0;
      solidVelocity[p][2]=0.0;
    }
#endif 

#if 0
    //set up rotating cylinder velcity
    const double angV = 1;
    for (int p=0; p<count; p++){
      double r = mag(solidPoint[p]-center);
      double angle = atan2(solidPoint[p][1]-center[1],solidPoint[p][0]-center[0]);
      solidVelocity[p][0] = -angV*r*sin(angle);
      solidVelocity[p][1] = angV*r*cos(angle);
      solidVelocity[p][2] = 0.0;
    }
#endif

    //set up temperature 
    const double InitT=300.0;
    for (int p=0; p<count; p++){
      solidTemperature[p]=InitT;
    }


    cout<<"count of particles is "<<count<<endl;
    //write out coordinate and velocity and particle type into file
    fp=fopen(file,"w");
    fprintf(fp,"%i\n",count);
    for(int p=0; p<count; p++){
      fprintf(fp, "%e\t%e\t%e\n", solidPoint[p][0],solidPoint[p][1],solidPoint[p][2]);
    } 
    for(int p=0; p<count; p++){
      fprintf(fp, "%e\t%e\t%e\n", solidVelocity[p][0],solidVelocity[p][1],solidVelocity[p][2]);
    }    
    for(int p=0; p<count; p++){
      fprintf(fp, "%i\n", type[p]);
    } 
    for(int p=0; p<count; p++){
      fprintf(fp, "%f\n", solidTemperature[p]);
    }
    fclose(fp);

}


const shared_ptr<Array<VecD3> > MPM::readCoordinates(const char *file)

{
    FILE *fp;
    int nMPM;
    double x=0, y=0, z=0;
   
    fp=fopen(file,"r");
    fscanf(fp,"%i\n",&nMPM);
    cout<<"number of particles is"<<nMPM<<endl;
    
    shared_ptr<Array<VecD3> > MPM_Points ( new Array<VecD3> (nMPM));
    //read in coordinate
    for(int i=0; i<nMPM; i++){
      fscanf(fp,"%lf\t%lf\t%lf\n", &x, &y, &z);
      (*MPM_Points)[i][0]=x;
      (*MPM_Points)[i][1]=y;
      (*MPM_Points)[i][2]=z;
    }
    fclose(fp);   
    return (MPM_Points);
}

const shared_ptr<Array<VecD3> > MPM::readVelocities(const char *file)

{
    FILE *fp;
    int nMPM;
    double vx=0, vy=0, vz=0;
    double x=0, y=0, z=0;
    fp=fopen(file,"r");
    fscanf(fp,"%i\n",&nMPM);
    
    shared_ptr<Array<VecD3> > MPM_Points ( new Array<VecD3> (nMPM));
    //read in cooridnate and skip
    for(int i=0; i<nMPM; i++){
      fscanf(fp,"%lf\t%lf\t%lf\n", &x, &y, &z);      
    }
    //read in velocity
    for(int i=0; i<nMPM; i++){
      fscanf(fp,"%lf\t%lf\t%lf\n", &vx, &vy, &vz);
      (*MPM_Points)[i][0]=vx;
      (*MPM_Points)[i][1]=vy;
      (*MPM_Points)[i][2]=vz;
    }
    fclose(fp);   
    return (MPM_Points);
}


const shared_ptr<Array<int> > MPM::readTypes(const char *file)

{
    FILE *fp;
    int nMPM;
    double x=0, y=0, z=0;
    int t=0;
    fp=fopen(file,"r");
    fscanf(fp,"%i\n",&nMPM);
    
    shared_ptr<Array<int> > MPM_Points ( new Array<int> (nMPM));
    //read in cooridnate and skip
    for(int i=0; i<nMPM; i++){
      fscanf(fp,"%lf\t%lf\t%lf\n", &x, &y, &z);      
    }
    //read in velocity and skip
    for(int i=0; i<nMPM; i++){
      fscanf(fp,"%lf\t%lf\t%lf\n", &x, &y, &z);      
    }
    //read in type
    for(int i=0; i<nMPM; i++){
      fscanf(fp,"%i\n", & t);
      (*MPM_Points)[i]=t;
    }
    fclose(fp);   
    return (MPM_Points);
}
const shared_ptr<Array<double> > MPM::readTemperatures(const char *file)

{
    FILE *fp;
    int nMPM;
    double x=0, y=0, z=0;
    int t=0;
    double temperature=0.0;
    fp=fopen(file,"r");
    fscanf(fp,"%i\n",&nMPM);
    
    shared_ptr<Array<double> > MPM_Points ( new Array<double> (nMPM));
    //read in cooridnate and skip
    for(int i=0; i<nMPM; i++){
      fscanf(fp,"%lf\t%lf\t%lf\n", &x, &y, &z);      
    }
    //read in velocity and skip
    for(int i=0; i<nMPM; i++){
      fscanf(fp,"%lf\t%lf\t%lf\n", &x, &y, &z);      
    }
    //read in type
    for(int i=0; i<nMPM; i++){
      fscanf(fp,"%i\n", & t);
    }
    //read in temperature
    for(int i=0; i<nMPM; i++){
      fscanf(fp,"%lf\n", & temperature);
      (*MPM_Points)[i]=temperature;
    }
    fclose(fp);   
    return (MPM_Points);
}

void MPM::Init(const shared_ptr<Array<VecD3> > coordinates,
	       const shared_ptr<Array<VecD3> > velocities,
	       const shared_ptr<Array<int> > types,
	       const shared_ptr<Array<double> > temperatures)
{

  const int n = (*coordinates).getLength();  //number of particles
  _particles.setCount(n);
  
  MPM::setCoordinates(coordinates);
  MPM::setVelocities(velocities);
  MPM::setTypes(types);
  MPM::setTemperatures(temperatures);

}


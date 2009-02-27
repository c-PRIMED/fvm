#include "MPM_Particles.h"

#include "StorageSite.h"
#include "CRConnectivity.h"


typedef MPM::VecD3 VecD3;

MPM::MPM():
  _particles(0),
  _coordinates()
{}

MPM::~MPM() { }


void MPM::setandwriteParticles(const char *file)
{
 //here, we want to build up a solid circle 
 //the number of solid points and coordiantes are written to file
    FILE *fp;
           
    VecD3 center;
    center[0]=0.5;
    center[1]=0.5;
    center[2]=0.0;

    int count=0;   
   
#if 1
    //set up particle cartesian coordinate
    int nX=100, nY=20, nZ=1;
    double gapX=0.5/nX, gapY=0.2/nY, gapZ=1.0/nZ;
    VecD3 solidPoint[nX*nY*nZ];
    VecD3 solidVelocity[nX*nY*nZ];

    for(int i=0; i<nX; i++){
      for(int j=0; j<nY; j++){
	for(int k=0; k<nZ; k++){	  
	  solidPoint[count][0]=i*gapX+0.25;
	  solidPoint[count][1]=j*gapY+0.4;
	  solidPoint[count][2]=k*gapZ;
	  count+=1;
	}
      }
    }

#endif

#if 0
    int nX=20, nY=200, nZ=1;
    double radius1=0., radius2=0.2;
    double gapR=(radius2-radius1)/nX, gapAngle=2*3.1415926/nY, gapZ=1.0/nZ;
    VecD3 solidPoint[nX*nY*nZ];
    VecD3 solidVelocity[nX*nY*nZ];

    //polar coordinate 
     for(int i=1; i<nX; i++){
      for(int j=0; j<nY; j++){
	for(int k=0; k<nZ; k++){	 
	  solidPoint[count][0]=(radius1+i*gapR)*(cos(j*gapAngle))+center[0];
	  solidPoint[count][1]=(radius1+i*gapR)*(sin(j*gapAngle))+center[1];
	  solidPoint[count][2]=k*gapZ+center[2];
	  count+=1;
	}
      }
    }

#endif

    //set up particle velocity
    for(int p=0; p<count; p++){
      solidVelocity[p][0]=0.0;
      solidVelocity[p][1]=0.0;
      solidVelocity[p][2]=0.0;
    }
    //write out coordinate and velocity into file
    fp=fopen(file,"w");
    fprintf(fp,"%i\n",count);
    for(int p=0; p<count; p++){
      fprintf(fp, "%lf\t%lf\t%lf\n", solidPoint[p][0],solidPoint[p][1],solidPoint[p][2]);
    } 
    for(int p=0; p<count; p++){
      fprintf(fp, "%lf\t%lf\t%lf\n", solidVelocity[p][0],solidVelocity[p][1],solidVelocity[p][2]);
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

void MPM::Init(const shared_ptr<Array<VecD3> > coordinates,
	       const shared_ptr<Array<VecD3> > velocities )
{

  const int n = (*coordinates).getLength();  //number of particles
  _particles.setCount(n);
  
  MPM::setCoordinates(coordinates);
  MPM::setVelocities(velocities);

}






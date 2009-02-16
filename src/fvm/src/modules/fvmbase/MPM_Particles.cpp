#include "MPM_Particles.h"

#include "StorageSite.h"
#include "CRConnectivity.h"


typedef MPM::VecD3 VecD3;

MPM::MPM():
  _particles(0),
  _coordinates()
{}

MPM::~MPM() { }


void MPM::setandwriteParticles(char *file)
{
 //here, we want to build up a solid circle 
 //the number of solid points and coordiantes are written to file
    FILE *fp;
   
    int nX=50, nY=50, nZ=1;
    double gapX=1.0/nX, gapY=1.0/nY, gapZ=1.0/nZ;
    double radius=0.2;
    VecD3 center;
    center[0]=0.5;
    center[1]=0.5;
    center[2]=0.0;

    int count=0;
    VecD3 temp;
    VecD3 solidPoint[nX*nY*nZ];

    VecD3 solidVelocity[nX*nY*nZ];
    //set up particle coordinate
    for(int i=0; i<nX; i++){
      for(int j=0; j<nY; j++){
	for(int k=0; k<nZ; k++){
	  temp[0]=i*gapX;
	  temp[1]=j*gapY;
	  temp[2]=k*gapZ;
	  VecD3 ds=temp-center;
	  if(mag2(ds) <= radius*radius){
	    solidPoint[count][0]=temp[0];
	    solidPoint[count][1]=temp[1];
	    solidPoint[count][2]=temp[2];
	   
	    count+=1;
	  }
	}
      }
    }
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


const shared_ptr<Array<VecD3> > MPM::readCoordinates(char *file)

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

const shared_ptr<Array<VecD3> > MPM::readVelocities(char *file)

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

void MPM::Init(shared_ptr<Array<VecD3> > coordinates,
	       shared_ptr<Array<VecD3> > velocities )
{

  const int n = (*coordinates).getLength();  //number of particles
  _particles.setCount(n);
  
  MPM::setCoordinates(coordinates);
  MPM::setVelocities(velocities);

}






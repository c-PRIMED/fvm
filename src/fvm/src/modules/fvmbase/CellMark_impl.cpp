// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "CellMark_impl.h"
#include "Mesh.h"
#include <iostream>
#include <string>


typedef Octree::Point Point;
typedef Octree::Bounds Bounds;

typedef Vector<double, 3> VecD3;
typedef Array<VecD3> VecD3Array;

void CellMark_Impl(Mesh& mesh, const GeomFields& geomFields, const string fileBase, 
		   Octree& O, MPM& solid, const int option)

{

    const StorageSite& cells = mesh.getCells();
    const int nCells = cells.getCount();          //need to include BC cells?
    const int nSelfCells = cells.getSelfCount();
    const VecD3Array& cellCentroid = dynamic_cast<const VecD3Array& > (geomFields.coordinate[cells]);
    const StorageSite& faces = mesh.getFaces();
    const int nFaces = faces.getCount();
    const CRConnectivity& faceCells = mesh.getAllFaceCells();
    const CRConnectivity& cellFaces = mesh.getCellFaces();
    const CRConnectivity& cellCells = mesh.getCellCells();
    const VecD3Array& faceArea =
          dynamic_cast<const VecD3Array&>(geomFields.area[faces]);
    const VecD3Array& faceCentroid = 
          dynamic_cast<const VecD3Array&> (geomFields.coordinate[faces]);
   
   
    const int writeOption = 0;
    /***************************************************************************************/
    //---build up Octree for cell centroid---//
    
    /*
    Point *points = new Point [nCells];

    for(int i=0; i<nCells; i++){
      points[i].coordinate=cellCentroid[i];
      points[i].cellIndex=i;
      points[i].code=0;
    }   
    
    //Octree build up parameters

    const int threshold=1;
    const int maxDepth=20;
    const int count=nCells;
    int currentDepth=0;
   
    //define a new Octree object

    //  Octree O; 

    //calculate entire domain bounds

    Bounds bounds=O.calcCubicBounds(points, count);

    //build Octree
    O.build(points, count, threshold, maxDepth, bounds, currentDepth);

#if 0
    //Report Octree    
    const string fileName1=fileBase+"Octree_report.dat";
    char* file1;
    file1=&fileName1[0];
    FILE *fp;
    fp=fopen(file1,"w");
    O.report(fp);
    fclose(fp);
#endif
    
    */
    //---so far, we have the Octree ready for search---//


   
    /**************************************************************************************/
    /*
    //---set up MPM particles---//
    string fileName2=fileBase+"MPMs.dat";
    char* file2;
    file2=&fileName2[0];
    

    // MPM solid;
    //initialize particle velocity and coordinate and write into file
    solid.setandwriteParticles(file2);

    //get coordinate
    const shared_ptr<VecD3Array> MPM_Coordinates = solid.readCoordinates(file2);

    //get velocity
    const shared_ptr<VecD3Array> MPM_Velocities = solid.readVelocities(file2);

    //get type
    const shared_ptr<Array<int> > MPM_Types = solid.readTypes(file2);

    //store all the information in MPM class
    solid.Init(MPM_Coordinates, MPM_Velocities, MPM_Types);
    */



    const StorageSite& particles = solid.getParticles();
    
    const int nMPM = particles.getCount();

    //cout<<"nMPM is "<<nMPM<<endl;
   
    const shared_ptr<VecD3Array>& MPM_Points = solid.getCoordinates();
    const shared_ptr<VecD3Array>& MPM_Vels = solid.getVelocities();
    const shared_ptr<Array<int> >& particleTypes = solid.getTypes();
    
    

    cout<<"count of cells is "<<nCells<<endl;
    /***************************************************************************************/
    //---Find out each solid point falls to which cells---//
    
    //a temporary array to store point to cell connectivity
    Array<int> MPM_PointstoCells (nMPM);
   
    //search options:
    //option 1: search the nearest cell to this particle and see if this particle falls to 
    //--->this cell; if not, search the neighbors of this cell, until find the cell which 
    //--->contains this particle
    //option 2: search the cells within a radius to this particle using Octree,
    //---> this radius is pre-estimated so that it must include the cell which contains 
    //---> the particle. Loop over the cells within this radius and find the cell containi//ng the particle
    //option 3: naive search over all cells to find out the cells within the radius, then
    //---> find out which cell contains the particle
    //option 4: naive search all the cells and find out which cell contains the particle
   
    if (option == 1){
      for(int p=0; p<nMPM; p++){
	const VecD3 MPM_point = (*MPM_Points)[p];
	MPM_PointstoCells[p]=-1;
	const int nearestCell = O.getNode(MPM_point);
	int flag=0;
	int inCellorNot=inCell(nearestCell, MPM_point, faceCells, cellFaces,
			       faceArea,faceCentroid);
	if (inCellorNot==1) {
	  flag = 1;
	  MPM_PointstoCells[p]=nearestCell;
	}
	int levelCount = 0;
	while (flag == 0 && levelCount <= 2){
	  levelCount ++;
	  const int nc = cellCells.getCount(nearestCell);
	  for (int c = 0; c < nc; c ++){
	    int cellCandidate = cellCells(nearestCell, c);
	    inCellorNot=inCell(cellCandidate, MPM_point, faceCells, cellFaces,
			       faceArea,faceCentroid);
	    if (inCellorNot==1) {
	      flag = 1;
	      MPM_PointstoCells[p]=cellCandidate;
	    }
	    // if (inCellorNot == 0){
	    //flag = 1;
	    //}
	  }	
	}
      }
    }
	


    if (option == 2){
      const double radius = 2.5;
      for(int p=0; p<nMPM; p++){
	const VecD3 MPM_point = (*MPM_Points)[p];
	vector<int> cellIndexList;
	O.getNodes(MPM_point, radius, cellIndexList);
	MPM_PointstoCells[p]=-1;
	for (int i=0; i< (int) cellIndexList.size(); i++)    {
	  int cellCandidate = cellIndexList[i];
	  const int inCellorNot=inCell(cellCandidate, MPM_point, faceCells, cellFaces,
				       faceArea,faceCentroid);
	  if (inCellorNot==1){
	    MPM_PointstoCells[p]=cellCandidate;
	  }
	}
      }
    }
    /*
    if (option == 3){
       double radius = 0.03;
       for(int p=0; p<nMPM; p++){
	 VecD3 MPM_point = (*MPM_Points)[p];
	 vector<int> cellIndexList=O.Naive_getNodes(MPM_point, count, points, radius);
	 MPM_PointstoCells[p]=-1;
	 for (int i=0; i< (int) cellIndexList.size(); i++)    {
	   int cellCandidate = cellIndexList[i];
	   const int inCellorNot=inCell(cellCandidate, MPM_point, faceCells, cellFaces,
				       faceArea,faceCentroid);
	   if (inCellorNot==1){
	     MPM_PointstoCells[p]=cellCandidate;
	   }
	 }
       }
    }
    */
    if (option == 4){
       for(int p=0; p<nMPM; p++){
	 VecD3 MPM_point = (*MPM_Points)[p];
	 MPM_PointstoCells[p]=-1;
	 for(int i=0; i< nCells; i++)    {
	    const int inCellorNot=inCell(i, MPM_point, faceCells, cellFaces,
				       faceArea,faceCentroid);
	    if (inCellorNot==1){
	      MPM_PointstoCells[p]=i;
	    }
	 }
       }
    }

	  
if (writeOption ==1)
  {
    FILE *fp;
    string fileName20=fileBase+"particles.dat";
    char* file20;
    file20=&fileName20[0];
    
    fp=fopen(file20,"w");
    for (int c=0; c<nMPM; c++){
      if((*particleTypes)[c]==1){
	fprintf(fp, "%i\t%e\t%e\t%e\n", c, (*MPM_Points)[c][0],(*MPM_Points)[c][1],(*MPM_Points)[c][2]);
      }
    }
  } 

    /************************************************************************************/
    //---create the CRconnectivity for solid and cells---//

    shared_ptr<CRConnectivity> particleCellsCR  =
         setParticleCells(particles, cells, MPM_PointstoCells);
    
    //store the CRconnectivity to mesh class

    mesh.setConnectivity( particles, cells, particleCellsCR);

    const CRConnectivity& particleCells = mesh.getConnectivity(particles, cells); 
    
if (writeOption ==1)
  {
    FILE *fp;
    string fileName14=fileBase+"particletocells.dat";
    char* file14;
    file14=&fileName14[0];
    
    fp=fopen(file14,"w");
  
    //test
    for(int c=0; c<nMPM; c++){
      //for each solid point, find out how many cells have connectivity to it
      int np = particleCells.getCount(c);
      //what they are
      for(int p=0; p<np; p++){
	fprintf(fp, "%i\t%i\n", c, particleCells(c, p));
      }
    } 
      fclose(fp);
  } 
    /************************************************************************************/
    //---create the CRconnectivity for cells and solid---//

    shared_ptr<CRConnectivity> cellParticlesCR = (*particleCellsCR).getTranspose();
    
    //store the CRconnectivity to mesh class

    mesh.setConnectivity(cells, particles, cellParticlesCR);

    const CRConnectivity& cellParticles = mesh.getConnectivity(cells, particles);


if (writeOption ==1)
  {   
    FILE *fp;
    string fileName13=fileBase+"celltoparticles.dat";
    char* file13;
    file13=&fileName13[0];
    
    fp=fopen(file13,"w");
    //test
    for(int c=0; c<nCells; c++){
      //for each solid point, find out how many cells have connectivity to it
      int np = cellParticles.getCount(c);
      //what they are
      for(int p=0; p<np; p++){
	//if (cellParticles(c,p) < 0)
	
	int pID = cellParticles(c,p);
	if ((*particleTypes)[pID] == 1)
	fprintf(fp,"%i\t%e\t%e\t%e\t%i\n", c, (*MPM_Points)[pID][0],(*MPM_Points)[pID][1],(*MPM_Points)[pID][2],(*particleTypes)[pID]);
      }
    } 
    fclose(fp);
  }
  
   

    /**************************mark cell**********************************/
    markCell ( mesh, nCells,nSelfCells, cellParticles, cellCells);

    //test mark cell
    
    if (writeOption ==1)
      reportCellMark (mesh, nCells, cellCentroid, fileBase);

 

    /***************************mark IBFaces********************/

    markIBFaces (mesh, nFaces, faceCells);
      
    const StorageSite&  ibFaces = mesh.getIBFaces();
   
    const Array<int>& ibFaceList = mesh.getIBFaceList();

if (writeOption ==1)
  {
    FILE *fp;
    string fileName15=fileBase+"ibfaces.dat";
    char* file15;
    file15=&fileName15[0];
    
    fp=fopen(file15,"w");
    //test
    for(int f=0; f<ibFaces.getCount(); f++){
      int fID = ibFaceList[f];
      fprintf(fp,"%i\t%e\t%e\t%e\n", fID, faceCentroid[fID][0],faceCentroid[fID][1],faceCentroid[fID][2]);
    }
    fclose(fp);
  }
    /*********check ibFaces***********************/
    checkIBFaces(ibFaceList, faceArea, faceCells, mesh);


    /************create ibFace to Particle connectivity****************/

    shared_ptr<CRConnectivity> ibFaceParticlesCR = 
      setibFaceParticles (mesh, ibFaces, ibFaceList, particles,faceCells, cellParticles, cellCells,  *particleTypes);

    //store the connectivity to mesh

    mesh.setConnectivity(ibFaces, particles, ibFaceParticlesCR);

    const CRConnectivity& ibFaceParticles = mesh.getConnectivity(ibFaces, particles);

   
    /************************create ibFace to Cells connectivity**********/

     shared_ptr<CRConnectivity> ibFaceCellsCR = 
       setibFaceCells (mesh, ibFaceList, ibFaces,  cells, faceCells, cellFaces, faceCentroid);
  
  //store the CRconnectivity to mesh class

     mesh.setConnectivity(ibFaces, cells, ibFaceCellsCR);   

     const CRConnectivity& ibFaceCells = mesh.getConnectivity(ibFaces, cells);

if (writeOption ==1)
  {   
    FILE *fp;
     string fileName11=fileBase+"ibfacetoparticle.dat";
     char* file11;
     file11=&fileName11[0];
   
     fp=fopen(file11,"w");
     //test
     for(int f=0; f<ibFaces.getCount(); f++){
       const int faceIndex = ibFaceList[f];
       //cout << faceCentroid[faceIndex] << endl;
      //for each ibface find out how many fluid cells have connectivity to it
      int nc = ibFaceParticles.getCount(f);
      //what they are
      for(int c=0; c<nc; c++){
	int pID = ibFaceParticles(f,c);
	//cout<<"face "<<faceIndex<<"  has  "<<ibFaceParticles(f, c)<<endl;
	//cout<<(*MPM_Points)[ibFaceParticles(f,c)]<<endl;
	fprintf(fp, "%i\t%i\t%e\t%e\t%e\t%i\n", faceIndex,pID, (*MPM_Points)[pID][0],(*MPM_Points)[pID][1],(*MPM_Points)[pID][2],(*particleTypes)[pID]);
	//fprintf(fp, "%i\t%i\t%lf\t%lf\t%lf\t%i\n", faceIndex,pID, (*MPM_Vels)[pID][0],(*MPM_Vels)[pID][1],(*MPM_Vels)[pID][2],(*particleTypes)[pID]);
      }
     } 
    

     string fileName12=fileBase+"ibfacetocell.dat";
     char* file12;
     file12=&fileName12[0];
    
     fp=fopen(file12,"w");
     //test
     for(int f=0; f<ibFaces.getCount(); f++){
       const int faceIndex = ibFaceList[f];
       //cout << faceCentroid[faceIndex] << endl;
      //for each ibface find out how many fluid cells have connectivity to it
      int nc = ibFaceCells.getCount(f);
      //what they are
      for(int c=0; c<nc; c++){
	//cout<<"face "<<faceIndex<<"  has  "<<ibFaceCells(f, c)<<endl;
	//cout<<(cellCentroid)[ibFaceCells(f,c)]<<endl;
	int cID = ibFaceCells(f,c);
	fprintf(fp, "%i\t%i\t%e\t%e\t%e\n",faceIndex,cID, (cellCentroid)[cID][0],(cellCentroid)[cID][1],(cellCentroid)[cID][2]);
      }
     }
     fclose(fp);
  }
    
 
}

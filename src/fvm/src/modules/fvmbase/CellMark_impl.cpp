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
    const int nCells = cells.getSelfCount();          //need to include BC cells?
    const VecD3Array& cellCentroid = dynamic_cast<const VecD3Array& > (geomFields.coordinate[cells]);
    const StorageSite& faces = mesh.getFaces();
    const CRConnectivity& faceCells = mesh.getAllFaceCells();
    const CRConnectivity& cellFaces = mesh.getCellFaces();
    const CRConnectivity& cellCells = mesh.getCellCells();
    const VecD3Array& faceArea =
          dynamic_cast<const VecD3Array&>(geomFields.area[faces]);
    const VecD3Array& faceCentroid = 
          dynamic_cast<const VecD3Array&> (geomFields.coordinate[faces]);
   
   

  
    /***************************************************************************************/
    //---build up Octree for cell centroid---//

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
    

    //---so far, we have the Octree ready for search---//


   
    /**************************************************************************************/
    //---set up MPM particles---//
    string fileName2=fileBase+"Output/MPM_ring.dat";
    char* file2;
    file2=&fileName2[0];
    

    // MPM solid;
    //initialize particle velocity and coordinate and write into file
    solid.setandwriteParticles(file2);

    //get coordinate
    const shared_ptr<VecD3Array> MPM_Coordinates = solid.readCoordinates(file2);

    //get velocity
    const shared_ptr<VecD3Array> MPM_Velocities = solid.readVelocities(file2);

    //store all the information in MPM class
    solid.Init(MPM_Coordinates, MPM_Velocities);
    
    const StorageSite& particles = solid.getParticles();
    
    const int nMPM = particles.getCount();

    //cout<<"nMPM is "<<nMPM<<endl;

    const shared_ptr<VecD3Array>& MPM_Points = solid.getCoordinates();
    


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
    //---> the particle. Loop over the cells within this radius and find the cell containing the particle
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

	while (flag == 0){
	  const int nc = cellCells.getCount(nearestCell);
	  for (int c = 0; c < nc; c ++){
	    int cellCandidate = cellCells(nearestCell, c);
	    inCellorNot=inCell(cellCandidate, MPM_point, faceCells, cellFaces,
			       faceArea,faceCentroid);
	    if (inCellorNot==1) {
	      flag = 1;
	      MPM_PointstoCells[p]=cellCandidate;
	    }
	    if (inCellorNot == 0){
	      flag = 1;
	    }
	  }	
	}
      }
    }
	


    if (option == 2){
      const double radius = 0.05;
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

	  
    /************************************************************************************/
    //---create the CRconnectivity for solid and cells---//

    shared_ptr<CRConnectivity> particleCellsCR  =
         setParticleCells(particles, cells, MPM_PointstoCells);
    
    //store the CRconnectivity to mesh class

    mesh.setConnectivity( particles, cells, particleCellsCR);

    //    const CRConnectivity& particleCells = mesh.getConnectivity(particles, cells); 
   
    /*
    //test
    for(int c=0; c<nMPM; c++){
      //for each solid point, find out how many cells have connectivity to it
      int np = particleCells.getCount(c);
      //what they are
      for(int p=0; p<np; p++){
	if (particleCells(c,p) < 0)
	cout<<"cell "<<c<<"  has  "<<particleCells(c,p)<<endl;
      }
    } 
    */
   
    /************************************************************************************/
    //---create the CRconnectivity for cells and solid---//

    shared_ptr<CRConnectivity> cellParticlesCR = (*particleCellsCR).getTranspose();
    
    //store the CRconnectivity to mesh class

    mesh.setConnectivity(cells, particles, cellParticlesCR);

    const CRConnectivity& cellParticles = mesh.getConnectivity(cells, particles);

    /*
    //test
    for(int c=0; c<nCells; c++){
      //for each solid point, find out how many cells have connectivity to it
      int np = cellParticles.getCount(c);
      //what they are
      for(int p=0; p<np; p++){
	if (cellParticles(c,p) < 0)
	cout<<"cell "<<c<<"  has  "<<cellParticles(c,p)<<endl;
      }
    } 
    */
  
   

    /**************************mark cell**********************************/
    markCell ( mesh, nCells, cellParticles, cellCells);

    //test mark cell
    
    reportCellMark (mesh, nCells, cellCentroid, fileBase);

 
#if 0
    /***************************mark IBFaces********************/

    markIBFaces (mesh, nCells, cellCells, cellFaces, faceCells);
   

    int count=0;
    //double check ibFaceCount
    for (f=0; f<nFaces; f++){
      int c0 = faceCells(f, 0);
      int c1 = faceCells(f, 1);
      int ibType0 = mesh.getIBTypeForCell(c0);
      int ibType1 = mesh.getIBTypeForCell(c1);
      if(!(ibType0==ibType1))
	count++;
    }



    const StorageSite&  ibFaces = mesh.getIBFaces();
   
    const Array<int>& ibFaceList = mesh.getIBFaceList();

    /************create ibFace to Particle connectivity****************/

    shared_ptr<CRConnectivity> ibFaceParticlesCR = 
      setibFaceParticles (mesh, ibFaces, ibFaceList, particles,faceCells, cellParticles);

    //store the connectivity to mesh

    mesh.setConnectivity(ibFaces, particles, ibFaceParticlesCR);

    const CRConnectivity& ibFaceParticles = mesh.getConnectivity(ibFaces, particles);

   
    /************************create ibFace to Cells connectivity**********/

     shared_ptr<CRConnectivity> ibFaceCellsCR = 
       setibFaceCells (mesh, ibFaceList, ibFaces,  cells, faceCells, O, faceCentroid);
  
  //store the CRconnectivity to mesh class

     mesh.setConnectivity(ibFaces, cells, ibFaceCellsCR);   

     //     const CRConnectivity& ibFaceCells = mesh.getConnectivity(ibFaces, cells);

     /*
     //test
     for(int f=0; f<ibFaces.getCount(); f++){
         //       const int faceIndex = ibFaceList[f];
       //cout << faceCentroid[faceIndex] << endl;
      //for each ibface find out how many fluid cells have connectivity to it
      int nc = ibFaceParticles.getCount(f);
      //what they are
      for(int c=0; c<nc; c++){
	//cout<<"face "<<faceIndex<<"  has  "<<ibFaceCells(f, c)<<endl;
	cout<<(*MPM_Points)[ibFaceParticles(f,c)]<<endl;
      }
    } 
    
     //test
     for(int f=0; f<ibFaces.getCount(); f++){
       const int faceIndex = ibFaceList[f];
       //cout << faceCentroid[faceIndex] << endl;
      //for each ibface find out how many fluid cells have connectivity to it
      int nc = ibFaceCells.getCount(f);
      //what they are
      for(int c=0; c<nc; c++){
	//cout<<"face "<<faceIndex<<"  has  "<<ibFaceCells(f, c)<<endl;
	cout<<(cellCentroid)[ibFaceCells(f,c)]<<endl;
      }
      } 
     */
    
#endif
   
 
}

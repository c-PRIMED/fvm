// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "CellMark.h"
typedef Vector<double, 3> VectorT3;
typedef Array<Vector<double, 3> > VectorT3Array;


int 
inCell(const int cellIndex, 
       const VectorT3& point, 
       const CRConnectivity& faceCells,
       const CRConnectivity& cellFaces,
       const VectorT3Array& faceArea,
       const VectorT3Array& faceCentroid)
{
 
  
  
  int faceNumber=cellFaces.getCount(cellIndex);
  int flag[faceNumber];
  int throwFlag=0;
  int sum=0;
  VectorT3 Af;
 
  //const Array<int>& cellFacesRow = cellFaces.getRow();
  //const Array<int>& cellFacesCol = cellFaces.getCol();

  //const Array<int>& faceCellsRow = faceCells.getRow();
  //const Array<int>& faceCellsCol = faceCells.getCol();
  
  // cout<<"cell index is"<<cellIndex<<endl;
 
  for (int nf=0; nf<faceNumber; nf++){
    const int f=cellFaces(cellIndex, nf);
   
    //first, judge c0 or c1 to define orientation
    const int c0=faceCells(f,0);
    //const int c1=faceCells(f,1);       
    if (cellIndex==c0){
     Af=(-faceArea[f]);
    }
    else {
     Af=(faceArea[f]);
    }  
  
    //calculate product    
    VectorT3 ds=point-faceCentroid[f];
    double product=dot(Af,ds);
    
    if(product > 0.0){ flag[nf]=1;}
    else if (product < 0.0) {flag[nf]=-1;}
    else {
      //cout<<cellIndex<<endl;
      throwFlag = 1;
    }
    sum+=flag[nf];
  }
  //particle is on face or vertex, throw it away
  if (throwFlag == 1){
    return (0);
  }
  //particle is in or out of cell
  else{
    if (sum==faceNumber){
      return (1);   
    }
    else return (-1);
  }
}





void markCell( Mesh& mesh, const int nCells, const int nSelfCells,
	       const CRConnectivity& cellParticles, const CRConnectivity& cellCells )
{

  //step 1: Mark the cells with solid particles as solid
  //and Mark the cells with no solid particles as fluid

   for(int c=0; c<nCells; c++){
     const int particleCount = cellParticles.getCount(c);
     //if no particle in cell, mark as fluid
     if (particleCount == 0) {
       mesh.setIBTypeForCell(c,Mesh::IBTYPE_FLUID);
     }
     //if has particle in cell, mark as solid
     else { mesh.setIBTypeForCell(c,Mesh::IBTYPE_SOLID); }
   }
   //step2: in solid cells, mark cells with no fluid neighbors as solid
   //and mark cells with at least one fluid neighbors as IB cells

    for (int c=0; c<nCells; c++){
      const int ibType = mesh.getIBTypeForCell(c);
      int flag;
      //search all solid cells
      if(ibType == Mesh::IBTYPE_SOLID){
	flag=1;  //true for solid cells
	const int ncNumber=cellCells.getCount(c);
	for(int nc=0; nc<ncNumber; nc++){
	  const int cellIndex=cellCells(c,nc);
	   if(mesh.getIBTypeForCell(cellIndex)==Mesh::IBTYPE_FLUID && cellIndex < nSelfCells){ 
	     //if(mesh.getIBTypeForCell(cellIndex)==Mesh::IBTYPE_FLUID){
	    flag=0;   
	  }
	}
	//if solid cell has at least one fluid cell neighbor, mark as IBM type
	if(flag==0) mesh.setIBTypeForCell(c,Mesh::IBTYPE_BOUNDARY);
      }
    }
   
}


void reportCellMark (const Mesh& mesh, const int nCells,
		     const VectorT3Array& cellCentroid,
		     const string fileBase)
{
    
    string fileName=fileBase+"CellMark.dat";
    char* file;
    file=&fileName[0];
    FILE *fp=fopen(file,"w");
    
    for(int c=0; c<nCells; c++){
      int ibType = mesh.getIBTypeForCell(c);
      fprintf(fp, "%i\t%i\n", c, ibType);
    } 
    fclose(fp);  

    string fileName1 = fileBase+"FluidCell.dat";
    string fileName2 = fileBase+"IBMCell.dat";
    string fileName3 = fileBase+"SolidCell.dat";
    char* file1;
    char* file2;
    char* file3;
    file1=&fileName1[0];
    file2=&fileName2[0];
    file3=&fileName3[0];
    FILE *fp1, *fp2, *fp3;
    fp1=fopen(file1,"w");
    fp2=fopen(file2,"w");
    fp3=fopen(file3,"w");
   
    for(int c=0; c<nCells; c++){
      int ibType = mesh.getIBTypeForCell(c);
      if(ibType == Mesh::IBTYPE_FLUID){
	fprintf(fp1, "%i\t%e\t%e\t%e\n", c, cellCentroid[c][0],cellCentroid[c][1],cellCentroid[c][2]);
      }
      else if(ibType==Mesh::IBTYPE_BOUNDARY){
	fprintf(fp2, "%i\t%e\t%e\t%e\n", c, cellCentroid[c][0],cellCentroid[c][1],cellCentroid[c][2]);
      }
      else if(ibType==Mesh::IBTYPE_SOLID){
	fprintf(fp3, "%i\t%e\t%e\t%e\n", c, cellCentroid[c][0],cellCentroid[c][1],cellCentroid[c][2]);
      }
    } 
    fclose(fp1);  
    fclose(fp2);  
    fclose(fp3);  
}


void markIBFaces(Mesh& mesh, const int nFaces, 		
		 const CRConnectivity& faceCells)
{
  //definition of ibFaces: the faces between IB cells and Fluid cells
  //first, count the number of ibFaces
    int ibFaceCount=0;
    for(int f=0; f<nFaces; f++){
      const int c0 = faceCells(f, 0);
      const int c1 = faceCells(f, 1);
      const int type0 = mesh.getIBTypeForCell(c0);
      const int type1 = mesh.getIBTypeForCell(c1);
      if(type0 == Mesh::IBTYPE_FLUID && type1 ==  Mesh::IBTYPE_BOUNDARY)
	ibFaceCount++;
      if(type1 == Mesh::IBTYPE_FLUID && type0 ==  Mesh::IBTYPE_BOUNDARY)
	ibFaceCount++;
    }
    cout<<"ibFaceCount is "<<ibFaceCount<<endl;
      
    //then, allocate an array for ibFace
    mesh.createIBFaceList(ibFaceCount);

    //insert the entries to ibface array
    ibFaceCount=0;
    for(int f=0; f<nFaces; f++){
      const int c0 = faceCells(f, 0);
      const int c1 = faceCells(f, 1);
      const int type0 = mesh.getIBTypeForCell(c0);
      const int type1 = mesh.getIBTypeForCell(c1);
      if(type0 == Mesh::IBTYPE_FLUID && type1 ==  Mesh::IBTYPE_BOUNDARY){
	mesh.addIBFace(ibFaceCount, f);
	ibFaceCount++;
      }
      if(type1 == Mesh::IBTYPE_FLUID && type0 ==  Mesh::IBTYPE_BOUNDARY){
	mesh.addIBFace(ibFaceCount, f);
	ibFaceCount++;
      }
    }


    //initialize storagesite ibFaces
    StorageSite&  ibFaces = mesh.getIBFaces();
    ibFaces.setCount(ibFaceCount);

}


void checkIBFaces(const Array<int> & ibFaceList,
		  const VectorT3Array& faceArea,
		  const CRConnectivity& faceCells,
		  const Mesh& mesh)
{

 // check if ibFaces form a closed curve //
   
  VectorT3 areaSum;
  areaSum[0] = 0.0;
  areaSum[1] = 0.0;
  areaSum[2] = 0.0;
  for(int i=0; i<ibFaceList.getLength();i++){
    const int fID = ibFaceList[i];
    const int C0 = faceCells(fID, 0);
    const int C1 = faceCells(fID, 1);
    if(mesh.getIBTypeForCell(C0)==Mesh::IBTYPE_FLUID
       && mesh.getIBTypeForCell(C1)==Mesh::IBTYPE_SOLID){
      areaSum += faceArea[fID];
    }
    else {
      areaSum += faceArea[fID];
    }
  }
  cout<<"sum of ibFace area is  "<<areaSum<<endl;
}


   

const shared_ptr<CRConnectivity> setibFaceParticles 
                          (const Mesh& mesh,
			   const StorageSite& ibFaces, 
			   const Array<int>& ibFaceList,
			   const StorageSite& particles,
			   const CRConnectivity& faceCells, 
			   const CRConnectivity& cellParticles,
			   const CRConnectivity& cellCells,
			   const Array<int>& particleType)
{

  //CR connectivity cellParticles includes all the particles located in each cell
  //here, we want to create ibFace to Particles in which only the surface particles are counted in
  //surface particles are identified by particle type 1
					       
  shared_ptr<CRConnectivity> ibFaceParticles (new CRConnectivity (ibFaces, particles));
  int maxcount = 0;
  int mincount = 1000;
  //initCount: new Array for row
  (*ibFaceParticles).initCount();
  const int rowSize = ibFaces.getCount();

  //specify the number of nonzeros for each row
  
  for(int p=0; p<rowSize; p++){
    const int faceIndex = ibFaceList [p];
    int C0 = faceCells(faceIndex, 0);
    int C1 = faceCells(faceIndex, 1);

    if (mesh.getIBTypeForCell(C1) == Mesh::IBTYPE_BOUNDARY)     
      {  C0 = C1; }
    //C0 is IBtype, C1 is fluid
   
    int nP = cellParticles.getCount(C0);
    int count=0;
    for(int n=0; n<nP; n++){
      int pID = cellParticles(C0, n);
      if(particleType[pID] == 1){
	count++;
      }
    }
#if 1
    //assuming only one fluid cell is used in interpolation, then at least three particles are needed to 
    //apply the linear least square method. If the current IB cell has less than three particles, then search 
    //neighbors for more particles
    if(count < 3){
      const int nbSize = cellCells.getCount(C0);
      for(int c=0; c < nbSize; c++){
	int cnb = cellCells(C0, c);
	nP = cellParticles.getCount(cnb);
	for (int n=0; n<nP; n++){
	  int pID = cellParticles(cnb, n);
	  if(particleType[pID] == 1){
	    count++;
	  }
	} 
      }
    }
#endif

    if(count>=maxcount) 	maxcount=count;   
    if(count<=mincount)	        mincount=count;
    
    (*ibFaceParticles).addCount(p, count);      
     
  }

  cout<<"max count of particles in IB Cells is "<<maxcount<<endl;
  cout<<"min count of particles in IB Cells is "<<mincount<<endl;
  //finishCount: allocate col array and reset row array
  //ready to get the entries for nonzeros
  (*ibFaceParticles).finishCount();

  //add in the entries for each row
  for(int p=0; p<rowSize; p++){
    const int faceIndex = ibFaceList [p];
    int C0 = faceCells(faceIndex, 0);
    int C1 = faceCells(faceIndex, 1);

    if (mesh.getIBTypeForCell(C1) == Mesh::IBTYPE_BOUNDARY)     
      {  C0 = C1; }
    //C0 is IBtype, C1 is fluid
    int count=0;
    int nP = cellParticles.getCount(C0);
    for(int n=0; n<nP; n++){
      int pID=cellParticles(C0,n);
      if(particleType[pID] == 1){
	count++;
	(*ibFaceParticles).add(p, pID);
      }
    }
#if 1  
    if(count < 3){
      const int nbSize = cellCells.getCount(C0);
      for(int c=0; c < nbSize; c++){
	int cnb = cellCells(C0, c);
	nP = cellParticles.getCount(cnb);
	for (int n=0; n<nP; n++){
	  int pID = cellParticles(cnb, n);
	  if(particleType[pID] == 1){
	    (*ibFaceParticles).add(p, pID);
	  }
	}
      }
    }
#endif

  }

  (*ibFaceParticles).finishAdd();
  
  return(ibFaceParticles);
}

const shared_ptr<CRConnectivity> setibFaceCells 
                          (const Mesh& mesh,
			   const Array<int>& ibFaceList,
			   const StorageSite& ibFaces, 
			   const StorageSite& cells,
			   const CRConnectivity& faceCells,
			   const CRConnectivity& cellFaces,
			   const VecD3Array& faceCentroid)
{

  shared_ptr<CRConnectivity> ibFaceCells (new CRConnectivity (ibFaces, cells));
  int maxcount=0;
  //initCount: new Array for row
  (*ibFaceCells).initCount();
  
  const int rowSize = ibFaces.getCount();

 
  //search level = 1, search only one fluid cell adjacent to IBface
  //search level = 2, search two levels of fluid cell neighborhood of the ibface
  const int searchLevel = 2 ;

  //specify the number of nonzeros for each row
 

  for(int p=0; p<rowSize; p++){
    const int IBfaceIndex = ibFaceList [p];
    int count=0;
    //find the fluid cells next to ibface
   
    int C0 = faceCells(IBfaceIndex, 0);      
    const int cellType = mesh.getIBTypeForCell(C0);
    if (cellType != Mesh::IBTYPE_FLUID){
      C0 = faceCells(IBfaceIndex,1);
    }
    count ++;
    if(searchLevel == 2){   
	
      const int nf = cellFaces.getCount(C0);
      for(int f=0; f<nf; f++){
	const int faceID = cellFaces(C0, f);
	if (faceID != IBfaceIndex){
	  const int CC0 = faceCells(faceID,0);
	  const int CC1 = faceCells(faceID,1);
	  if(CC0 != C0 && mesh.getIBTypeForCell(CC0) == Mesh::IBTYPE_FLUID){
	    count++;
	  }
	  if(CC1 != C0 && mesh.getIBTypeForCell(CC1) == Mesh::IBTYPE_FLUID){
	    count++;
	  }
	}
      }
    }
       
    (*ibFaceCells).addCount(p, count);   
    if (count>=maxcount)
      maxcount=count;
  }
    
  cout<<"max Cell neibhbors  "<<maxcount<<endl;
  
  //finishCount: allocate col array and reset row array
  //ready to get the entries for nonzeros
  (*ibFaceCells).finishCount();

  //add in the entries for each row
  for(int p=0; p<rowSize; p++){
    const int IBfaceIndex = ibFaceList [p];
    vector<int> cellIndexList;
    int C0 = faceCells(IBfaceIndex, 0);      
    const int cellType = mesh.getIBTypeForCell(C0);
    if (cellType != Mesh::IBTYPE_FLUID){
      C0 = faceCells(IBfaceIndex,1);
    }
    (*ibFaceCells).add(p, C0);

    if(searchLevel == 2){   
      const int nf = cellFaces.getCount(C0);
      for(int f=0; f<nf; f++){
	const int faceID = cellFaces(C0, f);
	if (faceID != IBfaceIndex){
	  const int CC0 = faceCells(faceID,0);
	  const int CC1 = faceCells(faceID,1);
	  if(CC0 != C0 && mesh.getIBTypeForCell(CC0) == Mesh::IBTYPE_FLUID){
	    (*ibFaceCells).add(p, CC0);
	  }
	  if(CC1 != C0 && mesh.getIBTypeForCell(CC1) == Mesh::IBTYPE_FLUID){
	    (*ibFaceCells).add(p, CC1);;
	  }
	}
      }
    }
  }
  (*ibFaceCells).finishAdd();
  
  return (ibFaceCells);
}


const  shared_ptr<CRConnectivity> setParticleCells
                          (const StorageSite& rowSite,
			   const StorageSite& colSite, 
			   const Array<int> & connectivity)
{
  const int rowSize = rowSite.getCount();
  //  const int colSize = colSite.getCount();

  shared_ptr<CRConnectivity> rowCol (new CRConnectivity (rowSite, colSite));

  //initCount: new Array for row
  (*rowCol).initCount();

  //specify the number of nonzeros for each row
  //here, each solid point only has connectivity with one cell
  //so for each row, it has only one nonzero
  for(int p=0; p<rowSize; p++){
    int value = connectivity[p];
    if (value != -1)
      (*rowCol).addCount(p, 1);
  }

  //finishCount: allocate col array and reset row array
  //ready to get the entries for nonzeros
  (*rowCol).finishCount();

  //add in the entries for each row
  for(int p=0; p<rowSize; p++){
    int value = connectivity[p];
    if (value != -1)
      (*rowCol).add(p, value);
  }

  (*rowCol).finishAdd();
  return(rowCol);

}

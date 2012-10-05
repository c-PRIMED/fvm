// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _VTKWRITER_H_
#define _VTKWRITER_H_

#include "atype.h"

#include "Mesh.h"
#include "ArrayWriter.h"
#include "Field.h"
#include "GeomFields.h"

// taken from VTK source
#define VTK_EMPTY_CELL        0
#define VTK_VERTEX            1
#define VTK_POLY_VERTEX       2
#define VTK_LINE              3
#define VTK_POLY_LINE         4
#define VTK_TRIANGLE          5
#define VTK_TRIANGLE_STRIP    6
#define VTK_POLYGON           7
#define VTK_PIXEL             8
#define VTK_QUAD              9
#define VTK_TETRA             10
#define VTK_VOXEL             11
#define VTK_HEXAHEDRON        12
#define VTK_WEDGE             13
#define VTK_PYRAMID           14
#define VTK_PENTAGONAL_PRISM  15
#define VTK_HEXAGONAL_PRISM   16
#define VTK_CONVEX_POINT_SET  41

template<class T>
class VTKWriter
{
public:
  typedef Array<T> TArray;
  typedef Array<int> IntArray;
  typedef Vector<T,3> VectorT3;
  typedef Array<Vector<T,3> > VectorT3Array;
  
  VTKWriter(const GeomFields& geomFields,
            const MeshList& meshes,
            const string fileName,
            const string& comment,
            const bool binary,
            const int atypeComponent,
            const bool surfaceOnly=false) :
    _geomFields(geomFields),
    _meshes(meshes),
    _fp(fopen(fileName.c_str(),"wb")),
    _binary(binary),
    _atypeComponent(atypeComponent),
    _surfaceOnly(surfaceOnly)
  {
    if (!_fp)
      throw CException("VTKWriter: cannot open file " + fileName +
                       "for writing");

    Field globalIndexField("gIndex");

    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];
        const StorageSite& nodes = mesh.getNodes();
        const int numNodes = nodes.getCount();
        shared_ptr<IntArray> gnPtr(new IntArray(numNodes));
        globalIndexField.addArray(nodes,gnPtr);
        *gnPtr = -1;
    }

    _gNodeCount = 0;
    _gCellCount = 0;
    int gCellNodeCount =0;
    
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];
        const StorageSite& nodes = mesh.getNodes();
        IntArray& gNodeIndex = dynamic_cast<IntArray&>(globalIndexField[nodes]);

        if (_surfaceOnly)
        {
            foreach(const FaceGroupPtr fgPtr, mesh.getAllFaceGroups())
            {
                const FaceGroup& fg = *fgPtr;
                if (fg.groupType!="interior")
                {
                    const StorageSite& faces = fg.site;
                    const int faceCount = faces.getCount();
                    const CRConnectivity& faceNodes = mesh.getFaceNodes(faces);

                    for(int f=0; f<faceCount; f++)
                    {
                        const int nFaceNodes = faceNodes.getCount(f);
                        for(int nn=0; nn<nFaceNodes; nn++)
                        {
                            const int node = faceNodes(f,nn);
                            if (gNodeIndex[node] == -1)
                              gNodeIndex[node] = _gNodeCount++;
                        }
                        _gCellCount++;
                        gCellNodeCount += nFaceNodes+1;
                        
                    }
                }
            }
        }
        else
        {
            const StorageSite& cells = mesh.getCells();
            const int numCells = cells.getSelfCount();
            const CRConnectivity& cellNodes = mesh.getCellNodes();
            const IntArray& ibType = dynamic_cast<const IntArray&>(_geomFields.ibType[cells]);
            
            for(int c=0; c<numCells; c++)
              if (ibType[c] == Mesh::IBTYPE_FLUID)
              {
                  const int nCellNodes = cellNodes.getCount(c);
                  for(int nn=0; nn<nCellNodes; nn++)
                  {
                      const int node = cellNodes(c,nn);
                      if (gNodeIndex[node] == -1)
                        gNodeIndex[node] = _gNodeCount++;
                  }
                  _gCellCount++;
                  gCellNodeCount += nCellNodes+1;
              }
        }
    }

    Array<Vector<double,3> > gNodeCoords(_gNodeCount);

    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];
        const StorageSite& nodes = mesh.getNodes();
        const int numNodes = nodes.getCount();

        const IntArray& gNodeIndex =
          dynamic_cast<const IntArray&>(globalIndexField[nodes]);

        const VectorT3Array& coords =
            dynamic_cast<const VectorT3Array&>(geomFields.coordinate[nodes]);

        for(int node=0; node<numNodes; node++)
        {
            const int gn = gNodeIndex[node];
            if (gn != -1)
            {
                for(int k=0; k<3; k++)
                  gNodeCoords[gn][k] =   coords[node][k];
            }
        }
    }
    
    fprintf(_fp,"# vtk DataFile Version 2.0\n");
    fprintf(_fp,"%s\n",comment.c_str());
    //    if (_surfaceOnly)
    //   fprintf(_fp,"ASCII\nDATASET POLYDATA\n");
    //else
      fprintf(_fp,"ASCII\nDATASET UNSTRUCTURED_GRID\n");
    
    fprintf(_fp,"POINTS %d double\n", _gNodeCount);

    for(int node=0; node<_gNodeCount; node++)
      fprintf(_fp,"%12.5le %12.5le %12.5le\n", gNodeCoords[node][0],
              gNodeCoords[node][1],gNodeCoords[node][2]);

    //if (_surfaceOnly)
    //  fprintf(_fp,"LINES");
    //else
      fprintf(_fp,"CELLS");
    
    fprintf(_fp," %d %d\n", _gCellCount, gCellNodeCount);

    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];
        const StorageSite& nodes = mesh.getNodes();
        const IntArray& gNodeIndex =
          dynamic_cast<const IntArray&>(globalIndexField[nodes]);

        if (_surfaceOnly)
        {
            foreach(const FaceGroupPtr fgPtr, mesh.getAllFaceGroups())
            {
                const FaceGroup& fg = *fgPtr;
                if (fg.groupType!="interior")
                {
                    const StorageSite& faces = fg.site;
                    const int faceCount = faces.getCount();
                    const CRConnectivity& faceNodes = mesh.getFaceNodes(faces);

                    for(int f=0; f<faceCount; f++)
                    {
                        const int nFaceNodes = faceNodes.getCount(f);
                        fprintf(_fp, "%d ", nFaceNodes);
                        for(int nn=0; nn<nFaceNodes; nn++)
                        {
                            const int node = faceNodes(f,nn);
                            fprintf(_fp, "%d ", gNodeIndex[node]);
                        }
                        fprintf(_fp,"\n");
                    }
                }
            }
        }
        else
        {
            
            const StorageSite& cells = mesh.getCells();
            const int numCells = cells.getSelfCount();
            const IntArray& ibType =
              dynamic_cast<const IntArray&>(_geomFields.ibType[cells]);
            
            const CRConnectivity& cellNodes = mesh.getCellNodes();
            
            for(int c=0; c<numCells; c++)
              if (ibType[c] == Mesh::IBTYPE_FLUID)
              {
                  const int nCellNodes = cellNodes.getCount(c);
                  fprintf(_fp, "%d ", nCellNodes);
                  for(int nn=0; nn<nCellNodes; nn++)
                  {
                      const int node = cellNodes(c,nn);
                      fprintf(_fp, "%d ", gNodeIndex[node]);
                  }
                  fprintf(_fp,"\n");
              }
        }
    }

    fprintf(_fp,"CELL_TYPES %d\n", _gCellCount);
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];
        const int dim = mesh.getDimension();
        if (_surfaceOnly)
        {
            foreach(const FaceGroupPtr fgPtr, mesh.getAllFaceGroups())
            {
                const FaceGroup& fg = *fgPtr;
                if (fg.groupType!="interior")
                {
                    const StorageSite& faces = fg.site;
                    const int faceCount = faces.getCount();
                    const CRConnectivity& faceNodes = mesh.getFaceNodes(faces);
                    
                    for(int f=0; f<faceCount; f++)
                    {
                        const int nFaceNodes = faceNodes.getCount(f);
                        int vtkCellType = VTK_POLYGON;
                        if (dim == 2)
                          vtkCellType = VTK_LINE;
                        else if (nFaceNodes == 3)
                          vtkCellType = VTK_TRIANGLE;
                        else if (nFaceNodes == 4)
                          vtkCellType = VTK_QUAD;
                        fprintf(_fp,"%d\n",vtkCellType);
                    }
                }
            }
        }
        else
        {
            const StorageSite& cells = mesh.getCells();
            const int numCells = cells.getSelfCount();
            const IntArray& ibType =
              dynamic_cast<const IntArray&>(_geomFields.ibType[cells]);
            
            const CRConnectivity& cellNodes = mesh.getCellNodes();
            
            for(int c=0; c<numCells; c++)
              if (ibType[c] == Mesh::IBTYPE_FLUID)
              {
                  const int nCellNodes = cellNodes.getCount(c);
                  int vtkCellType = VTK_CONVEX_POINT_SET;
                  if (dim == 2)
                  {
                      if (nCellNodes == 4)
                        vtkCellType = VTK_QUAD;
                      else if (nCellNodes == 3)
                        vtkCellType = VTK_TRIANGLE;
                      else
                        vtkCellType = VTK_POLYGON;
                  }
                  else
                  {
                      if (nCellNodes == 4)
                        vtkCellType = VTK_TETRA;
                      else if (nCellNodes == 8)
                        vtkCellType = VTK_HEXAHEDRON;
                      else if (nCellNodes == 5)
                        vtkCellType = VTK_PYRAMID;
                      else if (nCellNodes == 6)
                        vtkCellType = VTK_WEDGE;
                  }
                  fprintf(_fp,"%d\n",vtkCellType);
              }
        }
    }
    
    fprintf(_fp,"CELL_DATA %d\n", _gCellCount);

  }

  void init()
  {}
  

  void writeScalarField(const Field& field,
                        const string label)
  {
    const int numMeshes = _meshes.size();
    fprintf(_fp,"SCALARS %s double\n", label.c_str());
    fprintf(_fp,"LOOKUP_TABLE default\n");
    
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];
        if (_surfaceOnly)
        {
            const StorageSite& allFaces = mesh.getFaces();
            const TArray& a =
              dynamic_cast<const TArray&>(field[allFaces]);

            foreach(const FaceGroupPtr fgPtr, mesh.getAllFaceGroups())
            {
                const FaceGroup& fg = *fgPtr;
                if (fg.groupType!="interior")
                {
                    const StorageSite& faces = fg.site;
                    const int numFaces = faces.getSelfCount();
                    const int offset = faces.getOffset();
                    for(int f=0; f<numFaces; f++)
                    {
                        const int gf = f + offset;
                        fprintf(_fp,"%12.5le\n", (a[gf]));
                    }
                }
            }
        }
        else
        {
            const StorageSite& cells = mesh.getCells();
            const int numCells = cells.getSelfCount();
            const IntArray& ibType =
              dynamic_cast<const IntArray&>(_geomFields.ibType[cells]);
            
            const TArray& aCell = dynamic_cast<const TArray&>(field[cells]);
            for(int c=0; c<numCells; c++)
              if (ibType[c] == Mesh::IBTYPE_FLUID)
              {
                  fprintf(_fp,"%12.5le\n", (aCell[c]));
              }
        }
    }
  }
  
  void writeVectorField(const Field& field,const string label)
  {


    fprintf(_fp,"VECTORS %s double\n", label.c_str());
    const int numMeshes = _meshes.size();
    
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        if (_surfaceOnly)
        {
            const StorageSite& allFaces = mesh.getFaces();
            const VectorT3Array& v =
              dynamic_cast<const VectorT3Array&>(field[allFaces]);

            foreach(const FaceGroupPtr fgPtr, mesh.getAllFaceGroups())
            {
                const FaceGroup& fg = *fgPtr;
                if (fg.groupType!="interior")
                {
                    const StorageSite& faces = fg.site;
                    const int numFaces = faces.getSelfCount();
                    const int offset = faces.getOffset();
                    for(int f=0; f<numFaces; f++)
                    {
                        const int gf = f + offset;
                        fprintf(_fp,"%12.5le %12.5le %12.5le\n", (v[gf][0]),
                                (v[gf][1]),
                                (v[gf][2])
                                );
                    }
                }
            }
        }
        else
        {
            const StorageSite& cells = mesh.getCells();
            const int numCells = cells.getSelfCount();
            const IntArray& ibType =
              dynamic_cast<const IntArray&>(_geomFields.ibType[cells]);
            
            const VectorT3Array& aCell =
              dynamic_cast<const VectorT3Array&>(field[cells]);
            for(int c=0; c<numCells; c++)
              if (ibType[c] == Mesh::IBTYPE_FLUID)
              {
                  fprintf(_fp,"%12.5le %12.5le %12.5le\n", (aCell[c][0]),
                          (aCell[c][1]),
                          (aCell[c][2])
                          );
              }
        }
    }
  }

  void finish()
  {
    fclose(_fp);
  }

private:
  const GeomFields& _geomFields;
  const MeshList _meshes;
  FILE *_fp;
  const bool _binary;
  const int _atypeComponent;
  int _gNodeCount;
  int _gCellCount;
  const bool _surfaceOnly;
};
#endif

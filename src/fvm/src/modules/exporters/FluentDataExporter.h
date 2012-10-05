// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _FLUENTDATAEXPORTER_H_
#define _FLUENTDATAEXPORTER_H_

#include "atype.h"

#include "FluentReader.h"
#include "ArrayWriter.h"

template<class T>
class FluentDataExporter
{
public:
  typedef Array<T> TArray;
  
  FluentDataExporter(FluentReader& reader,
                     const string fileName,
                     const bool binary,
                     const int atypeComponent) :
    _reader(reader),
    _fp(fopen(fileName.c_str(),"wb")),
    _binary(binary),
    _atypeComponent(atypeComponent),
    _sectionID(_binary ? 3300 : 300)
  {
    if (!_fp)
      throw CException("FluentDataExporter: cannot open file " + fileName +
                       "for writing");
  }

  void init()
  {
    fprintf(_fp,"(4 (60 0 0 1 2 4 4 4 8 8 4))\n");
  }
  
  void closeSection()
  {    
    if (_binary)
      fprintf(_fp,"\nEnd of Binary Section   3300"); 
    fprintf(_fp,"))\n");
  }

  void writeScalarField(const Field& field,
                        const int fluentFieldId)
  {
    FaceZonesMap& faceZones = _reader.getFaceZones();
    CellZonesMap& cellZones = _reader.getCellZones();

    ScalarArrayWriter<T> writer(_binary,0,_atypeComponent);
    
    foreach(const CellZonesMap::value_type& pos, cellZones)
    {
        const FluentCellZone& cz = *pos.second;
        const Mesh& mesh = *cz.mesh;

        const StorageSite& cells = mesh.getCells();
        const StorageSite& faces = mesh.getFaces();

        if (field.hasArray(cells))
        {
            const TArray& aCell = dynamic_cast<const TArray&>(field[cells]);
            fprintf(_fp,"(%d (%d %d 1 0 1 %d %d)\n(",
                    _sectionID, fluentFieldId, cz.ID,
                    cz.iBeg+1,cz.iEnd+1);
            writer.write(_fp,aCell,0,cells.getSelfCount());
            closeSection();
            
            foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
            {
                const FaceGroup& fg = *fgPtr;
                const StorageSite& faces = fg.site;

                const FluentFaceZone& fz = *faceZones[fg.id];
                const CRConnectivity& faceCells = mesh.getFaceCells(faces);
                
                fprintf(_fp,"(%d (%d %d 1 0 1 %d %d)\n(",
                        _sectionID, fluentFieldId, fz.ID,
                        fz.iBeg+1,fz.iEnd+1);

                int cbeg = faceCells(0,1);
                
                writer.write(_fp,aCell,cbeg,faces.getSelfCount());
                closeSection();
            }
        }
        if (field.hasArray(faces))
        {
            const TArray& aFaces = dynamic_cast<const TArray&>(field[faces]);

            int faceOffset =0;
            foreach(int fzId, cz.interiorZoneIds)
            {
                const FluentFaceZone& fz = *faceZones[fzId];
                const int count = fz.iEnd-fz.iBeg+1;
                fprintf(_fp,"(%d (%d %d 1 0 1 %d %d)\n(",
                        _sectionID, fluentFieldId, fz.ID,
                        fz.iBeg+1,fz.iEnd+1);
                writer.write(_fp,aFaces,faceOffset,count);
                closeSection();
                faceOffset += count;
            }

            foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
            {
                const FaceGroup& fg = *fgPtr;
                const StorageSite& faces = fg.site;

                const FluentFaceZone& fz = *faceZones[fg.id];
                
                fprintf(_fp,"(%d (%d %d 1 0 1 %d %d)\n(",
                        _sectionID, fluentFieldId, fz.ID,
                        fz.iBeg+1,fz.iEnd+1);

                writer.write(_fp,aFaces,faces.getOffset(),faces.getSelfCount());
                closeSection();
            }
        }
                                                  
    }
  }
  
  void writeVectorField(const Field& field,const int fluentFieldId)
  {

    typedef Array<Vector<T,3> > VectorT3Array;

    FaceZonesMap& faceZones = _reader.getFaceZones();
    CellZonesMap& cellZones = _reader.getCellZones();

    for(int nd=0; nd<3; nd++)
    {
        VectorArrayWriter<T,3> writer(_binary,nd,_atypeComponent);
    
        foreach(const CellZonesMap::value_type& pos, cellZones)
        {
            const FluentCellZone& cz = *pos.second;
            const Mesh& mesh = *cz.mesh;

            const StorageSite& cells = mesh.getCells();

            if (field.hasArray(cells))
            {
                const VectorT3Array& aCell = dynamic_cast<const VectorT3Array&>(field[cells]);

                fprintf(_fp,"(%d (%d %d 1 0 1 %d %d)\n(",
                        _sectionID, fluentFieldId+nd, cz.ID,
                        cz.iBeg+1,cz.iEnd+1);
                writer.write(_fp,aCell,0,cells.getSelfCount());
                closeSection();
            
                foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
                {
                    const FaceGroup& fg = *fgPtr;
                    const StorageSite& faces = fg.site;

                    const FluentFaceZone& fz = *faceZones[fg.id];
                    const CRConnectivity& faceCells = mesh.getFaceCells(faces);
                
                    fprintf(_fp,"(%d (%d %d 1 0 1 %d %d)\n(",
                            _sectionID, fluentFieldId+nd, fz.ID,
                            fz.iBeg+1,fz.iEnd+1);

                    int cbeg = faceCells(0,1);
                
                    writer.write(_fp,aCell,cbeg,faces.getSelfCount());
                    closeSection();
                }
            }
        }
    }
    
  }

  void finish()
  {
    fclose(_fp);
  }

private:
  FluentReader& _reader;
  FILE *_fp;
  const bool _binary;
  const int _atypeComponent;
  const int _sectionID;
};
#endif

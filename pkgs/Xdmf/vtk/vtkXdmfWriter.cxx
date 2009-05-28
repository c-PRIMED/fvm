/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkXdmfWriter.cxx,v $
  Language:  C++
  Date:      $Date: 2009-03-13 20:01:09 $
  Version:   $Revision: 1.7 $


Copyright (c) 1993-2001 Ken Martin, Will Schroeder, Bill Lorensen  
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

 * Neither name of Ken Martin, Will Schroeder, or Bill Lorensen nor the names
   of any contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

 * Modified source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/
#include "vtkXdmfWriter.h"

#include "vtkGenericCell.h"
#include "vtkCell.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCellTypes.h"
#include "vtkCommand.h"
#include "vtkImageData.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"

#include "vtkRectilinearGrid.h"
#include "vtkStructuredGrid.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUniformGrid.h"

#include "vtkDataSetCollection.h"

#include "vtkSmartPointer.h"

#include "vtkCharArray.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkLongArray.h"
#include "vtkShortArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkUnsignedShortArray.h"
#include "vtkUnsignedIntArray.h"

#include "XdmfHDF.h"
#include "XdmfArray.h"
#include "XdmfAttribute.h"
#include "XdmfDOM.h"

#include "H5Spublic.h"

#include <vtkstd/map>
#include <vtksys/SystemTools.hxx>

#define vtkMAX(x, y) (((x)>(y))?(x):(y))

//----------------------------------------------------------------------------
  #define OUTPUTTEXT(a) vtkOutputWindowDisplayText(a);
  #define vtkGenericDebugMacro(a)  \
  { \
    vtkOStreamWrapper::EndlType endl; \
    vtkOStreamWrapper::UseEndl(endl); \
    vtkOStrStreamWrapper vtkmsg; \
    vtkmsg << " vtkXdmfWriter " a << "\n"; \
    OUTPUTTEXT(vtkmsg.str()); \
    vtkmsg.rdbuf()->freeze(0); \
  }
//---------------------------------------------------------------------------

struct vtkXdmfWriterInternal
{
  class CellType
    {
  public:
    CellType() : VTKType(0), NumPoints(0) {}
    CellType(const CellType& ct) : VTKType(ct.VTKType), NumPoints(ct.NumPoints) {}
    vtkIdType VTKType;
    vtkIdType NumPoints;
    vtkstd_bool operator<(const CellType& ct) const
      {
      return this->VTKType < ct.VTKType || this->VTKType == ct.VTKType && this->NumPoints < ct.NumPoints;
      }
    vtkstd_bool operator==(const CellType& ct) const
      {
      return this->VTKType == ct.VTKType && this->NumPoints == ct.NumPoints;
      }
    CellType& operator=(const CellType& ct)
      {
      this->VTKType = ct.VTKType;
      this->NumPoints = ct.NumPoints;
      return *this;
      }

    };
  typedef vtkstd::map<CellType, vtkSmartPointer<vtkIdList> > MapOfCellTypes;
  static void DetermineCellTypes(vtkPointSet *t, MapOfCellTypes& vec);
};


//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkXdmfWriter);
vtkCxxRevisionMacro(vtkXdmfWriter, "$Revision: 1.7 $");

//----------------------------------------------------------------------------
vtkXdmfWriter::vtkXdmfWriter()
{
  this->FileNameString = 0;
  this->HeavyDataSetNameString = 0;
  this->GridName = 0;
  this->DomainName = 0;
  this->CollectionName = 0;
  this->SetHeavyDataSetName( "XdmfData.h5" );
  this->SetGridName( "Unnamed" );
  this->CollectionType = XDMF_GRID_COLLECTION_UNSET;

  this->AllLight = 0;
  this->AllHeavy = 0;
  this->CurrIndent = 0;

  this->InputList = 0;

  this->HDF5ArrayName = 0;

  this->GridOnly = 0;
  this->InputsArePieces = 0;
  this->FullGridSize[0] = 0;
  this->FullGridSize[1] = 0;
  this->FullGridSize[2] = 0;
  this->TimeValue = -1E99;
  this->DocString = NULL;
  this->AppendGridsToDomain = 0;
  this->DOM = NULL;
}

//----------------------------------------------------------------------------
vtkXdmfWriter::~vtkXdmfWriter()
{
  this->SetHeavyDataSetName(0);
  this->SetFileNameString(0);
  if (this->InputList != NULL)
    {
    this->InputList->Delete();
    this->InputList = NULL;
    }
  this->SetHDF5ArrayName(0);
  this->SetDomainName(0);
  this->SetGridName(0);
  if (this->DocString)
  {
    delete []this->DocString;
  }
}

//----------------------------------------------------------------------------
void vtkXdmfWriter::SetFileName(const char* fname)
{
  if ( fname )
    {
    char* hname = new char [ strlen(fname) + 10 ]; // space for extension
    strcpy(hname, fname);
    size_t cc;
    for ( cc = strlen(hname); cc > 0; cc -- )
      {
      if ( hname[cc] == '.' )
        {
        break;
        }
      }
    if ( hname[cc] == '.' )
      {
      cc--;
      }
    if ( cc > 0 )
      {
      hname[cc+1] = 0;
      }
    strcat(hname, ".h5");
    this->SetHeavyDataSetName(hname);
    vtkDebugMacro(<< "Set Heavy Data Set Name: " << hname);
    }
  this->SetFileNameString(fname);
}

//----------------------------------------------------------------------------
const char* vtkXdmfWriter::GetFileName()
{
  return this->FileNameString;
}

//----------------------------------------------------------------------------
void vtkXdmfWriter::SetHeavyDataSetName( const char *name) 
{
  this->SetHeavyDataSetNameString(name);
  if ( name )
    {
    this->AllLight = 0;
    }
  else
    {
    this->AllLight = 1;
    this->AllHeavy = 0;
    }
  this->Modified();
}

//----------------------------------------------------------------------------
const char* vtkXdmfWriter::GetHeavyDataSetName()
{
  return this->HeavyDataSetNameString;
}

//----------------------------------------------------------------------------
int vtkXdmfWriter::ParseExistingFile(const char* filename)
{
  XdmfConstString CurrentFileName;
  vtkstd::string directory;

  if (!filename)
    {
    vtkErrorMacro("File name not set");
    return 1;
    }

  // First make sure the file exists.  This prevents an empty file
  // from being created on older compilers.
  if (!vtksys::SystemTools::FileExists(filename)) {
    vtkErrorMacro("Error opening file " << filename);
    return 1;
    }

  if ( !this->DOM )
    {
    this->DOM = new XdmfDOM();
    }

  //Tell the parser what the working directory is.
  directory = vtksys::SystemTools::GetFilenamePath(filename) + "/";
  if (directory == "")
    {
    directory = vtksys::SystemTools::GetCurrentWorkingDirectory() + "/";
    }
  directory = vtksys::SystemTools::ConvertToOutputPath(directory.c_str());
  this->DOM->SetWorkingDirectory(directory.c_str());

/*
  if ( !this->Internals->DataItem )
    {
    this->Internals->DataItem = new XdmfDataItem();
    this->Internals->DataItem->SetDOM(this->DOM);
    }
*/
  CurrentFileName = this->DOM->GetInputFileName();  
  // Only re-parse when necessary. Reparsing can be slow and can reset
  // parameters and stride that we read from file but then allow user to 
  // modify.
  this->DOM->SetInputFileName(filename);
  vtkDebugMacro(".!!............Preparing to Parse " << filename);
  this->DOM->Parse(filename);

  return 1;
}

//----------------------------------------------------------------------------
int vtkXdmfWriter::WriteHead( ostream& ost )
{
  ost << "<?xml version=\"1.0\" ?>" << "\n"
    << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" [" << "\n"
    << "<!ENTITY HeavyData \"" << this->HeavyDataSetNameString << "\">" <<"\n"
    << "]>" << "\n" << "\n" << "\n";
  this->Indent(ost);
  ost << "<Xdmf>";
  this->IncrementIndent();
  this->Indent(ost);

  return 1;
}

//----------------------------------------------------------------------------
int vtkXdmfWriter::WriteTail( ostream& ost )
{
  this->DecrementIndent();
  this->Indent(ost);
  ost << "</Xdmf>";
  this->Indent(ost);
  return 0;
}

//----------------------------------------------------------------------------
int vtkXdmfWriter::WriteCellArray( ostream& ost, vtkDataSet *ds, 
  const char* gridName, void* mapofcell, const void* celltype )
{
  vtkIdType PointsInPoly;
  vtkIdType i, j;

  vtkXdmfWriterInternal::CellType* ct = (vtkXdmfWriterInternal::CellType*)celltype;
  vtkXdmfWriterInternal::MapOfCellTypes* mc = (vtkXdmfWriterInternal::MapOfCellTypes*)mapofcell;
  PointsInPoly = ct->NumPoints;
  vtkIdList* il = (*mc)[*ct].GetPointer();

  ost << "<DataItem";
  this->IncrementIndent();
  this->Indent(ost);
  ost << " NumberType=\"Int\"";
  this->Indent(ost);
  ost << " Dimensions=\"" << il->GetNumberOfIds() << " " << PointsInPoly << "\"";
  this->Indent(ost);
  vtkIdList* cellPoints = vtkIdList::New();
  if( this->AllLight )
    {
    ost << " Format=\"XML\">";
    for( i = 0 ; i < il->GetNumberOfIds(); i++ )
      {
      this->Indent(ost);
      ds->GetCellPoints(il->GetId(i), cellPoints);
      if ( ct->VTKType == VTK_VOXEL )
        {
        // Hack for VTK_VOXEL
        ost << " " << cellPoints->GetId(0);
        ost << " " << cellPoints->GetId(1);
        ost << " " << cellPoints->GetId(3);
        ost << " " << cellPoints->GetId(2);
        ost << " " << cellPoints->GetId(4);
        ost << " " << cellPoints->GetId(5);
        ost << " " << cellPoints->GetId(7);
        ost << " " << cellPoints->GetId(6);
        }
      else if ( ct->VTKType == VTK_PIXEL )
        {
        // Hack for VTK_PIXEL
        ost << " " << cellPoints->GetId(0);
        ost << " " << cellPoints->GetId(1);
        ost << " " << cellPoints->GetId(3);
        ost << " " << cellPoints->GetId(2);
        }
      else
        {
        for( j = 0 ; j < PointsInPoly ; j++ )
          {
          ost << " " << cellPoints->GetId(j);
          }
        }
      }
    } 
  else 
    {
    // Create HDF File
    XdmfArray  Conns;
    XdmfHDF    H5;
    XdmfInt64  Dims[2];
    XdmfInt32  *Dp;

    const char* DataSetName = this->GenerateHDF5ArrayName(gridName, "Connections");
    ost << " Format=\"HDF\">";
    this->IncrementIndent();
    this->Indent(ost);
    ost << " " << DataSetName;
    this->DecrementIndent();
    Conns.SetNumberType( XDMF_INT32_TYPE );
    Dims[0] = il->GetNumberOfIds();
    Dims[1] = PointsInPoly;
    Conns.SetShape( 2, Dims );
    Dp = (XdmfInt32 *)Conns.GetDataPointer();
    for( i = 0 ; i < il->GetNumberOfIds(); i++ )
      {
      vtkIdType cid = il->GetId(i);
      ds->GetCellPoints(cid, cellPoints);
      if ( ct->VTKType == VTK_VOXEL )
        {
        // Hack for VTK_VOXEL
        *Dp++ = cellPoints->GetId(0);
        *Dp++ = cellPoints->GetId(1);
        *Dp++ = cellPoints->GetId(3);
        *Dp++ = cellPoints->GetId(2);
        *Dp++ = cellPoints->GetId(4);
        *Dp++ = cellPoints->GetId(5);
        *Dp++ = cellPoints->GetId(7);
        *Dp++ = cellPoints->GetId(6);
        }
      else if ( ct->VTKType == VTK_PIXEL )
        {
        // Hack for VTK_PIXEL
        *Dp++ = cellPoints->GetId(0);
        *Dp++ = cellPoints->GetId(1);
        *Dp++ = cellPoints->GetId(3);
        *Dp++ = cellPoints->GetId(2);
        }
      else
        {
        for( j = 0 ; j < PointsInPoly ; j++ )
          {
          *Dp++ = cellPoints->GetId(j);
          }
        }
      }
    H5.CopyType( &Conns );
    H5.CopyShape( &Conns );
    if( H5.Open( DataSetName, "rw" ) == XDMF_FAIL )
      {
      if( H5.CreateDataset( DataSetName ) == XDMF_FAIL ) 
        {
        vtkErrorMacro( "Can't Create Heavy Dataset " << DataSetName );
        return( -1 );
        }
      }
    H5.Write( &Conns );
    H5.Close();


    }
  cellPoints->Delete();
  this->DecrementIndent();
  this->Indent(ost);
  ost << "</DataItem>";
  return( il->GetNumberOfIds() );
}

//----------------------------------------------------------------------------
int vtkXdmfWriter::WritePoints( ostream& ost, vtkPoints *Points, vtkDataSet* dataSet,
  const char* gridName )
{
  int dims[3] = { -10, -10, -10 };
  return this->WriteVTKArray( ost, Points->GetData(), dataSet, 0, dims, "XYZ", 0, gridName, 
    this->AllLight );
}

//----------------------------------------------------------------------------
void vtkXdmfWriter::EndTopology( ostream& ost )
{
  this->DecrementIndent();
  this->Indent(ost);
  ost << "</Topology>";
}

//----------------------------------------------------------------------------
void vtkXdmfWriter::StartTopology( ostream& ost, const char* toptype, int rank, int *dims)
{
  ost << "<Topology ";
  this->IncrementIndent();
  this->Indent(ost);
  ost << " TopologyType=\"" << toptype << "\"";
  ost << " Dimensions=\"";
  int cc;
  for ( cc = rank-1; cc >= 0; cc -- )
    {
    if ( cc < rank - 1 )
      {
      ost << " ";
      }
    ost << dims[cc]; 
    }
  ost << "\">";
  //this->Indent(ost);
}
//----------------------------------------------------------------------------
void vtkXdmfWriter::StartTopology( ostream& ost, int cellType, vtkIdType numVert, vtkIdType numCells )
{
  ost << "<Topology ";
  this->IncrementIndent();
  switch( cellType ) 
    {
  case VTK_EMPTY_CELL :
    vtkDebugMacro(<< "Start Empty Cell");
  case VTK_VERTEX :
    vtkDebugMacro(<< "Start " <<  " VERTEX");
    ost << " TopologyType=\"POLYVERTEX\"";
    this->Indent(ost);
    break;
  case VTK_POLY_VERTEX :
    vtkDebugMacro(<< "Start " <<  " POLY_VERTEX");
    ost << " TopologyType=\"POLYVERTEX\"";
    this->Indent(ost);
    break;
  case VTK_LINE :
    vtkDebugMacro(<< "Start " <<  " LINE");
    ost << " TopologyType=\"POLYLINE\"";
    this->Indent(ost);
    ost << " NodesPerElement=\"" << numVert << "\"";
    this->Indent(ost);
    break;
  case VTK_POLY_LINE :
    vtkDebugMacro(<< "Start " <<  " POLY_LINE");
    ost << " TopologyType=\"POLYLINE\"";
    this->Indent(ost);
    ost << " NodesPerElement=\"" << numVert << "\"";
    this->Indent(ost);
    break;
  case VTK_TRIANGLE :
    vtkDebugMacro(<< "Start " <<  " TRIANGLE");
    ost << " TopologyType=\"TRIANGLE\"";
    this->Indent(ost);
    break;
  case VTK_TRIANGLE_STRIP :
    vtkDebugMacro(<< "Start " <<  " TRIANGLE_STRIP");
    ost << " TopologyType=\"TRIANGLE\"";
    this->Indent(ost);
    break;
  case VTK_POLYGON :
    vtkDebugMacro(<< "Start " <<  " POLYGON");
    ost << " TopologyType=\"POLYGON\"";
    this->Indent(ost);
    ost << " NodesPerElement=\"" << numVert << "\"";
    this->Indent(ost);
    break;
  case VTK_PIXEL :
    vtkDebugMacro(<< "Start " <<  " PIXEL");
    ost << " TopologyType=\"QUADRILATERAL\"";
    this->Indent(ost);
    break;
  case VTK_QUAD :
    vtkDebugMacro(<< "Start " <<  " QUAD");
    ost << " TopologyType=\"QUADRILATERAL\"";
    this->Indent(ost);
    break;
  case VTK_TETRA :
    vtkDebugMacro(<< "Start " <<  " TETRA");
    ost << " TopologyType=\"TETRAHEDRON\"";
    this->Indent(ost);
    break;
  case VTK_VOXEL :
    vtkDebugMacro(<< "Start " <<  " VOXEL");
    ost << " TopologyType=\"HEXAHEDRON\"";
    this->Indent(ost);
    break;
  case VTK_HEXAHEDRON :
    vtkDebugMacro(<< "Start " <<  " HEXAHEDRON");
    ost << " TopologyType=\"HEXAHEDRON\"";
    this->Indent(ost);
    break;
  case VTK_WEDGE :
    vtkDebugMacro(<< "Start " <<  " WEDGE");
    ost << " TopologyType=\"WEDGE\"";
    this->Indent(ost);
    break;
  case VTK_PYRAMID :
    vtkDebugMacro(<< "Start " <<  " PYRAMID");
    ost << " TopologyType=\"PYRAMID\"";
    this->Indent(ost);
    break;
  default :
    vtkErrorMacro("Unknown Topology Type");
    break;
    }
  ost << " Dimensions=\"" << numCells << "\">";
  this->Indent(ost);
}

//----------------------------------------------------------------------------
template<class AType, class NType>
vtkIdType vtkXdmfWriterWriteXMLScalar(vtkXdmfWriter* self, ostream& ost, 
  AType* array, vtkDataSet* dataSet, int* scaledExtent,
  const char* dataName, const char* arrayName, const char* gridName,
  const char* scalar_type, NType value,
  int allLight, int type,
  int dims[3],
  int cellData)
{
  if ( !array )
    {
    vtkErrorWithObjectMacro(self, "No array specified. Should be " 
      << scalar_type << " array");
    return -2;
    }
 
  int updateExtent[6];
  int extent[6];
  int useExtents = 0;
  int scaledDims[3] = { -1, -1, -1 };

  // vtkGenericDebugMacro(<< "Dataset: " << arrayName << "\n");
  int cc;
  if ( scaledExtent )
    {
    for ( cc = 0; cc < 3; cc ++ )
      {
      updateExtent[cc*2] = scaledExtent[cc*2];
      updateExtent[cc*2+1] = scaledExtent[cc*2+1];
      extent[cc*2] = scaledExtent[cc*2];
      extent[cc*2+1] = scaledExtent[cc*2+1];
      scaledDims[cc] = scaledExtent[cc*2+1] - scaledExtent[cc*2] + 1;
      }
    if ( array->GetNumberOfComponents() == 1 )
      {
      extent[0] = 0;
      extent[1] = array->GetNumberOfTuples()+1;
      }
    useExtents = 1;
    }
  else
    {
    switch ( dataSet->GetDataObjectType() )
      {
    case VTK_STRUCTURED_GRID:
        {
        vtkStructuredGrid* grid = vtkStructuredGrid::SafeDownCast(dataSet);
        grid->GetUpdateExtent(updateExtent);
        grid->GetExtent(extent);
        if ( cellData )
          {
          for ( cc = 0; cc < 3; cc ++ )
            {
            updateExtent[cc*2+1] -= 1;
            extent[cc*2+1] -= 1;
            updateExtent[cc*2+1] = vtkMAX(0, updateExtent[cc*2+1]);
            extent[cc*2+1] = vtkMAX(0, extent[cc*2+1]);
            }
          }
        useExtents = 1;
        }
      break;
    case VTK_RECTILINEAR_GRID:
        {
        vtkRectilinearGrid* grid = vtkRectilinearGrid::SafeDownCast(dataSet);
        grid->GetUpdateExtent(updateExtent);
        grid->GetExtent(extent);
        if ( cellData )
          {
          // vtkGenericDebugMacro(<< "Fix for cell arrays" << "\n");
          // vtkGenericDebugMacro(<< "Extent: " << extent[0] << " " << extent[1] << " " << extent[2] << " " << extent[3] << " " << extent[4] << " " << extent[5] << "\n");
          // vtkGenericDebugMacro(<< "updateExtent: " << updateExtent[0] << " " << updateExtent[1] << " " << updateExtent[2] << " " << updateExtent[3] << " " << updateExtent[4] << " " << updateExtent[5] << "\n");
          for ( cc = 0; cc < 3; cc ++ )
            {
            updateExtent[cc*2+1] -= 1;
            extent[cc*2+1] -= 1;
            updateExtent[cc*2+1] = vtkMAX(0, updateExtent[cc*2+1]);
            extent[cc*2+1] = vtkMAX(0, extent[cc*2+1]);
            }
          }
        useExtents = 1;
        }
      break;
    case VTK_UNIFORM_GRID:
        {
        vtkUniformGrid* grid = vtkUniformGrid::SafeDownCast(dataSet);
        grid->GetUpdateExtent(updateExtent);
        grid->GetExtent(extent);
        useExtents = 1;
        for ( cc = 0; cc < 6; cc ++ )
          {
          updateExtent[cc*2+1] -= cellData;
          extent[cc*2+1] -= cellData;
            updateExtent[cc*2+1] = vtkMAX(0, updateExtent[cc*2+1]);
            extent[cc*2+1] = vtkMAX(0, extent[cc*2+1]);
          }
        }
      break;
    case VTK_IMAGE_DATA:
    case VTK_STRUCTURED_POINTS:
        {
        vtkImageData* grid = vtkImageData::SafeDownCast(dataSet);
        grid->GetUpdateExtent(updateExtent);
        grid->GetExtent(extent);
        for ( cc = 0; cc < 3; cc ++ )
          {
          updateExtent[cc*2+1] -= cellData;
          extent[cc*2+1] -= cellData;
            updateExtent[cc*2+1] = vtkMAX(0, updateExtent[cc*2+1]);
            extent[cc*2+1] = vtkMAX(0, extent[cc*2+1]);
          }
        useExtents = 1;
        }
      break;
      }
    }

  if ( useExtents )
    {
    for ( cc = 0; cc < 3; cc ++ )
      {
      // vtkGenericDebugMacro(<< "Dims[" << cc << "]: " << dims[cc] << "\n");
      // cout << "Dims[" << cc << "]: " << dims[cc] << endl;
      }
    for ( cc = 0; cc < 3; cc ++ )
      {
      // vtkGenericDebugMacro(<< "UExt[" << cc << "]: " << updateExtent[cc*2] << " - " << updateExtent[cc*2+1] << "\n");
      // cout << "UExt[" << cc << "]: " << updateExtent[cc*2] << " - " << updateExtent[cc*2+1] << endl;
      }
    for ( cc = 0; cc < 3; cc ++ )
      {
      // vtkGenericDebugMacro(<< "RExt[" << cc << "]: " << extent[cc*2] << " - " << extent[cc*2+1] << "\n");
      // cout << "RExt[" << cc << "]: " << extent[cc*2] << " - " << extent[cc*2+1] << endl;
      }
    }

  if (!self->GetInputsArePieces() || (self->GetInputsArePieces() && self->CurrentInputNumber==0)) 
    {
    ost << "<DataItem";
    self->IncrementIndent();
    if ( dataName )
      {
      self->Indent(ost);
      ost << " Name=\"" << dataName << "\"";
      }
    self->Indent(ost);
    ost << " NumberType=\"" << scalar_type << "\"";
    self->Indent(ost);
    int precision = 1;
    switch ( type )
      {
      case XDMF_FLOAT64_TYPE: case XDMF_INT64_TYPE:
        precision = 8;
        break;
      case XDMF_FLOAT32_TYPE: case XDMF_INT32_TYPE: case XDMF_UINT32_TYPE:
        precision = 4;
        break;
      case XDMF_INT16_TYPE: case XDMF_UINT16_TYPE:
        precision = 2;
        break;
      }
    if ( precision > 1 )
      {
      ost << " Precision=\"" << precision << "\"";
      self->Indent(ost);
      }
    else if ( type == XDMF_FLOAT32_TYPE )
      {
      ost << " Precision=\"4\"";
      self->Indent(ost);
      }
    ost << " Dimensions=\"";
  
    // WritePoints sets dims to -10 as a signal 
    // But a single plane can set dims to -1
    if(dims[0] > -2){
        // cout << "Fixing Dims" << endl;
        dims[0] = vtkMAX(1, dims[0]);
        dims[1] = vtkMAX(1, dims[1]);
        dims[2] = vtkMAX(1, dims[2]);
        // cout << "New Dims = " << dims[0] << ", " << dims[1] << ", " << dims[2] << endl;
    }
    if ( dims[0] < 1 )
      {
      if ( scaledDims[0] < 1 )
        {
        ost << array->GetNumberOfTuples();
        }
      else
        {
        ost << scaledDims[0];
        }
      }
    else
      {
      if (self->GetInputsArePieces()) 
        {
        // vtkGenericDebugMacro(<< "Using Custom Disk Shape " << self->GetFullGridSize()[0] << " " << self->GetFullGridSize()[1] << " " << self->GetFullGridSize()[2] << "\n");
        ost << self->GetFullGridSize()[0] << " " << self->GetFullGridSize()[1] << " " << self->GetFullGridSize()[2];
        }
      else
        {
        ost << dims[2] << " " << dims[1] << " " << dims[0];
        }
      }
    if ( array->GetNumberOfComponents() > 1 )
      {
      ost << " " << array->GetNumberOfComponents();
      }
    ost << "\"";
    self->Indent(ost);
    }

  // if all light, piece writing is knackered
  NType val = value;
  if ( allLight )
    {
    ost << " Format=\"XML\">";
    vtkIdType jj, kk;
    vtkIdType xx, yy, zz;
    if ( useExtents )
      {
      int printOne = 0;
      // vtkGenericDebugMacro(<< "Use Extent" << "\n");
      jj = 0;
      vtkIdType idx = 0;
      vtkIdType size =
        (updateExtent[1]-updateExtent[0] + 1) *
        (updateExtent[3]-updateExtent[2] + 1) *
        (updateExtent[5]-updateExtent[4] + 1);
      // vtkGenericDebugMacro(<< "Size: " << size << " " << array->GetNumberOfComponents() << "\n");
      if ( size != array->GetNumberOfTuples() )
        {
        // Error
        vtkErrorWithObjectMacro(self, "Wrong number of tuples in the dataset. Found " << array->GetNumberOfTuples() << " expected: " << size);
        }

      for ( zz = extent[4]; zz <= extent[5]; ++ zz)
        {
        for ( yy = extent[2]; yy <= extent[3]; ++ yy)
          {
          for ( xx = extent[0]; xx <= extent[1]; ++ xx)
            {
            if (
              xx >= updateExtent[0] && xx <= updateExtent[1] &&
              yy >= updateExtent[2] && yy <= updateExtent[3] &&
              zz >= updateExtent[4] && zz <= updateExtent[5])
              {
              if ( jj % 3 == 0 )
                {
                self->Indent(ost);
                }
              for ( kk = 0; kk < array->GetNumberOfComponents(); kk ++ )
                {
                val = array->GetValue(idx * array->GetNumberOfComponents() + kk);
                ost << " " << val;
                printOne ++;
                }
              jj ++;
              }
            idx ++;
            }
          }
        }
      if ( !printOne )
        {
        abort();
        }
      }
    else
      {
      for ( jj = 0; jj < array->GetNumberOfTuples(); jj ++ )
        {
        if ( jj % 3 == 0 )
          {
          self->Indent(ost);
          }
        for ( kk = 0; kk < array->GetNumberOfComponents(); kk ++ )
          {
          val = array->GetValue(jj * array->GetNumberOfComponents() + kk);
          ost << " " << val;
          }
        }
      }
    }
  else
    {
    const char *DataSetName;
    XdmfHDF H5;
    XdmfArray Data;
    XdmfInt64 h5dims[5];
    int nh5dims = 0;

  // Jerry
  // return( array->GetNumberOfTuples() );

    DataSetName = self->GenerateHDF5ArrayName(gridName, arrayName);
    if (!self->GetInputsArePieces() || (self->GetInputsArePieces() && self->CurrentInputNumber==0)) 
      {
      ost << " Format=\"HDF\">";
      self->Indent(ost);
      ost << " " << DataSetName;
      // self->Indent(ost);
      }

    if ( dims[0] < 1 )
      {
      h5dims[0] = array->GetNumberOfTuples();
      h5dims[1] = array->GetNumberOfComponents();
      nh5dims = 2;
      }
    else
      {
      h5dims[0] = dims[2];
      h5dims[1] = dims[1];
      h5dims[2] = dims[0];
      h5dims[3] = array->GetNumberOfComponents();
      nh5dims = 4;
      if ( h5dims[3] <= 1 )
        {
        nh5dims = 3;
        }
      //vtkDebugMacro(<< "Use dims: " << h5dims[0] << " " << h5dims[1] << " " << h5dims[2]<< " " << h5dims[3] << " (" << nh5dims << ")" << "\n");
      }

    Data.SetNumberType(type);
    //vtkDebugMacro(<< "Data type: " << type << "\n");
    Data.SetShape(nh5dims, h5dims);

    /*  This is too slow
    vtkIdType jj, kk;
    vtkIdType pos = 0;
    for ( jj = 0; jj < array->GetNumberOfTuples(); jj ++ )
      {
      for ( kk = 0; kk < array->GetNumberOfComponents(); kk ++ )
        {
        Data.SetValue(pos, array->GetValue(pos));
        pos ++;
        }
      }
    */

    // vtkGenericDebugMacro(<< "Required: " << array->GetNumberOfTuples() << " * " << array->GetNumberOfComponents() << " =" << (array->GetNumberOfTuples() *  array->GetNumberOfComponents()) << "\n");

    if ( useExtents )
      {
      vtkIdType jj = 0;
      vtkIdType idx = 0;
      vtkIdType kk;
      vtkIdType xx, yy, zz;
      vtkIdType size =
        (updateExtent[1]-updateExtent[0] + 1) *
        (updateExtent[3]-updateExtent[2] + 1) *
        (updateExtent[5]-updateExtent[4] + 1);
      // vtkGenericDebugMacro(<< "Size: " << size << " " << array->GetNumberOfComponents() << "\n");
      if ( size != array->GetNumberOfTuples() )
        {
        // Error
        vtkErrorWithObjectMacro(self, "Wrong number of tuples in the dataset. Found " << array->GetNumberOfTuples() << " expected: " << size);
        }
      for ( zz = extent[4]; zz <= extent[5]; ++ zz)
        {
        for ( yy = extent[2]; yy <= extent[3]; ++ yy)
          {
          for ( xx = extent[0]; xx <= extent[1]; ++ xx)
            {
            if (
              xx >= updateExtent[0] && xx <= updateExtent[1] &&
              yy >= updateExtent[2] && yy <= updateExtent[3] &&
              zz >= updateExtent[4] && zz <= updateExtent[5])
              {
              for ( kk = 0; kk < array->GetNumberOfComponents(); kk ++ )
                {
                val = array->GetValue(idx * array->GetNumberOfComponents() + kk);
                Data.SetValue(jj++, val);
                }
              }
            idx ++;
            }
          }
        }
      }
    else
      {
      Data.SetValues(0, array->GetPointer(0), array->GetNumberOfTuples() * array->GetNumberOfComponents());
      }

    H5.CopyType( &Data );

    int src_npts=-1;

    if (self->GetInputsArePieces()) 
      {
      // Set the Full Array Shape
      // vtkGenericDebugMacro(<< "Using Custom Disk Shape " << "\n");
      h5dims[0] = self->GetFullGridSize()[2];
      h5dims[1] = self->GetFullGridSize()[1];
      h5dims[2] = self->GetFullGridSize()[0];
      H5.SetShape(nh5dims, h5dims);
      // Set a hyperslab to define only our piece
      XdmfInt64 start[4]  = {extent[4], extent[2], extent[0], 0};
      XdmfInt64 stride[4] = {1,1,1,1};
      XdmfInt64 count[4]  = {dims[2], dims[1], dims[0], h5dims[3]};
      H5.SelectHyperSlab(start, stride, count);
      src_npts = H5Sget_select_npoints( H5.GetDataSpace() );

      } 
    else
      {
      H5.CopyShape( &Data );
      }

    if( H5.Open( DataSetName, "rw" ) == XDMF_FAIL )
      {
      if( H5.CreateDataset( DataSetName ) == XDMF_FAIL ) 
        {
        vtkErrorWithObjectMacro(self, "Can't Create Heavy Dataset " << DataSetName);
        return( -1 );
        }
      }

    // Set the DataSpace again as it may have been reset in Open Call
    H5.CopyType( &Data );
    src_npts=-1;

    if (self->GetInputsArePieces()) 
      {
      // Set the Full Array Shape
      // vtkGenericDebugMacro(<< "Using Custom Disk Shape " << "\n");
      h5dims[0] = self->GetFullGridSize()[2];
      h5dims[1] = self->GetFullGridSize()[1];
      h5dims[2] = self->GetFullGridSize()[0];
      H5.SetShape(nh5dims, h5dims);
      // Set a hyperslab to define only our piece
      XdmfInt64 start[4]  = {
        static_cast<XdmfInt64>(extent[4]), 
        static_cast<XdmfInt64>(extent[2]), 
        static_cast<XdmfInt64>(extent[0]), 
        0};
      XdmfInt64 stride[4] = {1,1,1,1};
      XdmfInt64 count[4]  = {dims[2], dims[1], dims[0], h5dims[3]};
      H5.SelectHyperSlab(start, stride, count);
      src_npts = H5Sget_select_npoints( H5.GetDataSpace() );

      } 
    else
      {
      H5.CopyShape( &Data );
      }


    // vtkGenericDebugMacro(<< "Writing Disk Shape (offset) " << extent[4] << " " << extent[2] << " " << extent[0] << "\n");
    H5.Write( &Data );
    H5.Close();

    }

  if (!self->GetInputsArePieces() || (self->GetInputsArePieces() && self->CurrentInputNumber==0)) 
    {
    self->DecrementIndent();
    self->Indent( ost );
    ost << "</DataItem>";
   }
  return( array->GetNumberOfTuples() );
}

//----------------------------------------------------------------------------
int vtkXdmfWriter::WriteDataArray( ostream& ost, vtkDataArray* array, vtkDataSet* dataSet,
  int dims[3], const char* Name, const char* Center, int type, const char* gridName, int active, int cellData )
{
  const char* arrayName = Name;
  if ( array->GetName() )
    {
    arrayName = array->GetName();
    }

  if (!this->InputsArePieces || (this->InputsArePieces && this->CurrentInputNumber==0)) 
    {
    ost << "<Attribute";
    this->IncrementIndent();
    this->Indent(ost);
    if ( active )
      {
      ost << " Active=\"1\"";
      this->Indent(ost);
      }
    switch ( type )
      {
    case XDMF_ATTRIBUTE_TYPE_SCALAR:
      ost << " AttributeType=\"Scalar\"";
      break;
    case XDMF_ATTRIBUTE_TYPE_VECTOR:
      ost << " AttributeType=\"Vector\"";
      break;
    case XDMF_ATTRIBUTE_TYPE_TENSOR:
      ost << " AttributeType=\"Tensor\"";
      break;
    case XDMF_ATTRIBUTE_TYPE_MATRIX:
      ost << " AttributeType=\"Matrix\"";
      break;
    default:
      ost << " AttributeType=\"Unknown\"";
      }
    this->Indent(ost);
    ost << " Center=\"" << Center << "\"";
    this->Indent(ost);
    ost << " Name=\"" << arrayName << "\">";
    this->Indent(ost);
    }

  vtkIdType res = this->WriteVTKArray( ost, array, dataSet, 0, dims, arrayName, 0, gridName,
    this->AllLight, cellData );

  if (!this->InputsArePieces || (this->InputsArePieces && this->CurrentInputNumber==0)) 
    {
    this->DecrementIndent();
    this->Indent(ost);
    ost << "</Attribute>";
    this->Indent(ost);
    }

  return res;
}

//----------------------------------------------------------------------------
int vtkXdmfWriter::WriteVTKArray( ostream& ost, vtkDataArray* array, vtkDataSet* dataSet, 
  int *scaledExtent, int dims[3], const char* Name, const char* dataName, const char* gridName,
  int alllight, int cellData )
{
  vtkIdType res = -1;
  switch ( array->GetDataType() )
    {
  case VTK_UNSIGNED_CHAR:
    res = vtkXdmfWriterWriteXMLScalar(this, ost, vtkUnsignedCharArray::SafeDownCast(array),
      dataSet, scaledExtent, dataName, Name, gridName, "UChar", static_cast<short>(0), alllight,
      XDMF_UINT8_TYPE, dims, cellData);
    break;
  case VTK_UNSIGNED_SHORT:
    res = vtkXdmfWriterWriteXMLScalar(this, ost, vtkUnsignedShortArray::SafeDownCast(array),
      dataSet, scaledExtent, dataName, Name, gridName, "UInt", static_cast<short>(0), alllight,
      XDMF_UINT16_TYPE, dims, cellData);
    break;
  case VTK_UNSIGNED_INT:
    res = vtkXdmfWriterWriteXMLScalar(this, ost, vtkUnsignedIntArray::SafeDownCast(array),
      dataSet, scaledExtent, dataName, Name, gridName, "UInt", static_cast<int>(0), alllight,
      XDMF_UINT32_TYPE, dims, cellData);
    break;
  case VTK_CHAR:
    res = vtkXdmfWriterWriteXMLScalar(this, ost, vtkCharArray::SafeDownCast(array),
      dataSet, scaledExtent, dataName, Name, gridName, "Char", static_cast<short>(0), alllight,
      XDMF_INT8_TYPE, dims, cellData);
    break;
  case VTK_SHORT:
    res = vtkXdmfWriterWriteXMLScalar(this, ost, vtkShortArray::SafeDownCast(array),
      dataSet, scaledExtent, dataName, Name, gridName, "Int", static_cast<short>(0), alllight,
      XDMF_INT16_TYPE, dims, cellData);
    break;
  case VTK_INT:
    res = vtkXdmfWriterWriteXMLScalar(this, ost, vtkIntArray::SafeDownCast(array),
      dataSet, scaledExtent, dataName, Name, gridName, "Int", static_cast<int>(0), alllight,
      XDMF_INT32_TYPE, dims, cellData);
    break;
  case VTK_FLOAT:
    res = vtkXdmfWriterWriteXMLScalar(this, ost, vtkFloatArray::SafeDownCast(array),
      dataSet, scaledExtent, dataName, Name, gridName, "Float", static_cast<float>(0), alllight,
      XDMF_FLOAT32_TYPE, dims, cellData);
    break;
  case VTK_DOUBLE:
    res = vtkXdmfWriterWriteXMLScalar(this, ost, vtkDoubleArray::SafeDownCast(array),
      dataSet, scaledExtent, dataName, Name, gridName, "Float", static_cast<double>(0), alllight,
      XDMF_FLOAT64_TYPE, dims, cellData);
    break;
  case VTK_ID_TYPE:
    res = vtkXdmfWriterWriteXMLScalar(this, ost, vtkIdTypeArray::SafeDownCast(array),
      dataSet, scaledExtent, dataName, Name, gridName, "Int", static_cast<vtkIdType>(0), alllight,
      sizeof(vtkIdType) == sizeof(int) ? XDMF_INT32_TYPE : XDMF_INT64_TYPE, dims, cellData);
    break;
  default:
    vtkErrorMacro("Unknown scalar type: " << array->GetDataType());
    }
  if ( res == -2 )
    {
    vtkErrorMacro("Cannot convert array to specified type");
    }
  return res;
}

//----------------------------------------------------------------------------
void vtkXdmfWriter::WriteAttributes( ostream& ost, vtkDataSet* ds, const char* gridName )
{
  int extent[6];
  int cdims[3] = { -1, -1, -1 };
  int pdims[3] = { -1, -1, -1 };

  ds->GetUpdateExtent(extent);
  //vtkDebugMacro(<< "Extent: " << extent[0] << " "<< extent[1] << " "<< extent[2] << " "<< extent[3] << " "<< extent[4] << " " << extent[5] << "\n");
  if ( extent[1] >= extent[0] && extent[3] >= extent[2] && extent[5] >= extent[4] )
    {
    cdims[0] = pdims[0] = extent[1] - extent[0] +1;
    cdims[1] = pdims[1] = extent[3] - extent[2] +1;
    cdims[2] = pdims[2] = extent[5] - extent[4] +1;
    cdims[0] --;
    cdims[1] --;
    cdims[2] --;
    //vtkDebugMacro(<< "pDims: " << pdims[0] << " " << pdims[1] << " " << pdims[2] << "\n");
    //vtkDebugMacro(<< "CDims: " << cdims[0] << " " << cdims[1] << " " << cdims[2] << "\n");
    if ( cdims[0] <= 0 && cdims[1] <= 0 && cdims[2] <= 0 )
      {
      // Bogus dimensions.
      pdims[0] = pdims[1] = pdims[2] = -1;
      cdims[0] = cdims[1] = cdims[2] = -1;
      }
    }

  vtkCellData *CellData = ds->GetCellData();
  vtkPointData *PointData = ds->GetPointData();
  int cc;

  if( CellData )
    {
    for (cc=0; cc<CellData->GetNumberOfArrays(); cc++)
      {
      vtkDataArray* array = CellData->GetArray(cc);
      int type = XDMF_ATTRIBUTE_TYPE_NONE;
      if ( array == CellData->GetScalars()  || array->GetNumberOfComponents() == 1 )
        {
        type = XDMF_ATTRIBUTE_TYPE_SCALAR;
        }
      else if ( array == CellData->GetVectors() || array->GetNumberOfComponents() == 3 )
        {
        type = XDMF_ATTRIBUTE_TYPE_VECTOR;
        }
      else if ( array == CellData->GetTensors() || array->GetNumberOfComponents() == 6 )
        {
        type = XDMF_ATTRIBUTE_TYPE_TENSOR;
        }
      int active = 0;
      if ( array == CellData->GetScalars() ||
        array == CellData->GetVectors() ||
        array == CellData->GetTensors() )
        {
        active = 1;
        }
      char buffer[100];
      sprintf(buffer, "UnnamedCellArray%d", cc);
      this->WriteDataArray( ost, array, ds, cdims, buffer, "Cell", type, gridName, active, 1 );
      }
    }
  if( PointData )
    {
    for (cc=0; cc<PointData->GetNumberOfArrays(); cc++)
      {
      vtkDataArray* array = PointData->GetArray(cc);
      int type = XDMF_ATTRIBUTE_TYPE_NONE;
      if ( array == PointData->GetScalars() || array->GetNumberOfComponents() == 1 )
        {
        type = XDMF_ATTRIBUTE_TYPE_SCALAR;
        }
      else if ( array == PointData->GetVectors() || array->GetNumberOfComponents() == 3 )
        {
        type = XDMF_ATTRIBUTE_TYPE_VECTOR;
        }
      else if ( array == PointData->GetTensors() || array->GetNumberOfComponents() == 6 )
        {
        type = XDMF_ATTRIBUTE_TYPE_TENSOR;
        }
      int active = 0;
      if ( array == PointData->GetScalars() ||
        array == PointData->GetVectors() ||
        array == PointData->GetTensors() )
        {
        active = 1;
        }
      char buffer[100];
      sprintf(buffer, "UnnamedNodeArray%d", cc);
      this->WriteDataArray( ost, array, ds, pdims, buffer, "Node", type, gridName, active, 0 );
      }
    }
}

//----------------------------------------------------------------------------
int vtkXdmfWriter::WriteGrid( ostream& ost, const char* gridName, vtkDataSet* ds, 
  void* mapofcell, const void* celltype /* = 0 */)
{
  int type; 

  if( !ds) 
    {
    vtkErrorMacro("No Input Data Set");
    return( -1 );
    }
  // Only do XMl for the first piece if these are pieces of data
  if (!this->InputsArePieces || (this->InputsArePieces && this->CurrentInputNumber==0)) 
    {

    ost << "<Grid Name=\"" << gridName << "\"";
    if ( this->CollectionName )
      {
      ost << " Collection=\"" << this->CollectionName << "\"";
      }
    ost << ">";
    this->IncrementIndent();
    this->Indent(ost);
    type = ds->GetDataObjectType();
    if ( type == VTK_POLY_DATA || type == VTK_UNSTRUCTURED_GRID )
      {
      vtkPointSet *Polys = vtkPointSet::SafeDownCast(ds);
      const vtkXdmfWriterInternal::CellType* ct = (vtkXdmfWriterInternal::CellType*)celltype;
      vtkXdmfWriterInternal::MapOfCellTypes* mc = (vtkXdmfWriterInternal::MapOfCellTypes*)mapofcell;
      if ( ct == 0 )
        {
        ct = &mc->begin()->first;
        }
      this->StartTopology( ost, ct->VTKType, ct->NumPoints, (*mc)[*ct]->GetNumberOfIds() );
      this->WriteCellArray( ost, ds, gridName, mapofcell, ct );
      this->EndTopology( ost );
      this->Indent(ost);
      this->StartGeometry( ost, "XYZ" );
      this->WritePoints( ost, Polys->GetPoints(), ds, gridName);
      this->EndGeometry( ost );
      }
    else if ( type == VTK_STRUCTURED_POINTS || type == VTK_IMAGE_DATA || type == VTK_UNIFORM_GRID )
      {
      double Origin[3], Spacing[3];
      //int     Extent[6];
      int     updateExtent[6], wholeextent[6];
      vtkImageData *SGrid = static_cast<vtkImageData *>(ds);
      //SGrid->GetDimensions( Dims );
      SGrid->GetOrigin( Origin );
      SGrid->GetSpacing( Spacing );
      SGrid->GetWholeExtent(wholeextent);
      //SGrid->GetExtent( Extent );
      SGrid->GetUpdateExtent( updateExtent );
      if (this->InputsArePieces) 
      {
        this->StartTopology(ost, "3DCORECTMESH", 3, this->FullGridSize);
      }
      else {
        int Dims[3];
        Dims[0] = updateExtent[1] - updateExtent[0] + 1;
        Dims[1] = updateExtent[3] - updateExtent[2] + 1;
        Dims[2] = updateExtent[5] - updateExtent[4] + 1;
        this->StartTopology(ost, "3DCORECTMESH", 3, Dims);
      }
      this->EndTopology(ost);
      this->Indent(ost);
      this->StartGeometry( ost, "ORIGIN_DXDYDZ" );
      this->Indent(ost);

      int cc;
      for ( cc = 0; cc < 3; cc ++ )
        {
        Origin[cc] = Origin[cc] + Spacing[cc] * updateExtent[cc * 2];
        }

      // Origin
      ost << "<DataItem";
      this->IncrementIndent();
      this->Indent(ost);
      ost << " Name=\"Origin\"";
      this->Indent(ost);
      ost << " NumberType=\"Float\"";
      this->Indent(ost);
      ost << " Dimensions=\"3\"";
      this->Indent(ost);
      ost << " Format=\"XML\">";
      this->Indent(ost);
      ost << Origin[2] << " " << Origin[1] << " " << Origin[0];
      vtkDebugMacro(<< "-- Origin: " << Origin[2] << " " << Origin[1] << " " << Origin[0] << "\n");
      this->DecrementIndent();
      this->Indent(ost);
      ost << "</DataItem>";
      this->Indent(ost);

      // DX DY DZ
      ost << "<DataItem";
      this->IncrementIndent();
      this->Indent(ost);
      ost << " Name=\"Spacing\"";
      this->Indent(ost);
      ost << " NumberType=\"Float\"";
      this->Indent(ost);
      ost << " Dimensions=\"3\"";
      this->Indent(ost);
      ost << " Format=\"XML\">";
      this->Indent(ost);
      ost << Spacing[2] << " " << Spacing[1] << " " << Spacing[0];
      this->DecrementIndent();
      this->Indent(ost);
      ost << "</DataItem>";

      this->EndGeometry( ost );
      }
    else if ( type == VTK_STRUCTURED_GRID )
      {
      int     Dims[3];
      vtkStructuredGrid *SGrid = static_cast<vtkStructuredGrid *>(ds);
      SGrid->GetDimensions( Dims );
      this->StartTopology(ost, "3DSMESH", 3, Dims);
      this->EndTopology(ost);
      this->Indent(ost);
      this->StartGeometry(ost, "XYZ");
      this->WritePoints( ost, SGrid->GetPoints(), ds, gridName );
      this->EndGeometry(ost);
      }
    else if ( type == VTK_RECTILINEAR_GRID )
      {
      int     Dims[3];
      int     updateExtent[6];
      vtkDataArray  *Coord;
      vtkRectilinearGrid *RGrid = static_cast<vtkRectilinearGrid *>(ds);
      RGrid->GetDimensions( Dims );
      RGrid->GetUpdateExtent( updateExtent );
      Dims[0] = updateExtent[1] - updateExtent[0] + 1;
      Dims[1] = updateExtent[3] - updateExtent[2] + 1;
      Dims[2] = updateExtent[5] - updateExtent[4] + 1;
      this->StartTopology(ost, "3DRECTMESH", 3, Dims);
      this->EndTopology( ost );
      this->Indent(ost);
      this->StartGeometry(ost, "VXVYVZ");
      int dummydims[3];
      dummydims[0] = dummydims[1] = dummydims[2] = -1;
      // X Coordinated
      int scaledExtent[6];
      int extent[6];
      scaledExtent[2] = scaledExtent[4] = 0;
      scaledExtent[3] = scaledExtent[5] = 0;
      RGrid->GetExtent(extent);
      vtkDebugMacro(<< "Update extent: " << updateExtent[0] << " " << updateExtent[1] << " " << updateExtent[2] << " " << updateExtent[3] << " " << updateExtent[4]  << " " << updateExtent[5] << "\n");
      vtkDebugMacro(<< "Extent:        " << extent[0] << " " << extent[1] << " " << extent[2] << " " << extent[3] << " " << extent[4]  << " " << extent[5] << "\n");
      Coord = RGrid->GetXCoordinates();
//      Coord->Print(co ut);
      scaledExtent[0] = updateExtent[0] - extent[0];
      scaledExtent[1] = updateExtent[1] - extent[0];
      vtkDebugMacro(<< "Scaled Extent: " << scaledExtent[0] << " " << scaledExtent[1] << "\n");
      this->WriteVTKArray( ost, Coord, ds, scaledExtent, dummydims, "X Coordinates", "X", gridName, !this->AllHeavy);
      this->Indent(ost);
      // Y Coordinated
      Coord = RGrid->GetYCoordinates();
      scaledExtent[0] = updateExtent[2] - extent[2];
      scaledExtent[1] = updateExtent[3] - extent[2];
      this->WriteVTKArray( ost, Coord, ds, scaledExtent, dummydims, "Y Coordinates", "Y", gridName, !this->AllHeavy);
      this->Indent(ost);
      // Z Coordinated
      Coord = RGrid->GetZCoordinates();
      scaledExtent[0] = updateExtent[4] - extent[4];
      scaledExtent[1] = updateExtent[5] - extent[4];
      this->WriteVTKArray( ost, Coord, ds, scaledExtent, dummydims, "Z Coordinates", "Z", gridName, !this->AllHeavy);
      this->EndGeometry(ost);
      }
    this->Indent(ost);
    }

    if (!this->InputsArePieces || (this->InputsArePieces && this->CurrentInputNumber==0)) 
      {
      if (this->TimeValue!=-1E99) 
        {
        ost << "<Time Value=\"" << this->TimeValue << "\" />";
        this->Indent(ost);
        }
    }

  // do this for all pieces
  this->WriteAttributes(ost, ds, gridName);

  // Only do XMl for the first piece if these are pieces of data
  if (!this->InputsArePieces || (this->InputsArePieces && this->CurrentInputNumber==0)) 
    {
    this->DecrementIndent();
    this->Indent(ost);
    ost << "</Grid>";
    }

  return( 1 );
}

//----------------------------------------------------------------------------
void vtkXdmfWriter::Write()
{
  if ( this->AllLight && this->AllHeavy )
    {
    vtkErrorMacro("AllLight and AllHeavy are mutually exclusive. Please only select one.");
    return;
    }

  if ( !this->FileNameString )
    {
    vtkErrorMacro("No file name specified");
    return;
    }

  int numberOfInputs = this->GetNumberOfInputs();
  if ( numberOfInputs <= 0 )
    {
    vtkErrorMacro("No input or input of the wrong type");
    return;
    }

  bool exists = vtksys::SystemTools::FileExists(this->FileNameString);
  bool readHeader = false;
  if (this->AppendGridsToDomain && exists)
    {
/*
      this->ParseExistingFile(this->FileNameString);
      XdmfConstString test = this->DOM->Serialize(this->DOM->GetTree());
      if (test) {
        ofstream dummy("c:/xdmftest.xml");
        dummy << test;
      }
*/
      readHeader = this->ReadDocument(this->FileNameString);
    }

  ofstream ofs(this->FileNameString);
  if ( !ofs )
    {
    vtkErrorMacro("Cannot open file: " << this->FileNameString);
    return;
    }

  vtkDebugMacro(<< "Write to file: " << this->FileNameString << "\n");

  //
  // Write header or copy existing file contents
  //
  if (readHeader)
    {
    ofs << this->DocString;
    this->IncrementIndent();
    if (this->CollectionType != XDMF_GRID_COLLECTION_UNSET) 
      {
      this->IncrementIndent();
      }
    }
  else 
    {
    if ( !this->GridOnly )
      {
      this->WriteHead(ofs);
      ofs << "<Domain";
      if ( this->DomainName )
        {
        ofs << " Name=\"" << this->DomainName << "\"";
        }
      ofs << ">";
      }
    if (this->CollectionType != XDMF_GRID_COLLECTION_UNSET) 
      {
      this->IncrementIndent();
      this->Indent(ofs);
      ofs << "<Grid GridType=\"" << "Collection\"";
      ofs << " CollectionType=\"";
      switch (this->CollectionType) 
        {
        case XDMF_GRID_COLLECTION_TEMPORAL:
          ofs << "Temporal\">"; break;
        case XDMF_GRID_COLLECTION_SPATIAL:
          ofs << "Spatial\">"; break;
        case XDMF_GRID_COLLECTION_UNSET:
          ofs << "Unset\">"; break;
        }
      }
    }

  //
  // Write out grids
  //
  this->IncrementIndent();
  for (this->CurrentInputNumber=0; 
    this->CurrentInputNumber<numberOfInputs; 
    this->CurrentInputNumber++)
    {
    //
    vtkDataSet* ds = vtkDataSet::SafeDownCast(this->Inputs[this->CurrentInputNumber]);
    ds->Update();

    vtkstd::string arrayName;
    vtkDataArray* da = ds->GetFieldData()->GetArray("Name");
    if ( da )
      {
      vtkCharArray* nameArray = vtkCharArray::SafeDownCast(da);
      if ( nameArray )
        {
        arrayName = static_cast<char*>(nameArray->GetVoidPointer(0));
        }
      }
    if ( arrayName.empty() )
      {
      if ( this->GridName )
        {
        arrayName = this->GridName;
        }
      else
        {
        arrayName = "Unnamed";
        }
      }

    if ( numberOfInputs>1 && !this->InputsArePieces)
      {
      char buffer[100];
      sprintf(buffer, "%d", this->CurrentInputNumber);
      arrayName += buffer;
      }

    vtkXdmfWriterInternal::MapOfCellTypes cellTypes;
    vtkXdmfWriterInternal::DetermineCellTypes(vtkPointSet::SafeDownCast(ds), cellTypes);
    if ( cellTypes.size() > 1 )
      {
      vtkErrorMacro("Xdmf Writer only supports unstructured data of single cell type");
      continue;
      }
    if (!this->InputsArePieces || (this->InputsArePieces && this->CurrentInputNumber==0)) 
      {
      this->Indent(ofs);
      }
    if ( cellTypes.size() > 1 )
      {
      vtkXdmfWriterInternal::MapOfCellTypes::iterator it;
      int ct = 0;
      for ( it = cellTypes.begin(); it != cellTypes.end(); ++it )
        {
        ostrstream str;
        str << arrayName.c_str() << "_" << ct << ends;
        this->WriteGrid(ofs, str.str(), ds, &cellTypes, &(*it));
        str.rdbuf()->freeze(2);
        ct ++;
        }
      }
    else
      {
      this->WriteGrid(ofs, arrayName.c_str(), ds, &cellTypes, 0);
      }
    }

  this->DecrementIndent();
  if ( !this->GridOnly )
    {
    this->Indent( ofs );
    ofs << "</Domain>" << "\n";
    this->WriteTail(ofs);
    }
}
//----------------------------------------------------------------------------
void vtkXdmfWriter::CloseCollection()
{
  bool exists = vtksys::SystemTools::FileExists(this->FileNameString);
  if (this->AppendGridsToDomain && exists)
    {
    this->ReadDocument(this->FileNameString);
    ofstream ofs(this->FileNameString);
    if ( !ofs )
      {
      vtkErrorMacro("Cannot open file: " << this->FileNameString);
      return;
      }
    ofs << this->DocString;
    this->IncrementIndent();
    this->Indent( ofs );
    ofs << "</Grid>" << "\n";
    this->DecrementIndent();
    if ( !this->GridOnly )
      {
      this->Indent( ofs );
      ofs << "</Domain>" << "\n";
      this->WriteTail(ofs);
      }
    }
}
//----------------------------------------------------------------------------
void vtkXdmfWriter::SetInput(vtkDataSet* ds)
{
  this->SetNthInput(0, ds);
}

//----------------------------------------------------------------------------
// Add a dataset to the list of data to append.
void vtkXdmfWriter::AddInput(vtkDataObject *ds)
{
  this->vtkProcessObject::AddInput(ds);
}

//----------------------------------------------------------------------------
vtkDataObject *vtkXdmfWriter::GetInput(int idx)
{
  if (idx >= this->NumberOfInputs || idx < 0)
    {
    return NULL;
    }

  return (vtkDataObject *)(this->Inputs[idx]);
}

//----------------------------------------------------------------------------
// Remove a dataset from the list of data to append.
void vtkXdmfWriter::RemoveInput(vtkDataObject *ds)
{
  this->vtkProcessObject::RemoveInput(ds);
}

//----------------------------------------------------------------------------
vtkDataSetCollection *vtkXdmfWriter::GetInputList()
{
  int idx;

  if (this->InputList)
    {
    this->InputList->Delete();
    }
  this->InputList = vtkDataSetCollection::New();

  for (idx = 0; idx < this->NumberOfInputs; ++idx)
    {
    if (this->Inputs[idx] != NULL)
      {
      this->InputList->AddItem(static_cast<vtkDataSet*>((this->Inputs[idx])));
      }
    }  

  return this->InputList;
}

//----------------------------------------------------------------------------
void vtkXdmfWriter::StartGeometry( ostream& ost, const char* type )
{
  ost << "<Geometry GeometryType=\"" << type << "\">";
  this->IncrementIndent();
  this->Indent(ost);
}

//----------------------------------------------------------------------------
void vtkXdmfWriter::EndGeometry( ostream& ost )
{
  this->DecrementIndent();
  this->Indent(ost);
  ost << "</Geometry>";
}

//----------------------------------------------------------------------------
void vtkXdmfWriter::Indent(ostream& os)
{
  int cc;

  os << "\n";
  for ( cc = 0; cc < this->CurrIndent; cc ++ )
    {
    os << "  ";
    }
}

//----------------------------------------------------------------------------
const char* vtkXdmfWriter::GenerateHDF5ArrayName(const char* gridName, 
                                                 const char* array)
{
  if ( !this->HeavyDataSetNameString )
    {
    vtkErrorMacro("HeavyDataSetName is not yet specified");
    return 0;
    }
  size_t namelen = strlen(this->HeavyDataSetNameString) + strlen(array);
  if ( gridName )
    {
    namelen += strlen(gridName);
    }
  char *name = new char [ namelen + 10 ];
  // Should Use the ENTITY HeavyData
  if ( gridName )
    {
    sprintf(name, "%s:/%s/%s", this->HeavyDataSetNameString, gridName, array);
    //sprintf(name, "&HeavyData;:/%s/%s", gridName, array);
    }
  else
    {
    sprintf(name, "%s:/%s", this->HeavyDataSetNameString, array);
    //sprintf(name, "&HeavyData;:/%s", array);
    }
  this->SetHDF5ArrayName(name);
  delete [] name;
  return this->HDF5ArrayName;
}

//----------------------------------------------------------------------------
bool vtkXdmfWriter::ReadDocument(const char* filename)
{
  if (!vtksys::SystemTools::FileExists(filename))
  {
    return false;
  }
  ostrstream temp;
  ifstream infile(filename);
  char linebuffer[512];
  while (infile.good())
    {
    infile.getline(linebuffer,512); 
    if (vtksys::SystemTools::StringEndsWith(linebuffer, "</Domain>"))
      break;
    temp << linebuffer << "\n";
    }
  temp << ends;
  if (this->DocString)
  {
    delete []this->DocString;
  }
  this->DocString = strdup(temp.str());
  return true;
}
//----------------------------------------------------------------------------
void vtkXdmfWriter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
void vtkXdmfWriterInternal::DetermineCellTypes(vtkPointSet * t, vtkXdmfWriterInternal::MapOfCellTypes& vec)
{
  if ( !t )
    {
    return;
    }
  vtkIdType cc;
  vtkGenericCell* cell = vtkGenericCell::New();
  for ( cc = 0; cc < t->GetNumberOfCells(); cc ++ )
    {
    vtkXdmfWriterInternal::CellType ct;
    t->GetCell(cc, cell);
    ct.VTKType = cell->GetCellType();
    ct.NumPoints = cell->GetNumberOfPoints();
    vtkXdmfWriterInternal::MapOfCellTypes::iterator it = vec.find(ct);
    if ( it == vec.end() )
      {
      vtkIdList *l = vtkIdList::New();
      it = vec.insert(vtkXdmfWriterInternal::MapOfCellTypes::value_type(ct, 
          vtkSmartPointer<vtkIdList>(l))).first;
      l->Delete();
      }
    // it->second->InsertUniqueId(cc);;
    it->second->InsertNextId(cc);;
    }
  cell->Delete();
}

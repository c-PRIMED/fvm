    /*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkXdmfReader.cxx,v $
  Language:  C++
  Date:      $Date: 2009-04-06 13:16:03 $
  Version:   $Revision: 1.64 $


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
#include "vtkXdmfReader.h"

#include "vtkCallbackCommand.h"
#include "vtkDataArraySelection.h"
#include "vtkDataSet.h"
#include "vtkDoubleArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkRectilinearGrid.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStructuredGrid.h"
#include "vtkUnstructuredGrid.h"
#include "vtkXdmfDataArray.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkCellArray.h"
#include "vtkCharArray.h"
#include "vtkXMLParser.h"
#include "vtkImageData.h"
#include "vtkUniformGrid.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkCompositeDataPipeline.h"
#include "vtkMultiProcessController.h"
#include "vtkDataObjectTypes.h"
#include "vtkSmartPointer.h"

#include "XdmfArray.h"
#include "XdmfAttribute.h"
#include "XdmfDOM.h"
#include "XdmfDataDesc.h"
#include "XdmfDataItem.h"
#include "XdmfGrid.h"
#include "XdmfTopology.h"
#include "XdmfGeometry.h"
#include "XdmfTime.h"

#include <sys/stat.h>
#include <vtkstd/set>
#include <vtkstd/map>
#include <vtkstd/string>
#include <vtkstd/vector>
#include <vtksys/SystemTools.hxx>
#include <assert.h>
#include <functional>
#include <algorithm>

#define USE_IMAGE_DATA // otherwise uniformgrid

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkXdmfReader);
vtkCxxRevisionMacro(vtkXdmfReader, "$Revision: 1.64 $");

//----------------------------------------------------------------------------
vtkCxxSetObjectMacro(vtkXdmfReader,Controller,vtkMultiProcessController);

#if defined(_WIN32) && (defined(_MSC_VER) || defined(__BORLANDC__))
#  include <direct.h>
#  define GETCWD _getcwd
#else
#include <unistd.h>
#  define GETCWD getcwd
#endif

#define vtkMAX(x, y) (((x)>(y))?(x):(y))
#define vtkMIN(x, y) (((x)<(y))?(x):(y))

#define PRINT_EXTENT(x) "[" << (x)[0] << " " << (x)[1] << " " << (x)[2] << " " << (x)[3] << " " << (x)[4] << " " << (x)[5] << "]" 
//----------------------------------------------------------------------------
class vtkXdmfReaderGrid
{
public:
  //
  vtkXdmfReaderGrid() 
  {
      this->XMGrid = 0;
      this->DataDescription = 0;
      this->Enabled = 0;

      // this->Enabled = 1;
      this->vtkType = -1;
      this->Time = -1.0;
      this->TimeIndex = 0;
      this->isTemporal = 0; // 1 if Temporal 2 if Tempoal Collection
      this->isCollection = 0;
      this->isParallel = 0;

  }

  ~vtkXdmfReaderGrid()
  {
    delete this->XMGrid;
    vtkstd::vector<vtkXdmfReaderGrid*>::iterator it;
    for (it = Children.begin(); it != Children.end(); it++)
      {
      delete *it;
      }
  }

  vtkXdmfReaderGrid* GetChild(const char *childName);
  vtkXdmfReaderGrid* GetChild(int idx);
  vtkInformation* GetInformation()
    {
    if (!this->Information) 
      {
      this->Information = vtkSmartPointer<vtkInformation>::New();
      }
    return this->Information;
    }

  XdmfGrid       *XMGrid;
  XdmfDataDesc   *DataDescription;
  XdmfFloat64   Time;
  XdmfInt32     TimeIndex;
  XdmfInt32     isTemporal;
  XdmfInt32     isCollection;
  XdmfInt32     isParallel; // Read in parallel if possible
  vtkstd::string Name;
  int Enabled;
  int vtkType;

  vtkstd::vector<vtkXdmfReaderGrid*> Children;
  vtkSmartPointer<vtkInformation> Information;
};

//----------------------------------------------------------------------------
class vtkXdmfReaderMatchName : public vtkstd::binary_function<vtkXdmfReaderGrid *, const char *, bool>
{
  public:
    bool operator()(vtkXdmfReaderGrid *child, const char *name) const
    {
      return (child->Name==name);
    }
};
//----------------------------------------------------------------------------
vtkXdmfReaderGrid* vtkXdmfReaderGrid::GetChild(const char *gridName)
{
  vtkstd::vector<vtkXdmfReaderGrid *>::iterator it = 
    vtkstd::find_if(this->Children.begin(), this->Children.end(), 
    vtkstd::bind2nd(vtkXdmfReaderMatchName(),gridName));
  if (it!=this->Children.end()) 
  {
    return *it;
  }
  return NULL;
}

//----------------------------------------------------------------------------
vtkXdmfReaderGrid* vtkXdmfReaderGrid::GetChild(int idx)
{
  if ( idx<0 || idx>=(int)this->Children.size())
    {
    return NULL;
    }
  return this->Children[idx];
/*
  vtkstd::vector<vtkXdmfReaderGrid *>::iterator it;
  int cnt = 0;
  for ( it = this->Children.begin(); it != this->Children.end(); ++ it )
    {
    if ( cnt == idx )
      {
      return *it;
      }
    cnt ++;
    }
  return NULL;
*/
}

//============================================================================
class vtkXdmfReaderInternal
{
public:
  vtkXdmfReaderInternal()
  {
    this->DataItem = NULL;
    this->ArrayConverter = vtkXdmfDataArray::New();
    this->Data = NULL;
    this->DsmBuffer = NULL;
    this->InputString = 0;
    this->InputStringLength = 0;
    // this->ParallelLevel = 0;
    this->ParallelLevels.clear();
    this->LargestLevel = 0;
    this->UpdatePiece     = 0;
    this->UpdateNumPieces = 1;
    this->mostChildren = 0;
  }

  ~vtkXdmfReaderInternal()
  {
    if ( this->Data )
      {
      delete this->Data;
      this->Data = NULL;
      }    
    if ( this->DataItem )
      {
      delete this->DataItem;
      this->DataItem = 0;
      }
    this->ArrayConverter->Delete();
    this->ArrayConverter = NULL;
    delete [] this->InputString;
    this->InputString = 0;
  }


  // Returns the domain node for the domain with given name. If none such domain
  // exists, returns the 1st domain, if any.
  XdmfXmlNode GetDomain(const char* domainName)
    {
    if (domainName)
      {
      vtkstd::map<vtkstd::string, XdmfXmlNode>::iterator iter =
        this->DomainMap.find(domainName);
      if (iter != this->DomainMap.end())
        {
        return iter->second;
        }
      }
    if (this->DomainList.size() > 0)
      {
      return this->GetDomain(this->DomainList[0].c_str());
      }
    return 0;
    }

  vtkXdmfReaderGrid* GetGrid(const char* gridName);
  vtkXdmfReaderGrid* GetGrid(int idx);
  vtkXdmfReaderGrid *AddGrid(
    vtkXdmfReaderGrid *parent,
    const char *gridName);
  void DeleteChildren(vtkXdmfReaderGrid* parent);

  // Temporary method to update the list of arrays provided by the grid.
  int UpdateArrays(vtkXdmfReaderGrid *grid);
  int RequestGridInformation(vtkXdmfReaderGrid *grid, vtkInformation *destInfo);
  int FindParallelism(vtkXdmfReaderGrid *grid = 0);

  int RequestGridData(/*const char* currentGridName,*/
    vtkXdmfReaderGrid *grid,
    vtkDataObject *output,
    int timeIndex,
    int isSubBlock,
    double progressS, double progressE);
  
 
  vtkstd::vector<XdmfFloat64> TimeValues;
  vtkstd::vector<vtkstd::string> DomainList;
  vtkstd::map<vtkstd::string, XdmfXmlNode> DomainMap;
  vtkXdmfReaderGrid *Data;
  // vtkXdmfReaderGrid *ParallelLevel;
  vtkstd::vector<vtkXdmfReaderGrid*> ParallelLevels;
  vtkXdmfReaderGrid *LargestLevel;
  vtkXdmfReader* Reader;
  XdmfDataItem *DataItem;
  XdmfDsmBuffer *DsmBuffer;
  char *InputString;
  unsigned int InputStringLength;
  unsigned int mostChildren;

  // For converting arrays from XDMF to VTK format
  vtkXdmfDataArray *ArrayConverter;

  unsigned int UpdatePiece;
  unsigned int UpdateNumPieces;

};

//----------------------------------------------------------------------------
vtkXdmfReaderGrid* vtkXdmfReaderInternal::GetGrid(const char* gridName)
{
  if (!this->Data)
    {
    return NULL;
    }
  return this->Data->GetChild(gridName);
}

//----------------------------------------------------------------------------
vtkXdmfReaderGrid* vtkXdmfReaderInternal::GetGrid(int idx)
{
  if (!this->Data)
    {
    return NULL;
    }
  return this->Data->GetChild(idx);
}

//----------------------------------------------------------------------------
vtkXdmfReaderGrid* 
  vtkXdmfReaderInternal::AddGrid(
    vtkXdmfReaderGrid* parent,    
    const char *gridName)
{
  if ( !parent || !gridName )
    {
    return 0;
    }

  //finds or makes up an ActualGrid in the global list
//  vtkXdmfReaderGrid* grid = parent->GetChild(gridName);
//  if ( !grid )
//    {
    vtkXdmfReaderGrid* grid = new vtkXdmfReaderGrid; 
    grid->Name = gridName;
    parent->Children.push_back(grid);
//    }
  return grid;
}
void vtkXdmfReaderInternal::DeleteChildren(vtkXdmfReaderGrid* parent)
{
  if (!parent)
    {
    return;
    }

  vtkstd::vector<vtkXdmfReaderGrid*>::iterator it;
  for (it = parent->Children.begin(); it != parent->Children.end(); it++)
    {
    delete *it;
    }
  parent->Children.clear();
}

//============================================================================
class vtkXdmfReaderTester : public vtkXMLParser
{
public:
  vtkTypeMacro(vtkXdmfReaderTester, vtkXMLParser);
  static vtkXdmfReaderTester* New();
  int TestReadFile()
    {
      this->Valid = 0;
      if(!this->FileName)
        {
        return 0;
        }

      ifstream inFile(this->FileName);
      if(!inFile)
        {
        return 0;
        }

      this->SetStream(&inFile);
      this->Done = 0;

      this->Parse();

      if(this->Done && this->Valid )
        {
        return 1;
        }
      return 0;
    }
  void StartElement(const char* name, const char**)
    {
      this->Done = 1;
      if(strcmp(name, "Xdmf") == 0)
        {
        this->Valid = 1;
        }
    }

protected:
  vtkXdmfReaderTester()
    {
      this->Valid = 0;
      this->Done = 0;
    }

private:
  void ReportStrayAttribute(const char*, const char*, const char*) {}
  void ReportMissingAttribute(const char*, const char*) {}
  void ReportBadAttribute(const char*, const char*, const char*) {}
  void ReportUnknownElement(const char*) {}
  void ReportXmlParseError() {}

  int ParsingComplete() { return this->Done; }
  int Valid;
  int Done;
  vtkXdmfReaderTester(const vtkXdmfReaderTester&); // Not implemented
  void operator=(const vtkXdmfReaderTester&); // Not implemented
};

vtkStandardNewMacro(vtkXdmfReaderTester);

//============================================================================
vtkXdmfReader::vtkXdmfReader()
{  
  this->Internals = new vtkXdmfReaderInternal;
  this->Internals->Reader = this;

  this->DOM = 0;

  this->DomainName = 0;
  this->GridName = 0;
  this->NumberOfEnabledActualGrids=0;
  this->DsmBuffer = 0;
  
  for (int i = 0; i < 3; i ++ )
    {
    this->Stride[i] = 1;
    }

  this->Controller = 0;
  this->SetController(vtkMultiProcessController::GetGlobalController());

  // Setup the selection callback to modify this object when an array
  // selection is changed.
  this->PointDataArraySelection = vtkDataArraySelection::New();
  this->CellDataArraySelection = vtkDataArraySelection::New();
  this->SelectionObserver = vtkCallbackCommand::New();
  this->SelectionObserver->SetCallback(
    &vtkXdmfReader::SelectionModifiedCallback);
  this->SelectionObserver->SetClientData(this);
  this->PointDataArraySelection->AddObserver(vtkCommand::ModifiedEvent,
                                             this->SelectionObserver);
  this->CellDataArraySelection->AddObserver(vtkCommand::ModifiedEvent,
                                            this->SelectionObserver);

  this->OutputVTKType  = -1;
  this->OutputTemporal = 0;
  this->TimeStep       = 0;
  this->ActualTimeStep = 0;
  this->TimeStepRange[0] = 0;
  this->TimeStepRange[1] = 0;
  // this->DebugOn();
}

//----------------------------------------------------------------------------
vtkXdmfReader::~vtkXdmfReader()
{
  this->CellDataArraySelection->RemoveObserver(this->SelectionObserver);
  this->PointDataArraySelection->RemoveObserver(this->SelectionObserver);
  this->SelectionObserver->Delete();
  this->CellDataArraySelection->Delete();
  this->PointDataArraySelection->Delete();

  delete [] this->DomainName;
  this->DomainName = 0;
  
  delete this->Internals;

  if ( this->DOM )
    {
    delete this->DOM;
    }

  H5garbage_collect();  
  this->SetController(0);
}

void
vtkXdmfReader::SetDsmBuffer(void *Bufp){
    this->DsmBuffer = (XdmfDsmBuffer *)Bufp;
}

void *
vtkXdmfReader::GetDsmBuffer(){
    return(this->DsmBuffer);
}

//----------------------------------------------------------------------------
void vtkXdmfReader::SelectionModifiedCallback(vtkObject*, unsigned long,
                                              void* clientdata, void*)
{
  static_cast<vtkXdmfReader*>(clientdata)->Modified();
}

//----------------------------------------------------------------------------
int vtkXdmfReader::GetNumberOfPointArrays()
{
  return this->PointDataArraySelection->GetNumberOfArrays();
}

//----------------------------------------------------------------------------
int vtkXdmfReader::GetNumberOfCellArrays()
{
  return this->CellDataArraySelection->GetNumberOfArrays();
}

//----------------------------------------------------------------------------
const char* vtkXdmfReader::GetPointArrayName(int index)
{
  return this->PointDataArraySelection->GetArrayName(index);
}

//----------------------------------------------------------------------------
const char* vtkXdmfReader::GetCellArrayName(int index)
{
  return this->CellDataArraySelection->GetArrayName(index);
}

//----------------------------------------------------------------------------
int vtkXdmfReader::GetPointArrayStatus(const char* name)
{
  return this->PointDataArraySelection->ArrayIsEnabled(name);
}

//----------------------------------------------------------------------------
int vtkXdmfReader::GetCellArrayStatus(const char* name)
{
  return this->CellDataArraySelection->ArrayIsEnabled(name);
}

//----------------------------------------------------------------------------
void vtkXdmfReader::SetPointArrayStatus(const char* name, int status)
{
  vtkDebugMacro("Set point array \"" << name << "\" status to: " << status);
  if(status)
    {
    this->PointDataArraySelection->EnableArray(name);
    }
  else
    {
    this->PointDataArraySelection->DisableArray(name);
    }
}

//----------------------------------------------------------------------------
void vtkXdmfReader::SetCellArrayStatus(const char* name, int status)
{
  vtkDebugMacro("Set cell array \"" << name << "\" status to: " << status);
  if(status)
    {
    this->CellDataArraySelection->EnableArray(name);
    }
  else
    {
    this->CellDataArraySelection->DisableArray(name);
    }
}

//----------------------------------------------------------------------------
void vtkXdmfReader::EnableAllArrays()
{
  vtkDebugMacro("Enable all point and cell arrays");
  this->PointDataArraySelection->EnableAllArrays();
  this->CellDataArraySelection->EnableAllArrays();
}

//----------------------------------------------------------------------------
void vtkXdmfReader::DisableAllArrays()
{
  vtkDebugMacro("Disable all point and cell arrays");
  this->PointDataArraySelection->DisableAllArrays();
  this->CellDataArraySelection->DisableAllArrays();
}

//----------------------------------------------------------------------------
int vtkXdmfReader::GetNumberOfParameters()
{
  if(!this->DOM) 
    {
    return(-1);
    }
  // return(this->DOM->FindNumberOfParameters());
  return(0);
}

//----------------------------------------------------------------------------
int vtkXdmfReader::GetParameterType(const char *vtkNotUsed(parameterName))
{
  /*
    XdmfParameter *Param;
    if(!this->DOM) 
    {
    return(0);
    }
    Param = this->DOM->FindParameter(parameterName);
    if(Param) 
    {
    return(Param->GetParameterType());
    } 
    else 
    {
    return(0);
    }
  */
  return(0);
}
//----------------------------------------------------------------------------
const char *vtkXdmfReader::GetParameterTypeAsString(const char *vtkNotUsed(parameterName))
{
  /*
    if (this->GetParameterType(parameterName) == XDMF_PARAMETER_RANGE_TYPE) 
    {
    return("RANGE");
    } 
  */
  return("LIST");
}
//----------------------------------------------------------------------------
int vtkXdmfReader::GetParameterType(int vtkNotUsed(index))
{
  /*
    XdmfParameter *Param;
    if(!this->DOM) 
    {
    return(0);
    }
    Param = this->DOM->GetParameter(index);
    if(Param) 
    {
    return(Param->GetParameterType());
    } 
    else 
    {
    return(0);
    }
  */
  return(0);
}

//----------------------------------------------------------------------------
const char *vtkXdmfReader::GetParameterTypeAsString(int vtkNotUsed(index))
{
  /*
    if (this->GetParameterType(index) == XDMF_PARAMETER_RANGE_TYPE) 
    {
    return("RANGE");
    } 
  */
  return("LIST");
}

//----------------------------------------------------------------------------
int vtkXdmfReader::GetParameterRange(int vtkNotUsed(index), int vtkNotUsed(Shape)[3])
{
  /*
    XdmfParameter *Param;
    XdmfArray  *Parray;
    if(!this->DOM) 
    {
    return(0);
    }
    Param = this->DOM->GetParameter(index);
    if(Param) 
    {
    if( Param->GetParameterType() == XDMF_PARAMETER_RANGE_TYPE )
    {
    Parray = Param->GetArray();
    Shape[0] = Parray->GetValueAsInt64(0);
    Shape[1] = Parray->GetValueAsInt64(1);
    Shape[2] = Parray->GetValueAsInt64(2);
    } else {
    Shape[0] = 0;
    Shape[1] = 1;
    Shape[2] = Param->GetNumberOfElements();
    }
    return(Shape[2]);
    }
  */
  return(0);
}

//----------------------------------------------------------------------------
int vtkXdmfReader::GetParameterRange(const char *vtkNotUsed(parameterName), int vtkNotUsed(Shape)[3])
{
  /*
    XdmfParameter *Param;
    XdmfArray  *Parray;
    if(!this->DOM) 
    {
    return(0);
    }
    Param = this->DOM->FindParameter(parameterName);
    if(Param) 
    {
    if( Param->GetParameterType() == XDMF_PARAMETER_RANGE_TYPE )
    {
    Parray = Param->GetArray();
    Shape[0] = Parray->GetValueAsInt64(0);
    Shape[1] = Parray->GetValueAsInt64(1);
    Shape[2] = Parray->GetValueAsInt64(2);
    } 
    else 
    {
    Shape[0] = 0;
    Shape[1] = 1;
    Shape[2] = Param->GetNumberOfElements();
    }
    return(Shape[2]);
    }
  */
  return(0);
}

//----------------------------------------------------------------------------
const char *vtkXdmfReader::GetParameterRangeAsString(int vtkNotUsed(index))
{
  /*
    int Range[3];
    ostrstream StringOutput;
    if(this->GetParameterRange(index, Range) <= 0)
    {
    return(NULL);
    }
    StringOutput << ICE_64BIT_CAST Range[0] << " ";
    StringOutput << ICE_64BIT_CAST Range[1] << " ";
    StringOutput << ICE_64BIT_CAST Range[2] << ends;
    return(StringOutput.str());
  */
  return(0);
}

//----------------------------------------------------------------------------
const char *vtkXdmfReader::GetParameterRangeAsString(const char *vtkNotUsed(parameterName))
{
  /*
    int Range[3];
    ostrstream StringOutput;
    if (this->GetParameterRange(parameterName, Range) <= 0) 
    {
    return(NULL);
    }
    StringOutput << ICE_64BIT_CAST Range[0] << " ";
    StringOutput << ICE_64BIT_CAST Range[1] << " ";
    StringOutput << ICE_64BIT_CAST Range[2] << ends;
    return(StringOutput.str());
  */
  return(0);
}

//----------------------------------------------------------------------------
const char *vtkXdmfReader::GetParameterName(int vtkNotUsed(index))
{
  /*
    XdmfParameter *Param;
    if(!this->DOM) 
    {
    return(0);
    }
    Param = this->DOM->GetParameter(index);
    if(Param) 
    {
    return(Param->GetParameterName());
    } 
    else 
    {
    return(0);
    }
  */
  return(0);
}

//----------------------------------------------------------------------------
int vtkXdmfReader::SetParameterIndex(int vtkNotUsed(Index), int vtkNotUsed(CurrentIndex)) 
{
  /*
    XdmfParameter *Param;
    if(!this->DOM) 
    {
    return(0);
    }
    Param = this->DOM->GetParameter(Index);
    if(!Param) 
    {
    return(-1);
    }
    this->Modified();
    return(Param->SetCurrentIndex(CurrentIndex));
  */
  return(0);
}

//----------------------------------------------------------------------------
int vtkXdmfReader::GetParameterIndex(int vtkNotUsed(Index)) 
{
  /*
    XdmfParameter *Param;
    if(!this->DOM) 
    {
    return(0);
    }
    Param = this->DOM->GetParameter(Index);
    if(!Param) 
    {
    return(-1);
    }
    return(Param->GetCurrentIndex());
  */
  return(0);
}

//----------------------------------------------------------------------------
int vtkXdmfReader::SetParameterIndex(const char *vtkNotUsed(ParameterName),
                                     int vtkNotUsed(CurrentIndex)) 
{
  /*
    XdmfParameter *Param;
    int Status=-1;
    if(!this->DOM) 
    {
    return(0);
    }
    for(int i=0 ; i < this->DOM->FindNumberOfParameters() ;  i++)
    {
    Param = this->DOM->GetParameter(i);
    if(!Param) 
    {
    return(-1);
    }
    if(XDMF_WORD_CMP(Param->GetParameterName(), ParameterName))
    {
    Status = Param->SetCurrentIndex(CurrentIndex);
    this->Modified();
    if(Status <= 0 ) 
    {
    return(Status);
    }
    }
    }
    return(Status);
  */
  return(0);
}

//----------------------------------------------------------------------------
int vtkXdmfReader::GetParameterIndex(const char *vtkNotUsed(parameterName)) 
{
  /*
    XdmfParameter *Param;
    if(!this->DOM) 
    {
    return(0);
    }
    Param = this->DOM->FindParameter(parameterName);
    if(!Param) 
    {
    return(-1);
    }
    return(Param->GetCurrentIndex());
  */
  return(0);
}

//----------------------------------------------------------------------------
int vtkXdmfReader::GetParameterLength(int vtkNotUsed(index))
{
  /*
    XdmfParameter *Param;
    if(!this->DOM) 
    {
    return(0);
    }
    Param = this->DOM->GetParameter(index);
    if(Param) 
    {
    return(Param->GetNumberOfElements());
    } 
    else 
    {
    return(0);
    }
  */
  return(0);
}

//----------------------------------------------------------------------------
int vtkXdmfReader::GetParameterLength(const char *vtkNotUsed(parameterName))
{
  /*
    XdmfParameter *Param;
    if(!this->DOM) 
    {
    return(0);
    }
    Param = this->DOM->FindParameter(parameterName);
    if(Param) 
    {
    return(Param->GetNumberOfElements());
    } 
    else 
    {
    return(0);
    }
  */
  return 0;
}

//----------------------------------------------------------------------------
const char *vtkXdmfReader::GetParameterValue(const char *vtkNotUsed(parameterName)) 
{
  /*
    XdmfParameter *Param;
    if(!this->DOM) 
    {
    return(0);
    }
    Param = this->DOM->FindParameter(parameterName);
    if(!Param) 
    {
    return(0);
    }
    Param->Update();
    return(Param->GetParameterValue());
  */
  return 0;
}

//----------------------------------------------------------------------------
const char *vtkXdmfReader::GetParameterValue(int vtkNotUsed(index)) 
{
  /*
    XdmfParameter *Param;
    if(!this->DOM) 
    {
    return(0);
    }
    Param = this->DOM->GetParameter(index);
    if(!Param) 
    {
    return(0);
    }
    Param->Update();
    return(Param->GetParameterValue());
  */
  return 0;
}

//----------------------------------------------------------------------------
int vtkXdmfReader::GetNumberOfDomains()
{
  return this->Internals->DomainList.size();
}

//----------------------------------------------------------------------------
void vtkXdmfReader::SetDomainName(const char* domain)
{
  if ( this->DomainName == domain )
    {
    return;
    }
  if ( this->DomainName && domain && strcmp(this->DomainName, domain) == 0 )
    {
    return;
    }

  delete [] this->DomainName;
  this->DomainName = 0;
  if (domain)
    {
    this->DomainName = new char [ strlen(domain) + 1 ];
    strcpy(this->DomainName, domain);
    }

  // Since selected domain changed, we need to populate the grids available in
  // this domain.
  this->UpdateRootGrid();
}

//----------------------------------------------------------------------------
const char* vtkXdmfReader::GetDomainName(int idx)
{
  return this->Internals->DomainList[idx].c_str();
}

//----------------------------------------------------------------------------
void vtkXdmfReader::SetGridName(const char* grid)
{
  this->DisableAllGrids();
  this->EnableGrid(grid);
}

//----------------------------------------------------------------------------
int vtkXdmfReader::GetNumberOfGrids()
{
  if (this->Internals->Data)
    {
    return this->Internals->Data->Children.size();
    }
  return 0;
}

//----------------------------------------------------------------------------
const char* vtkXdmfReader::GetGridName(int idx)
{
  if ( idx<0 || idx>=(int)this->Internals->Data->Children.size())
    {
    return 0;
    }
  return this->Internals->Data->Children[idx]->Name.c_str();
/*
  vtkstd::vector<vtkXdmfReaderGrid *>::iterator it;
  int cnt = 0;
  for ( it = this->Internals->Data->Children.begin();
        it != this->Internals->Data->Children.end();
        ++ it )
    {
    if ( cnt == idx )
      {
      return (*it)->Name.c_str();
      }
    cnt ++;
    }
  return 0;
*/
}

//----------------------------------------------------------------------------
int vtkXdmfReader::GetGridIndex(const char* name)
{
  if ( !name )
    {
    return -1;
    }
  vtkstd::vector<vtkXdmfReaderGrid *>::iterator it;
  int cnt = 0;
  for ( it = this->Internals->Data->Children.begin();
        it != this->Internals->Data->Children.end();
        ++ it )
    {
    if ( (*it)->Name == name )
      {
      return cnt;
      }
    cnt ++;
    }

  return -1;
}

//----------------------------------------------------------------------------
void vtkXdmfReader::EnableGrid(const char* name)
{
  vtkDebugMacro("Enable grid " << name);
  vtkXdmfReaderGrid* grid = this->Internals->GetGrid(name);
  if ( !grid ) //work around an undo/redo crash
    {
    return;
    }
  if(!grid->Enabled)
    {
    ++this->NumberOfEnabledActualGrids;
    grid->Enabled = 1;
    this->Modified();
    //this->UpdateInformation();
    }
}

//----------------------------------------------------------------------------
void vtkXdmfReader::EnableGrid(int idx)
{
  vtkDebugMacro("Enable grid " << idx);
  vtkXdmfReaderGrid* grid = this->Internals->GetGrid(idx);
  if ( !grid )
    {
    return;
    }
  
  if(!grid->Enabled)
    {
    ++this->NumberOfEnabledActualGrids;
    grid->Enabled = 1;
    //this->PointDataArraySelection->RemoveAllArrays();
    //this->CellDataArraySelection->RemoveAllArrays();
    this->Modified();
    //this->UpdateInformation();
    }
}

//----------------------------------------------------------------------------
void vtkXdmfReader::EnableAllGrids()
{
  vtkDebugMacro("Enable all grids");
  vtkstd::vector<vtkXdmfReaderGrid *>::iterator it;
  int changed=0;
  for ( it = this->Internals->Data->Children.begin();
        it != this->Internals->Data->Children.end();
        ++ it )
    {
    if(!(*it)->Enabled)
      {
      (*it)->Enabled = 1;
      ++this->NumberOfEnabledActualGrids;
      changed=1;
      }
    }
  if(changed)
    {
    //this->PointDataArraySelection->RemoveAllArrays();
    //this->CellDataArraySelection->RemoveAllArrays();
    this->Modified();
    //this->UpdateInformation();
    }
}

//----------------------------------------------------------------------------
void vtkXdmfReader::DisableGrid(const char* name)
{
  vtkDebugMacro("Disable grid " << name);
  vtkXdmfReaderGrid* grid = this->Internals->GetGrid(name);
  if ( !grid )
    {
    return;
    }
  
  if(grid->Enabled)
    {
    grid->Enabled = 0;
    --this->NumberOfEnabledActualGrids;
    //this->PointDataArraySelection->RemoveAllArrays();
    //this->CellDataArraySelection->RemoveAllArrays();

    this->Modified();
    //this->UpdateInformation();
    }
}

//----------------------------------------------------------------------------
void vtkXdmfReader::DisableGrid(int idx)
{
  vtkDebugMacro("Disable grid " << idx);
  vtkXdmfReaderGrid* grid = this->Internals->GetGrid(idx);
  if ( !grid )
    {
    return;
    }
  
  if(grid->Enabled)
    {
    grid->Enabled = 0;
    --this->NumberOfEnabledActualGrids;
    //this->PointDataArraySelection->RemoveAllArrays();
    //this->CellDataArraySelection->RemoveAllArrays();
    this->Modified();
    //this->UpdateInformation();
    }
}

//----------------------------------------------------------------------------
void vtkXdmfReader::DisableAllGrids()
{
  vtkDebugMacro("Disable all grids");
  vtkstd::vector<vtkXdmfReaderGrid *>::iterator it;
  int changed=0;
  if (!this->Internals->Data) 
    {
    return;
    }
  for ( it = this->Internals->Data->Children.begin();
        it != this->Internals->Data->Children.end();
        ++ it )
    {
    if((*it)->Enabled)
      {
      (*it)->Enabled = 0;
      --this->NumberOfEnabledActualGrids;
      changed=1;
      }
    }
  if(changed)
    {
    //this->PointDataArraySelection->RemoveAllArrays();
    //this->CellDataArraySelection->RemoveAllArrays();
    this->Modified();
    //this->UpdateInformation();
    }
}

//----------------------------------------------------------------------------
void vtkXdmfReader::RemoveAllGrids()
{
  vtkDebugMacro("Remove all grids");
  vtkstd::vector<vtkXdmfReaderGrid *>::iterator it;
  if (this->Internals->Data)
    {
    delete this->Internals->Data;
    this->Internals->Data=NULL;
    }
  this->NumberOfEnabledActualGrids = 0;
  //this->PointDataArraySelection->RemoveAllArrays();
  //this->CellDataArraySelection->RemoveAllArrays();
  this->Modified();
  //this->UpdateInformation();
}

//----------------------------------------------------------------------------
int vtkXdmfReader::GetGridSetting(const char* name)
{
  vtkDebugMacro("GetGridSetting " << name);
  vtkXdmfReaderGrid* grid = this->Internals->GetGrid(name);
  if ( !grid )
    {
    return 0;
    }
  return grid->Enabled;
}

//----------------------------------------------------------------------------
int vtkXdmfReader::GetGridSetting(int idx)
{
  vtkDebugMacro("GetGridSetting " << idx);
  vtkXdmfReaderGrid* grid = this->Internals->GetGrid(idx);
  if ( !grid )
    {
    return 0;
    }
  return grid->Enabled;
}

//----------------------------------------------------------------------------
int vtkXdmfReader::UpdateDomains()
{
  // Populate the list of domains that the user can choose from.
  // And save a pointer that points to the one that the user chose.
  XdmfXmlNode domain = 0;
  int done = 0;
  //clear list
  this->Internals->DomainList.erase(this->Internals->DomainList.begin(),
                                    this->Internals->DomainList.end());
  this->Internals->DomainMap.clear();
  //populate list
  int cc;
  for ( cc = 0; !done; cc ++ )
    {
    ostrstream str1, str2;
    domain = this->DOM->FindElement("Domain", cc);
    if ( !domain )
      {
      break;
      }
    XdmfConstString domainName = this->DOM->Get( domain, "Name" );
    ostrstream str;
    if ( !domainName )
      {
      str << "Domain" << cc << ends;
      domainName = str.str();
      }
    this->Internals->DomainList.push_back(domainName);
    this->Internals->DomainMap[domainName] = domain;
    str.rdbuf()->freeze(0);
    }

  return 1;
}

//----------------------------------------------------------------------------
void vtkXdmfReader::UpdateRootGrid()
{
  cout << "UpdateRootGrid" << endl;
  if (!this->DomainName && this->GetNumberOfDomains() > 0)
    {
    vtkstd::string curDomain = this->GetDomainName(0);
    this->DomainName = new char [curDomain.length() + 1];
    strcpy(this->DomainName, curDomain.c_str());
    }
  XdmfXmlNode domain = this->Internals->GetDomain(this->DomainName);

  if ( !domain )
    {
    return;
    }

  vtkXdmfReaderGrid *ptr = this->Internals->Data;
  if (!ptr)
    {
    ptr = new vtkXdmfReaderGrid;
    ptr->Name = "DomainRoot";
    this->Internals->Data = ptr;
    }
  this->Internals->DeleteChildren(ptr);
  cout << "UpdateGrids Start" << endl;
  this->UpdateGrids(ptr, domain);
  cout << "UpdateGrids End" << endl;
  int nchildren = ptr->Children.size();
  this->OutputTemporal = 0;
  if (nchildren == 1)
    {
    this->OutputVTKType = this->Internals->GetGrid(0)->vtkType;    
    /*
    if (this->OutputVTKType==VTK_TEMPORAL_DATA_SET) 
      {
      this->OutputVTKType =  this->Internals->GetGrid(0)->GetChild(this->ActualTimeStep)->vtkType;
      this->OutputTemporal = 1;
      }
    */
    }
  else
    {
    this->OutputVTKType = VTK_MULTIBLOCK_DATA_SET;
    }

  // Repopulate the list of arrays.
  this->Internals->UpdateArrays(ptr);
}

//----------------------------------------------------------------------------
void vtkXdmfReader::UpdateGrids(vtkXdmfReaderGrid *parent, void *ParentNode)
{
  if (!parent || !ParentNode)
    {
    return;
    }
  

  int parentIsDomain = 0;
  XdmfXmlNode gridNode = 0;
  XdmfGrid  *xGrid;

  if(!parent->XMGrid){
      parentIsDomain = 1;
  }
  gridNode = this->DOM->FindElement("Grid", 0, (XdmfXmlNode)ParentNode);
  for (vtkIdType currentGrid = 0; gridNode; currentGrid ++ )
    {
    if ( !gridNode )
      {
      break;
      }

    // If Parent is Domain, create an XdmfGrid and UpdateInformation()
    // Otherwise since UpdateInformation() is recursive, there is
    // alread an XdmfGrid built internally.
    if(parentIsDomain){
        xGrid = new XdmfGrid;
        if(!xGrid){
            vtkErrorMacro("Can't create new XdmfGrid");
            break;
        }
        xGrid->SetName("DomainRoot");
        xGrid->SetDOM(this->DOM);
        xGrid->SetElement(gridNode);
        xGrid->UpdateInformation();
    }
    else
    {
        xGrid = parent->XMGrid->GetChild(currentGrid);
        if(!xGrid){
            vtkErrorMacro("parent->XMGrid is 0");
            break;
        }
        // Don't need to UpdateInformation(); it's recursive
    }
    
    // Make sure grid has a Name
    XdmfConstString gridName = xGrid->GetName();
    ostrstream str;
    if ( !gridName )
      {
      str << this->DOM->GetUniqueName("Grid") << ends;
      }
    else
      {
      str << gridName << ends;
      }
    gridName = str.str();
    
    vtkDebugMacro("Reading Light Data for " << gridName);
    // What Type of Grid
    XdmfInt32 gridType = xGrid->GetGridType();
    vtkXdmfReaderGrid *sub;

    if(gridType & XDMF_GRID_MASK)
    {
        // A collection or tree
      vtkDebugMacro(" Grid is a Collection/Tree ");
      sub = this->Internals->AddGrid(parent, gridName);
      sub->XMGrid = xGrid;
      sub->isCollection = 1;
      this->UpdateGrids(sub, gridNode);
      sub->vtkType = VTK_MULTIBLOCK_DATA_SET;
      // If grid is a collection, is it a temporal one?
      if(xGrid->GetCollectionType() == XDMF_GRID_COLLECTION_TEMPORAL)
        {
        vtkDebugMacro(" Grid is a TemporalDataSetCollection ");
        sub->isTemporal = 1;
        //
        // we will not output a TemporalDataset, but use this as a flag for later
        // sub->vtkType = VTK_TEMPORAL_DATA_SET;
        // THIS HAS CHANGED
        // now there is an isTemporal and isCollection flag.
        // Spatial collections and trees are multiblock, temporal 
        // collections use the output type of the child
        sub->vtkType = sub->GetChild(0)->vtkType;
        }
      }
    else
      {
      // It's a Uniform Grid
      sub = this->Internals->AddGrid(parent, gridName);
      if ( sub )
        {
        sub->XMGrid = xGrid;
        if ( !sub->XMGrid )
          {
            vtkErrorMacro("sub->XMGrid is 0");
            break;
          }
        if (sub->XMGrid->GetTopology()->GetClass() == XDMF_UNSTRUCTURED ) 
          {
          sub->vtkType = VTK_UNSTRUCTURED_GRID;
          } 
        else if (
          sub->XMGrid->GetTopology()->GetTopologyType() == XDMF_2DSMESH ||
          sub->XMGrid->GetTopology()->GetTopologyType() == XDMF_3DSMESH )
          { 
          sub->vtkType = VTK_STRUCTURED_GRID; 
          }
        else if (
          sub->XMGrid->GetTopology()->GetTopologyType() == XDMF_2DCORECTMESH ||
          sub->XMGrid->GetTopology()->GetTopologyType() == XDMF_3DCORECTMESH )
          {
#ifdef USE_IMAGE_DATA
          sub->vtkType = VTK_IMAGE_DATA;
#else
          sub->vtkType = VTK_UNIFORM_GRID;
#endif
          }
        else if (
          sub->XMGrid->GetTopology()->GetTopologyType() == XDMF_2DRECTMESH ||
          sub->XMGrid->GetTopology()->GetTopologyType() == XDMF_3DRECTMESH )
          {
          sub->vtkType = VTK_RECTILINEAR_GRID;
          }
        }
      }
    if( sub && (!parentIsDomain)) {
        sub->Enabled = 1;
    }
    str.rdbuf()->freeze(0);

    // Find the subgrids grids under this one
    gridNode = this->DOM->FindNextElement("Grid", gridNode);
    }
}

//----------------------------------------------------------------------------
void vtkXdmfReader::SetStride(int x, int y, int z)
{
  if ( x <= 0 || y <= 0 || z <= 0 )
    {
    vtkErrorMacro("Strides have to be greater than 0.");
    return;
    }
  vtkDebugMacro("Setting Stride to (" << x << "," << y << "," << z << ")");
  if ((this->Stride[0] != x)||(this->Stride[1] != y)||(this->Stride[2] != z))
    {
    this->Stride[0] = x;
    this->Stride[1] = y;
    this->Stride[2] = z;
    this->Modified();
    }
}

//----------------------------------------------------------------------------
const char * vtkXdmfReader::GetXdmfDOMHandle() 
{
  return( XdmfObjectToHandle( this->DOM ) );
}

//----------------------------------------------------------------------------
int vtkXdmfReader::CanReadFile(const char* fname)
{
  vtkXdmfReaderTester* tester = vtkXdmfReaderTester::New();
  tester->SetFileName(fname);
  int res = tester->TestReadFile();
  tester->Delete();
  return res;
}

//----------------------------------------------------------------------------

int vtkXdmfReader::FillOutputPortInformation(int,
                                             vtkInformation *info)
{ 
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataObject");
  return 1;
}

//----------------------------------------------------------------------------
void vtkXdmfReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "CellDataArraySelection: " << this->CellDataArraySelection 
     << endl;
  os << indent << "PointDataArraySelection: " << this->PointDataArraySelection 
     << endl;
  os << indent << "Domain: " << this->DomainName << endl;
  int cc;
  os << indent << "Grids:" << endl;
  for ( cc = 0; cc < this->GetNumberOfGrids(); ++ cc )
    {
    os << indent << " " << cc << ": " << this->GetGridName(cc) << " - "
       << (this->GetGridSetting(cc)?"enabled":"disabled") << endl;
    }
}


/////////////////////////////////////////////////////////////////////////////
//----------------------------------------------------------------------------
int vtkXdmfReader::ProcessRequest(vtkInformation* request,
                                  vtkInformationVector** inputVector,
                                  vtkInformationVector* outputVector)
{
  // create the output
  if(request->Has(vtkDemandDrivenPipeline::REQUEST_DATA_OBJECT()))
    {
    return this->RequestDataObject(outputVector);
    }
  return this->Superclass::ProcessRequest(request, inputVector, outputVector);
}

//----------------------------------------------------------------------------
bool vtkXdmfReader::ParseXML()
{
  cout << "Begin Parsing" << endl;
  // * Ensure that the required objects have been instantiated.
  if (!this->DOM)
    {
    this->DOM = new XdmfDOM();
    }
  if ( !this->Internals->DataItem )
    {
    this->Internals->DataItem = new XdmfDataItem();
    this->Internals->DataItem->SetDOM(this->DOM);
    }

  // * Check if the XML file needs to be re-read.
  bool modified = true;
  if (this->GetReadFromInputString())
    {
    const char* data=0;
    unsigned int data_length=0;
    if (this->InputArray)
      {
      vtkDebugMacro(<< "Reading from InputArray");
      data = this->InputArray->GetPointer(0);
      data_length = static_cast<unsigned int>(
        this->InputArray->GetNumberOfTuples()*
        this->InputArray->GetNumberOfComponents());
      }
    else if (this->InputString)
      {
      data = this->InputString;
      data_length = this->InputStringLength;
      }
    else
      {
      vtkErrorMacro("No input string specified");
      return false;
      }
    if (data && this->Internals->Data &&
      this->Internals->InputString &&
      data_length == this->Internals->InputStringLength &&
      STRNCASECMP(data, this->Internals->InputString, data_length) == 0)
      {
      modified = false;
      vtkDebugMacro("Input Text Unchanged ... skipping re-parse()");
      }
    else
      {
      vtkDebugMacro("Input Text Changed ...  re-parseing");
      delete this->Internals->Data;
      this->Internals->Data = 0;
      delete [] this->Internals->InputString;
      this->Internals->InputString = new char[data_length+1];
      memcpy(this->Internals->InputString, data, data_length);
      this->Internals->InputString[data_length] = 0; 
      this->Internals->InputStringLength = data_length;

      this->DOM->SetInputFileName(0);
      this->DOM->Parse(this->Internals->InputString);
      }
    }
  else
    {
    // Parse the file...
    if (!this->FileName )
      {
      vtkErrorMacro("File name not set");
      return false;
      }

    // First make sure the file exists.  This prevents an empty file
    // from being created on older compilers.
    if (!vtksys::SystemTools::FileExists(this->FileName))
      {
      vtkErrorMacro("Error opening file " << this->FileName);
      return false;
      }


    if (this->DOM->GetInputFileName() &&
      STRCASECMP(this->DOM->GetInputFileName(), this->FileName) == 0)
      {
      vtkDebugMacro("Filename Unchanged ... skipping re-parse()");
      modified = false;
      }
    else
      {
      vtkDebugMacro("Parsing file: " << this->FileName);

      //Tell the parser what the working directory is.
      vtkstd::string directory =
        vtksys::SystemTools::GetFilenamePath(this->FileName) + "/";
      if (directory == "")
        {
        directory = vtksys::SystemTools::GetCurrentWorkingDirectory() + "/";
        }
      //    directory = vtksys::SystemTools::ConvertToOutputPath(directory.c_str());
      this->DOM->SetWorkingDirectory(directory.c_str());
      this->DOM->SetInputFileName(this->FileName);
      this->DOM->Parse(this->FileName);
      }
    }
  cout << "Done Parsing" << endl;

  if (modified)
    {
    cout << "ParseXML" << endl;
    // Since the DOM was re-parsed we need to update the cache for domains
    // and re-populate the grids.
    this->UpdateDomains();
    this->UpdateRootGrid();
    }
  return true;
}

//----------------------------------------------------------------------------
int vtkXdmfReader::RequestDataObject(vtkInformationVector *outputVector)
{
  cout << "RequestDataObject: " << endl;
  if (!this->ParseXML())
    {
    return 0;
    }

  vtkDebugMacro("My output is a "
    << vtkDataObjectTypes::GetClassNameFromTypeId(this->OutputVTKType));

  //Look at the in memory structures and create an empty vtkDataObject of the
  //proper type for RequestData to fill in later, if needed.

  vtkDataObject *output= vtkDataObject::GetData(outputVector, 0);
  if (!output || output->GetDataObjectType() != this->OutputVTKType)
    {
    vtkInformation *outInfo=outputVector->GetInformationObject(0);
    output = vtkDataObjectTypes::NewDataObject(this->OutputVTKType);
    output->SetPipelineInformation(outInfo);
    outInfo->Set(vtkDataObject::DATA_EXTENT_TYPE(), output->GetExtentType());
    outInfo->Set(vtkDataObject::DATA_OBJECT(), output);
    output->Delete();
    }

//  this->GetPointDataArraySelection()->RemoveAllArrays();
//  this->GetCellDataArraySelection()->RemoveAllArrays();  
  return 1;
}

//-----------------------------------------------------------------------------
int vtkXdmfReaderInternal::FindParallelism(vtkXdmfReaderGrid *grid)
{
    int status = 0;
    int topLevel = 0;

    // Look for a level with enough parallelism
    // If none is found use the largest level (most children)
    if(!grid){
        this->ParallelLevels.clear();
        this->LargestLevel = 0;
        this->mostChildren = 0;
        grid = this->Data;
        topLevel = 1;
    }
    vtkstd::vector<vtkXdmfReaderGrid *>::iterator it;
        if((grid->Enabled) && (!grid->isTemporal) && (grid->Children.size() >= this->mostChildren)){
            this->mostChildren = grid->Children.size();
            this->LargestLevel = grid;
        }
    if((grid->Enabled) && (!grid->isTemporal) && (grid->Children.size() >= this->UpdateNumPieces)){
        // Found Enough Paralelism: Stop looking
        // this->ParallelLevel = grid;
        this->ParallelLevels.push_back(grid);
        grid->isParallel = 1;
        return 1;
    }
    // Process Children
    for ( it = grid->Children.begin();
        it != grid->Children.end();
        ++it )
    {
        {
        status = this->FindParallelism((*it));
        // if(status) return 1;
        }
    }
    if(topLevel && (this->ParallelLevels.size() == 0)){
        // Didn't find enough parallelism
        // Use largest
        if(this->LargestLevel){
            // this->ParallelLevel = this->LargestLevel;
            // this->ParallelLevels.push_back(this->LargestLevel);
            this->LargestLevel->isParallel = 1;
        }else{
            // cout << "No Parallelism, No LargestLevel, Are  Grids enabled ?" << endl;
        }
    }
    return 0;
}

void vtkXdmfReaderInternalUpdateArraysInternal(vtkXdmfReaderGrid* grid,
  vtkstd::set<vtkstd::string>& pointArrays, 
  vtkstd::set<vtkstd::string>& cellArrays)
{
  if (!grid)
    {
    return;
    }

  if (grid->Children.size()>0) 
    {
    //recurse inside composite datasets
    vtkstd::vector<vtkXdmfReaderGrid *>::iterator it;
    for (it = grid->Children.begin(); it != grid->Children.end(); ++it)
      {
      vtkXdmfReaderInternalUpdateArraysInternal((*it), pointArrays, cellArrays);
      }
    return;
    }

  XdmfGrid* xdmfGrid = grid->XMGrid;
  for (int kk=0; kk < xdmfGrid->GetNumberOfAttributes(); kk++)
    {
    XdmfAttribute *Attribute = xdmfGrid->GetAttribute(kk);
    const char *name = Attribute->GetName();
    if (name)
      {
      XdmfInt32 AttributeCenter = Attribute->GetAttributeCenter();
      if (AttributeCenter == XDMF_ATTRIBUTE_CENTER_GRID || 
        AttributeCenter == XDMF_ATTRIBUTE_CENTER_NODE)
        {
        pointArrays.insert(name);
        }
      else
        {
        cellArrays.insert(name);
        }
      }
    }
}
//-----------------------------------------------------------------------------
// Get list of arrays irrespective of their enable/disable states.
int vtkXdmfReaderInternal::UpdateArrays(vtkXdmfReaderGrid* grid)
{
  vtkDataArraySelection* pointDataArraySelection = 
    this->Reader->GetPointDataArraySelection();
  vtkDataArraySelection* cellDataArraySelection = 
    this->Reader->GetCellDataArraySelection();

  pointDataArraySelection->RemoveAllArrays();
  cellDataArraySelection->RemoveAllArrays();

  vtkstd::set<vtkstd::string> pointArrays;
  vtkstd::set<vtkstd::string> cellArrays;

  vtkXdmfReaderInternalUpdateArraysInternal(grid, pointArrays, cellArrays);
  vtkstd::set<vtkstd::string>::iterator iter;
  for (iter = pointArrays.begin(); iter != pointArrays.end();++iter)
    {
    pointDataArraySelection->AddArray(iter->c_str());
    }
  for (iter = cellArrays.begin(); iter != cellArrays.end();++iter)
    {
    cellDataArraySelection->AddArray(iter->c_str());
    }

  return 1;
}

//-----------------------------------------------------------------------------
int vtkXdmfReaderInternal::RequestGridInformation(
  vtkXdmfReaderGrid *grid, vtkInformation *destInfo)
{
  if (grid->Children.size()>0) 
    {
    //recurse inside composite datasets
    vtkDebugWithObjectMacro(this->Reader,
                            "Requesting info in collection or tree " 
                            << grid->Name.c_str());
    vtkstd::vector<vtkXdmfReaderGrid *>::iterator it;
    for ( it = grid->Children.begin();
          it != grid->Children.end();
          ++it )
      {
      //if ((*it)->Enabled ) //only top level grids can be disabled
        {
        this->RequestGridInformation((*it), destInfo);
        }
      }
    return 1;
    }
  vtkInformation *gridInfo = grid->GetInformation();

  XdmfInt32    Rank;
  XdmfInt64    Dimensions[ XDMF_MAX_DIMENSION ];
  XdmfInt64    EndExtent[ XDMF_MAX_DIMENSION ];
  //vtkDataArraySelection* pointDataArraySelection = 
  //  this->Reader->GetPointDataArraySelection();
  //vtkDataArraySelection* cellDataArraySelection = 
  //  this->Reader->GetCellDataArraySelection();
  int *readerStride = this->Reader->GetStride();
  
  XdmfGrid* xdmfGrid = grid->XMGrid;
  
  int kk;
  for( kk = 0 ; kk < xdmfGrid->GetNumberOfAttributes() ; kk++ )
    {
    /*
    XdmfAttribute       *Attribute;
    Attribute = xdmfGrid->GetAttribute( kk );
    const char *name = Attribute->GetName();
    if (name )
      {
      XdmfInt32 AttributeCenter = Attribute->GetAttributeCenter();
      if ( AttributeCenter == XDMF_ATTRIBUTE_CENTER_GRID || 
           AttributeCenter == XDMF_ATTRIBUTE_CENTER_NODE)
        {
        if ( !pointDataArraySelection->ArrayExists(name) )
          {
          pointDataArraySelection->AddArray(name);
          cout << "Add Array: " << name << endl;
          }
        }
      else
        {
        if ( !cellDataArraySelection->ArrayExists(name) )
          {
          cellDataArraySelection->AddArray(name);
          cout << "Add Cell Array: " << name << endl;
          }
        }
      }
      */
    }
  
  // Revised Initial Setup
  grid->DataDescription = xdmfGrid->GetTopology()->GetShapeDesc();
  Rank = grid->DataDescription->GetShape( Dimensions );
  int i;
  for(i = Rank ; i < XDMF_MAX_DIMENSION ; i++)
    {
    Dimensions[i] = 1;
    }
  // End Extent is Dim - 1
  EndExtent[0] = vtkMAX(0, Dimensions[0] - 1);
  EndExtent[1] = vtkMAX(0, Dimensions[1] - 1);
  EndExtent[2] = vtkMAX(0, Dimensions[2] - 1);
  // vtk Dims are i,j,k XDMF are k,j,i
  EndExtent[0] = vtkMAX(0, EndExtent[0]) / readerStride[2];
  EndExtent[1] = vtkMAX(0, EndExtent[1]) / readerStride[1];
  EndExtent[2] = vtkMAX(0, EndExtent[2]) / readerStride[0];
  vtkDebugWithObjectMacro(this->Reader, "EndExtents = " 
                          << (vtkIdType)EndExtent[0] << ", " 
                          << (vtkIdType)EndExtent[1] << ", " 
                          << (vtkIdType)EndExtent[2]);
  
  gridInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),
               0, EndExtent[2], 0, EndExtent[1], 0, EndExtent[0]);
  destInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),
               0, EndExtent[2], 0, EndExtent[1], 0, EndExtent[0]);
  
  vtkDebugWithObjectMacro(this->Reader, 
                          "Grid Type = " 
                          << xdmfGrid->GetTopology()->GetTopologyTypeAsString()
                          << " = " 
                          << xdmfGrid->GetTopology()->GetTopologyType());
  if( xdmfGrid->GetTopology()->GetClass() != XDMF_UNSTRUCTURED ) 
    {
    if( (xdmfGrid->GetTopology()->GetTopologyType() == XDMF_2DSMESH ) ||
        (xdmfGrid->GetTopology()->GetTopologyType() == XDMF_3DSMESH ) )
      {
      vtkDebugWithObjectMacro(this->Reader, 
                              "Setting Extents for vtkStructuredGrid");
      }
    else if ( xdmfGrid->GetTopology()->GetTopologyType() == XDMF_2DCORECTMESH|| 
              xdmfGrid->GetTopology()->GetTopologyType() == XDMF_3DCORECTMESH )
      {
      XdmfGeometry  *Geometry = xdmfGrid->GetGeometry();
      if ( Geometry->GetGeometryType() == XDMF_GEOMETRY_ORIGIN_DXDYDZ )
        { 
        // Update geometry so that origin and spacing are read
        Geometry->Update();
        XdmfFloat64 *origin = Geometry->GetOrigin();
        XdmfFloat64 *spacing = Geometry->GetDxDyDz();
        
        gridInfo->Set(vtkDataObject::ORIGIN(), 
                     origin[2], origin[1], origin[0]);
        gridInfo->Set(vtkDataObject::SPACING(), 
                     spacing[2], spacing[1], spacing[0]);
        destInfo->Set(vtkDataObject::ORIGIN(), 
                     origin[2], origin[1], origin[0]);
        destInfo->Set(vtkDataObject::SPACING(), 
                     spacing[2], spacing[1], spacing[0]);
        }
      }
    else  if ( xdmfGrid->GetTopology()->GetTopologyType() == XDMF_2DRECTMESH ||
               xdmfGrid->GetTopology()->GetTopologyType() == XDMF_3DRECTMESH )
      {
      vtkDebugWithObjectMacro(this->Reader, 
                              "Setting Extents for vtkRectilinearGrid");
      }
    else 
      {
      vtkErrorWithObjectMacro(this->Reader,"Unknown topology type: " 
                              << xdmfGrid->GetTopology()->GetTopologyType());
      }
    }

  return 1;
}
//-----------------------------------------------------------------------------
void vtkXdmfReader::AssignTimeIndex(vtkXdmfReaderGrid *ptr)
{
  if (!ptr)
    {
    return;
    }
    int actualTimeStep = vtkstd::find_if(
          this->Internals->TimeValues.begin(), 
          this->Internals->TimeValues.end(), 
          vtkstd::bind2nd( vtkstd::greater_equal<double>( ), ptr->Time)) 
          - this->Internals->TimeValues.begin();
    ptr->TimeIndex = actualTimeStep;
  int T = ptr->Children.size();
  for (int i=0; i<T; i++) {
    vtkXdmfReaderGrid *child = ptr->GetChild(i);
    this->AssignTimeIndex(child);
  }
}

void vtkXdmfReader::FindAllTimeValues(vtkXdmfReaderGrid *ptr)
{
    XdmfTime *time = 0;
  if (!ptr)
    {
    return;
    }
    if(ptr->XMGrid){
        time = ptr->XMGrid->GetTime();
    }
    if(time && (time->GetTimeType() == XDMF_TIME_UNSET)){
    }
    if(time && (time->GetTimeType() != XDMF_TIME_UNSET)){
        if(!ptr->isCollection) ptr->isTemporal = 1;
        this->OutputTemporal = 1;
        ptr->Time = time->GetValue();
        this->Internals->TimeValues.push_back(ptr->Time);
    }
  int T = ptr->Children.size();
  for (int i=0; i<T; i++) {
    vtkXdmfReaderGrid *child = ptr->GetChild(i);
    this->FindAllTimeValues(child);
  }
}
void vtkXdmfReader::FindTimeValues()
{
  vtkXdmfReaderGrid *ptr = this->Internals->Data;
  if (!ptr)
    {
    return;
    }
  this->Internals->TimeValues.clear();
  this->FindAllTimeValues(ptr);
  vtkstd::sort(this->Internals->TimeValues.begin(), this->Internals->TimeValues.end());
  vtkstd::vector<XdmfFloat64>::iterator  new_end = vtkstd::unique(this->Internals->TimeValues.begin(), this->Internals->TimeValues.end());
  this->Internals->TimeValues.erase(new_end, this->Internals->TimeValues.end());
  int T = this->Internals->TimeValues.size();
  for (int i=0; i<T; i++) {
  }
  this->AssignTimeIndex(this->Internals->Data);
  this->TimeStepRange[0] = 0;
  this->TimeStepRange[1] = this->Internals->TimeValues.size()-1;
}
//-----------------------------------------------------------------------------
int vtkXdmfReader::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  vtkDebugMacro("RequestInformation");
  
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  //say we can produce as many pieces are are desired
  outInfo->Set(
    vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(),-1);
  vtkXdmfReaderGrid *ptr = this->Internals->Data;
  if (!ptr)
    {
    return 1;
    }

  switch (this->OutputVTKType)
    {
    case VTK_POLY_DATA:
    case VTK_UNSTRUCTURED_GRID:
    case VTK_IMAGE_DATA:
    case VTK_UNIFORM_GRID:
    case VTK_RECTILINEAR_GRID:
    case VTK_STRUCTURED_GRID:
    {
    vtkXdmfReaderGrid *sptr = this->Internals->GetGrid(0);
//    sptr->Information = outInfo; // not sure if we need to do this now
    this->Internals->RequestGridInformation(sptr, outInfo);
    // release reference to meta info, we do not use it again so it is save for the vtkXdmfReaderGrid to get rid of it
    // NO! - don't release it here, otherwise the null pointer is checked later and a new information object
    // gets created - and all the fields are reset.
//    sptr->Information = 0;
    // and maybe we can now release the hold of it.
    }
    break;
    case VTK_MULTIBLOCK_DATA_SET:
    {
    vtkstd::vector<vtkXdmfReaderGrid *>::iterator it;
    for ( it = ptr->Children.begin();
          it != ptr->Children.end();
          ++it )
      {
      if ( (*it)->Enabled )
        {
        //composite data does not yet support meta info
        //so we do not yet put the results into outInfo
        this->Internals->RequestGridInformation((*it), outInfo);
        }
      }
    }
    break;
    default:
      vtkErrorMacro("Unknown type to generate information for");
      return 0;
    }

  this->ActualTimeStep = this->TimeStep;

  this->FindTimeValues();
  if (this->OutputTemporal) {
    if (this->Internals->TimeValues.size()>0) {
      outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), 
                    &this->Internals->TimeValues[0], 
                     this->Internals->TimeValues.size());
      double timeRange[2];
      timeRange[0] = this->Internals->TimeValues.front();
      timeRange[1] = this->Internals->TimeValues.back();
      outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);
      // Check which particular time was requested.
      if (!outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS())) 
      {
        double requestedTimeValue = this->Internals->TimeValues[this->TimeStep];
        // this should be the same, just checking for debug purposes
        this->ActualTimeStep = vtkstd::find_if(
          this->Internals->TimeValues.begin(), 
          this->Internals->TimeValues.end(), 
          vtkstd::bind2nd( vtkstd::greater_equal<double>( ), requestedTimeValue )) 
          - this->Internals->TimeValues.begin();
        vtkDebugMacro( << "Time step requested " << this->ActualTimeStep )
      } 
    }
  }

  return 1;
}

//----------------------------------------------------------------------------
int vtkXdmfReaderInternal::RequestGridData(
/*  const char* currentGridName,*/
  vtkXdmfReaderGrid *grid,
  vtkDataObject *output,
    int timeIndex,
  int isSubBlock,
  double progressS,
  double progressE)
{
  double totalProgFrac = progressE-progressS;
    vtkstd::vector<vtkXdmfReaderGrid *>::iterator it;

  // Handle Tempral Collection
  if(grid->isTemporal && grid->isCollection){
    for ( it = grid->Children.begin();
          it != grid->Children.end();
          ++it )
      {
          if((*it)->TimeIndex == timeIndex){
                int hasUpdateExtent, childhasUpdateExtent;
                int subBlock = isSubBlock;
                hasUpdateExtent = grid->GetInformation()->Has(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT());
                childhasUpdateExtent = (*it)->GetInformation()->Has(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT());
                if( hasUpdateExtent && (!childhasUpdateExtent)){
                    // cout << "Copying Update Extents" << endl;
                        (*it)->GetInformation()->CopyEntry(grid->GetInformation(), 
                            vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT());
                        (*it)->GetInformation()->CopyEntry(grid->GetInformation(), 
                            vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
                        (*it)->GetInformation()->CopyEntry(grid->GetInformation(), 
                            vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
                        (*it)->GetInformation()->CopyEntry(grid->GetInformation(), 
                            vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS());
                        subBlock = 0;
                }
                return(this->RequestGridData(
                (*it),
                output,
                timeIndex,
                subBlock,
                progressS,
                progressE
                ));
          }
      }

  }
  if((!grid->Enabled) || ((grid->isTemporal == 1) && (!grid->isCollection) && (timeIndex != grid->TimeIndex))){
      return 1;
  }
  if (grid->vtkType == VTK_MULTIBLOCK_DATA_SET)
    {
    //recurse inside composite datasets
    vtkDebugWithObjectMacro(this->Reader,
                            "Requesting data in collection or tree ");
    vtkMultiBlockDataSet *outMB = 
      vtkMultiBlockDataSet::SafeDownCast(output);
    if (!outMB)
      {
      vtkErrorWithObjectMacro(this->Reader,
                              "Wrong vtk output structure expected composite");
      return 0;
      }
    unsigned int outputGrid = 0;
    unsigned int IsParallel = 0;
    int nChildren = grid->Children.size();
    for ( it = grid->Children.begin();
          it != grid->Children.end();
          ++it )
      {
      double lprogressS = progressS 
        + totalProgFrac*(double)outputGrid/(double)nChildren;
      double lprogressE = progressS 
        + totalProgFrac*(double)(outputGrid+1)/(double)nChildren;

      IsParallel = grid->isParallel;
      vtkXdmfReaderGrid *child = (*it);
        // if((this->ParallelLevel != grid) ||
        //     ((this->ParallelLevel == grid) && ((outputGrid % this->UpdateNumPieces) == this->UpdatePiece)))
        // {
        if(!IsParallel ||
            (IsParallel && ((outputGrid % this->UpdateNumPieces) == this->UpdatePiece)))
        {
        vtkDataObject *soutput=
           vtkDataObjectTypes::NewDataObject(child->vtkType);
        if (soutput && soutput->IsA("vtkMultiBlockDataSet"))
          {
          outMB->SetBlock(outputGrid, vtkMultiBlockDataSet::SafeDownCast(soutput));
          }
        else
          {
          outMB->SetBlock(outputGrid, vtkDataSet::SafeDownCast(soutput));
          //outMB->SetMetaData(outputGrid, grid->GetInformation());
          }
        vtkDebugWithObjectMacro(
          this->Reader,
          "Recursively filling in ds " << outputGrid 
          << " a " << soutput->GetClassName() 
          << " from " << (*it)->Name.c_str());
        this->RequestGridData(
          /*(*it)->Name.c_str(),*/
          child,
          soutput,
          timeIndex,
          1,
          lprogressS,
          lprogressE
          );
        soutput->Delete();
        } else {
          outMB->SetBlock(outputGrid, static_cast<vtkDataSet*>(NULL));
          //outMB->SetMetaData(outputGrid, grid->GetInformation());
          }
        outputGrid++;
        }
    if (totalProgFrac>0.1)
      {
      this->Reader->UpdateProgress(progressE);
      }
    return 1;
    }

  vtkInformation* gridInfo = grid->GetInformation();
  int *readerStride = this->Reader->GetStride();
  
  vtkDataArraySelection* pointDataArraySelection = 
    this->Reader->GetPointDataArraySelection();
  vtkDataArraySelection* cellDataArraySelection = 
    this->Reader->GetCellDataArraySelection();
  
  // Handle single grid
  XdmfGrid* xdmfGrid = grid->XMGrid;
  XdmfDOM* xdmfDOM = xdmfGrid->GetDOM();
  
  if ( !grid->DataDescription )
    {
    grid->DataDescription = xdmfGrid->GetTopology()->GetShapeDesc();
    //continue;
    }
  
    // cout << "RequestGridData Node " << this->Reader->GetController()->GetLocalProcessId() << " Reading data for " << xdmfGrid->GetName() << endl;
  vtkDebugWithObjectMacro(this->Reader, 
                          "Reading Heavy Data for " << xdmfGrid->GetName());
#ifndef XDMF_NO_MPI
  xdmfGrid->SetDsmBuffer(this->DsmBuffer);
#endif
  xdmfGrid->Update();
  
  // True for all 3d datasets except unstructured grids
  int globalrank = 3;
  switch(xdmfGrid->GetTopology()->GetTopologyType())
    {
    case XDMF_2DSMESH: case XDMF_2DRECTMESH: case XDMF_2DCORECTMESH:
      globalrank = 2;
    }
  if ( xdmfGrid->GetTopology()->GetClass() == XDMF_UNSTRUCTURED )
    {
    globalrank = 1;
    }
  
  
  vtkIdType cc;
  XdmfXmlNode attrNode;
  XdmfXmlNode dataNode;
  XdmfInt64 start[4]  = { 0, 0, 0, 0 };
  XdmfInt64 stride[4] = { 1, 1, 1, 1 };
  XdmfInt64 count[4] = { 0, 0, 0, 0 };
  grid->DataDescription->GetShape(count);

  int upext[6];
  int whext[6];
  
  if( xdmfGrid->GetTopology()->GetClass() != XDMF_UNSTRUCTURED)
    {
    gridInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), whext);

    if(isSubBlock)
      {
      // the composite pipeline does not set update extent, so just take
      // the whole extent as update extent.
      gridInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), upext);
      }
    else if (gridInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT()))
      {
      gridInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), upext);
      }
    else 
      {
      memcpy(upext, whext, sizeof(int)*6);
      }
        
   /*
     cerr << "XDR(" << this << ") RD " 
     << "WholeExt" <<
     << whext[0] << ", " << whext[1] << " "
     << whext[2] << ", " << whext[3] << " "
     << whext[4] << ", " << whext[5] << endl;
     cerr << "XDR(" << this << ") RD " 
     << "UpdateExt" <<
     << upext[0] << ", " << upext[1] << " "
     << upext[2] << ", " << upext[3] << " "
     << upext[4] << ", " << upext[5] << endl;
   */
    start[2] = vtkMAX(0, upext[0]);
    start[1] = vtkMAX(0, upext[2]);
    start[0] = vtkMAX(0, upext[4]);
    
    count[2] = vtkMAX(1, upext[1] - upext[0]);
    count[1] = vtkMAX(1, upext[3] - upext[2]);
    count[0] = vtkMAX(1, upext[5] - upext[4]);
    }
  
  XdmfGeometry  *Geometry = xdmfGrid->GetGeometry();
  
  
  // Read Topology for Unstructured Grid
  if( xdmfGrid->GetTopology()->GetClass() == XDMF_UNSTRUCTURED ) 
    {
    vtkUnstructuredGrid  *vGrid = static_cast<vtkUnstructuredGrid *>(output);
    vtkCellArray                *verts;
    XdmfInt32           vType;
    XdmfInt32           NodesPerElement;
    vtkIdType           NumberOfElements;
    vtkIdType           i, j, index;    
    XdmfInt64           Length, *Connections;
    vtkIdType           *connections;
    int                 *cell_types, *ctp;
    
    vtkDebugWithObjectMacro(
      this->Reader,
      "Unstructured Topology is " 
      << xdmfGrid->GetTopology()->GetTopologyTypeAsString());
    switch ( xdmfGrid->GetTopology()->GetTopologyType() )
      {
      case  XDMF_POLYVERTEX :
        vType = VTK_POLY_VERTEX;
        break;
      case  XDMF_POLYLINE :
        vType = VTK_POLY_LINE;
        break;
      case  XDMF_POLYGON :
        vType = VTK_POLYGON;
        break;
      case  XDMF_TRI :
        vType = VTK_TRIANGLE;
        break;
      case  XDMF_QUAD :
        vType = VTK_QUAD;
        break;
      case  XDMF_TET :
        vType = VTK_TETRA;
        break;
      case  XDMF_PYRAMID :
        vType = VTK_PYRAMID;
        break;
      case  XDMF_WEDGE :
        vType = VTK_WEDGE;
        break;
      case  XDMF_HEX :
        vType = VTK_HEXAHEDRON;
        break;
      case  XDMF_EDGE_3 :
        vType = VTK_QUADRATIC_EDGE ;
        break;
      case  XDMF_TRI_6 :
        vType = VTK_QUADRATIC_TRIANGLE ;
        break;
      case  XDMF_QUAD_8 :
        vType = VTK_QUADRATIC_QUAD ;
        break;
      case  XDMF_TET_10 :
        vType = VTK_QUADRATIC_TETRA ;
        break;
      case  XDMF_PYRAMID_13 :
        vType = VTK_QUADRATIC_PYRAMID ;
        break;
      case  XDMF_WEDGE_15 :
        vType = VTK_QUADRATIC_WEDGE ;
        break;
      case  XDMF_HEX_20 :
        vType = VTK_QUADRATIC_HEXAHEDRON ;
        break;
      case XDMF_MIXED :
        vType = -1;
        break;
      default :
        XdmfErrorMessage("Unknown Topology Type = " 
                         << xdmfGrid->GetTopology()->GetTopologyType());
        return 1;
      }
    if( xdmfGrid->GetTopology()->GetTopologyType() != XDMF_MIXED)
      {
      NodesPerElement = xdmfGrid->GetTopology()->GetNodesPerElement();
      if ( xdmfGrid->GetTopology()->GetConnectivity()->GetRank() == 2 )
        {
        NodesPerElement = 
          xdmfGrid->GetTopology()->GetConnectivity()->GetDimension(1);
        }
    
      /* Create Cell Type Array */
      Length = xdmfGrid->GetTopology()->GetConnectivity()->GetNumberOfElements();
      Connections = new XdmfInt64[ Length ];
      xdmfGrid->GetTopology()->GetConnectivity()->GetValues(
        0, 
        Connections, 
        Length);
      
      NumberOfElements = 
        xdmfGrid->GetTopology()->GetShapeDesc()->GetNumberOfElements();
      ctp = cell_types = new int[ NumberOfElements ];
      
      /* Create Cell Array */
      verts = vtkCellArray::New();
      
      /* Get the pointer */
      connections = verts->WritePointer(
        NumberOfElements,
        NumberOfElements * ( 1 + NodesPerElement ));
      
      /* Connections : N p1 p2 ... pN */
      /* i.e. Triangles : 3 0 1 2    3 3 4 5   3 6 7 8 */
      index = 0;
      for( j = 0 ; j < NumberOfElements; j++ )
        {
        *ctp++ = vType;
        *connections++ = NodesPerElement;
        for( i = 0 ; i < NodesPerElement; i++ )
          {
          *connections++ = Connections[index++];
          }
        }
      } 
    else 
      {
      // Mixed Topology
      /* Create Cell Type Array */
      vtkIdTypeArray *IdArray;
      vtkIdType RealSize;
      
      Length = xdmfGrid->GetTopology()->GetConnectivity()->GetNumberOfElements();
      Connections = new XdmfInt64[ Length ];
      xdmfGrid->GetTopology()->GetConnectivity()->GetValues(
        0, 
        Connections, 
        Length );
      NumberOfElements = 
        xdmfGrid->GetTopology()->GetShapeDesc()->GetNumberOfElements();
      ctp = cell_types = new int[ NumberOfElements ];
      
      /* Create Cell Array */
      verts = vtkCellArray::New();
      
      /* Get the pointer. Make it Big enough ... too big for now */
      connections = verts->WritePointer(
        NumberOfElements,
        Length);
      //   Length);
      /* Connections : N p1 p2 ... pN */
      /* i.e. Triangles : 3 0 1 2    3 3 4 5   3 6 7 8 */
      index = 0;
      int sub = 0;
      for( j = 0 ; j < NumberOfElements; j++ )
        {
        switch ( Connections[index++] )
          {
          case  XDMF_POLYVERTEX :
            vType = VTK_POLY_VERTEX;
            NodesPerElement = Connections[index++];
            sub++;
            break;
          case  XDMF_POLYLINE :
            vType = VTK_POLY_LINE;
            NodesPerElement = Connections[index++];
            sub++;
            break;
          case  XDMF_POLYGON :
            vType = VTK_POLYGON;
            NodesPerElement = Connections[index++];
            sub++;
            break;
          case  XDMF_TRI :
            vType = VTK_TRIANGLE;
            NodesPerElement = 3;
            break;
          case  XDMF_QUAD :
            vType = VTK_QUAD;
            NodesPerElement = 4;
            break;
          case  XDMF_TET :
            vType = VTK_TETRA;
            NodesPerElement = 4;
            break;
          case  XDMF_PYRAMID :
            vType = VTK_PYRAMID;
            NodesPerElement = 5;
            break;
          case  XDMF_WEDGE :
            vType = VTK_WEDGE;
            NodesPerElement = 6;
            break;
          case  XDMF_HEX :
            vType = VTK_HEXAHEDRON;
            NodesPerElement = 8;
            break;
          case  XDMF_EDGE_3 :
            vType = VTK_QUADRATIC_EDGE ;
            NodesPerElement = 3;
            break;
          case  XDMF_TRI_6 :
            vType = VTK_QUADRATIC_TRIANGLE ;
            NodesPerElement = 6;
            break;
          case  XDMF_QUAD_8 :
            vType = VTK_QUADRATIC_QUAD ;
            NodesPerElement = 8;
            break;
          case  XDMF_TET_10 :
            vType = VTK_QUADRATIC_TETRA ;
            NodesPerElement = 10;
            break;
          case  XDMF_PYRAMID_13 :
            vType = VTK_QUADRATIC_PYRAMID ;
            NodesPerElement = 13;
            break;
          case  XDMF_WEDGE_15 :
            vType = VTK_QUADRATIC_WEDGE ;
            NodesPerElement = 15;
            break;
          case  XDMF_HEX_20 :
            vType = VTK_QUADRATIC_HEXAHEDRON ;
            NodesPerElement = 20;
            break;
          default :
            XdmfErrorMessage("Unknown Topology Type");
            return 1;
          }
        *ctp++ = vType;
        *connections++ = NodesPerElement;
        for( i = 0 ; i < NodesPerElement; i++ )
          {
          *connections++ = Connections[index++];
          }
        }
      // Resize the Array to the Proper Size
      IdArray = verts->GetData();
      RealSize = index - sub;
      vtkDebugWithObjectMacro(this->Reader, 
                              "Resizing to " << RealSize << " elements");
      IdArray->Resize(RealSize);
      }

    delete [] Connections;
    vGrid->SetCells(cell_types, verts);
    /* OK, because of reference counting */
    verts->Delete();
    delete [] cell_types;
    vGrid->Modified();
    }  // if( xdmfGrid->GetClass() == XDMF_UNSTRUCTURED ) 
  else if( xdmfGrid->GetTopology()->GetTopologyType() == XDMF_2DSMESH ||
           xdmfGrid->GetTopology()->GetTopologyType() == XDMF_3DSMESH )
    {
    vtkDebugWithObjectMacro(this->Reader, 
                            "Setting Extents for vtkStructuredGrid");
    vtkStructuredGrid  *vGrid = vtkStructuredGrid::SafeDownCast(output);
    vGrid->SetExtent(upext);    
    } 
  else if (xdmfGrid->GetTopology()->GetTopologyType() == XDMF_2DCORECTMESH ||
           xdmfGrid->GetTopology()->GetTopologyType() == XDMF_3DCORECTMESH ) 
    {
#ifdef USE_IMAGE_DATA
    vtkImageData *idata = vtkImageData::SafeDownCast(output);
#else
    vtkUniformGrid *idata = vtkUniformGrid::SafeDownCast(output);
#endif
    idata->SetExtent(upext);
    vtkDebugWithObjectMacro(this->Reader, <<"Setting Xdmf Extent : " 
      << upext[0] << " " << upext[1] << " " 
      << upext[2] << " " << upext[3] << " " 
      << upext[4] << " " << upext[5]);
    }
  else if ( xdmfGrid->GetTopology()->GetTopologyType() == XDMF_2DRECTMESH ||
            xdmfGrid->GetTopology()->GetTopologyType() == XDMF_3DRECTMESH )
    {
    vtkRectilinearGrid *vGrid = vtkRectilinearGrid::SafeDownCast(output);
    vGrid->SetExtent(upext);    
    }
  else
    {
    vtkErrorWithObjectMacro(this->Reader,"Do not understand topology type: " 
                            << xdmfGrid->GetTopology()->GetTopologyType());
    }
  // Read Geometry
  if( ( Geometry->GetGeometryType() == XDMF_GEOMETRY_X_Y_Z ) ||
      ( Geometry->GetGeometryType() == XDMF_GEOMETRY_XYZ ) ||
      ( Geometry->GetGeometryType() == XDMF_GEOMETRY_X_Y ) ||
      ( Geometry->GetGeometryType() == XDMF_GEOMETRY_XY ) )
    {
    XdmfInt64   Length;
    vtkPoints   *Points;
    vtkPointSet *Pointset = vtkPointSet::SafeDownCast(output);
    
    // Special flag, for structured data
    int structured_grid = 0;
    if ( vtkStructuredGrid::SafeDownCast(output) )
      {
      structured_grid = 1;
      }
    
    Points = Pointset->GetPoints();
    if( !Points )
      {
      vtkDebugWithObjectMacro(this->Reader,
                              "Creating vtkPoints" );
      Points = vtkPoints::New();
      Pointset->SetPoints( Points );
      // OK Because of Reference Counting
      Points->Delete();
      }
    
    if( Geometry->GetPoints() )
      {
      if( Points )
        {
        if ( Geometry->GetPoints()->GetNumberType() == XDMF_FLOAT32_TYPE )
          {
          if ( Points->GetData()->GetDataType() != VTK_FLOAT)
            {
            vtkFloatArray* da = vtkFloatArray::New();
            da->SetNumberOfComponents(3);
            Points->SetData(da);
            da->Delete();
            }
          }
        else // means == XDMF_FLOAT64_TYPE
          {
          if ( Points->GetData()->GetDataType() != VTK_DOUBLE )
            {
            vtkDoubleArray* da = vtkDoubleArray::New();
            da->SetNumberOfComponents(3);
            Points->SetData(da);
            da->Delete();
            }
          }
        
        Length = Geometry->GetPoints()->GetNumberOfElements();
        vtkDebugWithObjectMacro(this->Reader, 
                                "Setting Array of " << (int)Length << " = " 
                                << (int)Geometry->GetNumberOfPoints() 
                                << " Points");
        vtkIdType iskip[3] = { 0, 0, 0 };
        vtkIdType eskip[3] = { 0, 0, 0 };
        int strides_or_extents = 0;
        if ( structured_grid )
          {
          XdmfInt64     ii, jj, kk;
          vtkIdType numpoints = Geometry->GetNumberOfPoints();
          vtkIdType newnumpoints = ((upext[5] - upext[4] + 1) * (upext[3] - upext[2] + 1) * (upext[1] - upext[0] + 1));
          int cnt = 0;
          for (kk = upext[4]; kk <= upext[5]; kk ++ )
            {
            for ( jj = upext[2]; jj <= upext[3]; jj ++ )
              {
              for ( ii = upext[0]; ii <= upext[1]; ii ++ )
                {
                cnt ++;
                }
              }
            }
          newnumpoints = cnt;
          
          Points->SetNumberOfPoints(newnumpoints);
          vtkIdType dims[3];
          dims[0] = whext[1] - whext[0] + 1;
          dims[1] = whext[3] - whext[2] + 1;
          dims[2] = whext[5] - whext[4] + 1;
          iskip[0] = upext[0];
          iskip[1] = upext[2] * dims[0];
          iskip[2] = upext[4] * dims[0] * dims[1];
          eskip[0] = whext[1] - upext[1];
          eskip[1] = (whext[3] - upext[3]) * dims[0];
          eskip[2] = (whext[5] - upext[5]) * dims[0] * dims[1];
          if ( newnumpoints != numpoints )
            {
            strides_or_extents = 1;
            }
          }
        else
          {
          // Unstructured grid
          Points->SetNumberOfPoints( Geometry->GetNumberOfPoints() );
          int kk;
          for ( kk = 0; kk < 6; kk ++ )
            {
            upext[kk] = whext[kk];
            }
          }
        if( Geometry->GetPoints()->GetDataType() == XDMF_FLOAT32_TYPE && 
            !strides_or_extents) 
          {
          Geometry->GetPoints()->GetValues(
            0, 
            vtkFloatArray::SafeDownCast(Points)->GetPointer(0), Length );
          } 
        else if( Geometry->GetPoints()->GetDataType() == XDMF_FLOAT64_TYPE && 
                 !strides_or_extents) 
          {
          Geometry->GetPoints()->GetValues(
            0, 
            vtkDoubleArray::SafeDownCast(Points)->GetPointer(0), Length );
          } 
        else 
          {
          XdmfFloat64   *TmpPp;
          XdmfFloat64   *TmpPoints = new XdmfFloat64[ Length ];
          XdmfInt64     ii, jj, kk;
          
          Geometry->GetPoints()->GetValues( 0, TmpPoints, Length );
          TmpPp = TmpPoints;
          vtkIdType cnt = 0;
          
          if ( strides_or_extents )
            {
            XdmfInt64   Dimensions[3] = { 0, 0, 0 };
            XdmfTopology        *Topology = xdmfGrid->GetTopology();
            Topology->GetShapeDesc()->GetShape( Dimensions );
            cnt = 0;
            for (kk = whext[4]; kk < Dimensions[0]; kk ++ )
              {
              for ( jj = whext[2]; jj < Dimensions[1]; jj ++ )
                {
                for ( ii = whext[0]; ii < Dimensions[2]; ii ++ )
                  {
                  vtkIdType rii = ii / readerStride[0];
                  vtkIdType rjj = jj / readerStride[1];
                  vtkIdType rkk = kk / readerStride[2];
                  vtkIdType mii = ii % readerStride[0];
                  vtkIdType mjj = jj % readerStride[1];
                  vtkIdType mkk = kk % readerStride[2];
                  if ( 
                    rii >= upext[0] && rii <= upext[1] &&
                    rjj >= upext[2] && rjj <= upext[3] &&
                    rkk >= upext[4] && rkk <= upext[5] &&
                    mii == 0 && mjj == 0 && mkk == 0 )
                    {
                    // We are inside the extents
                    Points->SetPoint(cnt, TmpPp);
                    TmpPp += 3;
                    cnt ++;
                    }
                  else
                    {
                    TmpPp += 3;
                    }
                  }
                }
              }
            }
          else
            {
            cnt = 0;
            for( ii = 0 ; ii < Length / 3 ; ii++ )
              {
              Points->SetPoint(cnt, TmpPp);
              TmpPp += 3;
              cnt ++;
              }
            }
          delete [] TmpPoints;
          }
        Points->Modified();
        Pointset->Modified();
        } 
      else 
        {
        XdmfErrorMessage("Base Grid Has No Points");
        return 1;
        }
      } 
    else 
      {
      XdmfErrorMessage("No Points to Set");
      return 1;
      }
    if ( structured_grid )
      {
      stride[2] = readerStride[0];
      stride[1] = readerStride[1];
      stride[0] = readerStride[2];
      }
    }
  else if ( Geometry->GetGeometryType() == XDMF_GEOMETRY_ORIGIN_DXDYDZ )
    {
#ifdef USE_IMAGE_DATA
    vtkImageData *vGrid = vtkImageData::SafeDownCast(output);
#else
    vtkUniformGrid *vGrid = vtkUniformGrid::SafeDownCast(output);
#endif
    XdmfTopology *Topology = xdmfGrid->GetTopology();
    XdmfInt64   Dimensions[3] = { 0, 0, 0 };
    XdmfFloat64  *origin = Geometry->GetOrigin();
    XdmfFloat64 *spacing = Geometry->GetDxDyDz();
    Topology->GetShapeDesc()->GetShape( Dimensions );
    vGrid->SetDimensions(Dimensions[2], Dimensions[1], Dimensions[0]); //whole domain
    vGrid->SetExtent(upext); //subextent of domain to fill
    vGrid->SetOrigin(origin[2], origin[1], origin[0]); //xyz in space
    /*
    cerr << "XDR(" << this << ")::RGD Origin " 
         << origin[0] << "," << origin[1] << "," << origin[2] << endl;
    cerr << "XDR(" << this << ")::RGD Dimensions " 
         << Dimensions[0] << "," << Dimensions[1] << "," << Dimensions[2] << endl;
    cerr << "XDR(" << this << ")::RGD Spacing " 
         << spacing[0] << "," << spacing[1] << "," << spacing[2] << endl;
    cerr << "XDR(" << this << ")::RGD Extent " 
         << upext[0] << "," << upext[1] << " "
         << upext[2] << "," << upext[3] << " "
         << upext[4] << "," << upext[5] << endl;
    */
    stride[2] = readerStride[0];
    stride[1] = readerStride[1];
    stride[0] = readerStride[2];
    vGrid->SetSpacing(
      spacing[2]*readerStride[0], 
      spacing[1]*readerStride[1],
      spacing[0]*readerStride[2]);
    }
  else
    {
    // Special Rectilinear and CoRectilinear Geometries
    XdmfTopology        *Topology = xdmfGrid->GetTopology();
    vtkIdType Index;
    vtkRectilinearGrid *vGrid = vtkRectilinearGrid::SafeDownCast(output);
    if ( vGrid )
      {
      vtkDoubleArray      *XCoord, *YCoord, *ZCoord;
      XdmfFloat64 *Origin;
      XdmfInt64   Dimensions[3] = { 0, 0, 0 };
      
      // Make Sure Grid Has Coordinates
      Topology->GetShapeDesc()->GetShape( Dimensions );
      
      XCoord = vtkDoubleArray::New();
      vGrid->SetXCoordinates( XCoord );
      // OK Because of Reference Counting
      XCoord->Delete();   
      XCoord->SetNumberOfValues( count[2]+1 );
      YCoord = vtkDoubleArray::New();
      vGrid->SetYCoordinates( YCoord );
      // OK Because of Reference Counting
      YCoord->Delete();   
      YCoord->SetNumberOfValues( count[1]+1 );
      ZCoord = vtkDoubleArray::New();
      vGrid->SetZCoordinates( ZCoord );
      // OK Because of Reference Counting
      ZCoord->Delete();   
      ZCoord->SetNumberOfValues( count[0]+1 );

      // Build Vectors if nescessary
      if( Geometry->GetGeometryType() == XDMF_GEOMETRY_ORIGIN_DXDYDZ )
        {
        if( !Geometry->GetVectorX() )
          {
          Geometry->SetVectorX( new XdmfArray );
          Geometry->GetVectorX()->SetNumberType( XDMF_FLOAT32_TYPE );
          }
        if( !Geometry->GetVectorY() )
          {
          Geometry->SetVectorY( new XdmfArray );
          Geometry->GetVectorY()->SetNumberType( XDMF_FLOAT32_TYPE );
          }
        if( !Geometry->GetVectorZ() )
          {
          Geometry->SetVectorZ( new XdmfArray );
          Geometry->GetVectorZ()->SetNumberType( XDMF_FLOAT32_TYPE );
          }
        Geometry->GetVectorX()->SetNumberOfElements( Dimensions[2] );
        Geometry->GetVectorY()->SetNumberOfElements( Dimensions[1] );
        Geometry->GetVectorZ()->SetNumberOfElements( Dimensions[0] );
        Origin = Geometry->GetOrigin();
        Geometry->GetVectorX()->Generate(
          Origin[0],
          Origin[0] + ( Geometry->GetDx() * ( Dimensions[2] - 1 ) ) );
        Geometry->GetVectorY()->Generate(
          Origin[1],
          Origin[1] + ( Geometry->GetDy() * ( Dimensions[1] - 1 ) ) );
        Geometry->GetVectorZ()->Generate(
          Origin[2],
          Origin[2] + ( Geometry->GetDz() * ( Dimensions[0] - 1 ) ) );
        }
      int sstart[3];
      sstart[0] = start[0];
      sstart[1] = start[1];
      sstart[2] = start[2];
      
      vtkIdType cstart[3] = { 0, 0, 0 };
      vtkIdType cend[3];
      cstart[0] = vtkMAX(0, sstart[2]);
      cstart[1] = vtkMAX(0, sstart[1]);
      cstart[2] = vtkMAX(0, sstart[0]);
      
      cend[0] = start[2] + count[2]*readerStride[0]+1;
      cend[1] = start[1] + count[1]*readerStride[1]+1;
      cend[2] = start[0] + count[0]*readerStride[2]+1;
      
      vtkDebugWithObjectMacro(this->Reader,
                              "CStart: " << cstart[0] << ", " 
                              << cstart[1] << ", " 
                              << cstart[2] );
      vtkDebugWithObjectMacro(this->Reader,
                              "CEnd: " << cend[0] << ", " 
                              << cend[1] << ", " 
                              << cend[2] );
      
      // Set the Points
      for(Index=cstart[0], cc = 0 ; Index < cend[0] ; Index += readerStride[0])
        {
        XCoord->SetValue( cc++, 
                          Geometry->GetVectorX()->GetValueAsFloat64( Index ));
        } 
      for(Index = cstart[1], cc = 0 ; Index < cend[1]; Index+=readerStride[1] )
        {
        YCoord->SetValue( cc++ , 
                          Geometry->GetVectorY()->GetValueAsFloat64( Index ) );
        } 
      for(Index = cstart[2], cc = 0 ; Index < cend[2];Index += readerStride[2])
        {
        ZCoord->SetValue( cc++ , 
                          Geometry->GetVectorZ()->GetValueAsFloat64( Index ) );
        }
      
      stride[2] = readerStride[0];
      stride[1] = readerStride[1];
      stride[0] = readerStride[2];
      
      // vGrid->SetExtent(upext);    
      }
    }
  vtkDataSet *dataSet = static_cast<vtkDataSet *>(output);
  for ( cc = 0; cc < dataSet->GetPointData()->GetNumberOfArrays(); cc ++ )
    {
    dataSet->GetPointData()->RemoveArray(
      dataSet->GetPointData()->GetArrayName(cc));
    }
  int haveActive = 0;
  for( cc = 0 ; cc < xdmfGrid->GetNumberOfAttributes() ; cc++ )
    {
    XdmfInt32 AttributeCenter;
    XdmfInt32 AttributeType;
    int       Components;
    XdmfAttribute       *Attribute;
    
    Attribute = xdmfGrid->GetAttribute( cc );
    const char *name = Attribute->GetName();
    int status = 1;
    AttributeCenter = Attribute->GetAttributeCenter();
    AttributeType = Attribute->GetAttributeType();
    Components = 1;
    switch (AttributeType) 
      {
      case XDMF_ATTRIBUTE_TYPE_TENSOR :
        Components = 9;
        break;
      case XDMF_ATTRIBUTE_TYPE_VECTOR:
        Components = 3;
        break;
      default :
        Components = 1;
        break;
      }
    if (name )
      {
      if ( AttributeCenter == XDMF_ATTRIBUTE_CENTER_GRID || 
           AttributeCenter == XDMF_ATTRIBUTE_CENTER_NODE)
        {
        status = pointDataArraySelection->ArrayIsEnabled(name);
        }
      else
        {
        status = cellDataArraySelection->ArrayIsEnabled(name);
        }
      }
    if ( !status )
      {
      continue;
      }
    status = 1;
    vtkDebugWithObjectMacro(this->Reader,
                            "Array with name: " 
                            << name << " has status: " << status);
    // attrNode = this->DOM->FindElement("Attribute", cc);
    attrNode = Attribute->GetElement();
    // dataNode = this->DOM->FindElement("DataStructure", 0, attrNode);
    // Find the DataTransform or DataStructure below the <Attribute>
    dataNode = xdmfDOM->FindElement(NULL, 0, attrNode);
    this->DataItem->SetElement(dataNode);
    this->DataItem->UpdateInformation();
    grid->DataDescription = this->DataItem->GetDataDesc();
    
    if ( Attribute && status )
      {
      //Attribute->Update();
      // XdmfDataItem dataItem;
      XdmfArray *values;

      XdmfInt64 realcount[4] = { 0, 0, 0, 0 };
      realcount[0] = count[0];
      realcount[1] = count[1];
      realcount[2] = count[2];
      realcount[3] = count[3];
      if ( AttributeCenter == XDMF_ATTRIBUTE_CENTER_NODE ||
           AttributeCenter == XDMF_ATTRIBUTE_CENTER_GRID )
        {
        // Point count is 1 + cell extent
        realcount[0] ++;
        realcount[1] ++;
        realcount[2] ++;
        }
      /*
        XdmfArray *values = this->FormatMulti->ElementToArray( 
        dataNode, this->Internals->DataDescriptions[currentGrid] );
      */
      vtkDebugWithObjectMacro(this->Reader, 
                              "Topology class: " 
                              << xdmfGrid->GetTopology()->GetClassAsString());
#ifndef XDMF_NO_MPI
     this->DataItem->SetDsmBuffer(this->DsmBuffer);
#endif
      if(xdmfGrid->GetTopology()->GetClass() != XDMF_UNSTRUCTURED)
        {
        XdmfDataDesc* ds = grid->DataDescription;
        XdmfInt64 realdims[XDMF_MAX_DIMENSION];
        XdmfInt32 realrank = ds->GetShape(realdims);
        if ( realrank == 4 )
          {
          realcount[3] = realdims[3];
          }
        this->DataItem->GetDataDesc()->SelectHyperSlab(start, stride, realcount);
        vtkDebugWithObjectMacro(this->Reader,
                                "Dims = " << ds->GetShapeAsString()
                                << "Slab = " << ds->GetHyperSlabAsString());
        // Only works for the structured and rectilinear grid
        vtkDebugWithObjectMacro(this->Reader,
                                "Preparing to Read :" 
                                << xdmfDOM->Get(dataNode, "CData"));
        if (this->DataItem->Update()==XDMF_FAIL) 
        {
          vtkGenericWarningMacro(<<"Reading of HDF5 dataset failed");
        }
        values = this->DataItem->GetArray();
        }
      else 
        {
        this->DataItem->Update();
        values = this->DataItem->GetArray();
        }
      this->ArrayConverter->SetVtkArray( NULL );
      if ( values )
        {
        vtkDataArray* vtkValues = this->ArrayConverter->FromXdmfArray(
          values->GetTagName(), 1, globalrank, Components, 0);
        
        vtkDebugWithObjectMacro(this->Reader,
                                "Reading array: " << name );
        vtkValues->SetName(name);
        
        // Special Cases
        if( AttributeCenter == XDMF_ATTRIBUTE_CENTER_GRID ) 
          {
          // Implement XDMF_ATTRIBUTE_CENTER_GRID as PointData
          XdmfArray *tmpArray = new XdmfArray;
          
          vtkDebugWithObjectMacro(this->Reader,
                                  "Setting Grid Centered Values");
          tmpArray->CopyType( values );
          tmpArray->SetNumberOfElements( dataSet->GetNumberOfPoints() );
          tmpArray->Generate( values->GetValueAsFloat64(0), 
                              values->GetValueAsFloat64(0) );
          vtkValues->Delete();
          this->ArrayConverter->SetVtkArray( NULL );
          vtkValues=this->ArrayConverter->FromXdmfArray(tmpArray->GetTagName(), 1, 1, Components, 0);
          if( !name )
            {
            name = values->GetTagName();
            }
          vtkValues->SetName( name );
          delete tmpArray;
          AttributeCenter = XDMF_ATTRIBUTE_CENTER_NODE;
          }
        switch (AttributeCenter)
          {
          case XDMF_ATTRIBUTE_CENTER_NODE :
            dataSet->GetPointData()->RemoveArray(name);
            dataSet->GetPointData()->AddArray(vtkValues);
            if ( Attribute->GetActive() )
              {
              haveActive = 1;
              switch( AttributeType )
                {
                case XDMF_ATTRIBUTE_TYPE_SCALAR :
                  dataSet->GetPointData()->SetActiveScalars( name );
                  break;
                case XDMF_ATTRIBUTE_TYPE_VECTOR :
                  dataSet->GetPointData()->SetActiveVectors( name );
                  break;
                case XDMF_ATTRIBUTE_TYPE_TENSOR :
                  dataSet->GetPointData()->SetActiveTensors( name );
                  break;
                default :
                  break;
                }
              }
            break;
          case XDMF_ATTRIBUTE_CENTER_CELL :
            dataSet->GetCellData()->RemoveArray(name);
            dataSet->GetCellData()->AddArray(vtkValues);
            if ( Attribute->GetActive() )
              {
              haveActive = 1;
              switch( AttributeType )
                {
                case XDMF_ATTRIBUTE_TYPE_SCALAR :
                  dataSet->GetCellData()->SetActiveScalars( name );
                  break;
                case XDMF_ATTRIBUTE_TYPE_VECTOR :
                  dataSet->GetCellData()->SetActiveVectors( name );
                  break;
                case XDMF_ATTRIBUTE_TYPE_TENSOR :
                  dataSet->GetCellData()->SetActiveTensors( name );
                  break;
                default :
                  break;
                }
              }
            break;
          default : 
            vtkErrorWithObjectMacro(this->Reader,
                                    << "Can't Handle Values at " 
                                    <<Attribute->GetAttributeCenterAsString());
            break;
          }
        if ( vtkValues )
          {
          vtkValues->Delete();
          }
        if ( grid->DataDescription ) 
          {
          // delete grid->DataDescription;
          // grid->DataDescription = 0;
          }
        }
      }
    }
  if ( !haveActive )
    {
    vtkDataSetAttributes* fd = dataSet->GetPointData();
    for ( cc = 0; cc < fd->GetNumberOfArrays(); cc ++ )
      {
      vtkDataArray* ar = fd->GetArray(cc);
      switch ( ar->GetNumberOfComponents() )
        {
        case 1:
          fd->SetActiveScalars(ar->GetName());
          break;
        case 3:
          fd->SetActiveVectors(ar->GetName());
          break;
        case 6:
          fd->SetActiveTensors(ar->GetName());
          break;
        }
      }
    fd = dataSet->GetCellData();
    for ( cc = 0; cc < fd->GetNumberOfArrays(); cc ++ )
      {
      vtkDataArray* ar = fd->GetArray(cc);
      switch ( ar->GetNumberOfComponents() )
        {
        case 1:
          fd->SetActiveScalars(ar->GetName());
          break;
        case 3:
          fd->SetActiveVectors(ar->GetName());
          break;
        case 6:
          fd->SetActiveTensors(ar->GetName());
          break;
        }
      }
    }
  /*
  const char* name = currentGridName;
  vtkCharArray *nameArray = vtkCharArray::New();
  nameArray->SetName("Name");
  char *str = nameArray->WritePointer(0, strlen(name)+1);
  sprintf(str, "%s", name);
  output->GetFieldData()->AddArray(nameArray);
  nameArray->Delete();
  */

  // Release Duplicate Big Data 
  grid->XMGrid->Release();
  return 1;
}
//----------------------------------------------------------------------------
class WithinTolerance: public vtkstd::binary_function<double, double, bool>
{
public:
    result_type operator()(first_argument_type a, second_argument_type b) const
    {
      bool result = (fabs(a-b)<=(a*1E-6));
      return (result_type)result;
    }
};
//----------------------------------------------------------------------------
int vtkXdmfReader::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  if ( !this->GetReadFromInputString() && !this->FileName )
    {
    vtkErrorMacro("Not Reading from String and File name not set");
    return 0;
    }
  if ( !this->DOM )
    {
    return 0;
    }
  
  //long int starttime = this->GetMTime();

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkDataObject *outStructure = outInfo->Get(vtkDataObject::DATA_OBJECT());
  if (!outStructure)
    {
    vtkErrorMacro( << "No Output VTK structure");
    return VTK_ERROR;
    }
  vtkXdmfReaderGrid *ptr = this->Internals->Data;
  if (!ptr)
    {
    vtkErrorMacro( << "No Output XMReaderGrid");
    return 1;
    }

  this->Internals->UpdatePiece     = 0;
  this->Internals->UpdateNumPieces = 1;

  if (outInfo->Has(
        vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER()))
    {
    this->Internals->UpdatePiece = 
      outInfo->Get(
        vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
    }
  if (outInfo->Has(
        vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES()))
    {
    this->Internals->UpdateNumPieces =
      outInfo->Get(
        vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
    }
  this->Internals->DsmBuffer = this->DsmBuffer;

  vtkDebugMacro( << "UpdatePiece " << this->Internals->UpdatePiece << " : UpdateNumPieces " << this->Internals->UpdateNumPieces);

  //
  // Find the correct time step
  //
  this->ActualTimeStep = this->TimeStep;

  if (this->OutputTemporal) 
    {
    if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS()))
      {
      double requestedTimeValue = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS())[0];
      this->ActualTimeStep = vtkstd::find_if(
        this->Internals->TimeValues.begin(), 
        this->Internals->TimeValues.end(), 
        vtkstd::bind2nd( WithinTolerance( ), requestedTimeValue ))
        - this->Internals->TimeValues.begin();
      }
    // last check to stop dodgy crashes
    if (this->ActualTimeStep>=this->Internals->TimeValues.size()) {
      this->ActualTimeStep = 0;
    }
    // If the time step changes, the information is invalid
    vtkDebugMacro("Temporal Data Requested. Reset Information ");    
    vtkXdmfReaderGrid *sptr = this->Internals->GetGrid(0);
    sptr->Information = outInfo;
    // Set the time step we are producing in the output information
    outStructure->GetInformation()->Set(vtkDataObject::DATA_TIME_STEPS(), 
      &this->Internals->TimeValues[this->ActualTimeStep], 1);
    }
  this->Internals->FindParallelism();
  // else switch (this->OutputVTKType)
  switch (this->OutputVTKType)
    {
    case VTK_POLY_DATA:
    case VTK_UNSTRUCTURED_GRID:
    case VTK_IMAGE_DATA:
    case VTK_UNIFORM_GRID:
    case VTK_RECTILINEAR_GRID:
    case VTK_STRUCTURED_GRID:
    {
    ptr = this->Internals->GetGrid(0);
    //let the grid know what part it needs to fill in on this processor
    ptr->Information->CopyEntry(outInfo, 
       vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT());
    ptr->Information->CopyEntry(outInfo,
       vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
    ptr->Information->CopyEntry(outInfo,
       vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
    ptr->Information->CopyEntry(outInfo,
       vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS());

    vtkDebugMacro("Filling in atomic " 
                  << outStructure->GetClassName() 
                  << " with " << ptr->Name.c_str());    
    this->Internals->RequestGridData(ptr, outStructure, this->ActualTimeStep, 0,
                                           0.0, 1.0);
    }
    break;
    case VTK_MULTIBLOCK_DATA_SET:
    {
    vtkMultiBlockDataSet *outMB = 
      vtkMultiBlockDataSet::SafeDownCast(outStructure);
    if (!outMB)
      {
      vtkErrorMacro(<< "Wrong output type");
      return 0;
      }
    int outputGrid = 0;
    vtkstd::vector<vtkXdmfReaderGrid *>::iterator it;
    for ( it = ptr->Children.begin();
          it != ptr->Children.end();
          ++it )
      {
      double progressS = 
        (double)outputGrid / (double)this->NumberOfEnabledActualGrids;
      double progressE = 
        (double)(outputGrid+1) / (double)this->NumberOfEnabledActualGrids;

      vtkXdmfReaderGrid *child = (*it);      

      if ( (*it)->Enabled )
        {
        // if ((outputGrid%this->Internals->UpdateNumPieces) == this->Internals->UpdatePiece) 
        if(1)
          {
          vtkDataObject *output =
            vtkDataObjectTypes::NewDataObject(child->vtkType);
          if (output && output->IsA("vtkMultiBlockDataSet"))
            {
            outMB->SetBlock(outputGrid, vtkMultiBlockDataSet::SafeDownCast(output));
            }
          else
            {
            outMB->SetBlock(outputGrid, vtkDataSet::SafeDownCast(output));
            //outMB->SetMetaData(outputGrid, grid->GetInformation());
            }
          vtkDebugMacro(<< this->Internals->UpdatePiece << " filling in stripe " << outputGrid
                        << " with a " << output->GetClassName() 
                        << " from " << (*it)->Name.c_str());
          this->Internals->RequestGridData(
            /*(*it)->Name.c_str(),*/
            child,
            output,
            this->ActualTimeStep,
            1,
            progressS,
            progressE
            );
          output->Delete();
          }
        else
          {
          vtkDebugMacro(<< this->Internals->UpdatePiece << " ignoring stripe " << outputGrid 
                        << " from " << (*it)->Name.c_str());
          outMB->SetBlock(outputGrid, static_cast<vtkDataSet*>(NULL));
          //outMB->SetMetaData(outputGrid, grid->GetInformation());
          }
        this->UpdateProgress(progressE);
        outputGrid++;
        }
      }
    }
    break;
    default:
      vtkErrorMacro("Unknown type to generate information for");
      return 0;
    }

  //long int endtime = this->GetMTime();

  return 1;
}

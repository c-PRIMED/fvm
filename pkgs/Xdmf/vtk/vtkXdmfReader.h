/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkXdmfReader.h,v $
  Language:  C++
  Date:      $Date: 2009-04-06 13:16:06 $
  Version:   $Revision: 1.12 $

  Copyright (c) 1993-2002 Ken Martin, Will Schroeder, Bill Lorensen 
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkXdmfReader - read eXtensible Data Model and Format files
// .SECTION Description
// vtkXdmfReader is a source object that reads XDMF data.
// The output of this reader is a vtkMultiGroupDataSet with one group for
// every enabled grid in the domain.
// The superclass of this class, vtkDataReader, provides many methods for
// controlling the reading of the data file, see vtkDataReader for more
// information.
// .SECTION Caveats
// uses the XDMF API
// .SECTION See Also
// vtkDataReader

#ifndef __vtkXdmfReader_h
#define __vtkXdmfReader_h

#include "vtkDataReader.h"

class vtkDataObject;
class vtkDataArraySelection;
class vtkCallbackCommand;
class vtkMultiProcessController;
class vtkXdmfReaderInternal;
class vtkXdmfReaderGrid;

//BTX
class XdmfDsmBuffer;
class XdmfDOM;
//ETX

class VTK_EXPORT vtkXdmfReader : public vtkDataReader
{
public:
  static vtkXdmfReader* New();
  vtkTypeRevisionMacro(vtkXdmfReader, vtkDataReader);
  void PrintSelf(ostream& os, vtkIndent indent);

  // DOMAINS ///////////////////////////////////////////////////////////////
  // Description:
  // Get number of domains.
  int GetNumberOfDomains();

  // Description:
  // Get the name of domain at index.
  const char* GetDomainName(int idx);

  // Description:
  // Get/Set the current domain name. If none is set, the first domain will be
  // used. 
  virtual void SetDomainName(const char*);
  vtkGetStringMacro(DomainName);

  // GRIDS ///////////////////////////////////////////////////////////////////
  // Description:
  // Get number of grids in the current domain.
  int GetNumberOfGrids();

  // Description:
  // Get/Set the current grid name.
  void SetGridName(const char*);

  // Description:
  // Get the name of grid at index.
  const char* GetGridName(int idx);
  int GetGridIndex(const char* name);

  // Description:
  // Enable grids.
  void EnableGrid(const char* name);
  void EnableGrid(int idx);
  void EnableAllGrids();

  // Description:
  // Disable grids
  void DisableGrid(const char* name);
  void DisableGrid(int idx);
  void DisableAllGrids();
  void RemoveAllGrids(); // <<--FIXME: remove me.

  // Description:
  // Get current enable/disable of the grid
  int GetGridSetting(const char* name);
  int GetGridSetting(int idx);

  // ATTRIBUTES ///////////////////////////////////////////////////////////////
  // Description:
  // Get the data array selection tables used to configure which data
  // arrays are loaded by the reader.
  vtkGetObjectMacro(PointDataArraySelection, vtkDataArraySelection);
  vtkGetObjectMacro(CellDataArraySelection, vtkDataArraySelection);
  
  // Description:  
  // Get the number of point or cell arrays available in the input.
  int GetNumberOfPointArrays();
  int GetNumberOfCellArrays();
  
  // Description:
  // Get the name of the point or cell array with the given index in
  // the input.
  const char* GetPointArrayName(int index);
  const char* GetCellArrayName(int index);

  // Description:
  // Get/Set whether the point or cell array with the given name is to
  // be read.
  int GetPointArrayStatus(const char* name);
  int GetCellArrayStatus(const char* name);
  void SetPointArrayStatus(const char* name, int status);  
  void SetCellArrayStatus(const char* name, int status);  

  // Description:
  // Set whether the all point or cell arrays are to
  // be read.
  void EnableAllArrays();
  void DisableAllArrays();

  // PARAMETERS ///////////////////////////////////////////////////////////////
  // Description:
  // Get the number of Parameters
  int GetNumberOfParameters();

  // Description:
  // Get Parameter Type
  int GetParameterType(int index);
  int GetParameterType(const char *Name);
  const char *GetParameterTypeAsString(int index);
  const char *GetParameterTypeAsString(const char *Name);

  // Description:
  // Get start, stride, count
  int GetParameterRange(int index, int Shape[3]);
  int GetParameterRange(const char *Name, int Shape[3]);
  const char *GetParameterRangeAsString(int index);
  const char *GetParameterRangeAsString(const char *Name);

  // Description:
  // Get Parameter Name
  const char *GetParameterName(int index);

  // Description:
  // Set/Get Parameter Current Index
  int SetParameterIndex(const char *Name, int CurrentIndex); 
  int SetParameterIndex(int ParameterIndex, int CurrentIndex); 
  int GetParameterIndex(const char *Name);
  int GetParameterIndex(int index);

  // Description:
  // Get Length of Parameter
  int GetParameterLength(const char *Name);
  int GetParameterLength(int index);

  // Description:
  // Get the Current Value of the Parameter
  const char *GetParameterValue(int index);
  const char *GetParameterValue(const char *Name);

  // STRIDE ///////////////////////////////////////////////////////////////////
  // Description:
  // Set / get stride
  void SetStride(int x, int y, int z);
  void SetStride(int xyz[3])
    {
    this->SetStride(xyz[0], xyz[1], xyz[2]);
    }
  vtkGetVector3Macro(Stride, int);

  // MISCELANEOUS /////////////////////////////////////////////////////////////
  // Description:
  // Get the Low Level XdmfDOM
  const char *GetXdmfDOMHandle();

  // Description:
  // Get the Low Level XdmfGrid
  //Disable for now
  //const char *GetXdmfGridHandle(int idx);

  // Description:
  // Determine if the file can be readed with this reader.
  virtual int CanReadFile(const char* fname);
  
  // Description:
  // Set the controller used to coordinate parallel reading.
  void SetController(vtkMultiProcessController* controller);
  
  // Return the controller used to coordinate parallel reading. By default,
  // it is the global controller.
  vtkGetObjectMacro(Controller,vtkMultiProcessController);

  // Set DsmBubffer
  void SetDsmBuffer(void *Bufp);
  // Get DsmBubffer
  void *GetDsmBuffer();

  // Set the Timestep to be read. This is provided for compatibility
  // reasons only and should not be used. The correct way to
  // request time is using the UPDATE_TIME_STEPS information key
  // passed from downstream.
  vtkSetMacro(TimeStep, int);
  vtkGetMacro(TimeStep, int);

  // Description:
  // Save the range of valid timestep index values. This can be used by the PAraView GUI
  int TimeStepRange[2];
  vtkGetVector2Macro(TimeStepRange, int);

protected:
  vtkXdmfReader();
  ~vtkXdmfReader();

  // Description:
  // This methods parses the XML. Returns true on success. This method can be
  // called repeatedly. It has checks to ensure that the XML parsing is done
  // only if needed.
  bool ParseXML();

  virtual int ProcessRequest(vtkInformation *request,
                             vtkInformationVector **inputVector,
                             vtkInformationVector *outputVector);
  
  virtual int RequestDataObject(vtkInformationVector *outputVector);
  
  virtual int RequestData(vtkInformation *, vtkInformationVector **,
                          vtkInformationVector *);
  virtual int RequestInformation(vtkInformation *, vtkInformationVector **,
                                 vtkInformationVector *);
  virtual int FillOutputPortInformation(int port, vtkInformation *info);

  int  UpdateDomains();
  void UpdateRootGrid();
  void UpdateGrids(vtkXdmfReaderGrid *parent, void *ParentNode);
  void FindTimeValues();
  void FindAllTimeValues(vtkXdmfReaderGrid *ptr);
  void AssignTimeIndex(vtkXdmfReaderGrid *ptr);

  // Array selection helpers /////////////////////////////////////////////////
  static void SelectionModifiedCallback(vtkObject* caller, unsigned long eid,
                                        void* clientdata, void* calldata);

  vtkDataArraySelection* PointDataArraySelection;
  vtkDataArraySelection* CellDataArraySelection;

  vtkCallbackCommand* SelectionObserver;

  //
  vtkXdmfReaderInternal* Internals;
  XdmfDOM         *DOM;
  vtkMultiProcessController *Controller;

  char* DomainName;

  char* GridName;
  int NumberOfEnabledActualGrids;

  int Stride[3];

  int            GridsModified;
  int            OutputsInitialized;
  int            OutputVTKType;
  XdmfDsmBuffer *DsmBuffer;
  int            OutputTemporal;
  unsigned int   ActualTimeStep;
  int            TimeStep;
private:
  vtkXdmfReader(const vtkXdmfReader&); // Not implemented
  void operator=(const vtkXdmfReader&); // Not implemented  
};

#endif //__vtkXdmfReader_h

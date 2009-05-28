/*******************************************************************/
/*                               XDMF                              */
/*                   eXtensible Data Model and Format              */
/*                                                                 */
/*  Id : $Id: XdmfDsmH5Slice.cxx,v 1.1 2007-10-19 18:55:10 dave.demarle Exp $  */
/*  Date : $Date: 2007-10-19 18:55:10 $ */
/*  Version : $Revision: 1.1 $ */
/*                                                                 */
/*  Author:                                                        */
/*     Jerry A. Clarke                                             */
/*     clarke@arl.army.mil                                         */
/*     US Army Research Laboratory                                 */
/*     Aberdeen Proving Ground, MD                                 */
/*                                                                 */
/*     Copyright @ 2007 US Army Research Laboratory                */
/*     All Rights Reserved                                         */
/*     See Copyright.txt or http://www.arl.hpc.mil/ice for details */
/*                                                                 */
/*     This software is distributed WITHOUT ANY WARRANTY; without  */
/*     even the implied warranty of MERCHANTABILITY or FITNESS     */
/*     FOR A PARTICULAR PURPOSE.  See the above copyright notice   */
/*     for more information.                                       */
/*                                                                 */
/*******************************************************************/
#include "Xdmf.h"
#include "mpi.h"
#include "pthread.h"

// Node 0 Creates the Structure for the DataSet
// The DataSet is Size x JDIM x IDIM doubles
// Each Node then Writes 1 K Plane of Data
// The Last Node Reads the Corners of the DataSet and
// reads one J Plane Half way up
#define IDIM    15
#define JDIM    30
#define DataSetName "DSM:TestFile.h5:/TestDataSets/Values1"

// Fill an array with Values
// The Values corespond to the rank and position
void
FillArray(int rank, int size, XdmfArray *array){
    int i, j, cntr = 0;
    XdmfInt64   Dimensions[] = {1, JDIM, IDIM};
    XdmfFloat64 Value;

    array->SetShape(3, Dimensions);
    for(j=0;j<JDIM;j++){
        for(i=0;i<IDIM;i++){
            Value = rank + (j * .01) + (i * .0001);
            array->SetValue(cntr++, Value);
        }
    }
}

// MPI Barrier with a print
void
WaitForAll(int rank, int value, MPI_Comm comm){

    cout << "(" << rank << ") in barrier " << value << endl;
    MPI_Barrier(comm);
    cout << "(" << rank << ") out of barrier " << value << endl;
}

int
main( int argc, char **argv ) {

MPI_Comm        XdmfComm;
XdmfInt32       rank, size, status;
XdmfArray       Data;
XdmfDsmCommMpi  *MyComm = new XdmfDsmCommMpi;
XdmfDsmBuffer   *MyDsm = new XdmfDsmBuffer;
pthread_t       thread1;
XdmfInt32       Rank = 3;
XdmfInt64       Dimensions[] = { 1, JDIM, IDIM };
XdmfInt64       Start[] = {0, 0, 0};
XdmfInt64       Stride[] = {1, 1, 1};
XdmfInt64       Count[] = {1, JDIM, IDIM};
XdmfHDF         *ExternalDataSet;
XdmfArray       *InCoreCoordinates;
int              i, k, provided;
double          *DblPtr;



// Initialize MPI for Threads
MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);

// Initialize Communication Object
MyComm->Init();
// New Communicator for Xdmf Transactions
MyComm->DupComm( MPI_COMM_WORLD );

cout << "Hello from " << MyComm->GetId() << " of " << MyComm->GetTotalSize() << endl;
rank = MyComm->GetId();
size = Dimensions[0] = MyComm->GetTotalSize();

// Uniform Dsm : every node has a buffer of 1000000. Addresses are sequential
MyDsm->ConfigureUniform(MyComm, 1000000);

// Start another thread to handle DSM requests from other nodes
pthread_create(&thread1, NULL, &XdmfDsmBufferServiceThread, (void *)MyDsm);
// Wait for is to be ready
while(!MyDsm->GetThreadDsmReady()){
    // Spin
}
WaitForAll(rank, 1, MPI_COMM_WORLD);
// Begin HDF5 
ExternalDataSet = new XdmfHDF();
ExternalDataSet->SetDsmBuffer(MyDsm);
ExternalDataSet->SetNumberType(XDMF_FLOAT64_TYPE);
ExternalDataSet->SetShape( Rank, Dimensions );
Data.CopyType(ExternalDataSet);
FillArray(rank, size, &Data);

// Node 0 Creates the Structure for the File
// Writes Last Value to Make "File" the Correct size
if(rank == 0){
    XdmfInt64   BackCorner[] = {size-1,JDIM-1,IDIM-1};
    XdmfArray   Dummy;

    if(ExternalDataSet->Open( DataSetName, "rw" ) == XDMF_FAIL){
        // Shouldn't be necessary
        cout << "!! Node 0 Creating Dataset" << endl;
        ExternalDataSet->CreateDataset(DataSetName);
    }
    // Write Last Value
    ExternalDataSet->SelectCoordinates(1, BackCorner);
    Dummy.CopyType(ExternalDataSet);
    Dummy.SetNumberOfElements(1);
    ExternalDataSet->Write(&Dummy);
    ExternalDataSet->Close();
}
WaitForAll(rank, 2, MPI_COMM_WORLD);

// Everyone Writes One K Plane of Data
Start[0] = rank;
if(ExternalDataSet->Open( DataSetName, "rw" ) == XDMF_FAIL){
     cout << "(" << rank << ") !!! Error Dataset has not been created" << endl;
}
ExternalDataSet->SelectHyperSlab(Start, Stride, Count);
ExternalDataSet->Write(&Data);
ExternalDataSet->Close();

WaitForAll(rank, 3, MPI_COMM_WORLD);
// Last Node Reads Data
if(rank == (size - 1)){
    XdmfInt64   Corners[] = {0,0,0,  0,0,IDIM - 1,  size-1,JDIM-1,0, size-1,JDIM-1,IDIM-1};
    XdmfInt64   SliceStart[] = { 0,JDIM/2,0};
    XdmfInt64   SliceStride[] = { 1, 1, 1};
    XdmfInt64   SliceCount[] = {size, 1, IDIM};
    XdmfInt64   *SliceSize = SliceCount;

    InCoreCoordinates = new XdmfArray();
    InCoreCoordinates->CopyType(ExternalDataSet);
    // Read the corners
    InCoreCoordinates->SetNumberOfElements(4);
    ExternalDataSet->Open( DataSetName, "r" );
    ExternalDataSet->SelectCoordinates(4, Corners);
    ExternalDataSet->Read( InCoreCoordinates );
    cout << endl;
    cout << "4 of the Corners == ";
    DblPtr = (double *)InCoreCoordinates->GetDataPointer();
    for( k = 0 ; k < 4 ; k++ ){
          cout << " " << *DblPtr++;
    }
    cout << endl;
    // Read 1 J Plane
    InCoreCoordinates->SetShape(3, SliceSize);
    ExternalDataSet->SelectHyperSlab(SliceStart, SliceStride, SliceCount);
    ExternalDataSet->Read(InCoreCoordinates);
    cout << "Slice  == ";
    cout << endl;
    DblPtr = (double *)InCoreCoordinates->GetDataPointer();
    for( k = 0 ; k < size ; k++ ){
        cout << "    ";
        for(i=0;i<IDIM;i++){
          cout << " " << *DblPtr++;
        }
        cout << endl;
    }
    cout << endl;
    ExternalDataSet->Close();
    MyDsm->SendDone();
}

WaitForAll(rank, 4, MPI_COMM_WORLD);
pthread_join( thread1, NULL);

delete MyDsm;
delete MyComm;
MPI_Finalize();

return 0;
}




/*******************************************************************/
/*                               XDMF                              */
/*                   eXtensible Data Model and Format              */
/*                                                                 */
/*  Id : $Id: XdmfDsmH5.cxx,v 1.1 2007-10-19 18:55:10 dave.demarle Exp $  */
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

int DsmReady = 0;


void *DoServer(void *ptr){
    XdmfDsmBuffer   *MyDsm, *OrigDsm = (XdmfDsmBuffer *)ptr;
    XdmfInt32   op = 0;
    int rank;

    MyDsm = new XdmfDsmBuffer;
    MyDsm->Copy(OrigDsm);
    rank = MyDsm->GetComm()->GetId();
    cout << " Starting DSM Service on node " << MyDsm->GetComm()->GetId() << endl;
    DsmReady = 1;
    while(op != XDMF_DSM_OPCODE_DONE){
        op = 0;
        MyDsm->ServiceOnce(&op);
    }
    cout << " Ending DSM Service on node " << MyDsm->GetComm()->GetId() << endl;
}

void
WaitForAll(int rank, int value, MPI_Comm comm){

    cout << "(" << rank << ") in barrier " << value << endl;
    MPI_Barrier(comm);
    cout << "(" << rank << ") out of barrier " << value << endl;
}

int
main( int argc, char **argv ) {

MPI_Comm        XdmfComm;
XdmfInt32       rank, status;
XdmfArray       Data;
XdmfDsmCommMpi  *MyComm = new XdmfDsmCommMpi;
XdmfDsmBuffer   *MyDsm = new XdmfDsmBuffer;
pthread_t       thread1;
const char    *DataSetNameConst = "DSM:TestFile.h5:/TestDataSets/Values1";
XdmfInt32  Rank = 3;
XdmfInt64  Dimensions[] = { 10, 20, 30 };
XdmfHDF    *ExternalDataSet;
XdmfArray  *InCoreCoordinates;
char    *DataSetName;
int    i, k, provided;
double    *DataFromSomewhereElse;



// MPI_Init(&argc, &argv);
MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);


MyComm->Init();
// New Communicator for Xdmf Transactions
MyComm->DupComm( MPI_COMM_WORLD );

cout << "Hello from " << MyComm->GetId() << " of " << MyComm->GetTotalSize() << endl;
rank = MyComm->GetId();

// Uniform Dsm : every node has a buffer of 1000000. Addresses are sequential
MyDsm->ConfigureUniform(MyComm, 1000000);

pthread_create( &thread1, NULL, DoServer, MyDsm);
// pthread_create( &thread1, NULL, DoServer, &td);

while(!DsmReady){
    // Spin
}
cout << "Service Ready on " << rank << endl;
WaitForAll(rank, 1, MPI_COMM_WORLD);
// Begin HDF5 
ExternalDataSet = new XdmfHDF();
// ExternalDataSet->SetDebug(1);
ExternalDataSet->SetDsmBuffer(MyDsm);
Data.SetNumberType( XDMF_FLOAT64_TYPE );
Data.SetShape( Rank, Dimensions );
Data.Generate( 0, 5999 );
DataSetName = new char[ strlen(DataSetNameConst) + 1 ];
strcpy(DataSetName, DataSetNameConst);
if(rank == 0){
    ExternalDataSet->CopyType(&Data);
    ExternalDataSet->CopyShape(&Data);
    if(ExternalDataSet->Open( DataSetName, "rw" ) == XDMF_FAIL){
        ExternalDataSet->CreateDataset(DataSetName);
    }
    ExternalDataSet->Write(&Data);
    ExternalDataSet->Close();
}

WaitForAll(rank, 2, MPI_COMM_WORLD);
if(rank == 1){
    // cout << "(" << rank << ") Sleeping" << endl;
    // sleep(10);
    // cout << "(" << rank << ") Sending Done" << endl;
    InCoreCoordinates = new XdmfArray();
    InCoreCoordinates->SetNumberTypeFromString( "XDMF_FLOAT64_TYPE" );
    InCoreCoordinates->SetShapeFromString("4");
    ExternalDataSet->Open( DataSetName, "r" );
    ExternalDataSet->SelectCoordinatesFromString("0 0 0    0 0 29   9 19 0   9 19 29");
    ExternalDataSet->Read( InCoreCoordinates );
    ExternalDataSet->Close();
    cout << endl;
    cout << "4 of the Corners == ";
    DataFromSomewhereElse = (double *)InCoreCoordinates->GetDataPointer();
    for( k = 0 ; k < 4 ; k++ ){
          cout << " " << *DataFromSomewhereElse++;
    }
    cout << endl;
    MyDsm->SendDone();
}

WaitForAll(rank, 3, MPI_COMM_WORLD);
pthread_join( thread1, NULL);

delete MyDsm;
delete MyComm;
MPI_Finalize();

return 0;
}




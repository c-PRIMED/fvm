/*******************************************************************/
/*                               XDMF                              */
/*                   eXtensible Data Model and Format              */
/*                                                                 */
/*  Id : $Id: XdmfDsmThreads.cxx,v 1.1 2007-10-19 18:55:10 dave.demarle Exp $  */
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


void *DoServer(void *ptr){
    XdmfDsmBuffer   *MyDsm = (XdmfDsmBuffer *)ptr;
    XdmfInt32   op = 0;

    cout << " Starting DSM Service on node " << MyDsm->GetComm()->GetId() << endl;
    while(op != XDMF_DSM_OPCODE_DONE){
        op = 0;
        MyDsm->ServiceOnce(&op);
    }
    cout << " Ending DSM Service on node " << MyDsm->GetComm()->GetId() << endl;
}

int
main( int argc, char **argv ) {

XdmfInt32       rank, status;
XdmfArray       Data;
XdmfDsmCommMpi  *MyComm = new XdmfDsmCommMpi;
XdmfDsmBuffer   *MyDsm = new XdmfDsmBuffer;
pthread_t       thread1;

MPI_Init(&argc, &argv);

// New Communicator
MyComm->Init();
cout << "Hello from " << MyComm->GetId() << " of " << MyComm->GetTotalSize() << endl;
rank = MyComm->GetId();

// Uniform Dsm : every node has a buffer of 1000000. Addresses are sequential
MyDsm->ConfigureUniform(MyComm, 1000000);

MPI_Barrier(MPI_COMM_WORLD);
pthread_create( &thread1, NULL, DoServer, MyDsm);

Data.SetNumberType(XDMF_INT64_TYPE);
Data.SetNumberOfElements(100);
Data.Generate(rank + 10, rank + 110);

if(rank == 1){
    cout << "Before values = " << Data.GetValues(0, 10) << endl;
    cout << "Put #1" << endl;
    status = MyDsm->Put(100, sizeof(XdmfInt64), Data.GetDataPointer());
    if(status != XDMF_SUCCESS){
     cout << "Send Failed" << endl;
    }
    cout << "Put #2" << endl;
    status = MyDsm->Put(100, sizeof(XdmfInt64), Data.GetDataPointer());
    if(status != XDMF_SUCCESS){
     cout << "Send Failed" << endl;
    }
    cout << "Put #3" << endl;
    status = MyDsm->Put(100, sizeof(XdmfInt64), Data.GetDataPointer());
    if(status != XDMF_SUCCESS){
     cout << "Send Failed" << endl;
    }
}

MPI_Barrier(MPI_COMM_WORLD);
// Just write first value
status = MyDsm->Put(rank * sizeof(XdmfInt64), sizeof(XdmfInt64), Data.GetDataPointer());
if(status != XDMF_SUCCESS){
    cout << "Send Failed" << endl;
}

MPI_Barrier(MPI_COMM_WORLD);

if(rank == 1){
    // Get First 10 Values
    status = MyDsm->Get(0, 10 * sizeof(XdmfInt64), Data.GetDataPointer());
    if(status != XDMF_SUCCESS){
        cout << "Get Failed" << endl;
    }
    cout << "After values = " << Data.GetValues(0, 10) << endl;
    MyDsm->SendDone();
}

pthread_join( thread1, NULL);
MPI_Barrier(MPI_COMM_WORLD);

delete MyDsm;
delete MyComm;
MPI_Finalize();

return 0;
}




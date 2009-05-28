/*******************************************************************/
/*                               XDMF                              */
/*                   eXtensible Data Model and Format              */
/*                                                                 */
/*  Id : $Id: XdmfDsmExample.cxx,v 1.1 2007-10-19 18:55:10 dave.demarle Exp $  */
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

//using namespace std;

void *DoServer(void *ptr){
    XdmfDsmBuffer   *MyDsm = (XdmfDsmBuffer *)ptr;
    XdmfInt32   op = 0;

    while(op != XDMF_DSM_OPCODE_DONE){
        op = 0;
        MyDsm->ServiceUntilIdle(&op);
        if(op){
            cout << endl << "MyServer Data = " << MyDsm->GetStorage()->GetValues(0, 100) << endl;
        }
    }
}

int
main( int argc, char **argv ) {

int	            i, Data[256];
XdmfInt32       rank, status, who;
XdmfInt64       start, end;
XdmfDsmCommMpi  *MyComm = new XdmfDsmCommMpi;
XdmfDsmBuffer         *MyDsm = new XdmfDsmBuffer;

MPI_Init(&argc, &argv);

// New Communicator
MyComm->Init();
cout << "Hello from " << MyComm->GetId() << " of " << MyComm->GetTotalSize() << endl;
rank = MyComm->GetId();

// MyDsm->ConfigureUniform(MyComm, 1000000, 2, 50);
MyDsm->ConfigureUniform(MyComm, 1000000);

who = MyDsm->AddressToId(1500000);
MyDsm->GetAddressRangeForId(who, &start, &end);

// cout << "Address Range for " << who << " = " << start << " - " << end << endl;
// MyDsm->Put(0, 5000000, Data);
if(rank == 0){
    pthread_t thread1;

    cout << "Starting Thread" << endl;
    pthread_create( &thread1, NULL, DoServer, MyDsm);
    pthread_join( thread1, NULL);
    cout << "Thread Returned" << endl;
}
if(rank == 1000){
    XdmfArray   Data;
    XdmfInt32   op = 0;

    Data.SetNumberType(XDMF_INT64_TYPE);
    Data.SetNumberOfElements(100);
    Data.Generate(rank + 10, rank + 110);
    status = MyDsm->Put(0, Data.GetCoreLength(), Data.GetDataPointer());
    cout << endl << "MyServer Data = " << MyDsm->GetStorage()->GetValues(0, 100) << endl;
    while(op != XDMF_DSM_OPCODE_DONE){
        op = 0;
        MyDsm->ServiceUntilIdle(&op);
        if(op){
            cout << endl << "MyServer Data = " << MyDsm->GetStorage()->GetValues(0, 100) << endl;
        }
    }
    
    /*
    XdmfInt32   i, Opcode, Source;
    XdmfInt64   Address, Length;
    for(i=1;i<MyComm->GetTotalSize();i++){
    status = MyDsm->ReceiveCommandHeader(&Opcode, &Source, &Address, &Length);
    if(status = XDMF_SUCCESS){
        cout << "Receive From " << Source << " Address " << Address << " Length " << Length << endl;
    }else{
        cout << "Receive Failed" << endl;
    }
    }
    */
}else{
    XdmfArray   Data;

    if(rank == (MyComm->GetTotalSize() - 1)){
        sleep(1);
        MyDsm->SendDone();
    }else{
    Data.SetNumberType(XDMF_INT64_TYPE);
    Data.SetNumberOfElements(100);
    Data.Generate(rank + 10, rank + 110);
    // cout << endl << "Client Data = " << Data.GetValues(0, 110) << endl;

    status = MyDsm->Put(rank * 16, Data.GetCoreLength(), Data.GetDataPointer());
    /*
    XdmfInt32   Opcode, Dest = 0;
    XdmfInt64   Address = 0, Length = 100;
    Address = rank * 100;
    status = MyDsm->SendCommandHeader(Opcode, Dest, Address, Length);
    */
    if(status == XDMF_SUCCESS){
        cout << "Send Succeeded " << endl;
    }else{
        cout << "Send Failed" << endl;
    }
    }
}

MPI_Finalize();

return 0;
}




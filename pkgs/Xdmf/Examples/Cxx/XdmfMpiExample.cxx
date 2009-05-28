/*******************************************************************/
/*                               XDMF                              */
/*                   eXtensible Data Model and Format              */
/*                                                                 */
/*  Id : $Id: XdmfMpiExample.cxx,v 1.1 2007-10-19 18:55:10 dave.demarle Exp $  */
/*  Date : $Date: 2007-10-19 18:55:10 $ */
/*  Version : $Revision: 1.1 $ */
/*                                                                 */
/*  Author:                                                        */
/*     Jerry A. Clarke                                             */
/*     clarke@arl.army.mil                                         */
/*     US Army Research Laboratory                                 */
/*     Aberdeen Proving Ground, MD                                 */
/*                                                                 */
/*     Copyright @ 2002 US Army Research Laboratory                */
/*     All Rights Reserved                                         */
/*     See Copyright.txt or http://www.arl.hpc.mil/ice for details */
/*                                                                 */
/*     This software is distributed WITHOUT ANY WARRANTY; without  */
/*     even the implied warranty of MERCHANTABILITY or FITNESS     */
/*     FOR A PARTICULAR PURPOSE.  See the above copyright notice   */
/*     for more information.                                       */
/*                                                                 */
/*******************************************************************/
#include <Xdmf.h>

#include "mpi.h"

//using namespace std;

// Usage : XdmfFormatExample [ DataSetName ] 
int
main( int argc, char **argv ) {

XdmfHDF    *ExternalDataSet;
XdmfArray  *InCoreData;
XdmfArray  *InCoreSection;
XdmfInt32  Rank = 3;
// All Dimenions are "slowest changing first" : K,J,I
// and zero based
XdmfInt64  Dimensions[] = { 0, 20, 30 };
XdmfInt64  Start[] = { 0, 0, 0 };
XdmfInt64  Stride[] = { 1, 1, 1 };
XdmfInt64  Count[] = { 0, 20, 30 };
XdmfInt64  NumKPlanes = 10;

const char    *DataSetNameConst;
char    *DataSetName;
int    i, k;
double    *DataFromSomewhereElse;
int	size, rank;

if(argc > 1 ) {
  // i.e. NDGM:TestFile.h5:/TestDataSets/Values1
  DataSetNameConst = argv[1];
} else {
  // Domain:FileName:/HDF5Directory/..../HDF5DataSetName
  //  Domains : FILE, NDGM, GASS (Globus), CORE (malloc)
  DataSetNameConst = "FILE:TestFile.h5:/TestDataSets/Values1";
  }
DataSetName = new char[ strlen(DataSetNameConst) + 1 ];
strcpy(DataSetName, DataSetNameConst);

MPI_Init(&argc, &argv);
MPI_Comm_size(MPI_COMM_WORLD, &size);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);

cout << "Hello from Id " << rank << " of " << size << endl;
NumKPlanes = 2 * size;
Dimensions[0] = NumKPlanes;
Count[0] = 2;
Start[0] = rank * 2;

// Create Some Data
// XdmfArray and XdmfHDF ( and some others ) are derived
// classes from XdmfDataDesc. XdmfDataDesc has number type,
// number shape (rank and dimensions), and number selection.
// Selection is which "part" of the entire dataset we're
// dealing with. i.e WE could have a shape of [ 100, 200, 300 ]
// but just dealing with a block of [10, 20, 30 ] somewhere in
// the middle. By default Selction == Shape ... the whole thing.
InCoreData = new XdmfArray();
// InCoreData->SetGlobalDebug(rank);
InCoreData->SetNumberType( XDMF_FLOAT64_TYPE );
// InCoreData->SetShape( Rank, Dimensions );
InCoreData->SetShape( Rank, Count);
InCoreData->Generate( rank + 1, rank + 1);
// InCoreData->SelectHyperSlab(Start, Stride, Count);
// Convenience for :
// for( i = 0 ; i < InCoreData->GetNumberOfElements() ; i++ ){
//  InCoreData->SetValue( i, i );
//  }
//


// Create an external dataset if it doesn't exist
ExternalDataSet = new XdmfHDF();
// ExternalDataSet->SetDebug(rank);
// ExternalDataSet->SetCompression(9);
ExternalDataSet->SetNumberType( XDMF_FLOAT64_TYPE );
ExternalDataSet->SetShape( Rank, Dimensions );
// ExternalDataSet->CopyShape( InCoreData );
// In am MPI app the External Dataset would
// probable be much bigger that each node's InCore.
// i.e. ExternalDataSet->SetShape( 3, GloablDomainDimensions )
// but here we'll deal with just the InCore dataset size.
// ExternalDataSet->Open( DataSetName, "rw" );
if(ExternalDataSet->Open( DataSetName, "rw" ) == XDMF_FAIL){
	ExternalDataSet->CreateDataset(DataSetName);
	}
ExternalDataSet->SelectHyperSlab(Start, Stride, Count);
// InCoreData->SelectHyperSlab(Start, Stride, Count);

// Write the Data
ExternalDataSet->Write( InCoreData );
ExternalDataSet->Close();


// Read a J Plane down the middle
InCoreSection = new XdmfArray();
// Instead of allocating itself, use
// an external pointer : core dumps are
// your fault !!
DataFromSomewhereElse = new double[ NumKPlanes * 30 ];
InCoreSection->SetDataPointer( DataFromSomewhereElse );
InCoreSection->SetNumberTypeFromString( "XDMF_FLOAT64_TYPE" );
InCoreSection->SetNumberOfElements(NumKPlanes * 30);
ExternalDataSet->Open( DataSetName, "r" );
Count[0] = NumKPlanes;
Count[1] = 1;
Start[0] = 0;
ExternalDataSet->SelectHyperSlab(Start, Stride, Count);
ExternalDataSet->Read( InCoreSection );
ExternalDataSet->Close();

if(rank == 0){
cout << endl;
for( k = 0 ; k < NumKPlanes ; k++ ){
  cout << "k = " << k << ":";
  for( i = 0 ; i < 30 ; i++ ){
    cout << " " << *DataFromSomewhereElse++;
  }
  cout << endl;
}
}

delete [] DataSetName;

MPI_Finalize();

return 0;
}




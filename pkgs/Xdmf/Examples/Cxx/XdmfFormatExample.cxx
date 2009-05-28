/*******************************************************************/
/*                               XDMF                              */
/*                   eXtensible Data Model and Format              */
/*                                                                 */
/*  Id : $Id: XdmfFormatExample.cxx,v 1.1 2007-10-19 18:55:10 dave.demarle Exp $  */
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

//using namespace std;

// Usage : XdmfFormatExample [ DataSetName ] 
int
main( int argc, char **argv ) {

XdmfHDF    *ExternalDataSet;
XdmfArray  *InCoreData;
XdmfArray  *InCoreSection;
XdmfArray  *InCoreCoordinates;
XdmfInt32  Rank = 3;
// All Dimenions are "slowest changing first" : K,J,I
// and zero based
XdmfInt64  Dimensions[] = { 10, 20, 30 };

const char    *DataSetNameConst;
char    *DataSetName;
int    i, k;
double    *DataFromSomewhereElse;

if(argc > 1 ) {
  // i.e. NDGM:TestFile.h5:/TestDataSets/Values1
  DataSetNameConst = argv[1];
} else {
  // Domain:FileName:/HDF5Directory/..../HDF5DataSetName
  //  Domains : FILE, NDGM, GASS (Globus), CORE (malloc)
  DataSetNameConst = "SERIAL:TestFile.h5:/TestDataSets/Values1";
  }
DataSetName = new char[ strlen(DataSetNameConst) + 1 ];
strcpy(DataSetName, DataSetNameConst);

// Create Some Data
// XdmfArray and XdmfHDF ( and some others ) are derived
// classes from XdmfDataDesc. XdmfDataDesc has number type,
// number shape (rank and dimensions), and number selection.
// Selection is which "part" of the entire dataset we're
// dealing with. i.e WE could have a shape of [ 100, 200, 300 ]
// but just dealing with a block of [10, 20, 30 ] somewhere in
// the middle. By default Selction == Shape ... the whole thing.
InCoreData = new XdmfArray();
InCoreData->SetGlobalDebug(1);
InCoreData->SetNumberType( XDMF_FLOAT64_TYPE );
InCoreData->SetShape( Rank, Dimensions );
InCoreData->Generate( 0, 5999 );
// Convenience for :
// for( i = 0 ; i < InCoreData->GetNumberOfElements() ; i++ ){
//  InCoreData->SetValue( i, i );
//  }
//


// Create an external dataset if it doesn't exist
ExternalDataSet = new XdmfHDF();
ExternalDataSet->SetDebug(1);
// ExternalDataSet->SetCompression(9);
ExternalDataSet->CopyType( InCoreData );
ExternalDataSet->CopyShape( InCoreData );
// In am MPI app the External Dataset would
// probable be much bigger that each node's InCore.
// i.e. ExternalDataSet->SetShape( 3, GloablDomainDimensions )
// but here we'll deal with just the InCore dataset size.
if(ExternalDataSet->Open( DataSetName, "rw" ) == XDMF_FAIL){
        ExternalDataSet->CreateDataset(DataSetName);
        }

// Write the Data
ExternalDataSet->Write( InCoreData );
ExternalDataSet->Close();
// exit(1);

// Now Read in 4 values from the "corners"
InCoreCoordinates = new XdmfArray();
// Most Methods have a "FromString" convenience method.
// This makes wrapping for Python/Tcl/Java easier
InCoreCoordinates->SetNumberTypeFromString( "XDMF_FLOAT64_TYPE" );
InCoreCoordinates->SetShapeFromString("4");
// This is just the same dataset we wrote
// number type and shape get (re)set on Open
// ExternalDataSet->Open( DataSetName, "r" );
// Selections can be :
//  HyperSlab - Start, Stride, Count for Each Dimension
//  Coordinate - Parametric Coordinates
//  Function - Lex/Yacc stuff .... $0 * $1 / 1.2
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


InCoreSection = new XdmfArray();
// Instead of allocating itself, use
// an external pointer : core dumps are
// your fault !!
DataFromSomewhereElse = new double[ 150 ];
InCoreSection->SetDataPointer( DataFromSomewhereElse );
InCoreSection->SetNumberTypeFromString( "XDMF_FLOAT64_TYPE" );
// Make a 2D array to read back a section of the data
// We'll read in one "J" PLane
InCoreSection->SetShapeFromString("10 15");
ExternalDataSet->Open( DataSetName, "r" );
ExternalDataSet->SelectHyperSlabFromString("0 9 0", "1 1 2", "10 1 15");
// So from a 10x20x30 data set, start at 0,9,0. Stride 2 in Idim
// and read in 10 K x 1 J x ( 30 / 2 ) = 15 I ... a slice
ExternalDataSet->Read( InCoreSection );
ExternalDataSet->Close();

cout << endl;
for( k = 0 ; k < 10 ; k++ ){
  cout << "k = " << k << ":";
  for( i = 0 ; i < 15 ; i++ ){
    cout << " " << *DataFromSomewhereElse++;
  }
  cout << endl;
}

delete [] DataSetName;

return 0;
}




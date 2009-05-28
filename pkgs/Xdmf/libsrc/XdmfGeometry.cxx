/*******************************************************************/
/*                               XDMF                              */
/*                   eXtensible Data Model and Format              */
/*                                                                 */
/*  Id : $Id: XdmfGeometry.cxx,v 1.17 2009-01-26 21:15:21 clarke Exp $  */
/*  Date : $Date: 2009-01-26 21:15:21 $ */
/*  Version : $Revision: 1.17 $ */
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
#include "XdmfGeometry.h"

#include "XdmfTopology.h"
#include "XdmfDataItem.h"
#include "XdmfArray.h"
#include "XdmfDOM.h"
#include "XdmfHDF.h"

XdmfGeometry *GetXdmfGeometryHandle( void *Pointer ){
  //XdmfGeometry *tmp = (XdmfGeometry *)Pointer;
  return((XdmfGeometry *)Pointer);
  }

XdmfGeometry::XdmfGeometry() {
  this->SetElementName("Geometry");
  this->GeometryType = XDMF_GEOMETRY_NONE;
  this->Points = NULL;
  this->PointsAreMine = 1;
  this->VectorX = NULL;
  this->VectorY = NULL;
  this->VectorZ = NULL;
  this->SetOrigin( 0, 0, 0 );
  this->SetDxDyDz( 0, 0, 0 );
  }

XdmfGeometry::~XdmfGeometry() {
  if( this->PointsAreMine && this->Points )  delete this->Points;
  }

XdmfInt32
XdmfGeometry::Release()
{
  XdmfXmlNode node;
  XdmfInt32 Index = 0;
  XdmfXmlNode Node;

  if( this->PointsAreMine && this->Points ){
      delete this->Points;
      this->Points = NULL;
  }
  // this->NumberOfPoints = 0;
  Node = this->GetElement();
  node = this->DOM->FindDataElement(Index++, Node);
  // No Need to Release DataItems() since Data has been read 
  // and Stored in Internal Points
  return(XDMF_SUCCESS);
}
// Returns an existing DataItem or build a new one
XdmfDataItem *
XdmfGeometry::GetDataItem(XdmfInt32 Index, XdmfXmlNode Node){
    XdmfDataItem *di = NULL;
    XdmfXmlNode node;

    node = this->DOM->FindDataElement(Index, Node);
    if(node) {
        di = (XdmfDataItem *)this->GetCurrentXdmfElement(node);
    }
    if(!di){
        di = new XdmfDataItem;
        node = this->DOM->InsertNew(this->GetElement(), "DataItem");
        di->SetDOM(this->DOM);
        di->SetElement(node);
    }
    return(di);
}


XdmfInt32
XdmfGeometry::Build(){
    XdmfDataItem    *di = NULL;
    XdmfArray       *array;

    cout << "Building Geometry" << endl;
    if(XdmfElement::Build() != XDMF_SUCCESS) return(XDMF_FAIL);
    this->Set("GeometryType", this->GetGeometryTypeAsString());
    // Build Children from String , if provided
    if(this->BuildFromDataXml() == XDMF_SUCCESS) return(XDMF_SUCCESS);
    switch( this->GeometryType ){
      case XDMF_GEOMETRY_VXVYVZ:
            if(!this->VectorX || !this->VectorY || !this->VectorZ){
                XdmfErrorMessage("Vx Vy and Vz must be set");
                return(XDMF_FAIL);
            }
            // Vx
            di = this->GetDataItem(0, this->GetElement());
            di->SetArray(this->VectorX);
            if(this->VectorX->GetNumberOfElements() > 100) di->SetFormat(XDMF_FORMAT_HDF);
            di->Build();
            // Vy
            di = this->GetDataItem(1, this->GetElement());
            di->SetArray(this->VectorY);
            if(this->VectorY->GetNumberOfElements() > 100) di->SetFormat(XDMF_FORMAT_HDF);
            di->Build();
            // Vx
            di = this->GetDataItem(3, this->GetElement());
            di->SetArray(this->VectorZ);
            if(this->VectorZ->GetNumberOfElements() > 100) di->SetFormat(XDMF_FORMAT_HDF);
            di->Build();
        break;
      case XDMF_GEOMETRY_ORIGIN_DXDYDZ:
            // Origin
            di = this->GetDataItem(0, this->GetElement());
            di->SetFormat(XDMF_FORMAT_XML);
            array = di->GetArray();
            array->SetNumberOfElements(3);
            array->SetValues(0, this->Origin, 3);
            di->Build();
            // DxDyDz
            di = this->GetDataItem(1, this->GetElement());
            di->SetFormat(XDMF_FORMAT_XML);
            array = di->GetArray();
            array->SetNumberOfElements(3);
            array->SetValues(0, this->DxDyDz, 3);
            di->Build();
        break;
      default :
        if(this->Points){
            di = this->GetDataItem(0, this->GetElement());
            di->SetArray(this->Points);
            if(this->Points->GetNumberOfElements() > 100) di->SetFormat(XDMF_FORMAT_HDF);
            di->Build();
        }else{
            XdmfErrorMessage("XdmfGeometry->Points must be set for Geometry Type " << this->GetGeometryTypeAsString());
            return(XDMF_FAIL);
        }
        break;
    }
    return(XDMF_SUCCESS);
}

XdmfInt32
XdmfGeometry::Insert( XdmfElement *Child){
    if(Child && (
        XDMF_WORD_CMP(Child->GetElementName(), "DataItem") ||
        XDMF_WORD_CMP(Child->GetElementName(), "Information")
        )){
        return(XdmfElement::Insert(Child));
    }else{
        XdmfErrorMessage("Geometry can only Insert DataItem or Information elements");
    }
    return(XDMF_FAIL);
}

XdmfInt32
XdmfGeometry::SetOrigin( XdmfFloat64 X, XdmfFloat64 Y, XdmfFloat64 Z ){

this->Origin[0] = X;
this->Origin[1] = Y;
this->Origin[2] = Z;
return( XDMF_SUCCESS );
}

XdmfInt32
XdmfGeometry::SetOrigin( XdmfFloat64 *origin ){
return( this->SetOrigin( origin[0], origin[1], origin[2] ) );
}

XdmfInt32
XdmfGeometry::SetDxDyDz( XdmfFloat64 Dx, XdmfFloat64 Dy, XdmfFloat64 Dz ){
this->DxDyDz[0] = Dx;
this->DxDyDz[1] = Dy;
this->DxDyDz[2] = Dz;
return( XDMF_SUCCESS );
}


XdmfInt32
XdmfGeometry::SetDxDyDz( XdmfFloat64 *dxDyDz ){
return( this->SetDxDyDz( dxDyDz[0], dxDyDz[1], dxDyDz[2] ) );
}


XdmfArray *
XdmfGeometry::GetPoints(XdmfInt32 Create){
    if(!this->Points && Create){
        this->Points = new XdmfArray;
        this->PointsAreMine = 1;
    }
    return(this->Points);
}

XdmfInt32
XdmfGeometry::SetPoints( XdmfArray *points ){
    if(this->Points == points) return(XDMF_SUCCESS);
    if( this->PointsAreMine && this->Points ) delete this->Points;
    this->PointsAreMine = 0;
    this->Points = points;
    return( XDMF_SUCCESS );
    }

XdmfInt32
XdmfGeometry::SetGeometryTypeFromString( XdmfConstString geometryType ){

if( XDMF_WORD_CMP( geometryType, "X_Y_Z" ) ){
  this->GeometryType = XDMF_GEOMETRY_X_Y_Z;
  return(XDMF_SUCCESS);
  }
if( XDMF_WORD_CMP( geometryType, "X_Y" ) ){
  this->GeometryType = XDMF_GEOMETRY_X_Y;
  return(XDMF_SUCCESS);
  }
if( XDMF_WORD_CMP( geometryType, "XY" ) ){
  this->GeometryType = XDMF_GEOMETRY_XY;
  return(XDMF_SUCCESS);
  }
if( XDMF_WORD_CMP( geometryType, "XYZ" ) ){
  this->GeometryType = XDMF_GEOMETRY_XYZ;
  return(XDMF_SUCCESS);
  }
if( XDMF_WORD_CMP( geometryType, "ORIGIN_DXDYDZ" ) ){
  this->GeometryType = XDMF_GEOMETRY_ORIGIN_DXDYDZ;
  return(XDMF_SUCCESS);
  }
if( XDMF_WORD_CMP( geometryType, "VXVYVZ" ) ){
  this->GeometryType = XDMF_GEOMETRY_VXVYVZ;
  return(XDMF_SUCCESS);
  }
return( XDMF_FAIL );
}

XdmfString
XdmfGeometry::GetGeometryTypeAsString( void ){
static char Value[80];

switch( this->GeometryType ){
  case XDMF_GEOMETRY_VXVYVZ:
    strcpy( Value, "VXVYVZ" );
    break;
  case XDMF_GEOMETRY_ORIGIN_DXDYDZ:
    strcpy( Value, "ORIGIN_DXDYDZ" );
    break;
  case XDMF_GEOMETRY_X_Y_Z :
    strcpy( Value, "X_Y_Z" );
    break;
  case XDMF_GEOMETRY_X_Y :
    strcpy( Value, "X_Y" );
    break;
  case XDMF_GEOMETRY_XY :
    strcpy( Value, "XY" );
    break;
  default :
    strcpy( Value, "XYZ" );
    break;
  }
return( Value );
}

XdmfInt32
XdmfGeometry::UpdateInformation() {
XdmfConstString  Attribute;

if(XdmfElement::UpdateInformation() != XDMF_SUCCESS) return(XDMF_FAIL);
if( XDMF_WORD_CMP(this->GetElementType(), "Geometry") == 0){
    XdmfErrorMessage("Element type" << this->GetElementType() << " is not of type 'Geometry'");
    return(XDMF_FAIL);
}
Attribute = this->Get( "GeometryType" );
if(!Attribute){
    Attribute = this->Get( "Type" );
}
if( Attribute ){
  this->SetGeometryTypeFromString( Attribute );
} else {
  this->GeometryType = XDMF_GEOMETRY_XYZ;
}
if(!this->Name) this->SetName(GetUnique("Geometry_"));
return( XDMF_SUCCESS );
}

XdmfInt32
XdmfGeometry::Update() {


XdmfInt32  ArrayIndex;
XdmfInt64  Start[ XDMF_MAX_DIMENSION ];
XdmfInt64  Stride[ XDMF_MAX_DIMENSION ];
XdmfInt64  Count[ XDMF_MAX_DIMENSION ];
XdmfXmlNode     PointsElement;
XdmfArray       *points = NULL;
XdmfArray       *TmpArray;

if( this->GeometryType == XDMF_GEOMETRY_NONE ){
  if( this->UpdateInformation() == XDMF_FAIL ){
    XdmfErrorMessage("Can't Initialize From Element");
    return( XDMF_FAIL );
  }
}
if(XdmfElement::Update() != XDMF_SUCCESS) return(XDMF_FAIL);
ArrayIndex = 0;
if( ( this->GeometryType == XDMF_GEOMETRY_X_Y_Z ) ||
  ( this->GeometryType == XDMF_GEOMETRY_X_Y ) ||
  ( this->GeometryType == XDMF_GEOMETRY_XYZ ) ||
  ( this->GeometryType == XDMF_GEOMETRY_XY ) ){
 do {
  // Read the Data
  XdmfDebug("Reading Points ( SubElement #" << ArrayIndex + 1 << " )" );
  PointsElement = this->DOM->FindDataElement( ArrayIndex, Element );
  if( PointsElement ){
    XdmfDataItem PointsItem;
    if(PointsItem.SetDOM( this->DOM ) == XDMF_FAIL) return(XDMF_FAIL);
    if(PointsItem.SetElement(PointsElement, 0) == XDMF_FAIL) return(XDMF_FAIL);
    if(PointsItem.UpdateInformation() == XDMF_FAIL) return(XDMF_FAIL);
    PointsItem.SetDsmBuffer(this->DsmBuffer);
    if(PointsItem.Update() == XDMF_FAIL) return(XDMF_FAIL);
    TmpArray = PointsItem.GetArray();
    if( TmpArray ){
        if( !points ){
            switch( this->GeometryType ){
                case XDMF_GEOMETRY_X_Y_Z :
                    points = new XdmfArray;
                    points->CopyType( TmpArray );
                    points->SetNumberOfElements( TmpArray->GetNumberOfElements() * 3 );
                    break;
                case XDMF_GEOMETRY_XY :
                    points = new XdmfArray;
                    points->CopyType( TmpArray );
                    points->SetNumberOfElements( TmpArray->GetNumberOfElements() / 2 * 3 );
                    break;
                case XDMF_GEOMETRY_X_Y :
                    points = new XdmfArray;
                    points->CopyType( TmpArray );
                    points->SetNumberOfElements( TmpArray->GetNumberOfElements() * 3 );
                    break;
                default :
                    points = TmpArray;
                    // Assure DataItem Destructor does not delete XdmfArray
                    PointsItem.SetArrayIsMine(0);
                    break;
            }
        }
        // We Have made Points Rank = 1 if not XYZ
        switch( this->GeometryType ){
            case XDMF_GEOMETRY_X_Y_Z :
                    Start[0] = ArrayIndex;
                    Stride[0] = 3;
                    points->SelectHyperSlab( Start, Stride, NULL );
                    CopyArray( TmpArray, points);
                    this->NumberOfPoints = TmpArray->GetNumberOfElements();
                    break;
            case XDMF_GEOMETRY_XY :
                    Start[0] = TmpArray->GetNumberOfElements();
                    Start[1] = 3;
                    points->SetShape( 2 , Start );
                    Stride[0] = 1;
                    Stride[0] = 1;
                    Count[0] = TmpArray->GetNumberOfElements();
                    Count[1] = 2;
                    points->SelectHyperSlab( NULL, Stride, Count);
                    CopyArray( TmpArray, points);
                    this->NumberOfPoints = TmpArray->GetNumberOfElements() / 2 ;
                    break;
            case XDMF_GEOMETRY_X_Y :
                    Start[0] = ArrayIndex;
                    Stride[0] = 3;
                    points->SelectHyperSlab( Start, Stride, NULL );
                    CopyArray( TmpArray, points);
                    this->NumberOfPoints = TmpArray->GetNumberOfElements();
                    break;
            default :
                    // points = TmpArray so do nothing
                    this->NumberOfPoints = TmpArray->GetNumberOfElements() / 3;
                    break;
        }
    }
  } 
  ArrayIndex++;
 } while( ( ArrayIndex < 3 ) && ( PointsElement != NULL ) );
} else {
  if( this->GeometryType == XDMF_GEOMETRY_ORIGIN_DXDYDZ ) {
      XdmfDataItem PointsItem;
      PointsItem.SetDOM(this->DOM);
      PointsItem.SetDsmBuffer(this->DsmBuffer);
      XdmfDebug("Reading Origin and Dx, Dy, Dz" );
      PointsElement = this->DOM->FindDataElement(0, this->Element );
      if( PointsElement ){
        if(PointsItem.SetElement(PointsElement, 0) == XDMF_FAIL) return(XDMF_FAIL);
        if(PointsItem.UpdateInformation() == XDMF_FAIL) return(XDMF_FAIL);
        if(PointsItem.Update() == XDMF_FAIL) return(XDMF_FAIL);
        TmpArray = PointsItem.GetArray();
        if( TmpArray ){
            TmpArray->GetValues( 0, this->Origin, 3 );
        }
      PointsElement = this->DOM->FindDataElement(1, this->Element );
      if( PointsElement ){
        if(PointsItem.SetElement(PointsElement, 0) == XDMF_FAIL) return(XDMF_FAIL);
        if(PointsItem.UpdateInformation() == XDMF_FAIL) return(XDMF_FAIL);
        if(PointsItem.Update() == XDMF_FAIL) return(XDMF_FAIL);
        TmpArray = PointsItem.GetArray();
        if( TmpArray ){
          TmpArray->GetValues( 0, this->DxDyDz, 3 );
        }
      } else {
        XdmfErrorMessage("No Dx, Dy, Dz Specified");
        return( XDMF_FAIL );
      }
    } else {
      XdmfErrorMessage("No Origin Specified");
      return( XDMF_FAIL );
    }
  } else if( this->GeometryType == XDMF_GEOMETRY_VXVYVZ ) {
      XdmfDebug("Reading VectorX, VectorY, VectorZ" );
      PointsElement = this->DOM->FindDataElement(0, this->Element );
      if( PointsElement ){
        XdmfDataItem PointsItem;
        PointsItem.SetDOM(this->DOM);
        if(PointsItem.SetElement(PointsElement, 0) == XDMF_FAIL) return(XDMF_FAIL);
        if(PointsItem.UpdateInformation() == XDMF_FAIL) return(XDMF_FAIL);
        if(PointsItem.Update() == XDMF_FAIL) return(XDMF_FAIL);
        TmpArray = PointsItem.GetArray();
        if(!TmpArray){
            XdmfErrorMessage("Error Reading Points X Vector");
            return(XDMF_FAIL);
        }
        this->VectorX = TmpArray;
        PointsItem.SetArrayIsMine(0);
    } else {
      XdmfErrorMessage("No VectorX Specified");
      return( XDMF_FAIL );
      }
      PointsElement = this->DOM->FindDataElement(1, this->Element );
      if( PointsElement ){
        XdmfDataItem PointsItem;
        PointsItem.SetDOM(this->DOM);
        if(PointsItem.SetElement(PointsElement, 0) == XDMF_FAIL) return(XDMF_FAIL);
        if(PointsItem.UpdateInformation() == XDMF_FAIL) return(XDMF_FAIL);
        if(PointsItem.Update() == XDMF_FAIL) return(XDMF_FAIL);
        TmpArray = PointsItem.GetArray();
        if(!TmpArray){
            XdmfErrorMessage("Error Reading Points Y Vector");
            return(XDMF_FAIL);
        }
        this->VectorY = TmpArray;
        PointsItem.SetArrayIsMine(0);
    } else {
      XdmfErrorMessage("No VectorY Specified");
      return( XDMF_FAIL );
      }
      PointsElement = this->DOM->FindDataElement(2, this->Element );
      if( PointsElement ){
        XdmfDataItem PointsItem;
        PointsItem.SetDOM(this->DOM);
        if(PointsItem.SetElement(PointsElement, 0) == XDMF_FAIL) return(XDMF_FAIL);
        if(PointsItem.UpdateInformation() == XDMF_FAIL) return(XDMF_FAIL);
        if(PointsItem.Update() == XDMF_FAIL) return(XDMF_FAIL);
        TmpArray = PointsItem.GetArray();
        if(!TmpArray){
            XdmfErrorMessage("Error Reading Points Z Vector");
            return(XDMF_FAIL);
        }
        this->VectorZ = TmpArray;
        PointsItem.SetArrayIsMine(0);
    } else {
      XdmfErrorMessage("No VectorZ Specified");
      return( XDMF_FAIL );
      }
  }
}
if( points ){
    this->SetPoints( points );
    this->PointsAreMine = 1;
}
return( XDMF_SUCCESS );
}

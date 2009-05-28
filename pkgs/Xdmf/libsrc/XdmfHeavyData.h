/*******************************************************************/
/*                               XDMF                              */
/*                   eXtensible Data Model and Format              */
/*                                                                 */
/*  Id : $Id: XdmfHeavyData.h,v 1.2 2009-01-23 20:31:39 clarke Exp $  */
/*  Date : $Date: 2009-01-23 20:31:39 $ */
/*  Version : $Revision: 1.2 $ */
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
#ifndef __XdmfHeavyData_h
#define __XdmfHeavyData_h

#include "XdmfArray.h"

//! Container class for Heavy Data Access
/*!
This is an abstract convenience object for reading and writing
HeavyData Files. 
Datasets in HeavyDat are specified by :
\verbatim
  Domain:Filename:Pathname
where
  Domain = NDGM | FILE | CORE | GASS
    if Domain is not specified,
    FILE is assumed
  Filename = UNIX style Pathname of HeavyDat file
  Pathname = HeavyData Pathname inside HeavyData File
\endverbatim
*/
class XDMF_EXPORT XdmfHeavyData : public XdmfDataDesc {

public:
  XdmfHeavyData();
  ~XdmfHeavyData();

  XdmfConstString GetClassName() { return ( "XdmfHeavyData" ) ; };

//! Get the default NDGM Host for NDGM:File:/Dataset
        XdmfGetValueMacro(NdgmHost, XdmfConstString);
//! Set the default NDGM Host for NDGM:File:/Dataset
        void SetNdgmHost( XdmfConstString String ) { strcpy( this->NdgmHost, String ); }

//! Get the default Pathname for File:/Dataset
        XdmfGetValueMacro(WorkingDirectory, XdmfConstString);
//! Set the default Pathname for File:/Dataset
        void SetWorkingDirectory( XdmfConstString String );


//! Get the current domain
  XdmfGetValueMacro(Domain, XdmfConstString);
//! Set the current domain
  void SetDomain( XdmfConstString domain ) {
    strcpy( this->Domain, domain );
    } ;

//! Get the current filename
  XdmfGetValueMacro(FileName, XdmfConstString);
//! Set the current filename
  void SetFileName( XdmfConstString File );

//! Get the current HeavyData Dataset path
  XdmfGetValueMacro(Path, XdmfConstString);
//! Set the current HeavyData Dataset path
  void SetPath( XdmfConstString path ) {
    strcpy( this->Path, path );
    } ;

/*!
Get the current read/write access
values can be :
  "r"
  "w"
  "rw"
*/
  XdmfGetValueMacro(Access, XdmfConstString);
//! Set the access permissions
  void SetAccess( XdmfConstString access ) {
    strcpy( this->Access, access );
    } ;

protected:

  char    NdgmHost[XDMF_MAX_STRING_LENGTH];
  XdmfString WorkingDirectory;
  char    Access[XDMF_MAX_STRING_LENGTH];
  char    Domain[XDMF_MAX_STRING_LENGTH];
  XdmfString FileName;
  char    Path[XDMF_MAX_STRING_LENGTH];
};

/*
extern "C" {
extern XdmfString XdmfGetNdgmEntries( void );
extern void XdmfDeleteAllNdgmEntries( void );
extern XdmfInt64 XdmfAddNdgmEntry( XdmfString Name, XdmfInt64 Length );
  }
*/
#endif // __XdmfHeavyData_h

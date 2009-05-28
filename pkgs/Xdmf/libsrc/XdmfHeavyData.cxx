/*******************************************************************/
/*                               XDMF                              */
/*                   eXtensible Data Model and Format              */
/*                                                                 */
/*  Id : $Id: XdmfHeavyData.cxx,v 1.2 2009-01-23 20:31:39 clarke Exp $  */
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
#include "XdmfHeavyData.h"

XdmfHeavyData::XdmfHeavyData() {

  // Defaults
  this->SetDomain( "FILE" );
  this->FileName = 0;
  this->SetFileName( "XdmfHeavyData.dod" );
  this->SetPath( "/" );
  this->SetAccess( "r" );
  this->SetNdgmHost("");
  this->WorkingDirectory  = 0;
  this->SetWorkingDirectory("");

}

XdmfHeavyData::~XdmfHeavyData() {
  this->SetWorkingDirectory(0);
	this->SetFileName(0);
}

void XdmfHeavyData::SetWorkingDirectory( XdmfConstString String )
{
  if ( String == this->WorkingDirectory )
    {
    return;
    }
  if ( String && this->WorkingDirectory && strcmp(String, this->WorkingDirectory) == 0 )
    {
    return;
    }
  if ( this->WorkingDirectory )
    {
    delete [] this->WorkingDirectory;
    this->WorkingDirectory = 0;
    }
  if ( String )
    {
    this->WorkingDirectory = new char [ strlen(String) + 1 ];
    strcpy(this->WorkingDirectory, String);
    }
}

void XdmfHeavyData::SetFileName( XdmfConstString String )
{
  if ( String == this->FileName )
    {
    return;
    }
  if ( String && this->FileName && strcmp(String, this->FileName) == 0 )
    {
    return;
    }
  if ( this->FileName )
    {
    delete [] this->FileName;
    this->FileName = 0;
    }
  if ( String )
    {
    this->FileName = new char [ strlen(String) + 1 ];
    strcpy(this->FileName, String);
    }
}


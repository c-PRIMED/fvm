#!/bin/env python
#/*******************************************************************/
#/*                               XDMF                              */
#/*                   eXtensible Data Model and Format              */
#/*                                                                 */
#/*  Id : $Id: NdgmLs.py,v 1.2 2009-01-23 20:48:54 clarke Exp $  */
#/*  Date : $Date: 2009-01-23 20:48:54 $ */
#/*  Version : $Revision: 1.2 $ */
#/*                                                                 */
#/*  Author:                                                        */
#/*     Jerry A. Clarke                                             */
#/*     clarke@arl.army.mil                                         */
#/*     US Army Research Laboratory                                 */
#/*     Aberdeen Proving Ground, MD                                 */
#/*                                                                 */
#/*     Copyright @ 2002 US Army Research Laboratory                */
#/*     All Rights Reserved                                         */
#/*     See Copyright.txt or http://www.arl.hpc.mil/ice for details */
#/*                                                                 */
#/*     This software is distributed WITHOUT ANY WARRANTY; without  */
#/*     even the implied warranty of MERCHANTABILITY or FITNESS     */
#/*     FOR A PARTICULAR PURPOSE.  See the above copyright notice   */
#/*     for more information.                                       */
#/*                                                                 */
#/*******************************************************************/


import sys
import string
import Xdmf

class NdgmLs :
	def __init__( self, Hostname = None ) :
		self.BufferLength = 4000
		self.Ndgm = Xdmf.XdmfNDGM()
		self.Ndgm.SetModeToClient()
		self.entries = None
		if Hostname != None :
			self.Ndgm.SetNdgmHost( Hostname )
		
	def Ls( self ) :
		status = self.Ndgm.Open()
		if status <= 0 :
			print "Can't Connect to NDGM Server on " + self.Ndgm.GetNdgmHost()
			return None
		# Get the Length
		Length = self.Ndgm.GetTotalLength()
		# print 'NDGM on %s is %d bytes' % ( self.Ndgm.GetNdgmHost(), Length )
		self.RawString = string.split( Xdmf.XdmfGetNdgmEntries() )
		return self.RawString

	def Format( self, Raw = None) :
		if Raw == None :
			Raw = self.RawString
		NumberOfEntries = len( self.RawString ) / 4
		Index = 0
		self.entries = []
		for i in range( NumberOfEntries ) :
			# Skip -NDGM_ENTRY-
			Index += 1
			Name = self.RawString[ Index ]
			Index += 1
			Start = self.RawString[ Index ]
			Index += 1
			End = self.RawString[ Index ]
			Index += 1
			self.entries += ['%-20s %12d %12d' % (Name, int(Start), int(End))]
		return( self.entries )
			
		

if __name__ == '__main__' :
	argc = len( sys.argv )
	if argc > 1 :
		Hostname = sys.argv[ argc - 1 ]
	else :
		Hostname = None
	n  = NdgmLs( Hostname )
	l = n.Ls()
	if l == None :
		sys.exit( 1 )
	e = n.Format()
	NumberOfEntries = len( e )
	if NumberOfEntries < 1 :
		print '-1 : No Entries in NDGM'
	else :
		for i in range( NumberOfEntries ) :
			print '%2d : %s' % ( i, e[ i ])
		
		



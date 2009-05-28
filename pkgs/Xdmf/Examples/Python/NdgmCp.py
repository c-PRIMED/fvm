#!/bin/env python
#/*******************************************************************/
#/*                               XDMF                              */
#/*                   eXtensible Data Model and Format              */
#/*                                                                 */
#/*  Id : $Id: NdgmCp.py,v 1.2 2009-01-23 20:48:53 clarke Exp $  */
#/*  Date : $Date: 2009-01-23 20:48:53 $ */
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

import os
import sys
import string
import Xdmf
from NdgmLs import *


class NdgmCp :
	def __init__( self, From, To ) :
		self.Host = None
		self.From = From
		self.FromIsFile = 1
		self.To = To
		self.ToIsFile = 1
		l = string.split( From, ':' )
		if len(l) > 1 :
			self.FromIsFile = 0
			host, name = string.split( From, ':' )
			if len( host ) > 0 :
				if string.upper( host ) == 'NDGM' :
					host = None
				self.Host = host
			self.From = name
		l = string.split( To, ':' )
		if len(l) > 1 :
			self.ToIsFile = 0
			host, name = string.split( To, ':' )
			if len( host ) > 0 :
				if string.upper( host ) == 'NDGM' :
					host = None
				self.Host = host
			self.To = name
		if self.To == '.' :
			self.To = self.From
		self.Conn = NdgmLs( self.Host )
		self.Conn.Ls()
		self.Entries = self.Conn.Format()

	def Get( self ) :
		if self.FromIsFile :
			self.FromLength = os.path.getsize( self.From )
		else :
			for entry in self.Entries :
				entry = string.split( entry )
				if self.From == entry[0] :
					self.FromStart = int( entry[1] )
					self.FromEnd = int( entry[2] )
					self.FromLength = self.FromEnd - self.FromStart
		print 'Source Data is %d bytes' % self.FromLength

	def Put( self ) :
		if self.ToIsFile :
			if self.FromIsFile :
				Cmd = 'cp ' + self.From + ' ' + self.To
				print Cmd
				os.system( Cmd )
			else :
				From = str(self.FromStart) + ':' + str(self.FromEnd)
				Cmd = 'ice ndgm_cat ' + From + ' ' + self.To
				print Cmd
				os.system( Cmd )
		else :
			Found = 0
			for entry in self.Entries :
				entry = string.split( entry )
				if self.To == entry[0] :
					Found = 1
					self.ToStart = int( entry[1] )
					self.ToEnd = int( entry[2] )
					self.ToLength = self.ToEnd - self.ToStart
					break
			if Found :
				pass
			else :
				Xdmf.XdmfAddNdgmEntry( self.To, self.FromLength )
				self.Conn.Ls()
				self.Entries = self.Conn.Format()
				for entry in self.Entries :
					entry = string.split( entry )
					if self.To == entry[0] :
						Found = 1
						self.ToStart = int( entry[1] )
						self.ToEnd = int( entry[2] )
						self.ToLength = self.ToEnd - self.ToStart
						break
			if self.FromIsFile :
				From = self.From
			else :
				From = str(self.FromStart) + ':' + str(self.FromEnd)
			Cmd = 'ice ndgm_cat ' + From + ' ' + str( self.ToStart )
			print Cmd
			os.system( Cmd )

if __name__ == '__main__' :
	argc = len( sys.argv )
	n  = NdgmCp( sys.argv[ argc - 2 ] , sys.argv[ argc - 1 ] )
	n.Get()
	n.Put()
		
		



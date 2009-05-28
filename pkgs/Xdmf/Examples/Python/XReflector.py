#!/bin/env python

"""Reflection for XDMF and vtk"""

from Xdmf import *
from vtk import *
from libvtkXdmfPython import *
# from libVTKCommonPython import *
# from libVTKGraphicsPython import *
# from libVTKLocalPython import *

import string
import inspect


class XReflector :

	def __init__( self ) :
		self.Lookup = {}


	def RegisterObject(self, ObjectToRegister, Name ) :
		if (ObjectToRegister != None) and (Name != None) :
			self.Lookup[ Name ] = ObjectToRegister
			return( ObjectToRegister )
		return ( None )
	def CreateObject( self, Type, Name ) :

		try :
			NewObject = eval( Type + '()')
			self.RegisterObject( NewObject, Name )
			if Name != None :
				self.Lookup[ Name ] = NewObject
			return( NewObject )
		except :
			return( None )

	def CallMethod( self, ObjectToCall, Method, Args ) :
		# Methods = inspect.getmembers( XdmfDOM, inspect.ismethod)
		if hasattr( ObjectToCall, Method ) :
			RealArgs = ' '
			if Args != None :
				ArgsToTry = string.split( Args, ',' )
				# print 'ArgsToTry = ' + str(ArgsToTry)
				for Arg in ArgsToTry :
					RealArg = self.NameToObject( Arg )
					# print 'RealArg(' + Arg + ') = ' + str( RealArg )
					if RealArg == None :
						RealArgs += Arg + ','
					else :
						try :
							for item in RealArg :
								RealArgs += str(item) + ','
						except :
							RealArgs += 'self.NameToObject("' + Arg + '")'+ ','
					# print 'RealArgs = ' + RealArgs
			RealArgs = RealArgs[:-1]
			try :
				# print 'Calling ... ObjectToCall.' + Method + '(' + RealArgs + ')'
				Result = eval('ObjectToCall.' + Method + '(' + RealArgs + ')' )
				# print 'Back From Call'
				if Result == None :
					Result = 1
				return( Result )
			except :
				return( None )
		else :
			return None


	def NameToObject( self, Name ) :

		if self.Lookup.has_key( Name ) :
			return( self.Lookup[ Name ] )
		return( None )

	def ObjectToName ( self, Object ) :
		for Name in self.Lookup.keys() :
			if self.Lookup[ Name ] == Object :
				return( Name )
		return( None )

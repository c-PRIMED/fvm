#!/bin/env python

"""Reader for XDMF Grids"""

from Xdmf import *
from vtk import *
from libvtkXdmfPython import *

import string


class XdmfReader :

	def __init__( self, DOM, Node ) :
	
		self.DOM = DOM
		self.Node = Node
		self.GridDOM = None
		self.GridNode = None
		self.Reader = vtkXdmfReader()
#		self.Reader.DebugOn()
		# Set the Input XML
		FileName = DOM.Get( Node, "FileName" )
		if FileName == None :
			FileName = DOM.Get( Node, "Input" )
		if FileName == None :
			FileName = DOM.Get( Node, "InputFileName" )
		if FileName == None :
			FileName = DOM.Get( Node, "File" )
		if FileName == None :
			print 'No FileName set'
			return
		print 'SetInputFileName(' + FileName + ')'
		self.Reader.SetInputFileName( FileName )
		# Which Grid ?
		GridName = DOM.Get( Node, "Grid" )
		if GridName == None :
			GridName = DOM.Get( Node, "GridName" )
		if GridName :
			self.Reader.SetGridName(GridName)
		# Parse the XML
		self.Reader.Initialize()
		# Set Up Attributes
		Nparam = DOM.FindNumberOfElements( 'Parameter', Node )
		if Nparam > 0 :
#			print 'Setting %d Parameters' % Nparam
			DOMHandle = self.Reader.GetXdmfDOMHandle()
			RDOM = HandleToXdmfDOM(DOMHandle)
			for i in range( Nparam) :
				PNode = DOM.FindElement( 'Parameter', i, Node )
				Name = DOM.Get( PNode, "Name" )
				Value = DOM.Get( PNode, "CurrentIndex" )
				if not Value :
					Value = DOM.Get( PNode, "Value" )
				pnode = 1
				j = 0
				while pnode :
					pnode = RDOM.FindElement('Parameter', j, None)
					if pnode :
						pName = RDOM.Get(pnode, 'Name')
						if pName == Name :
							RDOM.Set(pnode, 'CurrentIndex', Value)
						j += 1
						pnode = RDOM.FindElement('Parameter', j, None)
		# Set Up Attributes
		Nattr = DOM.FindNumberOfElements( 'Attribute', Node )
		print 'Reading %d Attributes' % Nattr
		if Nattr > 0 :
			self.Reader.SetAllAttributeStatusOff()
			for i in range( Nattr ) :
				ANode = DOM.FindElement( 'Attribute', i, Node )
				Name = DOM.Get( ANode, "Name" )
				self.Reader.SetAttributeStatusOn( Name )

	def Update(self) :
		print 'Updating Grid'
		self.Reader.Update()
		return( self.Reader )

	def GetOutput(self) :
		return( self.Reader.GetOutput() )


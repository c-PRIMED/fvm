#!/bin/env python

"""XVis for Python"""

import sys
import getopt
import Xdmf

from XdmfReader import *
from XReflector import *

import string
import math


class XVis :

	def __init__( self ) :
		self.ParentObject = None
		self.Reflector = XReflector()
		self.loop = 1
		self.NdgmCheck = None
		self.Picker = None
		self.Mesa = 0
		self.LastX = 0.0
		self.LastY = 0.0
		self.LastZ = 0.0

	def PickMethod ( self ) :
		Picker = self.Picker
		PickMeasure = self.PickMeasure
		PickViewFile = self.PickViewFile
		if Picker.GetCellId < 0 :
			print 'No Cell Picked'
			return
		else :
			SelectionPoint = Picker.GetSelectionPoint()
			PickPosition = Picker.GetPickPosition()
			x = PickPosition[0]
			y = PickPosition[1]
			z = PickPosition[2]
			print 'X Y Z : %f %f %f' % ( x, y, z )
			dx = x - self.LastX
			dy = y - self.LastY
			dz = z - self.LastZ
			Length = math.sqrt( ( dx * dx ) + ( dy * dy ) + ( dz * dz))
			print 'Dx Dy Dz : %f %f %f' % ( dx, dy, dz)
			print 'Length : %f'  % Length
			self.LastX = x
			self.LastY = y
			self.LastZ = z
		if PickViewFile :
			print 'Writing View in ' + PickViewFile
			Window = Picker.GetRenderer().GetRenderWindow()
			WindowSize = Window.GetSize()
			# print 'Window Size = ' + str( WindowSize )
			WindowPosition = Window.GetPosition()
			# print 'Window Position = ' + str( WindowPosition )
			Camera = Picker.GetRenderer().GetActiveCamera()
			Position = Camera.GetPosition()
			# print 'Camera Position = ' + str( Position )
			FocalPoint = Camera.GetFocalPoint()
			# print 'Camera FocalPoint = ' + str( FocalPoint )
			ViewUp = Camera.GetViewUp()
			# print 'Camera ViewUp = ' + str( ViewUp )
			ViewFd = open( PickViewFile, "w" )
			ViewFd.write('<!--\n')	
			ViewFd.write('<Reference Target="RenWin" Method="SetPosition" Args="%d, %d"/>\n' % ( WindowPosition[0], WindowPosition[1]))
			ViewFd.write('<Reference Target="RenWin" Method="SetSize" Args="%d, %d"/>\n' % ( WindowSize[0], WindowSize[1]))
			ViewFd.write('<Reference Name="Camera" Target="Ren" Method="GetActiveCamera"/>\n')
			ViewFd.write('<Reference Target="Camera" Method="SetPosition" Args="%f, %f, %f"/>\n' % ( Position[0], Position[1], Position[2]))
			ViewFd.write('<Reference Target="Camera" Method="SetFocalPoint" Args="%f, %f, %f"/>\n' % ( FocalPoint[0], FocalPoint[1], FocalPoint[2]))
			ViewFd.write('<Reference Target="Camera" Method="SetViewUp" Args="%f, %f, %f"/>\n' % ( ViewUp[0], ViewUp[1], ViewUp[2]))
			ViewFd.write('<Reference Target="Ren" Method="ResetCameraClippingRange" />\n')
			ViewFd.write('-->\n')	
			ViewFd.write('<View>\n')
			ViewFd.write('<Window\n')
			ViewFd.write('\tSize="%d, %d"\n' % ( WindowSize[0], WindowSize[1]))
			ViewFd.write('\tPosition="%d, %d">\n' % ( WindowPosition[0], WindowPosition[1]))
			ViewFd.write('</Window>\n')
			ViewFd.write('<Camera\n')
			ViewFd.write('\tFocalPoint="%f, %f, %f"\n' % ( FocalPoint[0], FocalPoint[1], FocalPoint[2]))
			ViewFd.write('\tViewUp="%f, %f, %f"\n' % ( ViewUp[0], ViewUp[1], ViewUp[2]))
			ViewFd.write('\tPosition="%f, %f, %f">\n' % ( Position[0], Position[1], Position[2]))
			ViewFd.write('</Camera>\n')
			ViewFd.write('</View>\n')
			ViewFd.close()
	def Loop ( self, Node ) :
		DOM = self.DOM
		iRenName = DOM.Get( Node, "Interactor" )
		iRen = self.Reflector.NameToObject( iRenName )
		if iRen != None :
			iRen.Initialize()
		Continue = 1
		while Continue :
			if iRen != None :
				if hasattr( iRen, 'LoopOnce') :
					iRen.LoopOnce()
				if iRen.GetBreakLoopFlag() != 0 :
					Continue = 0
					break
			for i in range( DOM.GetNumberOfChildren( Node )) :
				ChildNode = DOM.GetChild( i, Node )
				ChildType = DOM.Get( ChildNode, "NodeType" )
				self.Visit( ChildNode )
	def Visit( self, Node ) :
		DOM = self.DOM
		Reflector = self.Reflector
		Type = DOM.Get( Node, "NodeType" )
		# print 'Visit ' + Type
		if Type == "Xvis" :
			return( 1 )
		elif Type == "Ndgm" :
			if self.NdgmCheck == None :
				self.NdgmCheck = Xdmf.XdmfNDGM()
				self.NdgmCheck.SetModeToClient()
				Status = self.NdgmCheck.Open()
				if Status <= 0 :
					print "Can't Connect to Ndgm"
					return None
			Ndgm = self.NdgmCheck
			BarrierNumber = DOM.Get( Node, "Barrier" )
			if BarrierNumber == None :
				BarrierNumber = 20
			else :
				BarrierNumber = int( BarrierNumber )
			CurrentValue = Ndgm.BarrierPoll( BarrierNumber )
			BarrierValue = DOM.Get( Node, "Value" )
			if BarrierValue == None :
				BarrierValue = CurrentValue
				DOM.Set(Node, "Value", str( CurrentValue )  )
			else :
				BarrierValue = int( BarrierValue )
			# print 'NDGM Current Value od %d = %d' % ( BarrierNumber, CurrentValue)
			if CurrentValue != BarrierValue :
				print 'Barrier Update #%d' % CurrentValue
				for i in range( DOM.GetNumberOfChildren( Node )) :
					ChildNode = DOM.GetChild( i, Node )
					ChildType = DOM.Get( ChildNode, "NodeType" )
					self.Visit( ChildNode )
				CurrentValue = Ndgm.BarrierPoll( BarrierNumber )
				# Avoid Multiple Update Locks
				DOM.Set(Node, "Value", str( CurrentValue ) )
			return( CurrentValue )
		elif Type == "Loop" :
			Status = self.Loop( Node )
			return( Status  )
		# elif Type == "Invoke" :
		# 	Name = DOM.Get(Node, "Target")
		# 	Method = DOM.Get(Node, "Method")
		# 	Args = DOM.Get(Node, "Args")
		# 	NewObject = Reflector.NameToObject( Name )
		# 	Status = Reflector.CallMethod( NewObject, Method, Args )
		# 	return( Status )
		elif Type == "Execute" :
			CData = string.strip(DOM.Get( Node, "CData" ))
			try :
				exec( CData )
				return( 1 )
			except :
				return( None )
		elif Type == "Print" :
			TargetName = DOM.Get(Node, "Target")
			NewLine = DOM.Get(Node, "NewLine" )
			TargetObject = Reflector.NameToObject( TargetName )
			if (NewLine != None ) and (string.upper( NewLine ) == "FALSE") :
				if TargetObject != None :
					print str( TargetObject ) ,
				else :
					# print "Target == None"
					print TargetName ,
			else :
				if TargetObject != None :
					print str( TargetObject )
				else :
					# print "Target == None"
					print TargetName
			return( TargetObject )
		elif (Type == "Reference") or (Type == "Invoke") :
			Name = DOM.Get( Node, "Name" )
			TargetName = DOM.Get(Node, "Target")
			Method = DOM.Get(Node, "Method")
			Args = DOM.Get(Node, "Args")
			if TargetName == None :
				TargetObject = self.ParentObject
			else :
				TargetObject = Reflector.NameToObject( TargetName )
			if Method == None :
				NewObject = TargetObject
			else :
				# print 'Calling Method ' + Method + ' on ' + str( TargetObject )
				NewObject = Reflector.CallMethod( TargetObject, Method, Args )
			# Name == None is OK
			Reflector.RegisterObject( NewObject, Name )
			return( NewObject )
		elif Type == "XdmfReader" :
			Name = DOM.Get(Node, "Name")
			NewObject = XdmfReader(DOM, Node)
			Reflector.RegisterObject( NewObject, Name )
			return( NewObject )
		else :
			Name = DOM.Get( Node, "Name" )
			if self.Mesa :
				if Type == 'PolyDataMapper' :
					Type = 'MesaPolyDataMapper'
				if Type == 'Actor' :
					Type = 'MesaActor'
				if Type == 'Renderer' :
					Type = 'MesaRenderer'
				if Type == 'RenderWindow' :
					Type = 'XMesaRenderWindow'
			# print 'Creating a ' + Type
			NewObject = Reflector.CreateObject( 'vtk' + Type, Name)
			if NewObject == None :
				NewObject = Reflector.CreateObject( Type, Name)
				if NewObject == None :
					# print "Can't Create " + Type
					return(None)
			if (Type == 'XdmfRenderWindowInteractor') or (Type == 'RenderWindowInteractor') :
				# Create a Picker
				self.Picker = vtkCellPicker()
				self.Picker.SetPickMethod( self.PickMethod )
				self.PickMeasure = DOM.Get( Node, "Measure" )
				self.PickViewFile = DOM.Get( Node, "ViewFile")
				NewObject.SetPicker( self.Picker )
		for i in range( DOM.GetNumberOfAttributes( Node )) :
			Attribute = DOM.GetAttribute( Node, i )
			if Attribute == 'Name' :
				pass
			else :
				Value = DOM.Get( Node, Attribute )
				# print 'Try ' + Attribute + ' = ' + Value
				Status = Reflector.CallMethod( NewObject, Attribute, Value)
				if Status == None :
					Status = Reflector.CallMethod( NewObject, 'Set' + Attribute, Value)
		return( NewObject )
				

	def Traverse( self, Node ) :

		DOM = self.DOM
		Reflector = self.Reflector
		Parent  = self.Visit( Node )
		self.ParentObject = Parent
		if Parent == None :
			return None	
		for i in range( DOM.GetNumberOfChildren( Node )) :
			ChildNode = DOM.GetChild( i, Node )
			ChildType = DOM.Get( ChildNode, "NodeType" )
			if ChildType == 'PolyDataMapper' :
				ChildType = 'Mapper'
			self.ParentObject = Parent
			Child = self.Traverse( ChildNode )
			# print 'Child = '  + str( Child )
			if ( Child != None ) and (ChildType != 'Reference')  :
				# Try Add
				# print '... ChildType = ' + str( ChildType )
				# print '.. Child = ' + str( Child )
				# print '... try ' + 'Add' + ChildType
				Attach = DOM.Get( ChildNode, "AttachAs" )
				if Attach == None :
					if hasattr( Parent, 'Add' + ChildType ) :
						# print 'Calling Add' + ChildType
						eval('Parent.Add' + ChildType + '( Child )')
					elif hasattr( Parent, 'Set' + ChildType ) :
						# print 'Calling Set' + ChildType
						eval('Parent.Set' + ChildType + '( Child )')
					elif hasattr( Child, 'GetOutput' ) :
						# print 'Calling SetInput'
						if hasattr( Parent, 'SetInput') :
							Parent.SetInput( Child.GetOutput())
				else :
					if hasattr( Child, 'GetOutput' ) :
						Attachment = Child.GetOutput()
					else :
						Attachment = Child
					if Attach == 'Actor2D' :
						eval('Parent.AddActor2D' + '( Attachment )')
					else :
						eval('Parent.Set' + Attach + '( Attachment )')
		return( Parent )

	def usage ( self, argv ) :
		print 'Usage ' + argv[0] + ' ' + \
			str( self.XVisOptions )

	XVisOptions = ["help", "loop", "mesa", "replace=", "xml=", "input=" ]

	def main( self, argv ) :
		argc = len( argv )
		try :
			# Only Accept Long Options
			opts, args = getopt.getopt(argv[1:],
					'', self.XVisOptions )

		except getopt.GetoptError:
			self.usage(argv)
			return(-1)
		# print 'opts = ' + str( opts )
		# print 'args = ' + str( args )

		if len(args) != 0 :
			self.usage(argv)
			return(-1)

		DOM = XdmfDOM()
		self.DOM = DOM
		self.defines = {}
		self.txt = None
		# self.InputFile = args[0]
		for o, a in opts:
			if o in ("-h", "--help"):
				self.usage(argv)
			if o == "--mesa" :
				self.Mesa = 1
			if o == "--loop" :
				self.loop = 1
			if o == "--xml" :
				self.DOM.Parse( a )
			if o == "--input" :
				a = os.path.abspath( a )
				dirname = os.path.dirname( a )
				fd = open(a, 'r')
				self.txt = fd.read()
				fd.close()
			if o == "--replace" : 
				Key, Value = string.split( a, ',')
				self.defines[ Key ] = Value
				
		# self.DOM.SetInputFileName( a )
		for key in self.defines.keys() :
			self.txt = string.replace( self.txt, key, self.defines[ key ] )
		# print 'txt = ' + self.txt
		self.DOM.Parse( self.txt )
		Node = self.DOM.FindElement( 'Xvis' )
		if Node == None :
			print "Can't Find XVis Node"
			return( -1 )
		self.Traverse( Node )
		return(0)

if __name__ == '__main__' :
	argv = sys.argv
	Xv = XVis()
	Xv.main( argv )

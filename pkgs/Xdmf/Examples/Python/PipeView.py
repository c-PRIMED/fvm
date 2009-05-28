#!/usr/bin/env python
#/*******************************************************************/
#/*                               XDMF                              */
#/*                   eXtensible Data Model and Format              */
#/*                                                                 */
#/*  Id : $Id: PipeView.py,v 1.2 2009-01-23 20:48:54 clarke Exp $  */
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

import getopt
import sys
import string


print 'Loading Xdmf'
import Xdmf
print 'Loading vtk'
from vtk import *
from libvtkXdmfPython import *
import vtkRenderWidget

import Tkinter
# from vtkPipeline import *
import vtkPipeline.vtkPipeline

class ViewAll:

	def __init__( self ):
		self.StartGrid = 0
		self.EndGrid = -1
		self.Attributes = []

	def FindGrids( self, FileName ) :
		DOM = Xdmf.XdmfDOM()
		DOM.SetInputFileName( FileName )
		DOM.Parse()
		self.NumberOfGrids = DOM.FindNumberOfElements( 'Grid' )

	def Done(self) :
		sys.exit()

	def SetUpPipeline( self, Ren ) :
		root = Tkinter.Tk()
		root.title("XDMF Viewer")
		wid = vtkRenderWidget.vtkTkRenderWidget (root, width=500, height=500)
		wid.pack (expand='true', fill='both')
		wid.bind ("<KeyPress-q>", 
			lambda e=None: self.Done() )
		renWin = wid.GetRenderWindow()
		renWin.AddRenderer(Ren)
		renWin.SetSize(500,500)
		renWin.Render ()

		pipe = vtkPipeline.vtkPipeline.vtkPipelineBrowser (root, renWin)
		pipe.browse ()
		root.mainloop ()

	def View ( self, FileName ):
		print 'Parsing ' + FileName

		if( self.EndGrid <= 0 ) :
			self.FindGrids( FileName )
			self.EndGrid = self.NumberOfGrids - 1
		
		Ren = vtkRenderer()
		for GridIndex in range( self.StartGrid, self.EndGrid + 1 ) :
			Reader = vtkXdmfReader()
			# The XML File is Input
			Reader.SetInputFileName( FileName )
			Reader.SetGridIndex( GridIndex )
			# Parse the XML but don't
			# yet read the Heavy Data (HDF5)
			# This is necessary so that vtk knows
			# the topology of the data since it
			# could be structured or unstructured
			Reader.Initialize()
			# Read XDMF Attributes
			Reader.SetAllAttributeStatusOff()
			Index = GridIndex - self.StartGrid
			if Index < len( self.Attributes ) :
				AttrIndex = self.Attributes[ Index ]
				print 'Setting Attribute %d On' % AttrIndex
				Reader.SetAttributeStatusOn( AttrIndex )
			
	
			GeometryFilter = vtkGeometryFilter()
			GeometryFilter.SetInput( Reader.GetOutput() )
			# Reader.Update() Triggers the Reading of HDF5
			#  this is triggered by the Viz Pipeline
			GeometryFilter.Update()


			# We can access the underlying XdmfGrid
			# if necessary
			Gridptr = Reader.GetXdmfGridHandle()
			Grid = Xdmf.HandleToXdmfGrid( Gridptr )
			print 'Grid %d has %d ' % (GridIndex, Grid.GetNumberOfElements()) + \
				Grid.GetTopologyTypeAsString()
			Nattr = Grid.GetNumberOfAttributes()
			print '<XdmfGrid> has %d Attributes' % Nattr
			for i in range( Nattr ) :
				Attribute = Grid.GetAttribute( i )
				print '\tAttribute #%d' % i
				print '\t\tName: ' + Attribute.GetName()
				print '\t\tCenter: ' + Attribute.GetAttributeCenterAsString()
				print '\t\tType: ' + Attribute.GetAttributeTypeAsString()
	
			# List All of the Available Arrays
			PointData = Reader.GetOutput().GetPointData()
			for i in range( PointData.GetNumberOfArrays() ) :
				Array = PointData.GetArray( i )
				print '\tArray #%d ' % i
				print '\t\tName: ' + Array.GetName()
				Min, Max = Array.GetRange()
				print '\t\tRange: %f -> %f' % (Min, Max)
	
	
	


			# Use the Third Attribute (Z Coordinate) for Color
			# GeometryFilter.GetInput().GetPointData().SetActiveScalars('Z Coordinate')

			Mapper = vtkPolyDataMapper()
			Mapper.SetInput( GeometryFilter.GetOutput() )
			ScalarRange = GeometryFilter.GetOutput().GetScalarRange()
			Mapper.SetScalarRange( ScalarRange )
			# Blue to Red
			Mapper.GetLookupTable().SetHueRange( .667, 0.0 )

			Actor = vtkActor()
			Actor.SetMapper( Mapper )
			# Actor.GetProperty().SetRepresentationToWireframe()
			Ren.AddActor( Actor )


		Ren.SetBackground(.2, .2, .2)

		self.SetUpPipeline(Ren)
		# RenWin = vtkRenderWindow()
		# RenWin.SetSize(500,500)
		# RenWin.AddRenderer(Ren)

		# iRen = vtkRenderWindowInteractor()
		# iRen.SetRenderWindow(RenWin)

		# Interact with the Mouse and Keyboard
		# iRen.Initialize()
		# iRen.Start()

	def usage (self) :
		print 'Options : --start=Grid# --end=Grid# --attribute=Index --attribute=Index ... File.xmf'
		sys.exit(0)

	def Options( self, opts ) :
		argc = len( opts )
		print '%d Args = ' % argc + str( opts )
		try :
			opts, args = getopt.getopt(opts,
					"aseh:",
					["help", "start=", "end=", "attribute=" ])

		except getopt.GetoptError:
			self.usage()
			sys.exit(2)
		print 'opts = ' + str( opts )
		print 'args = ' + str( args )
		output = None
		for o, a in opts:
			if o in ("-h", "--help"):
				self.usage()
				sys.exit()
			if o in ("-s", "--start"):
				self.StartGrid = int(a)
			if o in ("-e", "--end"):
				self.EndGrid = int(a)
			if o in ("-a", "--attribute"):
				print 'Appending ' + a
				self.Attributes.append( int(a) )
		print 'StartGrid %d' % self.StartGrid
		print 'EndGrid %d' % self.EndGrid
		print 'Attributes = ' + str( self.Attributes )
if __name__ == '__main__' :
	argc = len( sys.argv )
	viewer = ViewAll()
	viewer.Options( sys.argv[1:] )
	viewer.View( sys.argv[ argc - 1 ] )

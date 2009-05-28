#!/usr/bin/env python
#/*******************************************************************/
#/*                               XDMF                              */
#/*                   eXtensible Data Model and Format              */
#/*                                                                 */
#/*  Id : $Id: FromOBJ.py,v 1.2 2009-01-23 20:48:47 clarke Exp $  */
#/*  Date : $Date: 2009-01-23 20:48:47 $ */
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


print 'Loading vtk'
from libVTKCommonPython import *
from libVTKGraphicsPython import *
from libVTKImagingPython import *
from libVTKPatentedPython import *
from libVTKContribPython import *
from libVTKLocalPython import *

print 'Loading Xdmf'
import Xdmf


class FromOBJ :
	def __init__( self, FileName ) :
		self.Convert = 0
		self.FileName = FileName
		BaseList = string.split(FileName, '.')
		if len( BaseList ) == 1 :
			self.BaseName = BaseList[0]
		else :
			self.BaseName = string.join( BaseList[ : len( BaseList ) - 1 ] )
	def CreateXdmf( self ) :
		ObjReader = vtkOBJReader()
		ObjReader.SetFileName( self.FileName )
		print 'Reading ' + self.FileName
		ObjReader.Update()
		TriFilter = vtkTriangleFilter()
		TriFilter.SetInput( ObjReader.GetOutput() )
		TriFilter.Update()
		if self.Convert == 1 :
			Points = TriFilter.GetOutput().GetPoints()
			print 'Converting %d Points' % Points.GetNumberOfPoints()
			for i in range ( Points.GetNumberOfPoints() ) :
				x, y, z = Points.GetPoint( i )
				x = x * .0254
				y = y * .0254
				z = z * .0254
				Points.SetPoint( i, x, y, z )

		Merge = vtkCleanPolyData()
		Merge.SetTolerance(0)
		Merge.SetInput( TriFilter.GetOutput() )
		Normal = vtkPolyDataNormals()
		Normal.SetInput( Merge.GetOutput() )
		Normal.SetFeatureAngle(0)
		Normal.SplittingOff()
		Normal.ConsistencyOn()
		# Normal.DebugOn()
		Normal.Update()

		Writer = vtkXdmfDataSetWriter()
		Writer.SetInput( Normal.GetOutput() )
		Writer.SetHeavyDataSetName(self.BaseName + '.h5')
		Writer.WriteGrid()
		Writer.WriteAttributes()
		XML = Writer.GetXML()
		fd = open(self.BaseName + '.xml', "w" )
		fd.write("""<?xml version="1.0" ?>\n""")
		fd.write("""<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" [\n""")
		# fd.write("""<!ENTITY HeavyData "Zsu.h5">\n""")
		fd.write('<!ENTITY HeavyData "' + self.BaseName + '.h5">\n')
		fd.write("""]>\n""" )
		fd.write("""<Domain>\n<Grid>\n""")
		fd.write( XML )
		fd.write("""</Grid>\n</Domain>\n""")
		fd.close()

if __name__ == '__main__' :
	argc = len( sys.argv )
	FileName = sys.argv[ argc - 1 ]
	fobj = FromOBJ( FileName )
	if argc > 2 :
		print 'Converting From inches'
		fobj.Convert = 1
	fobj.CreateXdmf()

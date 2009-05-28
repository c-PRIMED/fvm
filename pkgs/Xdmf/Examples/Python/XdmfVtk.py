#!/usr/bin/env python

from vtk import *
from Xdmf import *
from libvtkXdmfPython import *
import string
import sys

platform = string.upper(sys.platform)

if (platform == 'WIN32') or (platform == 'CYGWIN') :
	class vtkXdmfRenderWindowInteractor ( vtkXdmfWin32RenderWindowInteractor ) :
			def __init__(self) :
				pass
else :
	class vtkXdmfRenderWindowInteractor ( vtkXdmfXRenderWindowInteractor ) :
			def __init__(self) :
				pass



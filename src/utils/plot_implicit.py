##/usr/bin/env python

from __future__ import division

import random
import struct
import sys

import numpy as np

#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.pyplot as plt

import vtk

SURF_PX     =  1   # Plane parallel to x-plane 
SURF_PY     =  2   # Plane parallel to y-plane 
SURF_PZ     =  3   # Plane parallel to z-plane 
SURF_PLANE  =  4   # Arbitrary plane
SURF_CYL_X  =  5   # Cylinder along x-axis
SURF_CYL_Y  =  6   # Cylinder along y-axis
SURF_CYL_Z  =  7   # Cylinder along z-axis
SURF_SPHERE =  8   # Sphere
SURF_BOX_X  =  9   # Box extending infinitely in x-direction
SURF_BOX_Y  = 10   # Box extending infinitely in y-direction
SURF_BOX_Z  = 11   # Box extending infinitely in z-direction
SURF_BOX    = 12   # Rectangular prism
SURF_GQ     = 13   # General quadratic surface

class Plot(object):

    def __init__(self):
        self.cells = {}
        self.cellLimits = {}
        self.points = []
        self.triangles = []

    def plot(self):


        print self.cells
        #return

        ren = vtk.vtkRenderer()
        renWin = vtk.vtkRenderWindow()
        renWin.AddRenderer(ren)
        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(renWin)
        ren.SetBackground(1, 1, 1)

        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName("output.vtp")

        for cellId,surfs in self.cells.iteritems():
            quads = []
            for surf in surfs:
                print surf
                quadric = vtk.vtkQuadric()
                a = [0 for ai in range(10)]
                if surf[0] == SURF_SPHERE:
                    x0,y0,z0,r = surf[2]
                    a[0] = 1      # x^2
                    a[1] = 1      # y^2
                    a[2] = 1      # z^2
                    a[3] = 0      # x*y
                    a[4] = 0      # y*z
                    a[5] = 0      # x*z
                    a[6] = -2*x0  # x
                    a[7] = -2*y0  # y
                    a[8] = -2*z0  # z
                    a[9] = x0**2 + y0**2 + z0**2 -r**2  # const
                elif surf[0] == SURF_PX:
                    x0 = surf[2]
                    a[0] = 0      # x^2
                    a[1] = 0      # y^2
                    a[2] = 0      # z^2
                    a[3] = 0      # x*y
                    a[4] = 0      # y*z
                    a[5] = 0      # x*z
                    a[6] = 1      # x
                    a[7] = 0      # y
                    a[8] = 0      # z
                    a[9] = -x0    # const
                elif surf[0] == SURF_PY:
                    y0 = surf[2]
                    a[0] = 0      # x^2
                    a[1] = 0      # y^2
                    a[2] = 0      # z^2
                    a[3] = 0      # x*y
                    a[4] = 0      # y*z
                    a[5] = 0      # x*z
                    a[6] = 0      # x
                    a[7] = 1      # y
                    a[8] = 0      # z
                    a[9] = -y0    # const
                elif surf[0] == SURF_PZ:
                    z0 = surf[2]
                    a[0] = 0      # x^2
                    a[1] = 0      # y^2
                    a[2] = 0      # z^2
                    a[3] = 0      # x*y
                    a[4] = 0      # y*z
                    a[5] = 0      # x*z
                    a[6] = 0      # x
                    a[7] = 0      # y
                    a[8] = 1      # z
                    a[9] = -z0    # const
                elif surf[0] == SURF_CYL_Z :
                    x0,y0,r = surf[2]
                    a[0] = 1      # x^2
                    a[1] = 1      # y^2
                    a[2] = 0      # z^2
                    a[3] = 0      # x*y
                    a[4] = 0      # y*z
                    a[5] = 0      # x*z
                    a[6] = -2*x0  # x
                    a[7] = -2*y0  # y
                    a[8] = 0      # z
                    a[9] = x0**2 + y0**2 - r**2    # const
                else: pass
                # adjust sign of coefficients based on specified sense
                s = surf[1]
                a = [ai*s for ai in a]
                quadric.SetCoefficients(*a)
                quads.append(quadric)

            combinedQuads = vtk.vtkImplicitBoolean()
            combinedQuads.SetOperationTypeToIntersection()
            for q in quads:
                combinedQuads.AddFunction(q)

            sample = vtk.vtkSampleFunction()
            sample.SetSampleDimensions(50, 50, 50)
            xmin,ymin,zmin,xmax,ymax,zmax = self.cellLimits[cellId]
            print self.cellLimits[cellId]
            sample.SetModelBounds(xmin,xmax,ymin,ymax,zmin,zmax)
            sample.SetImplicitFunction(combinedQuads)
            sample.ComputeNormalsOff()

            contour = vtk.vtkContourFilter()
            contour.SetInputConnection(sample.GetOutputPort())
            contour.SetValue(0, 0.0)

            contMapper = vtk.vtkPolyDataMapper()
            contMapper.SetInputConnection(contour.GetOutputPort())

            contActor = vtk.vtkActor()
            contActor.SetMapper(contMapper)

            contActor.GetProperty().SetColor(0,0,0)

            ren.AddActor(contActor)

            writer.SetInput(0,contour.GetOutput())

        writer.Write()

        iren.Initialize()
        renWin.Render()
        iren.Start()





    def load_file(self, filename):
        # Create binary reader for plot.out file
        plotFile = BinaryReader(filename)

        while True:
            try:
                cellId = plotFile.get_int()
                self.cells[cellId] = []
                self.cellLimits[cellId] = plotFile.get_double(6)
                n_surfs = plotFile.get_int()
                for s in range(n_surfs):
                    surfType = plotFile.get_int()
                    sign = plotFile.get_int()
                    n_coeffs = plotFile.get_int()
                    coeffs = plotFile.get_double(n_coeffs)
                    self.cells[cellId].append([surfType,sign,coeffs])
            except BinaryReaderError:
                break


class BinaryReader(object):

    def __init__(self, filename):
        """
        Initialize instance of Record object.
        """

        fh = open(filename, 'rb')
        self.data = fh.read()
        self.numBytes = len(self.data)

        self.reset()
        self.intSize    = struct.calcsize('i')
        self.longSize   = struct.calcsize('q')
        self.floatSize  = struct.calcsize('f')
        self.doubleSize = struct.calcsize('d')

    def get_data(self, n, typeCode, itemSize):
        """
        Returns one or more items of a specified type at the current
        position within the data list. If more than one item is read,
        the items are returned in a list.
        """
        if self.pos >= self.numBytes:
            raise BinaryReaderError("Already read all data from record")
        
        values = struct.unpack('{0}{1}'.format(n,typeCode), self.data[self.pos:self.pos+itemSize*n])
        self.pos += itemSize * n
        if n == 1:
            return values[0]
        else:
            return list(values)
        

    def get_int(self, n=1):
        """
        Returns one or more 4-byte integers.
        """
        return self.get_data(n,'i',self.intSize)
                             
    def get_long(self, n=1):
        """
        Returns one or more 8-byte integers.
        """
        return self.get_data(n,'q',self.longSize)
                             
    def get_float(self, n=1):
        """
        Returns one or more floats.
        """
        return self.get_data(n,'f',self.floatSize)
                             
    def get_double(self, n=1):
        """
        Returns one or more double
        """
        return self.get_data(n,'d',self.doubleSize)
                             
    def get_string(self, length, n=1):
        """
        Returns a string of a specified length starting at the current
        position in the data list.
        """

        if self.pos >= self.numBytes:
            raise BinaryReaderError("Already read all data from record")
        
        relevantData = self.data[self.pos:self.pos+length*n]
        (s,) = struct.unpack('{0}s'.format(length*n), relevantData)
        self.pos += length*n
        if n == 1:
            return s
        else:
            return [s[i*length:(i+1)*length] for i in range(n)]

    def reset(self):
        self.pos = 0


class BinaryReaderError(Exception):
    """Case class for all binary reader errors."""

    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return repr(self.msg)


if __name__ == '__main__':
    if len(sys.argv) >= 1:
        filename = sys.argv[1]

        p = Plot()

        # Load binary plot.out file
        print("Loading plotting file...")
        p.load_file(filename)

        # Display plot
        print("Generating plot...")
        p.plot()

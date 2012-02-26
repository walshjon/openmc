#!/usr/bin/env python

from __future__ import division

import random
import struct
import sys

import numpy as np



import vtk


class Plot(object):

    def __init__(self):
        self.cells = {}


    def plot(self):

        onlyPlotPoints = False
        onlyPlotPoints = True

        """ 
        if onlyPlotPoints:
            #mpl points - can be slow
            from mpl_toolkits.mplot3d import Axes3D
            import matplotlib.pyplot as plt
            fig = plt.figure()
            ax = Axes3D(fig)
            for cell,pts in self.cells.iteritems():
                ax.scatter(*zip(*pts))
            plt.show()
            return
        """

        ren = vtk.vtkRenderer()
        renWin = vtk.vtkRenderWindow()
        renWin.AddRenderer(ren)
        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(renWin)
        ren.SetBackground(1, 1, 1)

        #writer = vtk.vtkXMLUnstructuredGridWriter()
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName("output.vtp")
        writer.SetNumberOfPieces(len(self.cells.keys()))

        apd = vtk.vtkAppendPolyData()

        c = 0
        for cell,pts in self.cells.iteritems():
            print "Plotting {0} points...".format(len(pts))
            points = vtk.vtkPoints()
            ptVerts = vtk.vtkCellArray()
            cellData = vtk.vtkDoubleArray()
            cellData2 = vtk.vtkDoubleArray()
            cellData.SetName("openmc cell")
            cellData2.SetName("data array2")

            for i,(x,y,z) in enumerate(pts):
                points.InsertPoint(i,x, y, z)
                ptVerts.InsertNextCell(1)
                ptVerts.InsertCellPoint(i)
                cellData.InsertNextValue(cell)
                cellData2.InsertNextValue(i)

            pointPolyData = vtk.vtkPolyData()
            pointPolyData.SetPoints(points)
            pointPolyData.GetPointData().SetScalars(cellData)
            pointPolyData.GetPointData().AddArray(cellData2)
            pointPolyData.SetVerts(ptVerts)
            
            if onlyPlotPoints:
        
                ptMask = vtk.vtkMaskPoints()
                ptMask.SetInput(pointPolyData)

                sphere = vtk.vtkSphereSource()
                sphere.SetRadius(0.1)
                #sphere.SetResolution(6)

                # vtkGlyph3D takes two inputs: the input point set (SetInput) which can be
                # any vtkDataSet; and the glyph (SetSource) which must be a vtkPolyData.
                glyph = vtk.vtkGlyph3D()
                glyph.SetInputConnection(ptMask.GetOutputPort())
                glyph.SetSource(sphere.GetOutput())
                #glyph.SetVectorModeToUseNormal()
                #glyph.SetScaleModeToScaleByVector()
                glyph.SetScaleFactor(2.0)
                pointMapper = vtk.vtkPolyDataMapper()
                pointMapper.SetInputConnection(glyph.GetOutputPort())
                pointActor = vtk.vtkActor()
                pointActor.SetMapper(pointMapper)
                pointActor.GetProperty().SetColor(0.0, 0.79, 0.34)

                ren.AddActor(pointActor)
                #apd.AddInput(glyph.GetOutput())
                apd.AddInput(pointPolyData)


            else:
                # Construct the surface and create isosurface
                surf = vtk.vtkSurfaceReconstructionFilter()
                #surf.SetNeighborhoodSize(10)
                #surf.SetSampleSpacing(0.5)
                surf.SetInput(pointPolyData)

                cf = vtk.vtkContourFilter()
                cf.SetInputConnection(surf.GetOutputPort())
                #cf.SetInput(pointPolyData)
                cf.SetValue(0, 0.0)

                # Sometimes the contouring algorithm can create a volume whose gradient
                # vector and ordering of polygon (using the right hand rule) are
                # inconsistent. vtkReverseSense cures this problem.
                reverse = vtk.vtkReverseSense()
                reverse.SetInputConnection(cf.GetOutputPort())
                reverse.ReverseCellsOn()
                reverse.ReverseNormalsOn()

                map = vtk.vtkPolyDataMapper()
                map.SetInputConnection(reverse.GetOutputPort())
                map.ScalarVisibilityOff()

                surfaceActor = vtk.vtkActor()
                surfaceActor.SetMapper(map)
                surfaceActor.GetProperty().SetDiffuseColor(1.0000, 0.3882, 0.2784)
                surfaceActor.GetProperty().SetSpecularColor(1, 1, 1)
                surfaceActor.GetProperty().SetSpecular(.4)
                surfaceActor.GetProperty().SetSpecularPower(50)

                #delny = vtk.vtkDelaunay3D()
                #delny.SetInput(cf.GetOutput())
                #delny.SetTolerance(0.01)
                #delny.SetAlpha(0.2)
                #delny.BoundingTriangulationOff()

                ren.AddActor(surfaceActor)
                apd.AddInput(reverse.GetOutput())

                #writer.SetWritePiece(c)
                #c += 1
                #writer.SetInput(0,delny.GetOutput())
                #writer.SetInput(0,cf.GetOutput())

        writer.SetInput(0,apd.GetOutput())

        writer.Write()

        #iren.Initialize()
        #renWin.Render()
        #iren.Start()



    def load_file(self, filename):
        # Create binary reader for plot.out file
        plotFile = BinaryReader(filename)

        while True:
            try:
                cellId = plotFile.get_int(1)
                self.cells[cellId] = []
                for i in xrange(plotFile.get_int(1)):
                    self.cells[cellId].append(plotFile.get_double(3))
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

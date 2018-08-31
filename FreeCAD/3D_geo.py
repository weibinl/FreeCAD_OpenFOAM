import os
FREECADPATH='/usr/lib/freecad/lib'
import sys
sys.path.append(FREECADPATH)
import FreeCAD
import Units
import shutil
import Mesh
import MeshPart
import Part
import TemplateBuilder
import Draft
sys.path.insert(0,'/home/weibin/git/freecad_openfoam/FreeCAD_OpenFOAM/FreeCAD')
import CfdAnalysis
import CfdFluidMaterial
import CfdInitialiseFlowField
import CfdSolverFoam
import CfdPhysicsSelection
import CfdFluidBoundary
import CfdMesh
import tempfile
import CfdTools
import CfdMeshTools

# this will import the FreeCAD module

FREECAD_DOC_NAME='python_3D'
FREECAD_DOC_PATH='/home/weibin/FREECAD_DOC/'
FREECAD_DOC_EXTENSION='.fcstd' 
# extension to use when sae freecad file

WORKING_DOC=FreeCAD.newDocument(FREECAD_DOC_NAME)
FreeCAD.ActiveDocument=FreeCAD.getDocument(FREECAD_DOC_NAME)
post_diameter=10

LateralGap=10
DownstreamGap=10
TiltRatio=4 #lateral displacement vs lateral movement per bump, or # of rows required to displace one post in column
ColumnNumber=6 #one side, we built the symmetric model about the center bypass channel to reduce the required calculation time, also this is the number of channel, the post number should be ColumnNumber-1, count from the top
RowNumber=TiltRatio
MinimumBypassChannelWidth=25
LateralDisplacementPerRow=(DownstreamGap+post_diameter)/TiltRatio
DepthofChannel=10

points=[FreeCAD.Vector(post_diameter/2,0.0,0.0),FreeCAD.Vector(0.0,post_diameter/2,0.0),FreeCAD.Vector(-post_diameter/2,0.0,0.0),FreeCAD.Vector(0.0,-post_diameter/2,0.0)]
line = Draft.makeWire(points,closed=True,face=False,support=None)
Draft.autogroup(line)
f = FreeCAD.getDocument(FREECAD_DOC_NAME).addObject('Part::Extrusion', 'Extrude')
f = FreeCAD.getDocument(FREECAD_DOC_NAME).getObject('Extrude')
f.Base = FreeCAD.getDocument(FREECAD_DOC_NAME).getObject('DWire')
f.DirMode = "Normal"
f.DirLink = None
f.LengthFwd = 1.000000000000000
f.LengthRev = 0.000000000000000
f.Solid = True
f.Reversed = False
f.Symmetric = False
f.TaperAngle = 0.000000000000000
f.TaperAngleRev = 0.000000000000000
FreeCAD.ActiveDocument.recompute()
tilt_ratio=float(1.0/2.0)
Array = Draft.makeArray(FreeCAD.ActiveDocument.Extrude,FreeCAD.Vector(1,0,0),FreeCAD.Vector(0,1,0),ColumnNumber,RowNumber+1)
FreeCAD.getDocument(FREECAD_DOC_NAME).getObject("Array").IntervalX = (LateralGap+post_diameter, 0, 0)
FreeCAD.getDocument(FREECAD_DOC_NAME).getObject("Array").IntervalY = (LateralDisplacementPerRow, DownstreamGap+post_diameter, 0)
FreeCAD.ActiveDocument.recompute()

points=[FreeCAD.Vector(-3.0,-2.0,0.0),FreeCAD.Vector(13.0,-2.0,0.0),FreeCAD.Vector(13.0,10.0,0.0),FreeCAD.Vector(-3.0,10.0,0.0)]
line = Draft.makeWire(points,closed=True,face=False,support=None)
Draft.autogroup(line)
#globals()['Etra%s' % x]
i=1
Ex={}
Extrude1=FreeCAD.getDocument(FREECAD_DOC_NAME).addObject('Part::Extrusion', 'Extrude'+str(i))
Ex[1]=Extrude1
Ex[1] = FreeCAD.getDocument(FREECAD_DOC_NAME).getObject('Extrude'+str(i))
Ex[1].Base = FreeCAD.getDocument(FREECAD_DOC_NAME).getObject('DWire001')
Ex[1].DirMode = "Normal"
Ex[1].DirLink = None
Ex[1].LengthFwd = 2.000000000000000
Ex[1].LengthRev = 0.000000000000000
Ex[1].Solid = True
Ex[1].Reversed = False
Ex[1].Symmetric = False
Ex[1].TaperAngle = 0.000000000000000
Ex[1].TaperAngleRev = 0.000000000000000
FreeCAD.ActiveDocument.recompute()

FreeCAD.activeDocument().addObject("Part::Cut","Cut")
FreeCAD.activeDocument().Cut.Base = Ex[1]
FreeCAD.activeDocument().Cut.Tool = FreeCAD.activeDocument().Array
FreeCAD.ActiveDocument.recompute()


WORKING_DOC.saveAs(FREECAD_DOC_PATH+FREECAD_DOC_NAME+FREECAD_DOC_EXTENSION) 


'''
Gui.activateWorkbench("DraftWorkbench")
>>> 
>>> 
>>> 
>>> Gui.activateWorkbench("PartWorkbench")
>>> 
>>> 
>>> 
>>> 
>>> Gui.activateWorkbench("PartWorkbench")
>>> 
>>> 
>>> App.activeDocument().addObject("Part::Fuse","Fusion")
>>> App.activeDocument().Fusion.Base = App.activeDocument().Extrude001
>>> App.activeDocument().Fusion.Tool = App.activeDocument().Array
>>> Gui.activeDocument().hide("Extrude001")
>>> Gui.activeDocument().hide("Array")
>>> Gui.ActiveDocument.Fusion.ShapeColor=Gui.ActiveDocument.Extrude001.ShapeColor
>>> Gui.ActiveDocument.Fusion.DisplayMode=Gui.ActiveDocument.Extrude001.DisplayMode
>>> App.activeDocument().addObject("Part::Fuse","Fusion001")
>>> App.activeDocument().Fusion001.Base = App.activeDocument().Extrude001
>>> App.activeDocument().Fusion001.Tool = App.activeDocument().Array
>>> Gui.activeDocument().hide("Extrude001")
>>> Gui.activeDocument().hide("Array")
>>> Gui.ActiveDocument.Fusion001.ShapeColor=Gui.ActiveDocument.Extrude001.ShapeColor
>>> Gui.ActiveDocument.Fusion001.DisplayMode=Gui.ActiveDocument.Extrude001.DisplayMode
>>> Gui.SendMsgToActiveView("ViewFit")
>>> Gui.SendMsgToActiveView("ViewFit")
>>> Gui.ActiveDocument.setEdit('Fusion001',0)
>>> App.activeDocument().addObject("Part::Cut","Cut")
>>> App.activeDocument().Cut.Base = App.activeDocument().Extrude001
>>> App.activeDocument().Cut.Tool = App.activeDocument().Array
>>> Gui.activeDocument().hide("Extrude001")
>>> Gui.activeDocument().hide("Array")
>>> Gui.ActiveDocument.Cut.ShapeColor=Gui.ActiveDocument.Extrude001.ShapeColor
>>> Gui.ActiveDocument.Cut.DisplayMode=Gui.ActiveDocument.Extrude001.DisplayMode
'''
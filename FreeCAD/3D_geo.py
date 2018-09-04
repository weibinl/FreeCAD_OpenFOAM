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
PostDiameter=10
LateralGap=10
DownstreamGap=10
TiltRatio=6 #lateral displacement vs lateral movement per bump, or # of rows required to displace one post in column
ColumnNumber=4 #one side, we built the symmetric model about the center bypass channel to reduce the required calculation time, also this is the number of channel, the post number should be ColumnNumber-1, count from the top
RowNumber=TiltRatio
MinimumBypassChannelWidth=25
LateralDisplacementPerRow=(DownstreamGap+PostDiameter)/TiltRatio
DepthofChannel=10
MinimumBypassChannelWidth=25
XPosFirstPost=MinimumBypassChannelWidth/2+PostDiameter/2 #define the first post of the array
YPosFirstPost=PostDiameter/2
geom_counter=0

# generate simulation area
WidthofRegion=(ColumnNumber+1)*(LateralGap+PostDiameter)+MinimumBypassChannelWidth/2+PostDiameter/2
HeighofRegion=RowNumber*(DownstreamGap+PostDiameter)
points=[FreeCAD.Vector(0.0,0.0,0.0),FreeCAD.Vector(WidthofRegion,0.0,0.0),FreeCAD.Vector(WidthofRegion,HeighofRegion,0.0),FreeCAD.Vector(0.0,HeighofRegion,0.0)]
line = Draft.makeWire(points,closed=True,face=False,support=None)
Draft.autogroup(line)
Sim_area = FreeCAD.getDocument(FREECAD_DOC_NAME).addObject('Part::Extrusion', 'Sim_area')
Sim_area.Base = line
Sim_area.DirMode = "Normal"
Sim_area.DirLink = None
Sim_area.LengthFwd = 1.000000000000000
Sim_area.LengthRev = 0.000000000000000
Sim_area.Solid = True
Sim_area.Reversed = False
Sim_area.Symmetric = False
Sim_area.TaperAngle = 0.000000000000000
Sim_area.TaperAngleRev = 0.000000000000000
FreeCAD.ActiveDocument.recompute()

# main array
Main_Array={}

xpos=XPosFirstPost+LateralDisplacementPerRow
ypos=YPosFirstPost
points=[FreeCAD.Vector(xpos,ypos,0.0),FreeCAD.Vector(xpos+PostDiameter/2,ypos+PostDiameter/2,0.0),FreeCAD.Vector(xpos,ypos+PostDiameter,0.0),FreeCAD.Vector(xpos-PostDiameter/2,ypos+PostDiameter/2,0.0)]
line = Draft.makeWire(points,closed=True,face=False,support=None)
Draft.autogroup(line)
FreeCAD.ActiveDocument.recompute()

Array=Draft.makeArray(FreeCAD.ActiveDocument.DWire001,FreeCAD.Vector(1,0,0),FreeCAD.Vector(0,1,0),2,2)
Draft.autogroup(Array)
FreeCAD.ActiveDocument.recompute()
FreeCAD.getDocument(FREECAD_DOC_NAME).getObject("Array").IntervalX = (LateralGap+PostDiameter, 0, 0)
FreeCAD.getDocument(FREECAD_DOC_NAME).getObject("Array").IntervalY = (LateralDisplacementPerRow, DownstreamGap+PostDiameter, 0)
FreeCAD.getDocument(FREECAD_DOC_NAME).getObject("Array").NumberX = ColumnNumber
FreeCAD.getDocument(FREECAD_DOC_NAME).getObject("Array").NumberY = RowNumber-1

Main_Array[geom_counter+1] = FreeCAD.getDocument(FREECAD_DOC_NAME).addObject('Part::Extrusion', 'Extrude')
Main_Array[geom_counter+1].Base = Array
Main_Array[geom_counter+1].DirMode = "Normal"
Main_Array[geom_counter+1].DirLink = None
Main_Array[geom_counter+1].LengthFwd = 1.000000000000000
Main_Array[geom_counter+1].LengthRev = 0.000000000000000
Main_Array[geom_counter+1].Solid = True
Main_Array[geom_counter+1].Reversed = False
Main_Array[geom_counter+1].Symmetric = False
Main_Array[geom_counter+1].TaperAngle = 0.000000000000000
Main_Array[geom_counter+1].TaperAngleRev = 0.000000000000000
FreeCAD.ActiveDocument.recompute()
geom_counter=geom_counter+1
      
for CN in range(1,ColumnNumber+1,1):
    RN=RowNumber+1
    xpos=XPosFirstPost+(PostDiameter+LateralGap)*(CN-1)
    ypos=YPosFirstPost+(PostDiameter+DownstreamGap)*(RN-2)
    points=[FreeCAD.Vector(xpos,ypos,0.0),FreeCAD.Vector(xpos+PostDiameter/2,ypos+PostDiameter/2,0.0),FreeCAD.Vector(xpos,ypos+PostDiameter,0.0),FreeCAD.Vector(xpos-PostDiameter/2,ypos+PostDiameter/2,0.0)]
    line = Draft.makeWire(points,closed=True,face=False,support=None)
    Draft.autogroup(line)
    Main_Array[geom_counter+1] = FreeCAD.getDocument(FREECAD_DOC_NAME).addObject('Part::Extrusion', 'Extrude')
    Main_Array[geom_counter+1].Base = line
    Main_Array[geom_counter+1].DirMode = "Normal"
    Main_Array[geom_counter+1].DirLink = None
    Main_Array[geom_counter+1].LengthFwd = 1.000000000000000
    Main_Array[geom_counter+1].LengthRev = 0.000000000000000
    Main_Array[geom_counter+1].Solid = True
    Main_Array[geom_counter+1].Reversed = False
    Main_Array[geom_counter+1].Symmetric = False
    Main_Array[geom_counter+1].TaperAngle = 0.000000000000000
    Main_Array[geom_counter+1].TaperAngleRev = 0.000000000000000
    FreeCAD.ActiveDocument.recompute()
    geom_counter=geom_counter+1

#bypass channel
ByPassWidth={}
for i in range(1,RowNumber+1,1): 
    ByPassWidth[i]=(-0.0011*float(i)**2.0+0.2088*float(i)+18.403)/2.0

BypassWidthByOrder={}
for i in range(2,RowNumber+1,1): #diamond
   xpos=ByPassWidth[i]+PostDiameter/2.0
   ypos=YPosFirstPost+(PostDiameter+DownstreamGap)*(i-2)
   BypassWidthByOrder[i]=ByPassWidth[i]
   points=[FreeCAD.Vector(xpos,ypos,0.0),FreeCAD.Vector(xpos+PostDiameter/2,ypos+PostDiameter/2,0.0),FreeCAD.Vector(xpos,ypos+PostDiameter,0.0),FreeCAD.Vector(xpos-PostDiameter/2,ypos+PostDiameter/2,0.0)]
   line = Draft.makeWire(points,closed=True,face=False,support=None)
   Draft.autogroup(line)
   Main_Array[geom_counter+1] = FreeCAD.getDocument(FREECAD_DOC_NAME).addObject('Part::Extrusion', 'Extrude')
   Main_Array[geom_counter+1].Base = line
   Main_Array[geom_counter+1].DirMode = "Normal"
   Main_Array[geom_counter+1].DirLink = None
   Main_Array[geom_counter+1].LengthFwd = 1.000000000000000
   Main_Array[geom_counter+1].LengthRev = 0.000000000000000
   Main_Array[geom_counter+1].Solid = True
   Main_Array[geom_counter+1].Reversed = False
   Main_Array[geom_counter+1].Symmetric = False
   Main_Array[geom_counter+1].TaperAngle = 0.000000000000000
   Main_Array[geom_counter+1].TaperAngleRev = 0.000000000000000
   geom_counter=geom_counter+1
   FreeCAD.ActiveDocument.recompute()

for i in range(2,RowNumber+1,1): #rectangle
    xpos=ByPassWidth[i]+PostDiameter/2.0
    ypos=YPosFirstPost+(PostDiameter+DownstreamGap)*(i-2)
    width=XPosFirstPost+LateralDisplacementPerRow*(i-1)-ByPassWidth[i]-PostDiameter/2
    heigh=PostDiameter
    BypassWidthByOrder[i]=ByPassWidth[i]
    points=[FreeCAD.Vector(xpos,ypos,0.0),FreeCAD.Vector(xpos+width,ypos,0.0),FreeCAD.Vector(xpos+width,ypos+heigh,0.0),FreeCAD.Vector(xpos,ypos+heigh,0.0)]
    line = Draft.makeWire(points,closed=True,face=False,support=None)
    Draft.autogroup(line)
    Main_Array[geom_counter+1] = FreeCAD.getDocument(FREECAD_DOC_NAME).addObject('Part::Extrusion', 'Extrude')
    Main_Array[geom_counter+1].Base = line
    Main_Array[geom_counter+1].DirMode = "Normal"
    Main_Array[geom_counter+1].DirLink = None
    Main_Array[geom_counter+1].LengthFwd = 1.000000000000000
    Main_Array[geom_counter+1].LengthRev = 0.000000000000000
    Main_Array[geom_counter+1].Solid = True
    Main_Array[geom_counter+1].Reversed = False
    Main_Array[geom_counter+1].Symmetric = False
    Main_Array[geom_counter+1].TaperAngle = 0.000000000000000
    Main_Array[geom_counter+1].TaperAngleRev = 0.000000000000000
    geom_counter=geom_counter+1
    FreeCAD.ActiveDocument.recompute()

xpos=ByPassWidth[1]+PostDiameter/2
ypos=YPosFirstPost+(PostDiameter+DownstreamGap)*(RowNumber-1)
width=PostDiameter
heigh=PostDiameter
BypassWidthByOrder[RowNumber]=ByPassWidth[1]
BypassWidthByOrder[i]=ByPassWidth[i]
points=[FreeCAD.Vector(xpos,ypos,0.0),FreeCAD.Vector(xpos+width/2,ypos+heigh/2,0.0),FreeCAD.Vector(xpos,ypos+heigh,0.0),FreeCAD.Vector(xpos-width/2,ypos+heigh/2,0.0)]
line = Draft.makeWire(points,closed=True,face=False,support=None)
Draft.autogroup(line)
Main_Array[geom_counter+1] = FreeCAD.getDocument(FREECAD_DOC_NAME).addObject('Part::Extrusion', 'Extrude')
Main_Array[geom_counter+1].Base = line
Main_Array[geom_counter+1].DirMode = "Normal"
Main_Array[geom_counter+1].DirLink = None
Main_Array[geom_counter+1].LengthFwd = 1.000000000000000
Main_Array[geom_counter+1].LengthRev = 0.000000000000000
Main_Array[geom_counter+1].Solid = True
Main_Array[geom_counter+1].Reversed = False
Main_Array[geom_counter+1].Symmetric = False
Main_Array[geom_counter+1].TaperAngle = 0.000000000000000
Main_Array[geom_counter+1].TaperAngleRev = 0.000000000000000
geom_counter=geom_counter+1
FreeCAD.ActiveDocument.recompute()

xpos=ByPassWidth[1]+PostDiameter/2
ypos=YPosFirstPost+(PostDiameter+DownstreamGap)*(RowNumber-1)
width=XPosFirstPost-ByPassWidth[1]-PostDiameter/2
heigh=PostDiameter
BypassWidthByOrder[RowNumber]=ByPassWidth[1]
BypassWidthByOrder[i]=ByPassWidth[i]
points=[FreeCAD.Vector(xpos,ypos,0.0),FreeCAD.Vector(xpos+width,ypos,0.0),FreeCAD.Vector(xpos+width,ypos+heigh,0.0),FreeCAD.Vector(xpos,ypos+heigh,0.0)]
line = Draft.makeWire(points,closed=True,face=False,support=None)
Draft.autogroup(line)
Main_Array[geom_counter+1] = FreeCAD.getDocument(FREECAD_DOC_NAME).addObject('Part::Extrusion', 'Extrude')
Main_Array[geom_counter+1].Base = line
Main_Array[geom_counter+1].DirMode = "Normal"
Main_Array[geom_counter+1].DirLink = None
Main_Array[geom_counter+1].LengthFwd = 1.000000000000000
Main_Array[geom_counter+1].LengthRev = 0.000000000000000
Main_Array[geom_counter+1].Solid = True
Main_Array[geom_counter+1].Reversed = False
Main_Array[geom_counter+1].Symmetric = False
Main_Array[geom_counter+1].TaperAngle = 0.000000000000000
Main_Array[geom_counter+1].TaperAngleRev = 0.000000000000000
geom_counter=geom_counter+1
FreeCAD.ActiveDocument.recompute()



WallPostPosition={}
WallChannelWidth={}
for i in range(1,RowNumber+1,1):
   WallChannelWidth[i]=-0.0061*i**2+0.6*i

WallChannelWidthByOrder={}
for i in range(1,RowNumber/2,1): 
    xpos=WallChannelWidth[RowNumber/2-i]+PostDiameter+XPosFirstPost+(PostDiameter+LateralGap)*(ColumnNumber-1)+LateralDisplacementPerRow*i
    ypos=YPosFirstPost+(PostDiameter+DownstreamGap)*(i-1)
    width=PostDiameter
    heigh=PostDiameter
    BypassWidthByOrder[RowNumber]=ByPassWidth[1]
    BypassWidthByOrder[i]=ByPassWidth[i]
    points=[FreeCAD.Vector(xpos,ypos,0.0),FreeCAD.Vector(xpos+width/2,ypos+heigh/2,0.0),FreeCAD.Vector(xpos,ypos+heigh,0.0),FreeCAD.Vector(xpos-width/2,ypos+heigh/2,0.0)]
    line = Draft.makeWire(points,closed=True,face=False,support=None)
    Draft.autogroup(line)
    Main_Array[geom_counter+1] = FreeCAD.getDocument(FREECAD_DOC_NAME).addObject('Part::Extrusion', 'Extrude')
    Main_Array[geom_counter+1].Base = line
    Main_Array[geom_counter+1].DirMode = "Normal"
    Main_Array[geom_counter+1].DirLink = None
    Main_Array[geom_counter+1].LengthFwd = 1.000000000000000
    Main_Array[geom_counter+1].LengthRev = 0.000000000000000
    Main_Array[geom_counter+1].Solid = True
    Main_Array[geom_counter+1].Reversed = False
    Main_Array[geom_counter+1].Symmetric = False
    Main_Array[geom_counter+1].TaperAngle = 0.000000000000000
    Main_Array[geom_counter+1].TaperAngleRev = 0.000000000000000
    geom_counter=geom_counter+1
    FreeCAD.ActiveDocument.recompute()


for i in range(1,RowNumber/2,1): 
    xpos=WallChannelWidth[RowNumber/2-i]+PostDiameter+XPosFirstPost+(PostDiameter+LateralGap)*(ColumnNumber-1)+LateralDisplacementPerRow*i
    ypos=YPosFirstPost+(PostDiameter+DownstreamGap)*(i-1)-PostDiameter/2-DownstreamGap
    width=50
    heigh=PostDiameter*1.5+DownstreamGap*2
    BypassWidthByOrder[RowNumber]=ByPassWidth[1]
    BypassWidthByOrder[i]=ByPassWidth[i]
    points=[FreeCAD.Vector(xpos,ypos,0.0),FreeCAD.Vector(xpos+width,ypos,0.0),FreeCAD.Vector(xpos+width,ypos+heigh,0.0),FreeCAD.Vector(xpos,ypos+heigh,0.0)]
    line = Draft.makeWire(points,closed=True,face=False,support=None)
    Draft.autogroup(line)
    Main_Array[geom_counter+1] = FreeCAD.getDocument(FREECAD_DOC_NAME).addObject('Part::Extrusion', 'Extrude')
    Main_Array[geom_counter+1].Base = line
    Main_Array[geom_counter+1].DirMode = "Normal"
    Main_Array[geom_counter+1].DirLink = None
    Main_Array[geom_counter+1].LengthFwd = 1.000000000000000
    Main_Array[geom_counter+1].LengthRev = 0.000000000000000
    Main_Array[geom_counter+1].Solid = True
    Main_Array[geom_counter+1].Reversed = False
    Main_Array[geom_counter+1].Symmetric = False
    Main_Array[geom_counter+1].TaperAngle = 0.000000000000000
    Main_Array[geom_counter+1].TaperAngleRev = 0.000000000000000
    geom_counter=geom_counter+1
    FreeCAD.ActiveDocument.recompute()

for i in range (1,RowNumber/2+1,1):
    xpos=WallChannelWidth[RowNumber-i+1]+PostDiameter+XPosFirstPost+(PostDiameter+LateralGap)*(ColumnNumber-2)+LateralDisplacementPerRow*(i+RowNumber/2-1)
    ypos=YPosFirstPost+(PostDiameter+DownstreamGap)*(RowNumber/2-2+i)
    width=PostDiameter
    heigh=PostDiameter
    points=[FreeCAD.Vector(xpos,ypos,0.0),FreeCAD.Vector(xpos+width/2,ypos+heigh/2,0.0),FreeCAD.Vector(xpos,ypos+heigh,0.0),FreeCAD.Vector(xpos-width/2,ypos+heigh/2,0.0)]
    line = Draft.makeWire(points,closed=True,face=False,support=None)
    Draft.autogroup(line)
    WallPostPosition[i+RowNumber/2-1]=WallChannelWidth[RowNumber-i+1]+PostDiameter+XPosFirstPost+(PostDiameter+LateralGap)*(ColumnNumber-2)+LateralDisplacementPerRow*(i+RowNumber/2-1)
    WallChannelWidthByOrder[i+RowNumber/2-1]=WallChannelWidth[RowNumber-i+1]
    Main_Array[geom_counter+1] = FreeCAD.getDocument(FREECAD_DOC_NAME).addObject('Part::Extrusion', 'Extrude')
    Main_Array[geom_counter+1].Base = line
    Main_Array[geom_counter+1].DirMode = "Normal"
    Main_Array[geom_counter+1].DirLink = None
    Main_Array[geom_counter+1].LengthFwd = 1.000000000000000
    Main_Array[geom_counter+1].LengthRev = 0.000000000000000
    Main_Array[geom_counter+1].Solid = True
    Main_Array[geom_counter+1].Reversed = False
    Main_Array[geom_counter+1].Symmetric = False
    Main_Array[geom_counter+1].TaperAngle = 0.000000000000000
    Main_Array[geom_counter+1].TaperAngleRev = 0.000000000000000
    geom_counter=geom_counter+1
    FreeCAD.ActiveDocument.recompute()

for i in range(1,RowNumber/2+1,1):
    xpos=WallChannelWidth[RowNumber-i+1]+PostDiameter+XPosFirstPost+(PostDiameter+LateralGap)*(ColumnNumber-2)+LateralDisplacementPerRow*(i+RowNumber/2-1)
    ypos=YPosFirstPost+(PostDiameter+DownstreamGap)*(RowNumber/2-2+i)-PostDiameter/2
    width=50
    heigh=PostDiameter*1.5+DownstreamGap
    points=[FreeCAD.Vector(xpos,ypos,0.0),FreeCAD.Vector(xpos+width,ypos,0.0),FreeCAD.Vector(xpos+width,ypos+heigh,0.0),FreeCAD.Vector(xpos,ypos+heigh,0.0)]
    line = Draft.makeWire(points,closed=True,face=False,support=None)
    Draft.autogroup(line)
    WallPostPosition[i+RowNumber/2-1]=WallChannelWidth[RowNumber-i+1]+PostDiameter+XPosFirstPost+(PostDiameter+LateralGap)*(ColumnNumber-2)+LateralDisplacementPerRow*(i+RowNumber/2-1)
    WallChannelWidthByOrder[i+RowNumber/2-1]=WallChannelWidth[RowNumber-i+1]
    Main_Array[geom_counter+1] = FreeCAD.getDocument(FREECAD_DOC_NAME).addObject('Part::Extrusion', 'Extrude')
    Main_Array[geom_counter+1].Base = line
    Main_Array[geom_counter+1].DirMode = "Normal"
    Main_Array[geom_counter+1].DirLink = None
    Main_Array[geom_counter+1].LengthFwd = 1.000000000000000
    Main_Array[geom_counter+1].LengthRev = 0.000000000000000
    Main_Array[geom_counter+1].Solid = True
    Main_Array[geom_counter+1].Reversed = False
    Main_Array[geom_counter+1].Symmetric = False
    Main_Array[geom_counter+1].TaperAngle = 0.000000000000000
    Main_Array[geom_counter+1].TaperAngleRev = 0.000000000000000
    geom_counter=geom_counter+1
    FreeCAD.ActiveDocument.recompute()

xpos=WallChannelWidth[RowNumber/2]+PostDiameter+XPosFirstPost+(PostDiameter+LateralGap)*(ColumnNumber-1)
ypos=YPosFirstPost+(PostDiameter+DownstreamGap)*(RowNumber-1)
width=PostDiameter
heigh=PostDiameter
WallPostPosition[TiltRatio]=WallChannelWidth[RowNumber/2]+PostDiameter+XPosFirstPost+(PostDiameter+LateralGap)*(ColumnNumber-1)
WallChannelWidthByOrder[TiltRatio]=WallChannelWidth[RowNumber/2]
points=[FreeCAD.Vector(xpos,ypos,0.0),FreeCAD.Vector(xpos+width/2,ypos+heigh/2,0.0),FreeCAD.Vector(xpos,ypos+heigh,0.0),FreeCAD.Vector(xpos-width/2,ypos+heigh/2,0.0)]
line = Draft.makeWire(points,closed=True,face=False,support=None)
Draft.autogroup(line)
Main_Array[geom_counter+1] = FreeCAD.getDocument(FREECAD_DOC_NAME).addObject('Part::Extrusion', 'Extrude')
Main_Array[geom_counter+1].Base = line
Main_Array[geom_counter+1].DirMode = "Normal"
Main_Array[geom_counter+1].DirLink = None
Main_Array[geom_counter+1].LengthFwd = 1.000000000000000
Main_Array[geom_counter+1].LengthRev = 0.000000000000000
Main_Array[geom_counter+1].Solid = True
Main_Array[geom_counter+1].Reversed = False
Main_Array[geom_counter+1].Symmetric = False
Main_Array[geom_counter+1].TaperAngle = 0.000000000000000
Main_Array[geom_counter+1].TaperAngleRev = 0.000000000000000
geom_counter=geom_counter+1
FreeCAD.ActiveDocument.recompute()

xpos=WallChannelWidth[RowNumber/2]+PostDiameter+XPosFirstPost+(PostDiameter+LateralGap)*(ColumnNumber-1)
ypos=YPosFirstPost+(PostDiameter+DownstreamGap)*(RowNumber-1)-PostDiameter/2
width=50
heigh=PostDiameter*1.5+DownstreamGap
points=[FreeCAD.Vector(xpos,ypos,0.0),FreeCAD.Vector(xpos+width,ypos,0.0),FreeCAD.Vector(xpos+width,ypos+heigh,0.0),FreeCAD.Vector(xpos,ypos+heigh,0.0)]
line = Draft.makeWire(points,closed=True,face=False,support=None)
Draft.autogroup(line)
Main_Array[geom_counter+1] = FreeCAD.getDocument(FREECAD_DOC_NAME).addObject('Part::Extrusion', 'Extrude')
Main_Array[geom_counter+1].Base = line
Main_Array[geom_counter+1].DirMode = "Normal"
Main_Array[geom_counter+1].DirLink = None
Main_Array[geom_counter+1].LengthFwd = 1.000000000000000
Main_Array[geom_counter+1].LengthRev = 0.000000000000000
Main_Array[geom_counter+1].Solid = True
Main_Array[geom_counter+1].Reversed = False
Main_Array[geom_counter+1].Symmetric = False
Main_Array[geom_counter+1].TaperAngle = 0.000000000000000
Main_Array[geom_counter+1].TaperAngleRev = 0.000000000000000
geom_counter=geom_counter+1
FreeCAD.ActiveDocument.recompute()
'''
Cut={}
Cut[1]=Sim_area
print geom_counter
for i in range(1,geom_counter+1,1):
    Cut[i+1]=FreeCAD.activeDocument().addObject("Part::Cut","Cut")
    Cut[i+1].Base = Cut[i]
    Cut[i+1].Tool = Main_Array[i]
    FreeCAD.ActiveDocument.recompute()
    print i
'''
WORKING_DOC.saveAs(FREECAD_DOC_PATH+FREECAD_DOC_NAME+FREECAD_DOC_EXTENSION) 


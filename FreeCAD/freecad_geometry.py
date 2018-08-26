import os
FREECADPATH='/usr/lib/freecad/lib'
import sys
sys.path.append(FREECADPATH)
import FreeCAD

sys.path.insert(0,'/home/weibin/git/freecad_openfoam/FreeCAD_OpenFOAM/FreeCAD')
import CfdAnalysis
import CfdFluidMaterial
import CfdInitialiseFlowField
import CfdSolverFoam
import CfdPhysicsSelection
import CfdFluidBoundary
import CfdMesh
import CfdMeshTools
import Part
import tempfile
import CfdTools
# this will import the FreeCAD module

FREECAD_DOC_NAME='python_scrit_test'
FREECAD_DOC_PATH='/home/weibin/FREECAD_DOC/'
FREECAD_DOC_EXTENSION='.fcstd' 
# extension to use when sae freecad file

FreeCAD.newDocument(FREECAD_DOC_NAME)
FreeCAD.ActiveDocument=FreeCAD.getDocument(FREECAD_DOC_NAME)
CYLINDER_1=FreeCAD.ActiveDocument.addObject('Part::Cylinder','CYLINDER_1')
CYLINDER_1.Radius=5

#Gui.activateWorkbench("CfdOFWorkbench")

analysis = CfdAnalysis.makeCfdAnalysis('CfdAnalysis')

analysis.addObject(CfdFluidMaterial.makeCfdFluidMaterial('FluidProperties'))
analysis.addObject(CfdInitialiseFlowField.makeCfdInitialFlowField())
analysis.addObject(CfdSolverFoam.makeCfdSolverFoam())
analysis.addObject(CfdPhysicsSelection.makeCfdPhysicsSelection())

obj = FreeCAD.ActiveDocument.PhysicsModel
obj.Time = 'Steady'
obj.Phase = 'Single'
obj.Flow = 'Incompressible'
obj.Thermal = 'None'
obj.Turbulence = 'Laminar'
obj.gx = '0 mm/s^2'
obj.gy = '-9.8e+03 mm/s^2'
obj.gz = '0 mm/s^2'

mat = FreeCAD.ActiveDocument.FluidProperties
mat.Density = '1e+03 kg/m^3'
mat.DynamicViscosity = '0.001 kg/(m*s)'

# Values are converted to SI units and stored (eg. m/s)
init = FreeCAD.ActiveDocument.InitialiseFields.InitialVariables
init['PotentialFoam'] = True
init['UseInletUPValues'] = False
init['Ux'] = 0.0
init['Uy'] = 0.0
init['Uz'] = 0.0
init['Pressure'] = 0.0
init['alphas'] = {}
init['UseInletTemperatureValues'] = False
init['Temperature'] = 290.0
init['UseInletTurbulenceValues'] = False
init['omega'] = 0.994837673637
init['k'] = 0.01
init['Inlet'] = ''
FreeCAD.ActiveDocument.InitialiseFields.InitialVariables = init

# Values are converted to SI units and stored (eg. m/s)
analysis.addObject(CfdFluidBoundary.makeCfdFluidBoundary())
bc = FreeCAD.ActiveDocument.CfdFluidBoundary.BoundarySettings
bc['BoundaryType'] = 'inlet'
bc['BoundarySubtype'] = 'uniformVelocity'
bc['ThermalBoundaryType'] = 'fixedValue'
bc['VelocityIsCartesian'] = True
bc['Ux'] = 0.0
bc['Uy'] = 0.0
bc['Uz'] = -0.01
bc['VelocityMag'] = 0.0
bc['DirectionFace'] = ''
bc['ReverseNormal'] = True
bc['MassFlowRate'] = 0.0
bc['VolFlowRate'] = 0.0
bc['Pressure'] = 0.0
bc['SlipRatio'] = 0.0
bc['Temperature'] = 290.0
bc['HeatFlux'] = 0.0
bc['HeatTransferCoeff'] = 0.0
bc['TurbulenceInletSpecification'] = 'intensityAndLengthScale'
bc['TurbulentKineticEnergy'] = 0.01
bc['SpecificDissipationRate'] = 0.994837673637
bc['TurbulenceIntensity'] = 0.1
bc['TurbulenceLengthScale'] = 0.1
bc['PressureDropCoeff'] = 0.0
bc['ScreenWireDiameter'] = 0.0001
bc['ScreenSpacing'] = 0.0
bc['PorousBaffleMethod'] = 0
FreeCAD.ActiveDocument.CfdFluidBoundary.BoundarySettings = bc
FreeCAD.ActiveDocument.CfdFluidBoundary.Label = 'inlet'
FreeCAD.ActiveDocument.CfdFluidBoundary.References = []
FreeCAD.ActiveDocument.CfdFluidBoundary.References.append(('Cylinder', 'Face2'))
FreeCAD.ActiveDocument.recompute()

# Values are converted to SI units and stored (eg. m/s)
init = FreeCAD.ActiveDocument.InitialiseFields.InitialVariables
init['PotentialFoam'] = True
init['UseInletUPValues'] = False
init['Ux'] = 0.0
init['Uy'] = 0.0
init['Uz'] = 0.0
init['Pressure'] = 0.0
init['alphas'] = {}
init['UseInletTemperatureValues'] = False
init['Temperature'] = 290.0
init['UseInletTurbulenceValues'] = False
init['omega'] = 0.994837673637
init['k'] = 0.01
init['Inlet'] = ''
FreeCAD.ActiveDocument.InitialiseFields.InitialVariables = init
# Values are converted to SI units and stored (eg. m/s)
bc = FreeCAD.ActiveDocument.CfdFluidBoundary.BoundarySettings
bc['BoundaryType'] = 'inlet'
bc['BoundarySubtype'] = 'uniformVelocity'
bc['ThermalBoundaryType'] = 'fixedValue'
bc['VelocityIsCartesian'] = True
bc['Ux'] = 0.0
bc['Uy'] = 0.0
bc['Uz'] = -0.01
bc['VelocityMag'] = 0.0
bc['DirectionFace'] = ''
bc['ReverseNormal'] = True
bc['MassFlowRate'] = 0.0
bc['VolFlowRate'] = 0.0
bc['Pressure'] = 0.0
bc['SlipRatio'] = 0.0
bc['Temperature'] = 290.0
bc['HeatFlux'] = 0.0
bc['HeatTransferCoeff'] = 0.0
bc['TurbulenceInletSpecification'] = 'intensityAndLengthScale'
bc['TurbulentKineticEnergy'] = 0.01
bc['SpecificDissipationRate'] = 0.994837673637
bc['TurbulenceIntensity'] = 0.1
bc['TurbulenceLengthScale'] = 0.1
bc['PressureDropCoeff'] = 0.0
bc['ScreenWireDiameter'] = 0.0001
bc['ScreenSpacing'] = 0.0
bc['PorousBaffleMethod'] = 0
FreeCAD.ActiveDocument.CfdFluidBoundary.BoundarySettings = bc
FreeCAD.ActiveDocument.CfdFluidBoundary.Label = 'inlet'
FreeCAD.ActiveDocument.CfdFluidBoundary.References = []
FreeCAD.ActiveDocument.CfdFluidBoundary.References.append(('Cylinder', 'Face2'))
FreeCAD.ActiveDocument.recompute()

# Values are converted to SI units and stored (eg. m/s)
analysis.addObject(CfdFluidBoundary.makeCfdFluidBoundary())  ##remember to add this line whenever you add a new boundary condition
bc = FreeCAD.ActiveDocument.CfdFluidBoundary001.BoundarySettings
bc['BoundaryType'] = 'outlet'
bc['BoundarySubtype'] = 'staticPressure'
bc['ThermalBoundaryType'] = 'fixedValue'
bc['VelocityIsCartesian'] = True
bc['Ux'] = 0.0
bc['Uy'] = 0.0
bc['Uz'] = 0.0
bc['VelocityMag'] = 0.0
bc['DirectionFace'] = ''
bc['ReverseNormal'] = False
bc['MassFlowRate'] = 0.0
bc['VolFlowRate'] = 0.0
bc['Pressure'] = 0.0
bc['SlipRatio'] = 0.0
bc['Temperature'] = 290.0
bc['HeatFlux'] = 0.0
bc['HeatTransferCoeff'] = 0.0
bc['TurbulenceInletSpecification'] = 'intensityAndLengthScale'
bc['TurbulentKineticEnergy'] = 0.01
bc['SpecificDissipationRate'] = 0.994837673637
bc['TurbulenceIntensity'] = 0.1
bc['TurbulenceLengthScale'] = 0.1
bc['PressureDropCoeff'] = 0.0
bc['ScreenWireDiameter'] = 0.0001
bc['ScreenSpacing'] = 0.0
bc['PorousBaffleMethod'] = 0
FreeCAD.ActiveDocument.CfdFluidBoundary001.BoundarySettings = bc
FreeCAD.ActiveDocument.CfdFluidBoundary001.Label = 'outlet'
FreeCAD.ActiveDocument.CfdFluidBoundary001.References = []
FreeCAD.ActiveDocument.CfdFluidBoundary001.References.append(('Cylinder', 'Face3'))
FreeCAD.ActiveDocument.recompute()

# Values are converted to SI units and stored (eg. m/s)
analysis.addObject(CfdFluidBoundary.makeCfdFluidBoundary())  ##remember to add this line whenever you add a new boundary condition
bc = FreeCAD.ActiveDocument.CfdFluidBoundary002.BoundarySettings
bc['BoundaryType'] = 'wall'
bc['BoundarySubtype'] = 'fixed'
bc['ThermalBoundaryType'] = 'fixedValue'
bc['VelocityIsCartesian'] = True
bc['Ux'] = 0.0
bc['Uy'] = 0.0
bc['Uz'] = 0.0
bc['VelocityMag'] = 0.0
bc['DirectionFace'] = ''
bc['ReverseNormal'] = False
bc['MassFlowRate'] = 0.0
bc['VolFlowRate'] = 0.0
bc['Pressure'] = 0.0
bc['SlipRatio'] = 0.0
bc['Temperature'] = 290.0
bc['HeatFlux'] = 0.0
bc['HeatTransferCoeff'] = 0.0
bc['TurbulenceInletSpecification'] = 'intensityAndLengthScale'
bc['TurbulentKineticEnergy'] = 0.01
bc['SpecificDissipationRate'] = 0.994837673637
bc['TurbulenceIntensity'] = 0.1
bc['TurbulenceLengthScale'] = 0.1
bc['PressureDropCoeff'] = 0.0
bc['ScreenWireDiameter'] = 0.0001
bc['ScreenSpacing'] = 0.0
bc['PorousBaffleMethod'] = 0
FreeCAD.ActiveDocument.CfdFluidBoundary002.BoundarySettings = bc
FreeCAD.ActiveDocument.CfdFluidBoundary002.Label = 'wall'
FreeCAD.ActiveDocument.CfdFluidBoundary002.References = []
FreeCAD.ActiveDocument.CfdFluidBoundary002.References.append(('Cylinder', 'Face1'))
FreeCAD.ActiveDocument.recompute()

##start the mesh
Cylinder_Mesh=CfdMesh.makeCfdMesh('Cylinder_Mesh')
FreeCAD.ActiveDocument.ActiveObject.Part = FreeCAD.ActiveDocument.CYLINDER_1
analysis.addObject(FreeCAD.ActiveDocument.ActiveObject)

FreeCAD.ActiveDocument.Cylinder_Mesh.CharacteristicLengthMax = '0 mm'
FreeCAD.ActiveDocument.Cylinder_Mesh.MeshUtility = "snappyHexMesh"
FreeCAD.ActiveDocument.Cylinder_Mesh.ElementDimension = '3D'
FreeCAD.ActiveDocument.Cylinder_Mesh.CellsBetweenLevels = 3
FreeCAD.ActiveDocument.Cylinder_Mesh.EdgeRefinement = 0
FreeCAD.ActiveDocument.Cylinder_Mesh.PointInMesh = {'y': 0.0, 'x': 0.0, 'z': 0.0}

cart_mesh=CfdMeshTools.CfdMeshTools(Cylinder_Mesh)
cart_mesh.get_tmp_file_paths()
cart_mesh.setup_mesh_case_dir()
cart_mesh.get_region_data()
cart_mesh.write_mesh_case()
cart_mesh.write_part_file()

tmpdir = tempfile.gettempdir()
meshCaseDir = os.path.join(tmpdir, 'meshCase')
cmd = CfdTools.makeRunCommand('./Allmesh', meshCaseDir, source_env=False)
os.system(cmd[2])

## setup the simulation and write the file



FreeCAD.ActiveDocument.saveAs(FREECAD_DOC_PATH+FREECAD_DOC_NAME+FREECAD_DOC_EXTENSION)


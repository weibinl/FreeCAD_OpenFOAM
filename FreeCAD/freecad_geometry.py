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

FREECAD_DOC_NAME='python_scrit_test2'
FREECAD_DOC_PATH='/home/weibin/FREECAD_DOC/'
FREECAD_DOC_EXTENSION='.fcstd' 
# extension to use when sae freecad file

FreeCAD.newDocument(FREECAD_DOC_NAME)
FreeCAD.ActiveDocument=FreeCAD.getDocument(FREECAD_DOC_NAME)
CYLINDER_1=FreeCAD.ActiveDocument.addObject('Part::Cylinder','CYLINDER_1')
CYLINDER_1.Radius=2
CYLINDER_1.Height=20

#Gui.activateWorkbench("CfdOFWorkbench")

analysis = CfdAnalysis.makeCfdAnalysis('CfdAnalysis')

analysis.addObject(CfdFluidMaterial.makeCfdFluidMaterial('FluidProperties'))
analysis.addObject(CfdInitialiseFlowField.makeCfdInitialFlowField())
solver=analysis.addObject(CfdSolverFoam.makeCfdSolverFoam())
physics_model=obj = FreeCAD.ActiveDocument.addObject("App::FeaturePython", "PhysicsModel")
physics_model.addProperty("App::PropertyEnumeration", "Time", "Physics modelling","Resolve time dependence")
physics_model.Time = ['Steady', 'Transient']
physics_model.Time = 'Steady'

physics_model.addProperty("App::PropertyEnumeration", "Flow", "Physics modelling","Flow algorithm")
physics_model.Flow = ['Incompressible', 'Compressible', 'HighMachCompressible']
physics_model.Flow = 'Incompressible'

physics_model.addProperty("App::PropertyEnumeration", "Thermal", "Physics modelling","Thermal modelling")
physics_model.Thermal = ['None', 'Buoyancy', 'Energy']
physics_model.Thermal = 'None'

physics_model.addProperty("App::PropertyEnumeration", "Phase", "Physics modelling","Type of phases present")
physics_model.Phase = ['Single', 'FreeSurface']
physics_model.Phase = 'Single'

physics_model.addProperty("App::PropertyEnumeration", "Turbulence", "Physics modelling","Type of turbulence modelling")
physics_model.Turbulence = ['Inviscid', 'Laminar', 'RANS']
physics_model.Turbulence = 'Laminar'

physics_model.addProperty("App::PropertyEnumeration", "TurbulenceModel", "Physics modelling","Turbulence model")
physics_model.TurbulenceModel = ['kOmegaSST']

physics_model.addProperty("App::PropertyAcceleration", "gx", "Physics modelling","Gravitational acceleration vector (x component)")
physics_model.gx = '0 m/s^2'

physics_model.addProperty("App::PropertyAcceleration", "gy", "Physics modelling","Gravitational acceleration vector (y component)")
physics_model.gy = '-9.81 m/s^2'

physics_model.addProperty("App::PropertyAcceleration", "gz", "Physics modelling","Gravitational acceleration vector (z component)")
physics_model.gz = '0 m/s^2'

#set up material
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
bc['Uz'] = -1
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
FreeCAD.ActiveDocument.CfdFluidBoundary.References.append(('CYLINDER_1', 'Face2'))
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
FreeCAD.ActiveDocument.CfdFluidBoundary001.References.append(('CYLINDER_1', 'Face3'))
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
FreeCAD.ActiveDocument.CfdFluidBoundary002.References.append(('CYLINDER_1', 'Face1'))
FreeCAD.ActiveDocument.recompute()



##start the mesh
Cylinder_Mesh=CfdMesh.makeCfdMesh('Cylinder_Mesh')
FreeCAD.ActiveDocument.ActiveObject.Part = FreeCAD.ActiveDocument.CYLINDER_1
analysis.addObject(FreeCAD.ActiveDocument.ActiveObject)
FreeCAD.ActiveDocument.Cylinder_Mesh.MeshUtility = "cfMesh"
FreeCAD.ActiveDocument.Cylinder_Mesh.CharacteristicLengthMax = '0 mm'
FreeCAD.ActiveDocument.Cylinder_Mesh.ElementDimension = '3D'
FreeCAD.ActiveDocument.Cylinder_Mesh.CellsBetweenLevels = 3
FreeCAD.ActiveDocument.Cylinder_Mesh.EdgeRefinement = 0
FreeCAD.ActiveDocument.Cylinder_Mesh.PointInMesh = {'y': 0.0, 'x': 0.0, 'z': 0.0}


cart_mesh=CfdMeshTools.CfdMeshTools(Cylinder_Mesh)
cart_mesh.get_tmp_file_paths()
cart_mesh.setup_mesh_case_dir()
cart_mesh.get_region_data()
cart_mesh.write_mesh_case() ##remember to include /data file in the positon where you need to use defaultmesh
cart_mesh.write_part_file()

tmpdir = tempfile.gettempdir()
meshCaseDir = os.path.join(tmpdir,'meshCase')

cmd = CfdTools.makeRunCommand('./Allmesh', meshCaseDir, source_env=False)
os.system(cmd[2])

## setup solution 
wd=CfdTools.getTempWorkingDir()

import CfdCaseWriterFoam
solver_obj = CfdTools.getSolver(analysis)
case_folder = os.path.join(solver_obj.WorkingDir, solver_obj.InputCaseName)
cwd = os.curdir
# Perform initialisation here rather than __init__ in case of path changes
case_folder = os.path.expanduser(os.path.abspath(case_folder))
mesh_file_name = os.path.join(case_folder, solver_obj.InputCaseName, u".unv")
template_path = os.path.join(CfdTools.get_module_path(), "data", "defaults")
solverSettingsDict = CfdTools.getSolverSettings(solver_obj)
initial_conditions = CfdTools.getInitialConditions(analysis)
bc_group = CfdTools.getCfdBoundaryGroup(analysis)

def bafflesPresent():
    for b in bc_group:
        if b.BoundarySettings['BoundaryType'] == 'baffle':
            return True
    return False

initialisationZone_objs = CfdTools.getInitialisationZoneObjects(analysis)
zone_objs = CfdTools.getZoneObjects(analysis)
mesh_obj = CfdTools.getMesh(analysis)
material_objs = CfdTools.getMaterials(analysis)
porousZone_objs = CfdTools.getPorousZoneObjects(analysis)

def porousBafflesPresent():
    for b in bc_group:
        if b.BoundarySettings['BoundaryType'] == 'baffle' and \
            b.BoundarySettings['BoundarySubtype'] == 'porousBaffle':
            return True
    return False

phys_settings = dict(zip(physics_model.PropertiesList,(getattr(physics_model, prop) for prop in physics_model.PropertiesList)))
settings = {
            'physics': phys_settings,
            'fluidProperties': [],  # Order is important, so use a list
            'initialValues': init,
            'boundaries': dict((b.Label, b.BoundarySettings) for b in bc_group),
            'bafflesPresent': bafflesPresent(),
            'porousZones': {},
            'porousZonesPresent': False,
            'initialisationZones': {o.Label: o.initialisationZoneProperties for o in initialisationZone_objs},
            'initialisationZonesPresent': len(initialisationZone_objs) > 0,
            'zones': {o.Label: {'PartNameList': tuple(o.partNameList)} for o in zone_objs},
            'zonesPresent': len(zone_objs) > 0,
            'meshType': mesh_obj.Proxy.Type,
            'solver': solverSettingsDict,
            'system': {},
            'runChangeDictionary':False,
            'createPatchesFromSnappyBaffles': False
            }

def processSystemSettings():
    system_settings = settings['system']
    system_settings['FoamRuntime'] = CfdTools.getFoamRuntime()
    system_settings['CasePath'] = case_folder
    system_settings['TranslatedCasePath'] = CfdTools.translatePath(case_folder)
    system_settings['FoamPath'] = CfdTools.getFoamDir()
    system_settings['TranslatedFoamPath'] = CfdTools.translatePath(CfdTools.getFoamDir())
processSystemSettings()

def getSolverName():
        """ Solver name is selected based on selected physics. This should only be extended as additional physics are
        included. """
        solver = None
        if physics_model.Phase == 'Single':
            if len(material_objs) == 1:
                if physics_model.Flow == 'Incompressible':
                    if physics_model.Thermal == 'None':
                        if physics_model.Time == 'Transient':
                            solver = 'pimpleFoam'
                        else:
                            if porousZone_objs or porousBafflesPresent():
                                solver = 'porousSimpleFoam'
                            else:
                                solver = 'simpleFoam'
                    else:
                        raise RuntimeError("Only isothermal simulation currently supported for incompressible flow.")
                elif physics_model.Flow == 'HighMachCompressible':
                    solver = 'hisa'
                else:
                    raise RuntimeError(physics_model.Flow + " flow model currently not supported.")
            else:
                raise RuntimeError("Only one material object may be present for single phase simulation.")
        elif physics_model.Phase == 'FreeSurface':
            if physics_model.Time == 'Transient':
                if physics_model.Thermal == 'None':
                    if len(material_objs) == 2:
                        solver = 'interFoam'
                    elif len(material_objs) > 2:
                        solver = 'multiphaseInterFoam'
                    else:
                        raise RuntimeError("At least two material objects must be present for free surface simulation.")
                else:
                    raise RuntimeError("Only isothermal analysis is currently supported for free surface flow simulation.")
            else:
                raise RuntimeError("Only transient analysis is supported for free surface flow simulation.")
        else:
            raise RuntimeError(physics_model.Phase + " phase model currently not supported.")

        # Catch-all in case
        if solver is None:
            raise RuntimeError("No solver is supported to handle the selected physics with {} phases.".format(
                len(material_objs)))
        return solver

def processSolverSettings():
    solver_settings = settings['solver']
    if solver_settings['parallel']:
        if solver_settings['parallelCores'] < 2:
            solver_settings['parallelCores'] = 2
    solver_settings['solverName'] = getSolverName()

processSolverSettings()

def processFluidProperties():
        # self.material_obj currently stores everything as a string
        # Convert to (mostly) SI numbers for OpenFOAM
        
        for material_obj in material_objs:
            mp = {}
            mp['Name'] = material_obj.Label
            if 'Density' in material_obj.PropertiesList:
                mp['Density'] = Units.Quantity(material_obj.Density).getValueAs("kg/m^3").Value
            if 'DynamicViscosity' in material_obj.PropertiesList:
                if physics_model.Turbulence == 'Inviscid':
                    mp['DynamicViscosity'] = 0.0
                else:
                    mp['DynamicViscosity'] = Units.Quantity(material_obj.DynamicViscosity).getValueAs("kg/m/s").Value
                mp['KinematicViscosity'] = mp['DynamicViscosity']/mp['Density']
            if 'MolarMass' in material_obj.PropertiesList:
                # OpenFOAM uses kg/kmol
                mp['MolarMass'] = Units.Quantity(material_obj.MolarMass).getValueAs("kg/mol").Value*1000
            if 'Cp' in material_obj.PropertiesList:
                mp['Cp'] = Units.Quantity(material_obj.Cp).getValueAs("J/kg/K").Value
            if 'SutherlandConstant' in material_obj.PropertiesList:
                #mp['Cp'] = Units.Quantity(material_obj.SutherlandConstant).getValueAs("kg/m/s/K^0.5").Value
                # TODO workaround: Have to use wrong units as fractional units currently not supported
                mp['SutherlandConstant'] = Units.Quantity(material_obj.SutherlandConstant).getValueAs("kg/m/s").Value
            if 'SutherlandTemperature' in material_obj.PropertiesList:
                mp['SutherlandTemperature'] = Units.Quantity(material_obj.SutherlandTemperature).getValueAs("K").Value
            settings['fluidProperties'].append(mp)

processFluidProperties()

def processBoundaryConditions():
        """ Compute any quantities required before case build """
        for bc_name in settings['boundaries']:
            bc = settings['boundaries'][bc_name]
            if not bc['VelocityIsCartesian']:
                veloMag = bc['VelocityMag']
                face = bc['DirectionFace'].split(':')
                # See if entered face actually exists and is planar
                try:
                    selected_object = analysis.Document.getObject(face[0])
                    if hasattr(selected_object, "Shape"):
                        elt = selected_object.Shape.getElement(face[1])
                        if elt.ShapeType == 'Face' and CfdTools.is_planar(elt):
                            n = elt.normalAt(0.5, 0.5)
                            if bc['ReverseNormal']:
                               n = [-ni for ni in n]
                            velocity = [ni*veloMag for ni in n]
                            bc['Ux'] = velocity[0]
                            bc['Uy'] = velocity[1]
                            bc['Uz'] = velocity[2]
                        else:
                            raise RuntimeError
                    else:
                        raise RuntimeError
                except (SystemError, RuntimeError):
                    raise RuntimeError(bc['DirectionFace'] + " is not a valid, planar face.")
            if settings['solver']['solverName'] in ['simpleFoam', 'porousSimpleFoam', 'pimpleFoam']:
                bc['KinematicPressure'] = bc['Pressure']/settings['fluidProperties'][0]['Density']

            if bc['PorousBaffleMethod'] == 1:
                wireDiam = bc['ScreenWireDiameter']
                spacing = bc['ScreenSpacing']
                CD = 1.0  # Drag coeff of wire (Simmons - valid for Re > ~300)
                beta = (1-wireDiam/spacing)**2
                bc['PressureDropCoeff'] = CD*(1-beta)

            if settings['solver']['solverName'] in ['interFoam', 'multiphaseInterFoam']:
                # Make sure the first n-1 alpha values exist, and write the n-th one
                # consistently for multiphaseInterFoam
                sum_alpha = 0.0
                alphas_new = {}
                for i, m in enumerate(settings['fluidProperties']):
                    alpha_name = m['Name']
                    if i == len(settings['fluidProperties']) - 1:
                        if settings['solver']['solverName'] == 'multiphaseInterFoam':
                            alphas_new[alpha_name] = 1.0 - sum_alpha
                    else:
                        alpha = bc.get('alphas', {}).get(alpha_name, 0.0)
                        alphas_new[alpha_name] = alpha
                        sum_alpha += alpha
                bc['alphas'] = alphas_new

processBoundaryConditions()

def processInitialConditions():
        """ Do any required computations before case build. Boundary conditions must be processed first. """
        
        initial_values = settings['initialValues']
        if settings['solver']['solverName'] in ['simpleFoam', 'porousSimpleFoam', 'pimpleFoam']:
            mat_prop = settings['fluidProperties'][0]
            initial_values['KinematicPressure'] = initial_values['Pressure'] / mat_prop['Density']
        if settings['solver']['solverName'] in ['interFoam', 'multiphaseInterFoam']:
            # Make sure the first n-1 alpha values exist, and write the n-th one
            # consistently for multiphaseInterFoam
            sum_alpha = 0.0
            alphas_new = {}
            for i, m in enumerate(settings['fluidProperties']):
                alpha_name = m['Name']
                if i == len(settings['fluidProperties'])-1:
                    if settings['solver']['solverName'] == 'multiphaseInterFoam':
                        alphas_new[alpha_name] = 1.0-sum_alpha
                else:
                    alpha = initial_values.get('alphas', {}).get(alpha_name, 0.0)
                    alphas_new[alpha_name] = alpha
                    sum_alpha += alpha
            initial_values['alphas'] = alphas_new

        physics = settings['physics']

        #Find inlet if copying initial values from
        inlet_bc = None
        if initial_values['UseInletUPValues'] or \
           (physics['Thermal'] != 'None' and initial_values['UseInletThermalValues']) or \
           (physics['TurbulenceModel'] is not None and initial_values['UseInletTurbulenceValues']):
                first_inlet = None
                ninlets = 0
                for bc_name in settings['boundaries']:
                    bc = settings['boundaries'][bc_name]
                    if bc['BoundaryType'] in ['inlet', 'farField']:
                        ninlets = ninlets + 1
                        # Save first inlet in case match not found
                        if ninlets == 1:
                            first_inlet = bc
                        if initial_values['Inlet'] == bc_name:
                            inlet_bc = bc
                            break
                if inlet_bc is None:
                    if initial_values['Inlet']:
                        if ninlets == 1:
                            inlet_bc = first_inlet
                        else:
                            raise RuntimeError("Inlet {} not found to copy initial conditions from."
                                            .format(initial_values['Inlet']))
                    else:
                        inlet_bc = first_inlet
                if inlet_bc is None:
                    raise RuntimeError("No inlets found to copy initial conditions from.")

        # Copy velocity
        if initial_values['UseInletUPValues']:
            if inlet_bc['BoundaryType'] == 'farField':
                initial_values['Ux'] = inlet_bc['Ux']
                initial_values['Uy'] = inlet_bc['Uy']
                initial_values['Uz'] = inlet_bc['Uz']
                initial_values['Pressure'] = inlet_bc['Pressure']
            else:
                raise RuntimeError("Inlet type not appropriate to determine initial velocity and pressure.")

        if physics['Thermal'] == 'Energy':
            if inlet_bc['BoundaryType'] == 'inlet':
                if inlet_bc['ThermalBoundaryType'] == 'fixedValue':
                    initial_values['Temperature'] = inlet_bc['Temperature']
                else:
                    raise RuntimeError("Inlet type not appropriate to determine initial temperature")
            elif inlet_bc['BoundaryType'] == 'farField':
                initial_values['Temperature'] = inlet_bc['Temperature']
            else:
                raise RuntimeError("Inlet type not appropriate to determine initial temperature.")

        # Copy turbulence settings
        if physics['TurbulenceModel'] is not None:
            if initial_values['UseInletTurbulenceValues']:
                if inlet_bc['TurbulenceInletSpecification'] == 'TKEAndSpecDissipationRate':
                    initial_values['k'] = inlet_bc['TurbulentKineticEnergy']
                    initial_values['omega'] = inlet_bc['SpecificDissipationRate']
                elif inlet_bc['TurbulenceInletSpecification'] == 'intensityAndLengthScale':
                    if inlet_bc['BoundarySubtype'] == 'uniformVelocity' or \
                       inlet_bc['BoundarySubtype'] == 'characteristic':
                        Uin = (inlet_bc['Ux']**2 +
                               inlet_bc['Uy']**2 +
                               inlet_bc['Uz']**2)**0.5
                        I = inlet_bc['TurbulenceIntensity']
                        k = 3/2*(Uin*I)**2
                        Cmu = 0.09  # Standard turb model parameter
                        l = inlet_bc['TurbulenceLengthScale']
                        omega = k**0.5/(Cmu**0.25*l)
                        initial_values['k'] = k
                        initial_values['omega'] = omega
                    else:
                        raise RuntimeError(
                            "Inlet type currently unsupported for copying turbulence initial conditions.")
                else:
                    raise RuntimeError(
                        "Turbulence inlet specification currently unsupported for copying turbulence initial conditions")

processInitialConditions()

def clearCase(backup_path=None):
    """ Remove and recreate case directory, optionally backing up """
    output_path = case_folder
    if backup_path and os.path.isdir(output_path):
        shutil.move(output_path, backup_path)
    if os.path.isdir(output_path):
        shutil.rmtree(output_path)
    os.makedirs(output_path)  # mkdir -p

clearCase()

def exportZoneStlSurfaces():
        for zo in zone_objs:
            for i in range(len(zo.partNameList)):
                #shape = zo.shapeList[i].Shape
                path = os.path.join(solver_obj.WorkingDir,
                                    solver_obj.InputCaseName,
                                    "constant",
                                    "triSurface")
                if not os.path.exists(path):
                    os.makedirs(path)
                fname = os.path.join(path, zo.partNameList[i]+u".stl")
                sel_obj = analysis.Document.getObject(zo.partNameList[i])
                shape = sel_obj.Shape
                #meshStl = MeshPart.meshFromShape(shape, LinearDeflection = self.mesh_obj.STLLinearDeflection)
                meshStl = MeshPart.meshFromShape(shape, LinearDeflection = 0.1)
                meshStl.write(fname)
                print("Successfully wrote stl surface\n")

exportZoneStlSurfaces()

def processPorousZoneProperties():
        
        settings['porousZonesPresent'] = True
        porousZoneSettings = settings['porousZones']
        for po in porousZone_objs:
            pd = {'PartNameList': tuple(po.partNameList)}
            if po.porousZoneProperties['PorousCorrelation'] == 'DarcyForchheimer':
                pd['D'] = tuple(po.porousZoneProperties['D'])
                pd['F'] = tuple(po.porousZoneProperties['F'])
                pd['e1'] = tuple(po.porousZoneProperties['e1'])
                pd['e3'] = tuple(po.porousZoneProperties['e3'])
            elif po.porousZoneProperties['PorousCorrelation'] == 'Jakob':
                # Calculate effective Darcy-Forchheimer coefficients
                # This is for equilateral triangles arranged with the triangles pointing in BundleLayerNormal
                # direction (direction of greater spacing - sqrt(3)*triangleEdgeLength)
                pd['e1'] = tuple(po.porousZoneProperties['SpacingDirection'])  # OpenFOAM modifies to be orthog to e3
                pd['e3'] = tuple(po.porousZoneProperties['TubeAxis'])
                spacing = po.porousZoneProperties['TubeSpacing']
                d0 = po.porousZoneProperties['OuterDiameter']
                u0 = po.porousZoneProperties['VelocityEstimate']
                aspectRatio = po.porousZoneProperties['AspectRatio']
                kinVisc = settings['fluidProperties'][0]['KinematicViscosity']
                if kinVisc == 0.0:
                    raise ValueError("Viscosity must be set for Jakob correlation")
                if spacing < d0:
                    raise ValueError("Tube spacing may not be less than diameter")
                D = [0, 0, 0]
                F = [0, 0, 0]
                for (i, Sl, St) in [(0, aspectRatio*spacing, spacing), (1, spacing, aspectRatio*spacing)]:
                    C = 1.0/St*0.5*(1.0+0.47/(Sl/d0-1)**1.06)*(1.0/(1-d0/Sl))**(2.0-0.16)
                    D = C/d0*0.5*(u0*d0/kinVisc)**(1.0-0.16)
                    F = C*(u0*d0/kinVisc)**(-0.16)
                    D[i] = D
                    F[i] = F
                pd['D'] = tuple(D)
                pd['F'] = tuple(F)
                # Currently assuming zero drag parallel to tube bundle (3rd principal dirn)
            else:
                raise RuntimeError("Unrecognised method for porous baffle resistance")
            porousZoneSettings[po.Label] = pd

def processInitialisationZoneProperties():
        
        if settings['solver']['solverName'] in ['interFoam', 'multiphaseInterFoam']:
            # Make sure the first n-1 alpha values exist, and write the n-th one
            # consistently for multiphaseInterFoam
            for zone_name in settings['initialisationZones']:
                z = settings['initialisationZones'][zone_name]
                sum_alpha = 0.0
                if 'alphas' in z:
                    alphas_new = {}
                    for i, m in enumerate(settings['fluidProperties']):
                        alpha_name = m['Name']
                        if i == len(settings['fluidProperties'])-1:
                            if settings['solver']['solverName'] == 'multiphaseInterFoam':
                                alphas_new[alpha_name] = 1.0-sum_alpha
                        else:
                            alpha = z['alphas'].get(alpha_name, 0.0)
                            alphas_new[alpha_name] = alpha
                            sum_alpha += alpha
                    z['alphas'] = alphas_new

if porousZone_objs:
   processPorousZoneProperties()
processInitialisationZoneProperties()

TemplateBuilder.TemplateBuilder(case_folder, template_path, settings)

def floatEqual(a, b):
    """ Test whether a and b are equal within an absolute and relative tolerance """
    reltol = 10*sys.float_info.epsilon
    abstol = 1e-12  # Seems to be necessary on file read/write
    return abs(a-b) < abstol or abs(a - b) <= reltol*max(abs(a), abs(b))

def isSameGeometry(shape1, shape2):
    """ Copy of FemMeshTools.is_same_geometry, with fixes """
    # Check Area, CenterOfMass because non-planar shapes might not have more than one vertex defined
    same_Vertexes = 0
    # Bugfix: below was 1 - did not work for non-planar shapes
    if  len(shape1.Vertexes) > 0:
        # compare CenterOfMass
        # Bugfix: Precision seems to be lost on load/save
        if not floatEqual(shape1.CenterOfMass[0], shape2.CenterOfMass[0]) or \
                not floatEqual(shape1.CenterOfMass[1], shape2.CenterOfMass[1]) or \
                not floatEqual(shape1.CenterOfMass[2], shape2.CenterOfMass[2]):
            return False
        elif not floatEqual(shape1.Area, shape2.Area):
            return False
        else:
            # compare the Vertices
            for vs1 in shape1.Vertexes:
                for vs2 in shape2.Vertexes:
                    if floatEqual(vs1.X, vs2.X) and floatEqual(vs1.Y, vs2.Y) and floatEqual(vs1.Z, vs2.Z):
                        same_Vertexes += 1
                        # Bugfix: was 'continue' - caused false-negative with repeated vertices
                        break
            if same_Vertexes == len(shape1.Vertexes):
                return True
            else:
                return False

def matchFacesToTargetShape(ref_lists, shape):
    """ This function does a geometric matching of groups of faces much faster than doing face-by-face search
    :param ref_lists: List of lists of references - outer list is 'group' (e.g. boundary); refs are tuples
    :param shape: The shape to map to
    :return:  A list of tuples: (group index, reference) of matching refs for each face in shape
    """
    # Preserve original indices
    mesh_face_list = zip(shape.Faces, range(len(shape.Faces)))
    src_face_list = []
    for i, rl in enumerate(ref_lists):
        for br in rl:
            obj = FreeCAD.ActiveDocument.getObject(br[0])
            if not obj:
                raise RuntimeError("Referenced object '{}' not found - object may "
                                   "have been deleted".format(br[0]))
            try:
                bf = obj.Shape.getElement(br[1])
            except Part.OCCError:
                raise RuntimeError("Referenced face '{}:{}' not found - face may "
                                   "have been deleted".format(br[0], br[1]))
            src_face_list.append((bf, i, br))

    def compFn(x, y):
        if floatEqual(x, y):
            return 0
        elif x < y:
            return -1
        else:
            return 1

    # Sort boundary face list by centre of mass, x then y then z in case all in plane
    src_face_list.sort(cmp=compFn, key=lambda bf: bf[0].CenterOfMass.z)
    src_face_list.sort(cmp=compFn, key=lambda bf: bf[0].CenterOfMass.y)
    src_face_list.sort(cmp=compFn, key=lambda bf: bf[0].CenterOfMass.x)

    # Same sorting on mesh face list
    mesh_face_list.sort(cmp=compFn, key=lambda mf: mf[0].CenterOfMass.z)
    mesh_face_list.sort(cmp=compFn, key=lambda mf: mf[0].CenterOfMass.y)
    mesh_face_list.sort(cmp=compFn, key=lambda mf: mf[0].CenterOfMass.x)

    # Find faces with matching CofM
    i = 0
    j = 0
    j_match_start = 0
    matching = False
    candidate_mesh_faces = []
    for mf in mesh_face_list:
        candidate_mesh_faces.append([])
    while i < len(src_face_list) and j < len(mesh_face_list):
        bf = src_face_list[i][0]
        mf = mesh_face_list[j][0]
        if floatEqual(bf.CenterOfMass.x, mf.CenterOfMass.x):
            if floatEqual(bf.CenterOfMass.y, mf.CenterOfMass.y):
                if floatEqual(bf.CenterOfMass.z, mf.CenterOfMass.z):
                    candidate_mesh_faces[j].append((i, src_face_list[i][1], src_face_list[i][2]))
                    cmp = 0
                else:
                    cmp = (-1 if bf.CenterOfMass.z < mf.CenterOfMass.z else 1)
            else:
                cmp = (-1 if bf.CenterOfMass.y < mf.CenterOfMass.y else 1)
        else:
            cmp = (-1 if bf.CenterOfMass.x < mf.CenterOfMass.x else 1)
        if cmp == 0:
            if not matching:
                j_match_start = j
            j += 1
            matching = True
        elif cmp < 0:
            i += 1
            if matching:
                j = j_match_start
            matching = False
        elif cmp > 0:
            j += 1
            matching = False

    # Do comprehensive matching, and reallocate to original index
    successful_candidates = []
    for mf in mesh_face_list:
        successful_candidates.append([])
    for j in range(len(candidate_mesh_faces)):
        for k in range(len(candidate_mesh_faces[j])):
            i, nb, bref = candidate_mesh_faces[j][k]
            if isSameGeometry(src_face_list[i][0], mesh_face_list[j][0]):
                orig_idx = mesh_face_list[j][1]
                successful_candidates[orig_idx].append((nb, bref))

    return successful_candidates

def setupPatchNames():
        print ('Populating createPatchDict to update BC names')
        import CfdMeshTools
        # Init in case not meshed yet
        CfdMeshTools.CfdMeshTools(mesh_obj)
        
        settings['createPatches'] = {}
        
        mobj = mesh_obj

        # Make list of list of all boundary references for their corresponding boundary
        boundary_ref_lists = []
        for bc_id, bc_obj in enumerate(bc_group):
            boundary_ref_lists.append(bc_obj.References)

        # Match them up with faces in the meshed part
        matched_faces = matchFacesToTargetShape(boundary_ref_lists, mobj.Part.Shape)

        bc_lists = []
        for bc in bc_group:
            bc_lists.append([])
        for i in range(len(matched_faces)):
            if matched_faces[i]:
                nb, bref = matched_faces[i][0]
                bc_lists[nb].append(mobj.ShapeFaceNames[i])
                for k in range(len(matched_faces[i])-1):
                    nb2, bref2 = matched_faces[i][k+1]
                    

        for bc_id, bc_obj in enumerate(bc_group):
            bcDict = bc_obj.BoundarySettings
            bcType = bcDict["BoundaryType"]
            bcSubType = bcDict["BoundarySubtype"]
            patchType = CfdTools.getPatchType(bcType, bcSubType)
            settings['createPatches'][bc_obj.Label] = {
                'PatchNamesList': tuple(bc_lists[bc_id]),  # Tuple used so that case writer outputs as an array
                'PatchType': patchType
            }

        if mesh_obj.MeshRegionList:
            for regionObj in mesh_obj.MeshRegionList:
                if regionObj.Baffle:
                    settings['createPatchesFromSnappyBaffles'] = True

        if settings['createPatchesFromSnappyBaffles']:
            settings['createPatchesSnappyBaffles'] = {}
            # TODO Still need to include an error checker in the event that 
            # an internal baffle is created using snappy but is not linked up
            # with a baffle boundary condition (as in there is no baffle boundary condition which 
            # corresponds. Currently openfoam will throw a contextually
            # confusing error (only that the boundary does not exist). The primary difficulty with such a checker is 
            # that it is possible to define a boundary face as a baffle, which will be overridden
            # by the actual boundary name and therefore won't exist anymore. 
            for bc_id, bc_obj in enumerate(bc_group):
                bcDict = bc_obj.BoundarySettings
                bcType = bcDict["BoundaryType"]
                if bcType == "baffle":
                    tempBaffleList = []
                    tempBaffleListSlave = []
                    if mesh_obj.MeshRegionList:  # Can this if statement not be lumped with previous?
                        for regionObj in mesh_obj.MeshRegionList:
                            # print regionObj.Name
                            if regionObj.Baffle:
                                for sub in regionObj.References:
                                    # print sub[0].Name
                                    elems = sub[1]
                                    elt = FreeCAD.ActiveDocument.getObject(sub[0]).Shape.getElement(elems)
                                    if elt.ShapeType == 'Face':
                                        bcFacesList = bc_obj.Shape.Faces
                                        for bf in bcFacesList:
                                            isSameGeo = CfdTools.isSameGeometry(bf, elt)
                                            if isSameGeo:
                                                tempBaffleList.append(regionObj.Name+sub[0]+elems)
                                                tempBaffleListSlave.append(regionObj.Name+sub[0]+elems+"_slave")
                    settings['createPatchesSnappyBaffles'][bc_obj.Label] = {"PatchNamesList" : tuple(tempBaffleList),
                                                                            "PatchNamesListSlave" : tuple(tempBaffleListSlave)}

        # Add default faces
        flagName = False
        def_bc_list = []
        for i in range(len(matched_faces)):
            if not matched_faces[i]:
                def_bc_list.append(mobj.ShapeFaceNames[i])
                flagName = True
        if flagName:
            settings['createPatches']['defaultFaces'] = {
                'PatchNamesList': tuple(def_bc_list),
                'PatchType': "patch"
            }

setupPatchNames()
TemplateBuilder.TemplateBuilder(case_folder, template_path, settings)

def writeMesh():
        """ Convert or copy mesh files """
        if mesh_obj.Proxy.Type == "CfdMesh":
            import CfdMeshTools
            # Move Cartesian mesh files from temporary mesh directory to case directory
            cart_mesh = CfdMeshTools.CfdMeshTools(mesh_obj)
            cart_mesh = cart_mesh
            if mesh_obj.MeshUtility == "cfMesh":
                print("Writing Cartesian mesh\n")
                # cart_mesh.get_tmp_file_paths("cfMesh")  # Update tmp file locations
                cart_mesh.get_tmp_file_paths()  # Update tmp file locations
                CfdTools.copyFilesRec(cart_mesh.polyMeshDir, os.path.join(case_folder, 'constant', 'polyMesh'))
                CfdTools.copyFilesRec(cart_mesh.triSurfaceDir, os.path.join(case_folder, 'constant', 'triSurface'))
                # shutil.copy2(cart_mesh.temp_file_meshDict, os.path.join(self.case_folder,'system'))
                shutil.copy2(os.path.join(cart_mesh.meshCaseDir, 'system', 'meshDict'),
                             os.path.join(case_folder,'system'))
                shutil.copy2(os.path.join(cart_mesh.meshCaseDir,'Allmesh'),case_folder)
                shutil.copy2(os.path.join(cart_mesh.meshCaseDir,'log.cartesianMesh'),case_folder)
                shutil.copy2(os.path.join(cart_mesh.meshCaseDir,'log.surfaceFeatureEdges'),case_folder)

            elif mesh_obj.MeshUtility == "snappyHexMesh":
                print("Writing snappyHexMesh generated Cartesian mesh\n")
                # cart_mesh.get_tmp_file_paths("snappyHexMesh")  # Update tmp file locations
                cart_mesh.get_tmp_file_paths()  # Update tmp file locations
                CfdTools.copyFilesRec(cart_mesh.polyMeshDir, os.path.join(case_folder,'constant','polyMesh'))
                CfdTools.copyFilesRec(cart_mesh.triSurfaceDir, os.path.join(case_folder,'constant','triSurface'))
                # shutil.copy2(cart_mesh.temp_file_blockMeshDict, os.path.join(self.case_folder,'system'))
                # shutil.copy2(cart_mesh.temp_file_snappyMeshDict, os.path.join(self.case_folder,'system'))
                # shutil.copy2(cart_mesh.temp_file_surfaceFeatureExtractDict, os.path.join(self.case_folder,'system'))
                shutil.copy2(os.path.join(cart_mesh.meshCaseDir, 'system', 'blockMeshDict'),
                             os.path.join(case_folder,'system'))
                shutil.copy2(os.path.join(cart_mesh.meshCaseDir, 'system', 'snappyHexMeshDict'),
                             os.path.join(case_folder,'system'))
                shutil.copy2(os.path.join(cart_mesh.meshCaseDir, 'system', 'surfaceFeatureExtractDict'),
                             os.path.join(case_folder,'system'))
                shutil.copy2(os.path.join(cart_mesh.meshCaseDir, 'Allmesh'), case_folder)
                shutil.copy2(os.path.join(cart_mesh.meshCaseDir, 'log.blockMesh'), case_folder)
                shutil.copy2(os.path.join(cart_mesh.meshCaseDir, 'log.surfaceFeatureExtract'), case_folder)
                shutil.copy2(os.path.join(cart_mesh.meshCaseDir, 'log.snappyHexMesh'), case_folder)

            elif mesh_obj.MeshUtility == "gmsh":
                print("Writing gmsh generated mesh\n")
                cart_mesh.get_tmp_file_paths()  # Update tmp file locations
                CfdTools.copyFilesRec(cart_mesh.polyMeshDir, os.path.join(case_folder,'constant','polyMesh'))
                CfdTools.copyFilesRec(cart_mesh.triSurfaceDir, os.path.join(case_folder,'constant','gmsh'))
                shutil.copy2(os.path.join(cart_mesh.meshCaseDir, 'Allmesh'), case_folder)
                shutil.copy2(os.path.join(cart_mesh.meshCaseDir, 'log.gmshToFoam'), case_folder)

            if mesh_obj.ElementDimension == '2D':
                shutil.copy2(os.path.join(os.path.join(cart_mesh.meshCaseDir,'system'),'extrudeMeshDict'), 
                             os.path.join(case_folder, 'system'))
                shutil.copy2(os.path.join(cart_mesh.meshCaseDir, 'log.extrudeMesh'), case_folder)
        else:
            raise RuntimeError("Unrecognised mesh type")
writeMesh()
fname = os.path.join(case_folder, "Allrun")
import stat
s = os.stat(fname)
os.chmod(fname, s.st_mode | stat.S_IEXEC)
CfdTools.movePolyMesh(case_folder)
#move the case to desired location
os.system("cd /tmp")
os.system("cp -r /tmp/case/ /home/weibin/freecad_create_openfoam/")
#start the simulation
case="/home/weibin/freecad_create_openfoam/case"
from PyFoam . Applications . Decomposer import Decomposer
from PyFoam . Applications . Runner import Runner
from PyFoam . Applications . PlotRunner import PlotRunner
#Runner (args=["--progress","Decomposer","proc",8,"-case",case])
os.system("cd /home/weibin/freecad_create_openfoam/case && bash Allrun")

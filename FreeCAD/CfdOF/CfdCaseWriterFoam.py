# ***************************************************************************
# *                                                                         *
# *   Copyright (c) 2015 - Qingfeng Xia <qingfeng.xia eng ox ac uk>         *
# *   Copyright (c) 2017 - Alfred Bogaers (CSIR) <abogaers@csir.co.za>      *
# *   Copyright (c) 2017 - Johan Heyns (CSIR) <jheyns@csir.co.za>           *
# *   Copyright (c) 2017 - Oliver Oxtoby (CSIR) <ooxtoby@csir.co.za>        *
# *                                                                         *
# *   This program is free software; you can redistribute it and/or modify  *
# *   it under the terms of the GNU Lesser General Public License (LGPL)    *
# *   as published by the Free Software Foundation; either version 2 of     *
# *   the License, or (at your option) any later version.                   *
# *   for detail see the LICENCE text file.                                 *
# *                                                                         *
# *   This program is distributed in the hope that it will be useful,       *
# *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
# *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
# *   GNU Library General Public License for more details.                  *
# *                                                                         *
# *   You should have received a copy of the GNU Library General Public     *
# *   License along with this program; if not, write to the Free Software   *
# *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  *
# *   USA                                                                   *
# *                                                                         *
# ***************************************************************************

import FreeCAD
import CfdTools
from CfdTools import cfdMessage
import os
import os.path
import shutil
from PySide import QtCore
from PySide.QtCore import QRunnable, QObject
import Units
import TemplateBuilder


class CfdCaseWriterFoam:
    def __init__(self, analysis_obj):
        self.analysis_obj = analysis_obj
        self.solver_obj = CfdTools.getSolver(analysis_obj)
        self.physics_model, isPresent = CfdTools.getPhysicsModel(analysis_obj)
        self.mesh_obj = CfdTools.getMesh(analysis_obj)
        self.material_objs = CfdTools.getMaterials(analysis_obj)
        self.bc_group = CfdTools.getCfdBoundaryGroup(analysis_obj)
        self.initial_conditions, isPresent = CfdTools.getInitialConditions(analysis_obj)
        self.porousZone_objs = CfdTools.getPorousZoneObjects(analysis_obj)
        self.initialisationZone_objs = CfdTools.getInitialisationZoneObjects(analysis_obj)
        self.zone_objs = CfdTools.getZoneObjects(analysis_obj)
        self.mesh_generated = False

    def writeCase(self):
        """ writeCase() will collect case settings, and finally build a runnable case. """
        cfdMessage("Start to write case to folder {}\n".format(self.solver_obj.WorkingDir))
        cwd = os.curdir
        if not os.path.exists(self.solver_obj.WorkingDir):
            raise IOError("Path " + self.solver_obj.WorkingDir + " does not exist.")

        # Perform initialisation here rather than __init__ in case of path changes
        self.case_folder = os.path.join(self.solver_obj.WorkingDir, self.solver_obj.InputCaseName)
        self.case_folder = os.path.expanduser(os.path.abspath(self.case_folder))
        self.mesh_file_name = os.path.join(self.case_folder, self.solver_obj.InputCaseName, u".unv")

        self.template_path = os.path.join(CfdTools.get_module_path(), "data", "defaults")

        solverSettingsDict = CfdTools.getSolverSettings(self.solver_obj)

        # Collect settings into single dictionary
        if not self.mesh_obj:
            raise RuntimeError("No mesh object found in analysis")
        phys_settings = dict(zip(self.physics_model.PropertiesList,
                                 (getattr(self.physics_model, prop) for prop in self.physics_model.PropertiesList)))
        self.settings = {
            'physics': phys_settings,
            'fluidProperties': [],  # Order is important, so use a list
            'initialValues': self.initial_conditions,
            'boundaries': dict((b.Label, b.BoundarySettings) for b in self.bc_group),
            'bafflesPresent': self.bafflesPresent(),
            'porousZones': {},
            'porousZonesPresent': False,
            'initialisationZones': {o.Label: o.initialisationZoneProperties for o in self.initialisationZone_objs},
            'initialisationZonesPresent': len(self.initialisationZone_objs) > 0,
            'zones': {o.Label: {'PartNameList': tuple(o.partNameList)} for o in self.zone_objs},
            'zonesPresent': len(self.zone_objs) > 0,
            'meshType': self.mesh_obj.Proxy.Type,
            'solver': solverSettingsDict,
            'system': {},
            'runChangeDictionary':False
            }

        self.processSystemSettings()
        self.processSolverSettings()
        self.processFluidProperties()
        self.processBoundaryConditions()
        self.processInitialConditions()
        self.clearCase()

        self.exportZoneStlSurfaces()
        if self.porousZone_objs:
            self.processPorousZoneProperties()
        self.processInitialisationZoneProperties()

        self.settings['createPatchesFromSnappyBaffles'] = False
        cfdMessage("Matching boundary conditions ...\n")
        self.setupPatchNames()

        TemplateBuilder.TemplateBuilder(self.case_folder, self.template_path, self.settings)
        self.writeMesh()

        # Update Allrun permission - will fail silently on Windows
        fname = os.path.join(self.case_folder, "Allrun")
        import stat
        s = os.stat(fname)
        os.chmod(fname, s.st_mode | stat.S_IEXEC)

        # Move mesh files, after being edited, to polyMesh.org
        CfdTools.movePolyMesh(self.case_folder)

        cfdMessage("Successfully wrote {} case to folder {}\n".format(
                   self.solver_obj.SolverName, self.solver_obj.WorkingDir))
        return True

    def getSolverName(self):
        """ Solver name is selected based on selected physics. This should only be extended as additional physics are
        included. """
        solver = None
        if self.physics_model.Phase == 'Single':
            if len(self.material_objs) == 1:
                if self.physics_model.Flow == 'Incompressible':
                    if self.physics_model.Thermal == 'None':
                        if self.physics_model.Time == 'Transient':
                            solver = 'pimpleFoam'
                        else:
                            if self.porousZone_objs or self.porousBafflesPresent():
                                solver = 'porousSimpleFoam'
                            else:
                                solver = 'simpleFoam'
                    else:
                        raise RuntimeError("Only isothermal simulation currently supported for incompressible flow.")
                elif self.physics_model.Flow == 'HighMachCompressible':
                    solver = 'hisa'
                else:
                    raise RuntimeError(self.physics_model.Flow + " flow model currently not supported.")
            else:
                raise RuntimeError("Only one material object may be present for single phase simulation.")
        elif self.physics_model.Phase == 'FreeSurface':
            if self.physics_model.Time == 'Transient':
                if self.physics_model.Thermal == 'None':
                    if len(self.material_objs) == 2:
                        solver = 'interFoam'
                    elif len(self.material_objs) > 2:
                        solver = 'multiphaseInterFoam'
                    else:
                        raise RuntimeError("At least two material objects must be present for free surface simulation.")
                else:
                    raise RuntimeError("Only isothermal analysis is currently supported for free surface flow simulation.")
            else:
                raise RuntimeError("Only transient analysis is supported for free surface flow simulation.")
        else:
            raise RuntimeError(self.physics_model.Phase + " phase model currently not supported.")

        # Catch-all in case
        if solver is None:
            raise RuntimeError("No solver is supported to handle the selected physics with {} phases.".format(
                len(self.material_objs)))
        return solver

    def processSolverSettings(self):
        solver_settings = self.settings['solver']
        if solver_settings['parallel']:
            if solver_settings['parallelCores'] < 2:
                solver_settings['parallelCores'] = 2
        solver_settings['solverName'] = self.getSolverName()

    def processSystemSettings(self):
        system_settings = self.settings['system']
        system_settings['FoamRuntime'] = CfdTools.getFoamRuntime()
        system_settings['CasePath'] = self.case_folder
        system_settings['TranslatedCasePath'] = CfdTools.translatePath(self.case_folder)
        system_settings['FoamPath'] = CfdTools.getFoamDir()
        system_settings['TranslatedFoamPath'] = CfdTools.translatePath(CfdTools.getFoamDir())

    def clearCase(self, backup_path=None):
        """ Remove and recreate case directory, optionally backing up """
        output_path = self.case_folder
        if backup_path and os.path.isdir(output_path):
            shutil.move(output_path, backup_path)
        if os.path.isdir(output_path):
            shutil.rmtree(output_path)
        os.makedirs(output_path)  # mkdir -p

    # Mesh

    def writeMesh(self):
        """ Convert or copy mesh files """
        if self.mesh_obj.Proxy.Type == "CfdMesh":
            import CfdMeshTools
            # Move Cartesian mesh files from temporary mesh directory to case directory
            self.cart_mesh = CfdMeshTools.CfdMeshTools(self.mesh_obj)
            cart_mesh = self.cart_mesh
            if self.mesh_obj.MeshUtility == "cfMesh":
                print("Writing Cartesian mesh\n")
                # cart_mesh.get_tmp_file_paths("cfMesh")  # Update tmp file locations
                cart_mesh.get_tmp_file_paths()  # Update tmp file locations
                CfdTools.copyFilesRec(cart_mesh.polyMeshDir, os.path.join(self.case_folder, 'constant', 'polyMesh'))
                CfdTools.copyFilesRec(cart_mesh.triSurfaceDir, os.path.join(self.case_folder, 'constant', 'triSurface'))
                # shutil.copy2(cart_mesh.temp_file_meshDict, os.path.join(self.case_folder,'system'))
                shutil.copy2(os.path.join(cart_mesh.meshCaseDir, 'system', 'meshDict'),
                             os.path.join(self.case_folder,'system'))
                shutil.copy2(os.path.join(cart_mesh.meshCaseDir,'Allmesh'),self.case_folder)
                shutil.copy2(os.path.join(cart_mesh.meshCaseDir,'log.cartesianMesh'),self.case_folder)
                shutil.copy2(os.path.join(cart_mesh.meshCaseDir,'log.surfaceFeatureEdges'),self.case_folder)

            elif self.mesh_obj.MeshUtility == "snappyHexMesh":
                print("Writing snappyHexMesh generated Cartesian mesh\n")
                # cart_mesh.get_tmp_file_paths("snappyHexMesh")  # Update tmp file locations
                cart_mesh.get_tmp_file_paths()  # Update tmp file locations
                CfdTools.copyFilesRec(cart_mesh.polyMeshDir, os.path.join(self.case_folder,'constant','polyMesh'))
                CfdTools.copyFilesRec(cart_mesh.triSurfaceDir, os.path.join(self.case_folder,'constant','triSurface'))
                # shutil.copy2(cart_mesh.temp_file_blockMeshDict, os.path.join(self.case_folder,'system'))
                # shutil.copy2(cart_mesh.temp_file_snappyMeshDict, os.path.join(self.case_folder,'system'))
                # shutil.copy2(cart_mesh.temp_file_surfaceFeatureExtractDict, os.path.join(self.case_folder,'system'))
                shutil.copy2(os.path.join(cart_mesh.meshCaseDir, 'system', 'blockMeshDict'),
                             os.path.join(self.case_folder,'system'))
                shutil.copy2(os.path.join(cart_mesh.meshCaseDir, 'system', 'snappyHexMeshDict'),
                             os.path.join(self.case_folder,'system'))
                shutil.copy2(os.path.join(cart_mesh.meshCaseDir, 'system', 'surfaceFeatureExtractDict'),
                             os.path.join(self.case_folder,'system'))
                shutil.copy2(os.path.join(cart_mesh.meshCaseDir, 'Allmesh'), self.case_folder)
                shutil.copy2(os.path.join(cart_mesh.meshCaseDir, 'log.blockMesh'), self.case_folder)
                shutil.copy2(os.path.join(cart_mesh.meshCaseDir, 'log.surfaceFeatureExtract'), self.case_folder)
                shutil.copy2(os.path.join(cart_mesh.meshCaseDir, 'log.snappyHexMesh'), self.case_folder)

            elif self.mesh_obj.MeshUtility == "gmsh":
                print("Writing gmsh generated mesh\n")
                cart_mesh.get_tmp_file_paths()  # Update tmp file locations
                CfdTools.copyFilesRec(cart_mesh.polyMeshDir, os.path.join(self.case_folder,'constant','polyMesh'))
                CfdTools.copyFilesRec(cart_mesh.triSurfaceDir, os.path.join(self.case_folder,'constant','gmsh'))
                shutil.copy2(os.path.join(cart_mesh.meshCaseDir, 'Allmesh'), self.case_folder)
                shutil.copy2(os.path.join(cart_mesh.meshCaseDir, 'log.gmshToFoam'), self.case_folder)

            if self.mesh_obj.ElementDimension == '2D':
                shutil.copy2(os.path.join(os.path.join(cart_mesh.meshCaseDir,'system'),'extrudeMeshDict'), 
                             os.path.join(self.case_folder, 'system'))
                shutil.copy2(os.path.join(cart_mesh.meshCaseDir, 'log.extrudeMesh'), self.case_folder)
        else:
            raise RuntimeError("Unrecognised mesh type")

    def setupMesh(self, updated_mesh_path, scale):
        if os.path.exists(updated_mesh_path):
            CfdTools.convertMesh(self.case_folder, updated_mesh_path, scale)

    def processFluidProperties(self):
        # self.material_obj currently stores everything as a string
        # Convert to (mostly) SI numbers for OpenFOAM
        settings = self.settings
        for material_obj in self.material_objs:
            mp = {}
            mp['Name'] = material_obj.Label
            if 'Density' in material_obj.PropertiesList:
                mp['Density'] = Units.Quantity(material_obj.Density).getValueAs("kg/m^3").Value
            if 'DynamicViscosity' in material_obj.PropertiesList:
                if self.physics_model.Turbulence == 'Inviscid':
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

    def processBoundaryConditions(self):
        """ Compute any quantities required before case build """
        settings = self.settings
        for bc_name in settings['boundaries']:
            bc = settings['boundaries'][bc_name]
            if not bc['VelocityIsCartesian']:
                veloMag = bc['VelocityMag']
                face = bc['DirectionFace'].split(':')
                # See if entered face actually exists and is planar
                try:
                    selected_object = self.analysis_obj.Document.getObject(face[0])
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

    def processInitialConditions(self):
        """ Do any required computations before case build. Boundary conditions must be processed first. """
        settings = self.settings
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

    # Zones

    def exportZoneStlSurfaces(self):
        for zo in self.zone_objs:
            import Mesh
            for i in range(len(zo.partNameList)):
                #shape = zo.shapeList[i].Shape
                path = os.path.join(self.solver_obj.WorkingDir,
                                    self.solver_obj.InputCaseName,
                                    "constant",
                                    "triSurface")
                if not os.path.exists(path):
                    os.makedirs(path)
                fname = os.path.join(path, zo.partNameList[i]+u".stl")
                import MeshPart
                sel_obj = self.analysis_obj.Document.getObject(zo.partNameList[i])
                shape = sel_obj.Shape
                #meshStl = MeshPart.meshFromShape(shape, LinearDeflection = self.mesh_obj.STLLinearDeflection)
                meshStl = MeshPart.meshFromShape(shape, LinearDeflection = 0.1)
                meshStl.write(fname)
                print("Successfully wrote stl surface\n")

    def processPorousZoneProperties(self):
        settings = self.settings
        settings['porousZonesPresent'] = True
        porousZoneSettings = settings['porousZones']
        for po in self.porousZone_objs:
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
                kinVisc = self.settings['fluidProperties']['KinematicViscosity']
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

    def processInitialisationZoneProperties(self):
        settings = self.settings
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

    def bafflesPresent(self):
        for b in self.bc_group:
            if b.BoundarySettings['BoundaryType'] == 'baffle':
                return True
        return False

    def porousBafflesPresent(self):
        for b in self.bc_group:
            if b.BoundarySettings['BoundaryType'] == 'baffle' and \
               b.BoundarySettings['BoundarySubtype'] == 'porousBaffle':
                return True
        return False

    def setupPatchNames(self):
        print ('Populating createPatchDict to update BC names')
        import CfdMeshTools
        # Init in case not meshed yet
        CfdMeshTools.CfdMeshTools(self.mesh_obj)
        settings = self.settings
        settings['createPatches'] = {}
        bc_group = self.bc_group
        mobj = self.mesh_obj

        # Make list of list of all boundary references for their corresponding boundary
        boundary_ref_lists = []
        for bc_id, bc_obj in enumerate(bc_group):
            boundary_ref_lists.append(bc_obj.References)

        # Match them up with faces in the meshed part
        matched_faces = CfdTools.matchFacesToTargetShape(boundary_ref_lists, mobj.Part.Shape)

        bc_lists = []
        for bc in bc_group:
            bc_lists.append([])
        for i in range(len(matched_faces)):
            if matched_faces[i]:
                nb, bref = matched_faces[i][0]
                bc_lists[nb].append(mobj.ShapeFaceNames[i])
                for k in range(len(matched_faces[i])-1):
                    nb2, bref2 = matched_faces[i][k+1]
                    if nb2 != nb:
                        cfdMessage(
                            "Boundary '{}' reference {}:{} also assigned as "
                            "boundary '{}' reference {}:{}\n".format(
                                bc_group[nb].Label, bref[0], bref[1], bc_group[nb2].Label, bref2[0], bref2[1]))

        for bc_id, bc_obj in enumerate(bc_group):
            bcDict = bc_obj.BoundarySettings
            bcType = bcDict["BoundaryType"]
            bcSubType = bcDict["BoundarySubtype"]
            patchType = CfdTools.getPatchType(bcType, bcSubType)
            settings['createPatches'][bc_obj.Label] = {
                'PatchNamesList': tuple(bc_lists[bc_id]),  # Tuple used so that case writer outputs as an array
                'PatchType': patchType
            }

        if self.mesh_obj.MeshRegionList:
            for regionObj in self.mesh_obj.MeshRegionList:
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
                    if self.mesh_obj.MeshRegionList:  # Can this if statement not be lumped with previous?
                        for regionObj in self.mesh_obj.MeshRegionList:
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

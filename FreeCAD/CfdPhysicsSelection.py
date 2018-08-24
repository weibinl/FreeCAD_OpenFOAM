import os
import os.path
import FreeCAD
import CfdTools
from femcommands.manager import CommandManager

def makeCfdPhysicsSelection(name="PhysicsModel"):
    obj = FreeCAD.ActiveDocument.addObject("App::FeaturePython", name)
    _CfdPhysicsModel(obj)
    return obj


class _CfdPhysicsModel:
    """ The CFD Physics Model """
    def __init__(self, obj):
        obj.Proxy = self
        self.Type = "PhysicsModel"
        self.initProperties(obj)

    def initProperties(self, obj):
        # obj.supportedProperties()
        # ['App::PropertyBool', 'App::PropertyBoolList', 'App::PropertyFloat', 'App::PropertyFloatList',
        #  'App::PropertyFloatConstraint', 'App::PropertyPrecision', 'App::PropertyQuantity',
        #  'App::PropertyQuantityConstraint', 'App::PropertyAngle', 'App::PropertyDistance', 'App::PropertyLength',
        #  'App::PropertyArea', 'App::PropertyVolume', 'App::PropertySpeed', 'App::PropertyAcceleration',
        #  'App::PropertyForce', 'App::PropertyPressure', 'App::PropertyInteger', 'App::PropertyIntegerConstraint',
        #  'App::PropertyPercent', 'App::PropertyEnumeration', 'App::PropertyIntegerList', 'App::PropertyIntegerSet',
        #  'App::PropertyMap', 'App::PropertyString', 'App::PropertyUUID', 'App::PropertyFont',
        #  'App::PropertyStringList', 'App::PropertyLink', 'App::PropertyLinkChild', 'App::PropertyLinkGlobal',
        #  'App::PropertyLinkSub', 'App::PropertyLinkSubChild', 'App::PropertyLinkSubGlobal', 'App::PropertyLinkList',
        #  'App::PropertyLinkListChild', 'App::PropertyLinkListGlobal', 'App::PropertyLinkSubList',
        #  'App::PropertyLinkSubListChild', 'App::PropertyLinkSubListGlobal', 'App::PropertyMatrix',
        #  'App::PropertyVector', 'App::PropertyVectorDistance', 'App::PropertyPosition', 'App::PropertyDirection',
        #  'App::PropertyVectorList', 'App::PropertyPlacement', 'App::PropertyPlacementList',
        #  'App::PropertyPlacementLink', 'App::PropertyColor', 'App::PropertyColorList', 'App::PropertyMaterial',
        #  'App::PropertyMaterialList', 'App::PropertyPath', 'App::PropertyFile', 'App::PropertyFileIncluded',
        #  'App::PropertyPythonObject', 'App::PropertyExpressionEngine', 'Part::PropertyPartShape',
        #  'Part::PropertyGeometryList', 'Part::PropertyShapeHistory', 'Part::PropertyFilletEdges',
        #  'Fem::PropertyFemMesh', 'Fem::PropertyPostDataObject']

        if 'Time' not in obj.PropertiesList:
            obj.addProperty("App::PropertyEnumeration", "Time", "Physics modelling",
                            "Resolve time dependence")
            obj.Time = ['Steady', 'Transient']
            obj.Time = 'Steady'

        if 'Flow' not in obj.PropertiesList:
            obj.addProperty("App::PropertyEnumeration", "Flow", "Physics modelling",
                            "Flow algorithm")
            obj.Flow = ['Incompressible', 'Compressible', 'HighMachCompressible']
            obj.Flow = 'Incompressible'

        if 'Thermal' not in obj.PropertiesList:
            obj.addProperty("App::PropertyEnumeration", "Thermal", "Physics modelling",
                            "Thermal modelling")
            obj.Thermal = ['None', 'Buoyancy', 'Energy']
            obj.Thermal = 'None'

        if 'Phase' not in obj.PropertiesList:
            obj.addProperty("App::PropertyEnumeration", "Phase", "Physics modelling",
                            "Type of phases present")
            obj.Phase = ['Single', 'FreeSurface']
            obj.Phase = 'Single'

        if 'Turbulence' not in obj.PropertiesList:
            obj.addProperty("App::PropertyEnumeration", "Turbulence", "Physics modelling",
                        "Type of turbulence modelling")
            obj.Turbulence = ['Inviscid', 'Laminar', 'RANS']
            obj.Turbulence = 'Laminar'

        if 'TurbulenceModel' not in obj.PropertiesList:
            obj.addProperty("App::PropertyEnumeration", "TurbulenceModel", "Physics modelling",
                            "Turbulence model")
            obj.TurbulenceModel = ['kOmegaSST']

        if 'gx' not in obj.PropertiesList:
            obj.addProperty("App::PropertyAcceleration", "gx", "Physics modelling",
                            "Gravitational acceleration vector (x component)")
            obj.gx = '0 m/s^2'
        if 'gy' not in obj.PropertiesList:
            obj.addProperty("App::PropertyAcceleration", "gy", "Physics modelling",
                            "Gravitational acceleration vector (y component)")
            obj.gy = '-9.81 m/s^2'
        if 'gz' not in obj.PropertiesList:
            obj.addProperty("App::PropertyAcceleration", "gz", "Physics modelling",
                            "Gravitational acceleration vector (z component)")
            obj.gz = '0 m/s^2'

    def execute(self, obj):
        return




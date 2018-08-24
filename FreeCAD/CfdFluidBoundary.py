import FreeCAD
import _CfdFluidBoundary


def makeCfdFluidBoundary(name="CfdFluidBoundary"):
    ''' makeCfdFluidBoundary([name]): Creates a fluid boundary condition'''
    obj = FreeCAD.ActiveDocument.addObject("Part::FeaturePython", name)
    _CfdFluidBoundary._CfdFluidBoundary(obj)
    
    return obj


import FreeCAD
import Part


class PartFeature:
    "Part containing CfdFluidBoundary faces"
    def __init__(self, obj):
        obj.Proxy = self


class _CfdFluidBoundary(PartFeature):
    "The CfdFluidBoundary object"
    def __init__(self, obj):
        PartFeature.__init__(self, obj)

        obj.Proxy = self
        self.Type = "CfdFluidBoundary"
        obj.addProperty("App::PropertyPythonObject", "References")
        obj.addProperty("App::PropertyPythonObject", "BoundarySettings")

        # Default settings
        obj.References = []
        obj.BoundarySettings = {'BoundaryType': '',
                                'BoundarySubtype': '',
                                'VelocityIsCartesian': True,
                                'Ux': 0.0,
                                'Uy': 0.0,
                                'Uz': 0.0,
                                'VelocityMag': 0.0,
                                'DirectionFace': '',
                                'ReverseNormal': False,
                                'Pressure': 0.0,
                                'SlipRatio': 0.0,
                                'VolFlowRate': 0.0,
                                'MassFlowRate': 0.0,
                                'PorousBaffleMethod': 0,
                                'PressureDropCoeff': 0.0,
                                'ScreenWireDiameter': 0.0,
                                'ScreenSpacing': 0.0,
                                'ThermalBoundaryType': 'zeroGradient',
                                'Temperature': 293,
                                'HeatFlux': 0,
                                'HeatTransferCoeff': 0,
                                'TurbulenceInletSpecification': 'intensityAndLengthScale',
                                'TurbulentKineticEnergy': 0.01,
                                'SpecificDissipationRate': 1,
                                'TurbulenceIntensity': 0.1,
                                'TurbulenceLengthScale': 0.1,
                                'alphas': {}}

    def execute(self, obj):
        ''' Create compound part at recompute. '''
        docName = str(obj.Document.Name)
        doc = FreeCAD.getDocument(docName)
        listOfFaces = []
        for i in range(len(obj.References)):
            ref = obj.References[i]
            selection_object = doc.getObject(ref[0])
            if selection_object is not None:  # May have been deleted
                try:
                    listOfFaces.append(selection_object.Shape.getElement(ref[1]))
                except Part.OCCError:  # Face may have been deleted
                    pass
        if len(listOfFaces) > 0:
            obj.Shape = Part.makeCompound(listOfFaces)
        else:
            obj.Shape = Part.Shape()
        
        return

    def __getstate__(self):
        return None

    def __setstate__(self, state):
        return None



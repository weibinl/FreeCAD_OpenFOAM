# ***************************************************************************
# *                                                                         *
# *   Copyright (c) 2013-2015 - Juergen Riegel <FreeCAD@juergen-riegel.net> *
# *   Copyright (c) 2017 - Alfred Bogaers <abogaers@csir.co.za>             *
# *   Copyright (c) 2017 - Oliver Oxtoby <ooxtoby@csir.co.za>               *
# *   Copyright (c) 2017 - Johan Heyns <jheyns@csir.co.za>                  *
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

__title__ = "_TaskPanelCfdInitialiseInternalFlowField"
__author__ = ""
__url__ = "http://www.freecadweb.org"


import FreeCAD
import os
import sys
import os.path
import CfdTools
from CfdTools import inputCheckAndStore, setInputFieldQuantity
import Units

if FreeCAD.GuiUp:
    import FreeCADGui
    from PySide import QtCore
    from PySide import QtCore
    from PySide import QtGui
    from PySide.QtCore import Qt
    from PySide.QtGui import QApplication
    import FemGui


class _TaskPanelCfdInitialiseInternalFlowField:
    '''The editmode TaskPanel for InitialVariables objects'''
    def __init__(self, obj, physics_model, boundaries, material_objs):
        # Update object with latest properties for backward-compatibility
        obj.Proxy.initProperties(obj)

        FreeCADGui.Selection.clearSelection()
        self.obj = obj
        self.physicsModel = physics_model
        self.boundaries = boundaries
        self.material_objs = material_objs
        self.InitialVariables = self.obj.InitialVariables.copy()

        self.form = FreeCADGui.PySideUic.loadUi(os.path.join(os.path.dirname(__file__),
                                                             "TaskPanelCfdInitialiseInternalField.ui"))

        self.form.basicPropertiesFrame.setVisible(False)
        self.form.radioButtonPotentialFlow.toggled.connect(self.radioUPChanged)
        self.form.radioButtonUseInletValuesUP.toggled.connect(self.radioUPChanged)
        self.form.turbulencePropertiesFrame.setVisible(False)

        self.form.Ux.valueChanged.connect(self.UxChanged)
        self.form.Uy.valueChanged.connect(self.UyChanged)
        self.form.Uz.valueChanged.connect(self.UzChanged)
        self.form.pressure.valueChanged.connect(self.PChanged)

        self.form.comboFluid.currentIndexChanged.connect(self.comboFluidChanged)
        self.form.inputVolumeFraction.valueChanged.connect(self.inputVolumeFractionChanged)

        self.form.checkUseInletValuesThermal.stateChanged.connect(self.checkUseInletTemperatureValueChanged)
        self.form.inputTemperature.valueChanged.connect(self.inputTemperatureChanged)

        self.form.checkUseInletValues.stateChanged.connect(self.checkUseInletValuesChanged)
        self.form.comboInlets.currentIndexChanged.connect(self.comboInletsChanged)
        self.form.inputk.valueChanged.connect(self.inputkChanged)
        self.form.inputOmega.valueChanged.connect(self.inputOmegaChanged)

        self.populateUiBasedOnPhysics()

    def populateUiBasedOnPhysics(self):
        self.form.radioButtonPotentialFlow.setToolTip(
            "Initialise the velocity field using an incompressible, potential "
            "or irrotational flow assumption.")
        potential_foam = self.InitialVariables['PotentialFoam']
        use_inlet_UP = self.InitialVariables['UseInletUPValues']
        if potential_foam:
            self.form.radioButtonPotentialFlow.toggle()
        elif use_inlet_UP:
            self.form.radioButtonUseInletValuesUP.toggle()
        else:
            self.form.radioButtonSpecifyValues.toggle()

        setInputFieldQuantity(self.form.Ux, str(self.InitialVariables['Ux'])+"m/s")
        setInputFieldQuantity(self.form.Uy, str(self.InitialVariables['Uy'])+"m/s")
        setInputFieldQuantity(self.form.Uz, str(self.InitialVariables['Uz'])+"m/s")
        setInputFieldQuantity(self.form.pressure, str(self.InitialVariables['Pressure'])+"kg/m/s^2")

        # Add volume fraction fields
        if self.physicsModel.Phase != 'Single':
            self.form.volumeFractionsFrame.setVisible(True)
            mat_names = []
            for m in self.material_objs:
                mat_names.append(m.Label)
            self.form.comboFluid.clear()
            self.form.comboFluid.addItems(mat_names[:-1])
        else:
            self.form.volumeFractionsFrame.setVisible(False)
            self.form.comboFluid.clear()

        if self.physicsModel.Turbulence in ['RANS', 'LES']:
            self.form.turbulencePropertiesFrame.setVisible(True)
        else:
            self.form.turbulencePropertiesFrame.setVisible(False)

        use_inlet_turb = self.InitialVariables['UseInletTurbulenceValues']
        self.form.checkUseInletValues.setChecked(use_inlet_turb)
        setInputFieldQuantity(self.form.inputk, str(self.InitialVariables.get('k'))+"m^2/s^2")
        setInputFieldQuantity(self.form.inputOmega, str(self.InitialVariables.get('omega'))+"rad/s")

        use_inlet_temp = self.InitialVariables['UseInletTemperatureValues']
        self.form.checkUseInletValuesThermal.setChecked(use_inlet_temp)
        setInputFieldQuantity(self.form.inputTemperature, str(self.InitialVariables['Temperature'])+"K")

        # Add any inlets to the list
        for b in self.boundaries:
            if b.BoundarySettings['BoundaryType'] in ['inlet', 'farField']:
                self.form.comboInlets.addItems([b.Label])
        self.form.comboInlets.setCurrentIndex(self.form.comboInlets.findText(self.InitialVariables['Inlet']))

        self.updateUi()

    def updateUi(self):
        potential_foam = bool(self.InitialVariables['PotentialFoam'])
        use_inlet_UP = bool(self.InitialVariables['UseInletUPValues'])
        use_inlet_turb = bool(self.InitialVariables['UseInletTurbulenceValues'])
        use_inlet_temp = bool(self.InitialVariables['UseInletTemperatureValues'])
        self.form.basicPropertiesFrame.setVisible(not (potential_foam or use_inlet_UP))
        if self.physicsModel.Thermal == 'Energy':
            self.form.energyFrame.setVisible(not use_inlet_temp)
            self.form.bouyancyFrame.setVisible(False)
        elif self.physicsModel.Thermal == 'Buoyancy':
            self.form.energyFrame.setVisible(False)
            self.form.bouyancyFrame.setVisible(True)
        else:
            self.form.thermalPropertiesFrame.setVisible(False)
        self.form.frameInlets.setVisible(
            (use_inlet_UP or use_inlet_turb or use_inlet_temp) and self.form.comboInlets.count() > 1)
        self.form.kEpsilonFrame.setVisible(False)
        self.form.kOmegaSSTFrame.setVisible(False)
        self.form.SpalartAlmerasFrame.setVisible(False)
        if self.physicsModel.TurbulenceModel == 'kOmegaSST':
            self.form.kOmegaSSTFrame.setVisible(not use_inlet_turb)

    def radioUPChanged(self):
        self.InitialVariables['PotentialFoam'] = bool(self.form.radioButtonPotentialFlow.isChecked())
        self.InitialVariables['UseInletUPValues'] = bool(self.form.radioButtonUseInletValuesUP.isChecked())
        self.updateUi()

    def UxChanged(self, value):
        inputCheckAndStore(value, "m/s", self.InitialVariables, 'Ux')

    def UyChanged(self, value):
        inputCheckAndStore(value, "m/s", self.InitialVariables, 'Uy')

    def UzChanged(self, value):
        inputCheckAndStore(value, "m/s", self.InitialVariables, 'Uz')

    def PChanged(self, value):
        inputCheckAndStore(value, "kg/m/s^2", self.InitialVariables, 'Pressure')

    def getMaterialName(self, index):
        return self.material_objs[index].Label

    def comboFluidChanged(self, index):
        if 'alphas' not in self.InitialVariables:
            self.InitialVariables['alphas'] = {}
        setInputFieldQuantity(self.form.inputVolumeFraction,
                              str(self.InitialVariables['alphas'].get(self.getMaterialName(index), 0.0)))

    def inputVolumeFractionChanged(self, value):
        inputCheckAndStore(value, "m/m", self.InitialVariables['alphas'], self.form.comboFluid.currentText())

    def inputTemperatureChanged(self, value):
        inputCheckAndStore(value, "K", self.InitialVariables, 'Temperature')

    def checkUseInletTemperatureValueChanged(self, checked):
        self.InitialVariables['UseInletTemperatureValues'] = bool(checked)
        self.updateUi()

    def checkUseInletValuesChanged(self, checked):
        self.InitialVariables['UseInletTurbulenceValues'] = bool(checked)
        self.updateUi()

    def comboInletsChanged(self, index):
        self.InitialVariables['Inlet'] = self.form.comboInlets.currentText()

    def inputkChanged(self, text):
        inputCheckAndStore(text, "m^2/s^2", self.InitialVariables, 'k')

    def inputOmegaChanged(self, text):
        inputCheckAndStore(text, "rad/s", self.InitialVariables, 'omega')

    def accept(self):
        doc = FreeCADGui.getDocument(self.obj.Document)
        doc.resetEdit()

        FreeCADGui.doCommand("\n# Values are converted to SI units and stored (eg. m/s)")
        FreeCADGui.doCommand("init = FreeCAD.ActiveDocument.{}.InitialVariables".format(self.obj.Name))
        FreeCADGui.doCommand("init['PotentialFoam'] = {}".format(self.InitialVariables['PotentialFoam']))
        FreeCADGui.doCommand("init['UseInletUPValues'] = {}".format(self.InitialVariables['UseInletUPValues']))
        FreeCADGui.doCommand("init['Ux'] = {}".format(self.InitialVariables['Ux']))
        FreeCADGui.doCommand("init['Uy'] = {}".format(self.InitialVariables['Uy']))
        FreeCADGui.doCommand("init['Uz'] = {}".format(self.InitialVariables['Uz']))
        FreeCADGui.doCommand("init['Pressure'] = {}".format(self.InitialVariables['Pressure']))
        FreeCADGui.doCommand("init['alphas'] = {}")
        if len(self.material_objs) > 1:
            for i in range(len(self.material_objs)-1):
                alphaName = self.getMaterialName(i)
                FreeCADGui.doCommand("init['alphas']['{}'] = {}".format(
                    alphaName, self.InitialVariables['alphas'].get(alphaName, 0.0)))
        FreeCADGui.doCommand("init['UseInletTemperatureValues'] "
                             "= {}".format(self.InitialVariables['UseInletTemperatureValues']))
        FreeCADGui.doCommand("init['Temperature'] "
                             "= {}".format(self.InitialVariables['Temperature']))
        FreeCADGui.doCommand("init['UseInletTurbulenceValues'] "
                             "= {}".format(self.InitialVariables['UseInletTurbulenceValues']))
        FreeCADGui.doCommand("init['omega'] = {}".format(self.InitialVariables['omega']))
        FreeCADGui.doCommand("init['k'] = {}".format(self.InitialVariables['k']))
        FreeCADGui.doCommand("init['Inlet'] = '{}'".format(self.InitialVariables['Inlet']))
        FreeCADGui.doCommand("FreeCAD.ActiveDocument.{}.InitialVariables = init".format(self.obj.Name))

    def reject(self):
        doc = FreeCADGui.getDocument(self.obj.Document)
        doc.resetEdit()

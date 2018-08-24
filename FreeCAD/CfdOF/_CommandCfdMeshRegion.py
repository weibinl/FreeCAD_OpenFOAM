# ***************************************************************************
# *                                                                         *
# *   Copyright (c) 2016 - Bernd Hahnebach <bernd@bimstatik.org>            *
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
import platform
try:
    from femcommands.manager import CommandManager
except ImportError:  # Backward compatibility
    from PyGui.FemCommands import FemCommands as CommandManager
import FreeCADGui
from PySide import QtCore
import os
import CfdTools


class _CommandMeshRegion(CommandManager):
    """ The Cfd_MeshRegion command definition """
    def __init__(self):
        super(_CommandMeshRegion, self).__init__()
        icon_path = os.path.join(CfdTools.get_module_path(), "Gui", "Resources", "icons", "meshRegion.png")
        self.resources = {'Pixmap': icon_path,
                          'MenuText': QtCore.QT_TRANSLATE_NOOP("Cfd_MeshRegion", "Mesh region"),
                          'Accel': "M, R",
                          'ToolTip': QtCore.QT_TRANSLATE_NOOP("Cfd_MeshRegion", "Creates a mesh region")}
        self.is_active = 'with_femmesh'

    def Activated(self):
        FreeCAD.ActiveDocument.openTransaction("Create MeshRegion")
        FreeCADGui.addModule("CfdMeshRegion")
        sel = FreeCADGui.Selection.getSelection()
        if len(sel) == 1:
            sobj = sel[0]
            if len(sel) == 1 \
                    and hasattr(sobj, "Proxy") \
                    and sobj.Proxy.Type == "CfdMesh":
                FreeCADGui.doCommand("CfdMeshRegion.makeCfdMeshRegion(App.ActiveDocument." + sobj.Name + ")")
                FreeCADGui.ActiveDocument.setEdit(FreeCAD.ActiveDocument.ActiveObject.Name)

        FreeCADGui.Selection.clearSelection()

FreeCADGui.addCommand('Cfd_MeshRegion', _CommandMeshRegion())

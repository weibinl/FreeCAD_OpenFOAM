from freecad_import import import_freecad_path
FREECADPATH='/usr/lib/freecad/lib'
FreeCAD=import_freecad_path(FREECADPATH)
# this will import the FreeCAD module

FREECAD_DOC_NAME='python_scrit_test'
FREECAD_DOC_PATH='/home/weibin/FREECAD_DOC/'
FREECAD_DOC_EXTENSION='.fcstd' 
# extension to use when sae freecad file

WORKING_DOC=FreeCAD.newDocument(FREECAD_DOC_NAME)
CYLINDER_1=WORKING_DOC.addObject('Part::Cylinder','CYLINDER_1')
CYLINDER_1.Radius=4
WORKING_DOC.saveAs(FREECAD_DOC_PATH+FREECAD_DOC_NAME+FREECAD_DOC_EXTENSION)
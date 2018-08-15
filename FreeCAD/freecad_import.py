def import_freecad_path(FREECAD_PAHT):
    import sys
    sys.path.append(FREECAD_PAHT)
    import FreeCAD
    return FreeCAD
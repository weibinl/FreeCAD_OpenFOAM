# Computational fluid dynamics (CFD) workbench for FreeCAD

This workbench aims to help users set up and run CFD analyses. It guides the user in selecting the relevant physics, 
specifying the material properties, generating a mesh, assigning boundary conditions and setting the solver settings
before running the simulation. Where possible, best practices are included to improve the stability of the solvers.

![screenshot](https://forum.freecadweb.org/download/file.php?id=35618)

#### File format compatibility

As the workbench is still in its early development phase and evolving rapidly, new developments may sometimes
change the saved object format. This may mean that some objects do not load correctly
from previously saved files, and must be re-generated. For better 
forward-compatibilty, it is advised that you save the Python script used to generate an analysis
(right-click in Python console and choose 'Save history as'). 

## Features

### Current:

* Incompressible, laminar flow (simpleFoam, pimpleFoam).
* Basic multiphase capability (interFoam, multiphaseInterFoam)
* Basic material data base.
* Flow initialisation with a potential solver.
* Tetrahedral meshing using GMSH including multiple region meshing using FEM workbench functionality.
* Post processing using paraview.
* Porous regions and porous baffles.
* Runs on Windows 7-10
* Unit testing
* Cut-cell Cartesian meshing with boundary layers (cfMesh).
* Cut-cell Cartesian meshing with porous media (snappyHexMesh).
* Extension to turbulent flow using RANS (k-w SST).
* New case builder using a modular template structure
* Macro scripting

### Planned:

* Conjugate heat transfer.

### Platforms supported

#### Linux: 

Any system on which FreeCAD and the prerequisites listed below can be installed. The following have been tested:
* Ubuntu 16.04 
* Fedora 24-25

#### Windows:

* Windows 7 (tested)
* Windows 8 (not yet tested)
* Windows 10 (tested)

Note: 64-bit version is required

#### MacOSX:

Not tested, but a POSIX system. Possible to install and run OpenFOAM. 
      
=============================================
  
## Getting started

### Prerequisites

The CFD workbench depends on the following external software, some of
which can be automatically installed (see below for instructions).

- [Latest release version of FreeCAD (0.17)](https://github.com/FreeCAD/FreeCAD/releases/tag/0.17)
 or [latest development version (0.18 prerelease; requires git commit 13528 or later)](https://github.com/FreeCAD/FreeCAD/releases)  
- [OpenFOAM (version 5.x currently tested; most other recent versions should work.)](http://openfoam.org/download/)  
- [Paraview](http://www.paraview.org/)  
- [GMSH (version 2.13 or later)](http://gmsh.info/)  
- [cfMesh (version 1.1.2 updated to compile with OpenFOAM v5.x)](https://sourceforge.net/projects/cfmesh-cfdof/)

### Setting up CFD workbench

#### Windows

The latest FreeCAD build can be obtained from
https://www.freecadweb.org/wiki/Download and the latest
CFD workbench can be installed into it using the Addon manager:

* After running the installer or extracting the .7z archive to a directory <FreeCAD-directory>,
run FreeCAD in place (<FreeCAD-directory\bin\FreeCAD.exe). 
* Select Tools | Addon manager ...
* Select CfdOF in the list of workbenches, and click "Install/update"
* Restart FreeCAD
* For installation of dependencies, see below.

Note: The CFD workbench can be updated at any time through the Addon manager.

##### Dependency installation

Dependencies can be checked and installed
conveniently from the CFD Preferences panel in FreeCAD.
In the FreeCAD window, select Edit | Preferences ... and
choose "CFD". 

The OpenFOAM installation is via the [blueCFD-Core](http://bluecfd.github.io/Core/Downloads/) package (version 2017-2),
with which Paraview comes bundled. This can be installed
manually using the above link, or by clicking the relevant
button in the Preferences panel described above.

Set the OpenFOAM install directory in the preferences
panel to \<blueCFD install directory\>\OpenFOAM-5.x
 (It will be automatically detected in the default install
location.)

Likewise, cfMesh can be installed from the 
Preferences panel. cfMesh is automatically built from source 
inside the OpenFOAM environment if installed from the 
Preferences panel. Note that this is a lengthy process.

Choosing the "Check dependencies" option will verify that all 
prerequisites have been successfully installed.

#### Linux

The latest release or development version of FreeCAD can be obtained from 
https://github.com/FreeCAD/FreeCAD/releases (.AppImage files).
The .AppImage containers may cause library version conflicts
when running OpenFOAM from within FreeCAD. 
The [Ubuntu PPA daily build](https://www.freecadweb.org/wiki/Install_on_Unix)
packages are an alternative binary option. Otherwise, FreeCAD can be built 
from the source code at 
https://github.com/FreeCAD/FreeCAD\. 

The latest CFD workbench can be installed into FreeCAD using the Addon manager:

* Run FreeCAD 
* Select Tools | Addon manager ...
* Select CfdOF in the list of workbenches, and click "Install/update"
* Restart FreeCAD
* For installation of dependencies, see below.


##### Dependency installation

Dependencies can be checked and some of them installed
conveniently from the CFD Preferences panel in FreeCAD.
In the FreeCAD window, select Edit | Preferences ... and
choose "CFD".

However, in Linux, manual installation is required for 
[OpenFOAM](http://openfoam.org/download/),
[Paraview](http://www.paraview.org/) and
[GMSH](http://gmsh.info/), which should be
installed using your distribution's package manager
or the links above.

Set the OpenFOAM install directory in the preferences
panel - typical install locations are /home/user/OpenFOAM/OpenFOAM-5.x
or /opt/openfoam5 (It will be automatically detected in common default install
locations.)

cfMesh can be installed using the Preferences panel described above,
and will be downloaded and built from the source
code inside your OpenFOAM installation if you have
not already done so yourself. Note that this is a lengthy process.

Choosing the "Check dependencies" option will verify that all 
prerequisites have been successfully installed.


## Feedback

### Submitting Bugs

Please discuss issues on the [CfdOF dedicated FreeCAD forum](https://forum.freecadweb.org/viewforum.php?f=37).

## Development

It is asked that developers should only add functionality or code that is working and can be tested. Dead code, even
portions included for possible future functionality, reduces function clarity and increases the maintenance overhead. 
Our philosophy is 'Do the basics well' and therefore robust operation takes precedence over extended functionality.

### Testing

Unit testing is currently under development. Where possible, it is asked that all new functionality should be included
in the unit test framework.


### Style guide

For consistency please follow [PEP8](https://www.python.org/dev/peps/pep-0008/)
1. Use 4 spaces per indentation level (spaces are preferred over tabs).
2. Limit all lines to a maximum of 120 characters.
3. Break lines before binary operators.
4. Blank lines 
    
    - Surround top-level function and class definitions with two lines.

    - Definitions inside a class are surrounded by a single line.
    
5. Imports should usually be on separate lines.
6. Comments

    - Docstrings always use """triple double-quotes"""
    
    - Block comment starts with a # and a single space and are indented to the same level as that code
    
    - Use inline comments sparingly. They are on the same line as a statement and should be separated by at least two
 spaces from the statement. 

7. Avoid trailing whitespaces
8. Naming convention

    - ClassNames (Camel)
    - variable_names_without_capitals (Underscore)
    - CONSTANTS_USE_CAPITALS (Uppercase)
    - functions_without_capitals (underscore, preferred as it follows PEP8)
    - functionsWithoutCapitals (Camel instead of underscore is accepted as it is widely used within FreeCAD)
    - __class_attribute (Double leading underscore)


## Lead developers

Oliver Oxtoby (CSIR, 2016-2018) <oliveroxtoby@gmail.com>  
Johan Heyns (CSIR, 2016-2018) <jaheyns@gmail.com>  
Alfred Bogaers (CSIR, 2016-2018) <abogaers@csir.co.za>    

We would like to thank Eskom for their financial contribution as well as the following people for their contribution to the code:

Qingfeng Xia (2015) - Original framework;   
Klaus Sembritzki (2017) - Multiphase


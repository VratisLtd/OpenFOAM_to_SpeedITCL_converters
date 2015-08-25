OpenFOAM_to_SpeedITCL_converters
================================


# Introduction

SpeedIT FLOWCL™ is a general-purpose finite volume-based Computational Fluid 
Dynamics (CFD) flow solver that fully runs on GPU with OpenCL. It accelerates 
pressure-velocity solution procedure for Navier-Stokes equation, namely SIMPLE 
and PISO. SpeedIT FLOW is easy to use by CFD engineers with all levels of 
expertise. It can be invoked easily as a standalone application. 

You can get more information on our [website](www.vratis.com)

This software is an OpenFOAM to SpeedIT FlowCL converter. It allows for a 
conversion of ready OpenFOAM test cases to SpeedIT Flow format. 

OpenFOAM to SpeedIT Flow converters can be downloaded from the 
[GitHub](https://github.com/VratisLtd/OpenFOAM\_to\_SpeedITCL_converters)

# Requirements

* GCC 4.7.3 
* CMAKE ver. 2.8.11 and newer,
* [Boost 1.55](http://www.boost.org/),
* [OpenFOAM ver. 2.0.1](http://www.openfoam.org/archive/2.0.1/download/)

# Installation

Set up the OpenFOAM environment by sourcing one of the configuration files. 
Usualy $WM\_PROJECT\_DIR/etc/bashrc. For more information please check 
[OpenFOAM installation website](http://www.openfoam.org/download/source.php), 
section Setting Environment Variables.

The configuration procedure is based on CMAKE system.
Create a new directory for the converter outside the OpenFOAM\_to\_SpeedIT\_converters 
source directory. Navigate to this directory and execute _cmake_ followed by a 
path to the downloaded converters code directory

$ cmake /path\_to\_downloaded\_converters\_code\_directory/

As a result a a set of directories and files will be created (including 
_Makefile_). Now it is enough to execute _make_ and the OpenFoam2SpeedIT file 
will be created.

If you want to make changes to the project you can open it with QtCreator. 
Select File -> Open file or project, navigate to the created directory and 
choose cmake_install.cmake file.

# Usage

You can execute the converter either by navigating to the test case directory 
and executing the _OpenFoam2SpeedIT_ there or from outside the case directory 
with the --case option:

$ OpenFoam2SpeedIT --case /path\_to\_test\_case/

As a result a _sitfcl_ directory will be created. In this directory all the data 
needed for the launch of the SpeedIT Flow solvers is kept.

Converter have additional options:

  **--help**                produce help message

  **--case arg** (=./)      path to case

  **--out arg** (=./sitfcl) output path

 **--binary arg** (=0)     use binary mode

  **--precision arg** (=16) number of digits in non-binary mode


# Contact

For technical questions please contact us at <support@vratis.com>.

For general inquires please contact us at <info@vratis.com>.

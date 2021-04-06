# s3m-dev
---------

This is the **source code of S3M**, a spatially explicit and hydrology-oriented cryospheric model that successfully reconstructs seasonal snow and glacier evolution through time and that can be natively coupled with distributed hydrologic models. 

Model physics include precipitation-phase partitioning, snow and glacier energy and mass balances, snow rheology and hydraulics, and a data-assimilation protocol. 

S3M includes an explicit representation of the spatial patterns of snow liquid-water content, an hybrid approach to snowmelt that decouples the radiation- and temperature-driven contributions, the implementation of the Delta-h parametrization for distributed ice-thickness change, and the inclusion of a distributed debris-driven melt factor. 

Components
**********
All codes and datasets are freely available, and users can obtain them from our GitHub repository [1_].

Prerequisites
*************

In order to use S3M, users are strongly recommended to check whether the following characteristics, libraries and packages are available and correctly installed on their machine.

Usually, S3M is installed on **Linux Debian/Ubuntu 64bit** environment and all libraries, packages and applications must be compilled and/or installed accordingly.

All codes, subroutines and scripts are developed using **Fortran** (version 2003 and greater) programming language [2_].

Our GitHub repository with details on how to install all needed libraries and environments is called **fp-envs** [3_] and the users can find it in our repository hosted by GitHub [1_].

Fortran libraries
-----------------

S3M needs the netcdf4 library to read input provided by other preprocessing tools, as well as to write output for external applications (such as Panoply, cdo, ncview ...).
In order to set and compile netcdf4 library, a bash script named **setup_fp_env_system.sh** is provided [3_]. 
This script will automatically download **zlib** [4_], **hdf5** [5_] and **netcdf4** [6_] libraries from online repositories; after downloading source-compressed archives, the script automatically creates a folder in “$HOME/user/fp_libs_system/” where all libraries will be compiled and installed. During installation, an environment file named “fp_env_system” is created for saving LD_LIBRARY_PATH (for native code libraries) and PATH (for executables) references of installed libraries.





References
**********
| [1_] CIMA Hydrology and Hydraulics GitHub Repository
| [2_] Fortran programming language
| [3_] FloodPROOFS virtual environment tools
| [4_] ZLIB compression library
| [5_] HDF5 data software library 
| [6_] NetCDF4 data software library 
| [7_] Hydrological Model Continuum codes
| [8_] NetBeans Apache IDE 
| [9_] GDB 

.. _1: https://github.com/c-hydro
.. _2: https://en.wikipedia.org/wiki/Fortran
.. _3: https://github.com/c-hydro/fp-env
.. _4: https://zlib.net/
.. _5: https://www.hdfgroup.org/solutions/hdf5/
.. _6: https://www.unidata.ucar.edu/
.. _7: https://github.com/c-hydro/hmc-dev
.. _8: https://netbeans.apache.org/
.. _9: https://www.gnu.org/software/gdb/
.. _license: LICENSE.rst
.. _changelog: CHANGELOG.rst
.. _authors: AUTHORS.rst

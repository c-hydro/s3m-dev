S3M, CIMA's cryospheric model
============================

Welcome to the GitHub repository of **S3M**, a spatially explicit and hydrology-oriented cryospheric model that successfully reconstructs seasonal snow and glacier evolution through time and that can be natively coupled with distributed hydrologic models. 

Model physics include precipitation-phase partitioning, snow and glacier energy and mass balances, snow rheology and hydraulics, and a data-assimilation protocol. 

S3M includes an explicit representation of the spatial patterns of snow liquid-water content, an hybrid approach to snowmelt that decouples the radiation- and temperature-driven contributions, the implementation of the Delta-h parametrization for distributed ice-thickness change, and the inclusion of a distributed debris-driven melt factor. 

Components
-----------------
All codes and datasets are freely available, and users can obtain them from our GitHub repository [1_].

Prerequisites
-----------------

In order to use S3M, users are strongly recommended to check whether the following characteristics, libraries and packages are available and correctly installed on their machine.

Usually, S3M is installed on **Linux Debian/Ubuntu 64bit** environment and all libraries, packages and applications must be compilled and/or installed accordingly.

All codes, subroutines and scripts are developed using **Fortran** (version 2003 and greater) programming language [2_].

Our GitHub repository with details on how to install all needed libraries and environments is called **fp-envs** [3_] and the users can find it in our repository hosted by GitHub [1_].

Fortran libraries
-----------------

S3M needs the netcdf4 library to read input provided by other preprocessing tools, as well as to write output for external applications (such as Panoply, cdo, ncview ...).
In order to set and compile netcdf4 library, a bash script named **setup_fp_env_system.sh** is provided [3_]. 
This script will automatically download **zlib** [4_], **hdf5** [5_] and **netcdf4** [6_] libraries from online repositories; after downloading source-compressed archives, the script automatically creates a folder in “$HOME/user/fp_libs_system/” where all libraries will be compiled and installed. During installation, an environment file named “fp_env_system” is created for saving LD_LIBRARY_PATH (for native code libraries) and PATH (for executables) references of installed libraries.

S3M libraries
-------------
After preparing all necessary libraries and environmental settings, the source files of S3M must be compiled to run properly [7_]. Usually, source files are compiled using **GNU compilers** (such as gcc and gfortran), which have to be pre-installed on user’s machine. To compile and optimize S3M codes, a bash file called **configure.sh** is provided [7_]. Using this setup file, user will have to answer some questions to guide the machine through compiling S3M.
Usually, options are set as follows:

    • set compiler type [1] for using GNU/GFortran compiler;
    • set optimization option [2] for using production mode; 
    • set profiler option [2] for skipping profiling, which is used to control model performances;
    • set NetCDF4 library [1] for using NetCDF4 input and output files format. 
       
The user will then be given an opportunity to set the path of the NetCDF4 library, for example in the form “$HOME/user/fp_libs_system/nc4”. By pressing the ENTER key instead of specifying the path, the script will try to locate this library automatically. 

Finally, the user will be allowed to choose a name for the .x compiled model. Here again, by pressing the ENTER key instead of specifying the name, the script will assign this name automatically. 

Netbeans configuration
----------------------
In order to further develop S3M, users can set their own devoloping environment using a IDE, such as NetBeans Apache [8_]. Such environments are particularly helpful with code development and debugging, also allowing one to analyze codes and memory usage. Here, the configuration of a debug workspace in Apache NetBeans IDE will be presented.

First, packages for C, C++, and Fortran programming languages must be installed in Apache NetBeans; to complete this step, users have to install the package related to C/C++ language in NetBeans following these instructions:

  1) Tools --> Plugins --> Settings --> Tick "NetBeans 8.2 Plugin Portal" --> Close 

and then reboot the Apache NetBeans IDE.
Second, users should create a New Project following this path: 

  2) File --> New Project --> Categories :: Samples :: C/C++ --> Projects :: Fortran Hello World Application --> Next --> Choose Name --> Close

After creating this folder project, users should first remove all default source code under the ‘’Source Files’’ directory on the left menu, where the name of the project is visible, and second move all required source files of S3M into the system physical folder associated to this project. The first step must be performed from Netbeans, while the second step must be performed externally from NetBeans using standard copy-paste. Third, source-code files must now be imported in the NetBeans project using the left menu where the name of the project is visible, right-clicking on the project, and choosing:

  3) Add existing item ...

Performing this action, a form to select files will be opened. After selecting all required files from the project folder, all source files will now be available into the NetBeans project under the ‘’Source Files’’ directory. Note that copy-pasting source files into the project folder ensures that NetBeans will correctly link the project to the source code located in the same physical folder of the NetBeans project. Despite this step not being compulsory, NetBeans does not physically create files into the project folder if one added items directly from another folder within Netbeans, meaning the user may end up inadvertently modifying source code from other folders, including other branches or versions of S3M. 

Next steps cover the configuration of dependencies in the project. Particularly, the main task is to link the NetCDF library to the project.
For configuring the NetCDF4 in Apache NetBeans IDE, users must first right-click on the project name on the left menu where the name of the project is visible and then follow: 

  4) Properties --> Linker --> Libraries,
  
then click on the […] icon, select ‘’Add Library…’’ and find the following files in /path_to_netcdf/

  netcdff.a and netcdff.so 
  
Please select files with the "double f" for fortran libraries! Press OK.

Second, right-click on the project name on the left menu where the name of the project is visible and then follow: 

  5) Properties --> Linker --> Additional Options
     
Here, the user should include the following instructions, where ‘’ path_to_netcdf’’ must be replaced with the actual path to the NetCDF library in user’s machine:

  -I/path_to_netcdf/include/ 
  -L/path_to_netcdf/lib/ 
  -lnetcdff -lnetcdff   
      
Please note that the first line starts with a capital i, for ‘’Include’’!  Press OK.

Third, right-click on the project name on the left menu where the name of the project is visible and then follow: 

  6) Properties --> Fortran Compiler --> Additional Options

      -I/path_to_netcdf/include/ 
      -L/path_to_netcdf/lib/ 
      -lnetcdff -lnetcdff  
      
Here again, ‘’ path_to_netcdf’’ must be replaced with the actual path to the NetCDF library in user’s machine. Press OK. 

Fourth, right-click on the project name on the left menu where the name of the project is visible and then follow: 

  7) Properties --> Fortran Compiler --> Additional Options
  
      gfortran: -cpp -DLIB_NC
      ifort: -fpp -DLIB_NC  

Press OK. 

Fifth, right-click on the project name on the left menu where the name of the project is visible and then follow:

  8) Properties --> Run --> Environment --> NewValue
  
      Name: LD_LIBRARY_PATH 
      Value: $LD_LIBRARY_PATH:/path_to_necdf/lib/
      
Here again, ‘’ path_to_netcdf’’ must be replaced with the actual path to the NetCDF library in user’s machine. Press OK.

Once the NetCDF4 are linked, it will be possible to compile each source file in NetBeans using the F9 key. Note that the bash file called **configure.sh** [7_] specifies the correct order for compiling source files. 

After performing all these steps, the users have to set the debug command to run S3M using, for instance, a namelist file of a study case.  In order to do so, right-click on the project name on the left menu where the name of the project is visible and then follow::  
  
  9) Properties --> Debug --> Debug Command 
  	
  "${OUTPUT_PATH}" domainname.info.txt

Where domainname.info.txt is the namelist file with all model settings (see an example on this GitHub repository). 

After setting the environment and all needed options for running the model, the users will be able to achieve a deeper understanding of S3M by, e.g., using breakpoints and all the features available in the **gdb** debugging library [9_].

Contribute and Guidelines
-----------------

We are happy if you want to contribute. Please raise an issue explaining what is missing or if you find a bug. We will also gladly accept pull requests in our master branch for new features or bug fixes.

Authors
-----------------

All authors involved in the library development for S3M are reported in this authors_ file.

License
-----------------

By accessing or using S3M, code, data or documentation, you agree to be bound by the S3M license available. See the license_ for details. 

Changelog
-----------------

All notable changes and bug fixing fir this project will be documented in this changelog_ file.
    
    
References
-----------------
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
.. _3: https://github.com/c-hydro/fp-envs
.. _4: https://zlib.net/
.. _5: https://www.hdfgroup.org/solutions/hdf5/
.. _6: https://www.unidata.ucar.edu/
.. _7: https://github.com/c-hydro/s3m-dev
.. _8: https://netbeans.apache.org/
.. _9: https://www.gnu.org/software/gdb/
.. _license: LICENSE.rst
.. _changelog: CHANGELOG.rst
.. _authors: AUTHORS.rst

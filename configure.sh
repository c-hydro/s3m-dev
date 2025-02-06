#!/bin/bash -e

#-----------------------------------------------------------------------------------------
# Script option(s)
Script="S3M Library Builder"
Version="1.2.0"
Date='2025/02/06'

# Other option(s)
Archive_Default="s3m-5.3.4.tar.gz"
# Default folder
Lib_Dir_Deps_Default="$HOME/fp_libs_s3m"
# Executables folder
Lib_Dir_Exec_Default="$Lib_Dir_Deps_Default/s3m"
# Compilation Mode
Lib_Building_Default=false
# Executable name
Exec_Default='S3M_Model_V5_3_4.x'
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Start - Script
echo "----------------------------------------------------------------"
echo "$Script - Version $Version "
echo "Script to set, compile and build S3M model"
echo "----------------------------------------------------------------"
echo ""
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Parse script argument(s)
echo "----------------------------------------------------------------"
echo "Parse argument(s) ... "

# Get arguments number and values
Args_N=$#
Args_Values=$@

echo ""
echo " => Script arguments number: $Args_N"
echo " => Script arguments values: $Args_Values"
echo ""
echo " => Script arguments 1 - Archive Name [string: filename]-> $1"
echo " => Script arguments 2 - Directory of dependencies [string: path] -> $2"
echo " => Script arguments 3 - Directory of S3M executable [string: path] -> $3"
echo " => Script arguments 4 - Compilation Mode [boolean: {true, false}] -> $4"
echo ""

# Set argument(s)
if [ $# -eq 0 ]; then
	Archive=$Archive_Default
	Lib_Dir_Deps=$Lib_Dir_Deps_Default
	Lib_Dir_Exec=$Lib_Dir_Exec_Default
	Lib_Building_Automatic=$Lib_Building_Default
	echo " => Script arguments - SET [None] DEFAULT [1,2,3,4]"
elif [ $# -eq 1 ]; then
	Archive=$1
	Lib_Dir_Deps=$Lib_Dir_Deps_Default
	Lib_Dir_Exec=$Lib_Dir_Exec_Default
	Lib_Building_Automatic=$Lib_Building_Default
	echo " => Script arguments - SET [1] DEFAULT [2,3,4]"
elif [ $# -eq 2 ]; then
	Archive=$1
	Lib_Dir_Deps=$2
	Lib_Dir_Exec="$Lib_Dir_Deps/s3m"
	Lib_Building_Automatic=$Lib_Building_Default
	echo " => Script arguments - SET [1,2] DEFAULT [3,4]"
elif [ $# -eq 3 ]; then
	Archive=$1
	Lib_Dir_Deps=$2
	Lib_Dir_Exec=$3
	Lib_Building_Automatic=$Lib_Building_Default
	echo " => Script arguments - SET [1,2,3] DEFAULT [4]"
elif [ $# -eq 4 ]; then
	Archive=$1
	Lib_Dir_Deps=$2
	Lib_Dir_Exec=$3
	Lib_Building_Automatic=$4
	echo " => Script arguments - SET [1,2,3,4] DEFAULT [None]"
else
	echo " => Script arguments - FATAL ERROR"
	exit 1
fi

# Message about compilation mode
if $Lib_Building_Automatic ; then 
	echo " => Script compilation mode: AUTOMATIC"
else
	echo " => Script compilation mode: MANUAL"
fi

echo ""
echo "Parse Argument(s)  ... OK!"
echo "----------------------------------------------------------------"
echo ""
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Automatically checking for netcdf library
echo "----------------------------------------------------------------"
echo "Step 0 - Configure ==> Detection of NetCDF4 path ... "
echo ""

if (find $Lib_Dir_Deps -type d  -path '*nc4'); then
    NC_Dir_Default=$(find $Lib_Dir_Deps -type d  -path '*nc4')
    echo "NetCDF4 complete library path set using automatic detection [$NC_Dir_Default]"
else
    NC_Dir_Default="/$HOME/fp_libs_system/nc4/"
    echo "NetCDF4 complete library path set using a DEFAULT path [$NC_Dir_Default]"
fi

echo ""
echo "Step 0 - Configure ==> Detection of NetCDF4 ... OK!"
echo "----------------------------------------------------------------"
echo ""
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Menu to set compiler type
echo "----------------------------------------------------------------"
echo "Step 1 - Configure ==> Set compiler type ... "
echo ""
if $Lib_Building_Automatic ; then  
	echo " ==> GNU/GFortran Compiler selected in automatic mode"
	Comp_Name="GNU/GFortran"
	Comp_Exec="gfortran"
	Comp_Obj="-c -g -O2 -cpp "
else
	PS3=' ==> Please enter your choice: '
	Comp_Opts=( " GNU/GFortran  INTEL/Fortran  Quit" )
	select Opt in $Comp_Opts; do
		if [ "$Opt" = "Quit" ]; then
			echo Quit
			exit
		elif [ "$Opt" = "GNU/GFortran" ]; then
			Comp_Name="GNU/GFortran"
			Comp_Exec="gfortran"
			Comp_Obj="-c -g -O2 -cpp "
			break
		elif [ "$Opt" = "INTEL/Fortran" ]; then
			Comp_Name="INTEL/Fortran"
			Comp_Exec="ifort"
			Comp_Obj="-c -g -O2 -fpp "
			break
		else
			echo 'Bad Option!'
		fi
	done
fi

echo " ==> Compiler Name: " $Comp_Name "; Compiler Exec: " $Comp_Exec "; Comp Obj: " $Comp_Obj
echo ""
echo "Step 1 - Configure ==> Set compiler type ... OK!"
echo "----------------------------------------------------------------"
echo ""
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Menu to set optimization option
echo "----------------------------------------------------------------"
echo "Step 2 - Configure ==> Set optimization option ... "
echo ""

if $Lib_Building_Automatic ; then  
	echo " ==> Production optimization for GNU/GFortran compiler selected in automatic mode"
	Optim_Opt="Production"
	Optim_Exec='-O3 -march=native -Ofast -funroll-loops -fimplicit-none  -Wall  -Wline-truncation  -fwhole-file  -std=f2008 -fall-intrinsics '
else
	PS3=' ==> Please enter your choice: '
	Optim_Opts=( " Debug  Production  Quit" )
	select Opt in $Optim_Opts; do
		if [ "$Opt" = "Quit" ]; then

			echo Quit
			exit

		elif [ "$Opt" = "Debug" ]; then
			Optim_Opt="Debug"
			
			if [ "$Comp_Name" = "GNU/GFortran" ]; then
				Optim_Exec='-O2 -g3 -ggdb -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -std=f2008  -pedantic  -fbacktrace  -fall-intrinsics '
			elif [ "$Comp_Name" = "INTEL/Fortran" ]; then
				Optim_Exec='-O2 -static -static-intel -assume byterecl -align dcommons -fast '
			fi
			break

		elif [ "$Opt" = "Production" ]; then
			Optim_Opt="Production"

			if [ "$Comp_Name" = "GNU/GFortran" ]; then
				Optim_Exec='-O3 -march=native -Ofast -funroll-loops -fimplicit-none  -Wall  -Wline-truncation  -fwhole-file  -std=f2008 -fall-intrinsics '
			elif [ "$Comp_Name" = "INTEL/Fortran" ]; then
				Optim_Exec='-O2 -static -static-intel -assume byterecl -align dcommons -fast '
			fi
			break

		else
			echo 'Bad Option!'
		fi
	done
fi

echo " ==> Optimization Option: " $Optim_Opt "; Optimization Exec: " $Optim_Exec
echo ""
echo "Step 2 - Configure ==> Set optimization option ... OK!"
echo "----------------------------------------------------------------"
echo ""
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Menu to set profiler option
echo "----------------------------------------------------------------"
echo "Step 3 - Configure ==> Set profiler option ... "
echo ""

if $Lib_Building_Automatic ; then 
	echo " ==> Profiler option selected in automatic mode"
	Prof_Opt=""
else
	PS3=' ==> Please enter your choice: '
	Prof_Opts=( " Yes  No  Quit" )
	select Opt in $Prof_Opts; do
		if [ "$Opt" = "Quit" ]; then
			echo Quit
			exit
		elif [ "$Opt" = "Yes" ]; then
			Prof_Opt="-pg"
			break
		elif [ "$Opt" = "No" ]; then
			Prof_Opt=""
			break
		else
			echo 'Bad Option!'
		fi
	done
fi

echo " ==> Profiler Option: " $Prof_Opt
echo "" 
echo "Step 3 - Configure ==> Set profiler option ... OK!"
echo "----------------------------------------------------------------"
echo ""
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Menu to set netCDF4 library
echo "----------------------------------------------------------------"
echo "Step 4 - Configure ==> Set NetCDF4 library ... "
echo ""

if $Lib_Building_Automatic ; then 
	echo " ==> NetCDF4 library for GNU/GFortran compiler selected in automatic mode"

	NC_Dir=$NC_Dir_Default

	NC_Inc=$NC_Dir/include \
	NC_Lib=$NC_Dir/lib \
	NC_Libs="-lnetcdff -lnetcdf"
	NC_Opt="-DLIB_NC"

	Comp_Obj=$Comp_Obj$NC_Opt 
else
	PS3=' ==> Please enter your choice: '
	NC_Opts=( " Yes  No  Quit" )
	select Opt in $NC_Opts; do
		if [ "$Opt" = "Quit" ]; then
			echo Quit
			exit
		elif [ "$Opt" = "Yes" ]; then
			read -r -p ' ==> Please enter NetCDF4 complete library path: ' NC_Dir
			
			if [ -z "$NC_Dir" ]; then
				NC_Dir=$NC_Dir_Default
	  			echo " ==> NetCDF4 complete library path set using DEFAULT path!"
			else
				echo " ==> NetCDF4 complete library path set by USER"
			fi

			NC_Inc=$NC_Dir/include \
			NC_Lib=$NC_Dir/lib \
			NC_Libs="-lnetcdff -lnetcdf"
			NC_Opt="-DLIB_NC"

			Comp_Obj=$Comp_Obj$NC_Opt 	

			break
		elif [ "$Opt" = "No" ]; then
			NC_Dir=''
			NC_Inc=""
			NC_Lib=""
			NC_Libs=""
			NC_Opt=""

			Comp_Obj=$Comp_Obj$NC_Opt 

			break
		else
			echo 'Bad Option!'
		fi
	done
fi

echo " ==> NetCDF4 path: " $NC_Dir "; NetCDF4 Option: " $NC_Opt "; NetCDF4 Comp Obj: " $Comp_Obj
echo ""
echo "Step 4 - Configure ==> Set NetCDF4 library ... OK!"
echo "----------------------------------------------------------------"
echo ""
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Set executable S3M name
echo "----------------------------------------------------------------"
echo "Step 5 - Configure ==> Set S3M executable name ... "
echo ""

if $Lib_Building_Automatic ; then 
	echo " ==> S3M executable name selected in automatic mode"
	Exec=$Exec_Default
else
	read -r -p ' ==> Please enter S3M executable name: ' Exec
	if [ -z "$Exec" ]; then
		Exec=$Exec_Default
	  	echo " ==> S3M executable name [Default]: " $Exec 
	else
		echo " ==> S3M executable name [User]: " $Exec 
	fi
fi

echo ""
echo "Step 5 - Configure ==> Set S3M executable name ... OK"
echo "----------------------------------------------------------------"
echo ""
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Set stack size (kbytes, -s) unlimited
echo "----------------------------------------------------------------"
echo "Step 1 - Compile ==> Set stack size to unlimited ... "
echo ""

ulimit -s unlimited
ulimit -a

echo ""
echo "Step 1 - Compile ==> Set stack size to unlimited ... OK!"
echo "----------------------------------------------------------------"
echo ""
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Expand archive file
echo "----------------------------------------------------------------"
echo "Step 2 - Compile ==> Expand archive file ... "
echo ""

# Get current folder
Current_Dir=${PWD}

if [ -e $Archive ]; then

    echo " ==> tar xvfz $Archive -C $Archive_Dir"
    tar xvfz $Archive -C $Archive_Dir --strip-components=1
    
    echo " ==> Files in compressed format -- Unzipped!"
    Archive_Dir=$Current_Dir/temp
    
    if [ -d "$Archive_Dir" ]; then
        rm -R $Archive_Dir
    fi
    
    mkdir $Archive_Dir
    
else
    echo " ==> Files in uncompressed format -- Skipped!"
    Archive_Dir=$Current_Dir
fi

echo "" 
echo "Step 2 - Compile ==> Expand archive file ... DONE!"
echo "----------------------------------------------------------------"
echo ""
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Build module(s) and subroutine(s)
echo "----------------------------------------------------------------"
echo "Step 3 - Compile ==> Build module(s) and subroutine(s) ... "
echo ""

cd $Archive_Dir

$Comp_Exec $Comp_Obj gnufor2.f90 $Prof_Opt
$Comp_Exec $Comp_Obj S3M_Module_Tools_Debug.f90 $Prof_Opt
$Comp_Exec $Comp_Obj S3M_Module_Args.f90 $Prof_Opt
$Comp_Exec $Comp_Obj S3M_Module_Namelist.f90 $Prof_Opt
$Comp_Exec $Comp_Obj S3M_Module_Tools_Interp.f90 $Prof_Opt
$Comp_Exec $Comp_Obj S3M_Module_Tools_Generic.f90 $Prof_Opt
$Comp_Exec $Comp_Obj S3M_Module_Tools_IO.f90 $Prof_Opt -I${NC_Inc} -L${NC_Lib} ${NC_Libs} 
$Comp_Exec $Comp_Obj S3M_Module_Tools_Time.f90 $Prof_Opt
$Comp_Exec $Comp_Obj S3M_Module_Vars_Loader.f90 $Prof_Opt
$Comp_Exec $Comp_Obj S3M_Module_Vars_Manager.f90 $Prof_Opt
$Comp_Exec $Comp_Obj S3M_Module_Phys_Snow_Apps_PhasePart.f90 $Prof_Opt
$Comp_Exec $Comp_Obj S3M_Module_Phys_Snow_Apps_Glaciers.f90 $Prof_Opt
$Comp_Exec $Comp_Obj S3M_Module_Phys_Snow_Apps_Melting.f90 $Prof_Opt
$Comp_Exec $Comp_Obj S3M_Module_Phys_Snow_Apps_Hydraulics.f90 $Prof_Opt
$Comp_Exec $Comp_Obj S3M_Module_Phys_Snow_Apps_Density.f90 $Prof_Opt
$Comp_Exec $Comp_Obj S3M_Module_Phys_Snow_Apps_Assimilation.f90 $Prof_Opt
$Comp_Exec $Comp_Obj S3M_Module_Phys_Snow.f90 $Prof_Opt
$Comp_Exec $Comp_Obj S3M_Module_Phys.f90 $Prof_Opt
$Comp_Exec $Comp_Obj S3M_Module_Data_AssSWE_Gridded.f90 $Prof_Opt -I${NC_Inc} -L${NC_Lib} ${NC_Libs}
$Comp_Exec $Comp_Obj S3M_Module_Data_Forcing_Gridded.f90 $Prof_Opt -I${NC_Inc} -L${NC_Lib} ${NC_Libs}
$Comp_Exec $Comp_Obj S3M_Module_Data_Updating_Gridded.f90 $Prof_Opt -I${NC_Inc} -L${NC_Lib} ${NC_Libs}
$Comp_Exec $Comp_Obj S3M_Module_Data_Output_Gridded.f90 $Prof_Opt -I${NC_Inc} -L${NC_Lib} ${NC_Libs}
$Comp_Exec $Comp_Obj S3M_Module_Data_Restart_Gridded.f90 $Prof_Opt -I${NC_Inc} -L${NC_Lib} ${NC_Libs}
$Comp_Exec $Comp_Obj S3M_Module_Data_Static_Gridded.f90 $Prof_Opt -I${NC_Inc} -L${NC_Lib} ${NC_Libs}
$Comp_Exec $Comp_Obj S3M_Module_Info_Gridded.f90 $Prof_Opt -I${NC_Inc} -L${NC_Lib} ${NC_Libs}
$Comp_Exec $Comp_Obj S3M_Module_Info_Time.f90 $Prof_Opt -I${NC_Inc} -L${NC_Lib} ${NC_Libs}

echo ""
echo "Step 3 - Compile ==> Build module(s) and subroutine(s) ... DONE!"
echo "----------------------------------------------------------------"
echo ""
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Link object(s) and create executable(s)
echo "----------------------------------------------------------------"
echo "Step 4 - Compile ==> Link object(s) and create S3M model executable ... "
echo ""

$Comp_Exec $Optim_Exec gnufor2.o S3M_Module_Tools_Debug.o S3M_Module_Args.o S3M_Module_Namelist.o S3M_Module_Tools_Interp.o S3M_Module_Tools_Generic.o S3M_Module_Tools_IO.o S3M_Module_Tools_Time.o S3M_Module_Vars_Loader.o S3M_Module_Vars_Manager.o S3M_Module_Phys_Snow_Apps_PhasePart.o S3M_Module_Phys_Snow_Apps_Glaciers.o S3M_Module_Phys_Snow_Apps_Melting.o S3M_Module_Phys_Snow_Apps_Hydraulics.o S3M_Module_Phys_Snow_Apps_Density.o S3M_Module_Phys_Snow_Apps_Assimilation.o S3M_Module_Phys_Snow.o S3M_Module_Phys.o S3M_Module_Data_AssSWE_Gridded.o S3M_Module_Data_Forcing_Gridded.o S3M_Module_Data_Updating_Gridded.o S3M_Module_Data_Output_Gridded.o S3M_Module_Data_Restart_Gridded.o S3M_Module_Data_Static_Gridded.o S3M_Module_Info_Gridded.o S3M_Module_Info_Time.o S3M_Main.f90 -o $Exec $Prof_Opt -I${NC_Inc} -L${NC_Lib} ${NC_Libs}

echo ""
echo "Step 4 - Compile ==> Link object(s) and create S3M model executable ... DONE!"
echo "----------------------------------------------------------------"
echo ""
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Change option(s) if S3M executable
echo "----------------------------------------------------------------"
echo "Step 5 - Compile ==> Change option(s) of S3M model executable ... "
echo ""

if [ -d "$Lib_Dir_Exec" ]; then
    rm -R $Lib_Dir_Exec
fi
mkdir -p $Lib_Dir_Exec

if [ -e $Exec ]; then
    chmod +x $Exec
    cp -r $Archive_Dir/$Exec $Lib_Dir_Exec/$Exec
    echo " ==> $Exec copied in library folder ... DONE!"
else
    echo " ==> $Exec copied in library folder ... FAILED!"
fi

echo ""
echo "Step 5 - Compile ==> Change option(s) of S3M model executable ... DONE!"
echo "----------------------------------------------------------------"
echo ""
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# End - Script
echo "----------------------------------------------------------------"
echo "$Script - Version $Version "
echo "Script to set, compile and build S3M model"
echo "COMPLETED - Bye, Bye"
echo "----------------------------------------------------------------"
echo ""
#-----------------------------------------------------------------------------------------








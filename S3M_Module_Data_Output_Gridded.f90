!------------------------------------------------------------------------------------------    
! File:   S3M_Module_Data_Output_Gridded.f90
! Author(s): Fabio Delogu, Francesco Silvestro, Simone Gabellani, Francesco Avanzi
!
! Created on May 6, 2015, 4:36 PM
! Last update on November 16, 2020 04:45 PM
!
! Module to save output results (only NC is supported!!!).
!------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------
! Module Header
module S3M_Module_Data_Output_Gridded
    
    !------------------------------------------------------------------------------------------
    ! External module(s) for all subroutine in this module
#ifdef LIB_NC
    use netcdf
#endif
    
    use S3M_Module_Namelist,        only:   oS3M_Namelist
    use S3M_Module_Vars_Loader,     only:   oS3M_Vars
    
#ifdef LIB_NC
    use S3M_Module_Tools_IO,        only:   S3M_Tools_IO_Put2d_Binary_INT, &
                                            S3M_Tools_IO_Put2d_Binary_DBL, &
                                            S3M_Tools_IO_Put3d_Binary_INT, &
                                            S3M_Tools_IO_Put3d_Binary_DBL, &
                                            S3M_Tools_IO_Put2d_NC, &
                                            S3M_Tools_IO_Put3d_NC, &
                                            S3M_Tools_IO_PutTime_DBL_NC, &
                                            S3M_Tools_IO_PutTime_STR_NC, &
                                            check
#else
    use S3M_Module_Tools_IO,        only:   S3M_Tools_IO_Put2d_Binary_INT, &
                                            S3M_Tools_IO_Put2d_Binary_DBL, &
                                            S3M_Tools_IO_Put3d_Binary_INT, &
                                            S3M_Tools_IO_Put3d_Binary_DBL
                                            
#endif
                                                                         
    use S3M_Module_Tools_Debug
    use S3M_Module_Tools_Generic,   only:   S3M_Tools_Generic_ReplaceText, & 
                                            S3M_Tools_Generic_CreateFolder, &
                                            S3M_Tools_Generic_ZipFile, &
                                            S3M_Tools_Generic_RemoveFile, &
                                            transpose3Dvar
    
    use S3M_Module_Tools_Time,      only:   S3M_Tools_Time_DateDiff
    
    use gnufor2 
    
    ! Implicit none for all subroutines in this module
    implicit none
    !------------------------------------------------------------------------------------------
    
contains
    
    !------------------------------------------------------------------------------------------
    ! Subroutine to manage output gridded data
    subroutine S3M_Data_Output_Gridded_Cpl( iID, sTime, &
                                            iRowsStart, iRowsEnd, iColsStart, iColsEnd, &
                                            iTime, iDaySteps)
        
        !------------------------------------------------------------------------------------------
        ! Variable(s)
        integer(kind = 4)           :: iID
        integer(kind = 4)           :: iRows, iCols
        integer(kind = 4)           :: iRowsStart, iRowsEnd, iColsStart, iColsEnd
        integer(kind = 4)           :: iTime, iDaySteps
        
        integer(kind = 4)           :: iFlagIceMassBalance, iFlagOutputMode
        
        character(len = 19)         :: sTime
        character(len = 700)        :: sPathData_Output
        character(len = 700)        :: sCommandCreateFolder
        character(len = 2)          :: sWYstart
        
        real(kind = 4), dimension(iRowsEnd - iRowsStart + 1, iColsEnd - iColsStart + 1)     ::  a2dVarDEM
        
        real(kind = 4), dimension(iRowsEnd - iRowsStart + 1, iColsEnd - iColsStart + 1)     ::  a2dVarLat, a2dVarLon
        integer(kind = 4), dimension(iRowsEnd - iRowsStart + 1, iColsEnd - iColsStart + 1)  ::  a2iVarAgeS
        real(kind = 4), dimension(iRowsEnd - iRowsStart + 1, iColsEnd - iColsStart + 1)     ::  a2dVarRainFall, a2dVarSnowFall, &
                                                                                                a2dVarPrecip, &
                                                                                                a2dVarSnowFallDayCum, &
                                                                                                a2dVarSWE_D, &
                                                                                                a2dVarRho_D, &
                                                                                                a2dVarSWE_W, a2dVarTheta_W, &
                                                                                                a2dVarSWE, a2dVarRhoS, a2dVarH_S, &
                                                                                                a2dVarRhoS0, & 
                                                                                                a2dVarSnowMask, a2dVarMeltingS, &
                                                                                                a2dVarMeltingSDayCum, &
                                                                                                a2dVarRefreezingS, a2dVarAlbedoS, &
                                                                                                a2dVarMeltingG, &
                                                                                                a2dVarMeltingGCumWY, &
                                                                                                a2dVarChangeThickness, &
                                                                                                a2dVarOutflow, & 
                                                                                                a2dVarReff, &
                                                                                                a2dVarIceThick, &
                                                                                                a2dVarTaC_MeanDays10
                                                                                                
        real(kind = 4), dimension(iRowsEnd - iRowsStart + 1, iColsEnd - iColsStart + 1, iDaySteps)   :: a3dVarTaC_1Days

        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Initialize variable(s)
        a2dVarDEM = -9999; a2dVarLat = -9999; a2dVarLon = -9999 
        
        a2dVarRainFall = -9999; a2dVarSnowFall = -9999; a2dVarPrecip = -9999; a2dVarSnowFallDayCum = -9999;
        a2dVarSWE_D = -9999; a2dVarRho_D = -9999; a2dVarSWE_W = -9999; a2dVarTheta_W = -9999
        a2dVarSWE = -9999; a2dVarRhoS = -9999; a2dVarH_S = -9999; a2dVarRhoS0 = -9999;
        a2dVarSnowMask = -9999; a2dVarMeltingS = -9999; a2dVarMeltingSDayCum = -9999; a2dVarRefreezingS = -9999;
        a2dVarAlbedoS = -9999; a2dVarMeltingG = -9999; a2dVarMeltingGCumWY = -9999.0; a2iVarAgeS = 0;
        a2dVarTaC_MeanDays10 = -9999.0; a2dVarChangeThickness = -9999.0;
        a2dVarOutflow = -9999; a2dVarReff = -9999; 
        a2dVarIceThick = -9999; 
        
        sCommandCreateFolder = ""
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Defining iRows and iCols (output data)
        iRows = iRowsEnd - iRowsStart + 1
        iCols = iColsEnd - iColsStart + 1
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Get global information
        sPathData_Output = oS3M_Namelist(iID)%sPathData_Output_Gridded
        sCommandCreateFolder = oS3M_Namelist(iID)%sCommandCreateFolder
        iDaySteps = oS3M_Namelist(iID)%iDaySteps
        iFlagIceMassBalance = oS3M_Namelist(iID)%iFlagIceMassBalance
        sWYstart = oS3M_Namelist(iID)%sWYstart 
        iFlagOutputMode = oS3M_Namelist(iID)%iFlagOutputMode

        ! Info start
        call mprintf(.true., iINFO_Extra, ' Data :: Output gridded ... ' )
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Replace general path with specific time feature(s)
        call S3M_Tools_Generic_ReplaceText(sPathData_Output, '$yyyy', sTime(1:4))
        call S3M_Tools_Generic_ReplaceText(sPathData_Output, '$mm', sTime(6:7))
        call S3M_Tools_Generic_ReplaceText(sPathData_Output, '$dd', sTime(9:10))
        call S3M_Tools_Generic_ReplaceText(sPathData_Output, '$HH', sTime(12:13))
        call S3M_Tools_Generic_ReplaceText(sPathData_Output, '$MM', sTime(15:16))
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Create output folder
        call S3M_Tools_Generic_CreateFolder(sCommandCreateFolder, sPathData_Output, .true.)
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Get variable(s) from global workspace (terrain)
        a2dVarDEM = oS3M_Vars(iID)%a2dDem
        a2dVarLat = oS3M_Vars(iID)%a2dLat; 
        a2dVarLon = oS3M_Vars(iID)%a2dLon; 
        
        ! Get variable(s) from global workspace (snow physics)
        !SNOWFALL AND RAINFALL
        a2dVarRainFall = oS3M_Vars(iID)%a2dRainFall           
        a2dVarSnowFall = oS3M_Vars(iID)%a2dSnowFall
        a2dVarPrecip = oS3M_Vars(iID)%a2dPrecip
        a2dVarSnowFallDayCum = oS3M_Vars(iID)%a2dSnowFallDayCum
        
        !DRY WET SNOW
        a2dVarSWE_D = oS3M_Vars(iID)%a2dSWE_D
        a2dVarRho_D =  oS3M_Vars(iID)%a2dRho_D        
        a2dVarSWE_W = oS3M_Vars(iID)%a2dSWE_W
        a2dVarTheta_W = oS3M_Vars(iID)%a2dTheta_W           
        
        !BULK SNOW
        a2dVarSWE = oS3M_Vars(iID)%a2dSWE
        a2dVarRhoS = oS3M_Vars(iID)%a2dRhoS
        a2dVarRhoS0 = oS3M_Vars(iID)%a2dRhoS0        
        a2dVarH_S = oS3M_Vars(iID)%a2dH_S        
        a2dVarSnowMask = oS3M_Vars(iID)%a2dMaskS
        
        !SNOWMELT AND REFREEZING
        a2dVarMeltingS = oS3M_Vars(iID)%a2dMelting
        a2dVarMeltingSDayCum = oS3M_Vars(iID)%a2dMeltingDayCum
        a2dVarRefreezingS = oS3M_Vars(iID)%a2dRefreezingS
        a2dVarMeltingG = oS3M_Vars(iID)%a2dMeltingG
        a2dVarAlbedoS = oS3M_Vars(iID)%a2dAlbedo_Snow 
        a2iVarAgeS = oS3M_Vars(iID)%a2iAge
        
        !SNOW HYDRAULICS
        a2dVarOutflow = oS3M_Vars(iID)%a2dOutflow        
        a2dVarReff = oS3M_Vars(iID)%a2dReff        
        
        !GLACIERS
        a2dVarIceThick = oS3M_Vars(iID)%a2dIceThick
        a2dVarMeltingGCumWY = oS3M_Vars(iID)%a2dMeltingGCumWY
        a2dVarChangeThickness = oS3M_Vars(iID)%a2dChangeThickness
        
        !MISCELLANEOUS
        a3dVarTaC_1Days = oS3M_Vars(iID)%a3dTaC_Days1
        a2dVarTaC_MeanDays10 = oS3M_Vars(iID)%a2dTaC_MeanDays10
     
        !------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        ! Subroutine for writing sequential netCDF output data 
        call S3M_Data_Output_Gridded_NC(iID, &
                                sPathData_Output, &
                                iRows, iCols, &
                                iTime, sTime, iDaySteps, &
                                iFlagIceMassBalance, iFlagOutputMode, &
                                a2dVarRainFall, a2dVarSnowFall, a2dVarPrecip, a2dVarSnowFallDayCum, &
                                a2dVarSWE_D, a2dVarRho_D, a2dVarSWE_W, a2dVarTheta_W, &
                                a2dVarSWE, a2dVarRhoS, a2dVarH_S, a2dVarRhoS0, a2dVarSnowMask, &
                                a2dVarMeltingS, a2dVarMeltingSDayCum, a2dVarRefreezingS, a2dVarAlbedoS, &
                                a2dVarMeltingG, a2dVarMeltingGCumWY, a2dVarChangeThickness, a2iVarAgeS, &
                                a2dVarTaC_MeanDays10, &
                                a2dVarOutflow, a2dVarReff, &
                                a2dVarIceThick, &
                                a2dVarLat, a2dVarLon, a3dVarTaC_1Days, &
                                sWYstart)                    
        !------------------------------------------------------------------------------------------ 
            
        !------------------------------------------------------------------------------------------
        ! Info end
        call mprintf(.true., iINFO_Extra, ' Data :: Output gridded ... OK' )
        !------------------------------------------------------------------------------------------
        
    end subroutine S3M_Data_Output_Gridded_Cpl
    !------------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------------
    ! Subroutine to write netCDF gridded data output
#ifdef LIB_NC
    subroutine S3M_Data_Output_Gridded_NC(iID,  &
                                          sPathData_Output, &
                                          iRows, iCols, &
                                          iTime, sTime, iDaySteps, &
                                          iFlagIceMassBalance, iFlagOutputMode, &
                                          a2dVarRainFall, a2dVarSnowFall, a2dVarPrecip, a2dVarSnowFallDayCum, &
                                          a2dVarSWE_D, a2dVarRho_D, a2dVarSWE_W, a2dVarTheta_W, &
                                          a2dVarSWE, a2dVarRhoS, a2dVarH_S, a2dVarRhoS0, a2dVarSnowMask, &
                                          a2dVarMeltingS, a2dVarMeltingSDayCum, a2dVarRefreezingS, a2dVarAlbedoS, &
                                          a2dVarMeltingG, a2dVarMeltingGCumWY, a2dVarChangeThickness, a2iVarAgeS, &
                                          a2dVarTaC_MeanDays10, &
                                          a2dVarOutflow, a2dVarReff, &
                                          a2dVarIceThick, &
                                          a2dVarLat, a2dVarLon, a3dVarTaC_1Days, &
                                          sWYstart)
                                      
        !------------------------------------------------------------------------------------------
        ! Variable(s)
        integer(kind = 4)                       :: iID                  
           
        character(len = 700), intent(in)        :: sPathData_Output
        character(len = 700)                    :: sFileNameData_Output
        character(len = 700)                    :: sCommandZipFile
        character(len = 256)                    :: sVarName, sVarNameLong
        character(len = 256)                    :: sVarGridMap, sVarDescription, sVarCoords
        character(len = 256)                    :: sVarUnits
        character(len = 2)                      :: sWYstart
        integer(kind = 4), intent(in)           :: iRows, iCols
        
        integer(kind = 4)                       :: iTime, iDaySteps, iFlagIceMassBalance, iTimeStep, iFlagOutputMode
        character(len = 19)                     :: sTime, sTimeSave, sTimeRef
        
        real(kind = 4)                          :: dVarMissingValue
        real(kind = 4)                          :: dVarXLLCorner, dVarYLLCorner
        real(kind = 4)                          :: dVarCellSizeX, dVarCellSizeY
        
        integer(kind = 4)       :: iErr
        integer(kind = 4)       :: iFileID
        integer(kind = 4)       :: iID_Dim_Rows, iID_Dim_Cols, iID_Dim_Time_1Day, iID_Dim_Time_DBL, iID_Dim_Time_STR
        
        !SNOWFALL AND RAINFALL
        real(kind = 4), dimension(iRows, iCols), intent(in)         :: a2dVarRainFall, a2dVarSnowFall, a2dVarPrecip, &
                                                                       a2dVarSnowFallDayCum
        
        !DRY WET SNOW
        real(kind = 4), dimension(iRows, iCols), intent(in)         :: a2dVarSWE_D, a2dVarRho_D, a2dVarSWE_W, &
                                                                       a2dVarTheta_W
        
        !BULK SNOW
        real(kind = 4), dimension(iRows, iCols), intent(in)         :: a2dVarSWE, a2dVarRhoS, a2dVarRhoS0, a2dVarH_S, &
                                                                       a2dVarSnowMask
                                                                       
        !SNOWMELT AND REFREEZING
        real(kind = 4), dimension(iRows, iCols), intent(in)         :: a2dVarMeltingS, a2dVarMeltingSDayCum, &
                                                                       a2dVarRefreezingS, a2dVarMeltingG, a2dVarAlbedoS, &
                                                                       a2dVarTaC_MeanDays10
        integer(kind = 4), dimension(iRows, iCols), intent(in)      :: a2iVarAgeS

        !SNOW HYDRAULICS
        real(kind = 4), dimension(iRows, iCols), intent(in)         :: a2dVarOutflow, a2dVarReff 
        
        !GLACIERS
        real(kind = 4), dimension(iRows, iCols), intent(in)         :: a2dVarIceThick, a2dVarMeltingGCumWY, &
                                                                       a2dVarChangeThickness

        !MISCELLANEOUS
        real(kind = 4), dimension(iRows, iCols), intent(in)         :: a2dVarLat
        real(kind = 4), dimension(iRows, iCols), intent(in)         :: a2dVarLon
        real(kind = 4), dimension(iRows, iCols, iDaySteps)          :: a3dVarTaC_1Days
        
        real(kind = 8)          :: dTimeHours_Diff, dTimeSeconds_Diff
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Initialize variable(s)
        sCommandZipFile = ""
        sTimeRef = '1970-01-01_00:00:00'
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Get information from global workspace
        sCommandZipFile = oS3M_Namelist(iID)%sCommandZipFile
        
        dVarXLLCorner = oS3M_Namelist(iID)%dXLLCornerL
        dVarYLLCorner = oS3M_Namelist(iID)%dYLLCornerL
        dVarCellSizeX = oS3M_Namelist(iID)%dXCellSizeL
        dVarCellSizeY = oS3M_Namelist(iID)%dYCellSizeL
        
        ! Info start       
        call mprintf(.true., iINFO_Extra, ' Data :: Output gridded :: NetCDF ... ' )
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Define elapsed hours from reference date
        call S3M_Tools_Time_DateDiff(sTime, sTimeRef, dTimeSeconds_Diff, dTimeHours_Diff)
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Filename output (example: S3M_201404300000.nc)
        sFileNameData_Output = trim(sPathData_Output)//"S3M_"// &
        sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
        sTime(12:13)//sTime(15:16)// &
        ".nc"

        ! Info netCDF filename
        call mprintf(.true., iINFO_Verbose, ' Save filename (result gridded): '//trim(sFileNameData_Output)//' ... ')
        !------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        ! Create netcdf file
        call check( nf90_create(trim(sFileNameData_Output), NF90_NETCDF4, iFileID) )
	
        ! Dimension(s)
        call check( nf90_def_dim(iFileID, "time", 1, iID_Dim_Time_DBL) )
        call check( nf90_def_dim(iFileID, "time_str_length", 19, iID_Dim_Time_STR) )
        call check( nf90_def_dim(iFileID, "time_1day", iDaySteps, iID_Dim_Time_1Day) )
        call check( nf90_def_dim(iFileID, "Y", iRows, iID_Dim_Rows) )
        call check( nf90_def_dim(iFileID, "X", iCols, iID_Dim_Cols) )
        !------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        ! Global attribute(s)
        sTimeSave(1:len_trim(sTime)) = sTime
        call check( nf90_put_att(iFileID, NF90_GLOBAL, "time_coverage_end", sTimeSave) )
        call check( nf90_put_att(iFileID, NF90_GLOBAL, "xllcorner", dVarXLLCorner) )
        call check( nf90_put_att(iFileID, NF90_GLOBAL, "yllcorner", dVarYLLCorner) )
        call check( nf90_put_att(iFileID, NF90_GLOBAL, "xcellsize", dVarCellSizeX) )
        call check( nf90_put_att(iFileID, NF90_GLOBAL, "ycellsize", dVarCellSizeY) )
        call check( nf90_put_att(iFileID, NF90_GLOBAL, "nrows", iRows) )        
        call check( nf90_put_att(iFileID, NF90_GLOBAL, "ncols", iCols) )
        !------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        ! Definition mode OFF - Data mode ON
        call check( nf90_enddef(iFileID))
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Writing time variable(s) in netcdf output file
        ! TIMES
        sVarName = 'times'; sVarNameLong = 'times definition of output datasets'; 
        sVarDescription = 'times';
        sVarUnits = 'time units'; 
        sVarCoords = ''; iTimeStep = 1
        call S3M_Tools_IO_PutTime_STR_NC(iFileID, iID_Dim_Time_STR, & 
                             sVarName, sVarNameLong, sVarDescription, &
                             sVarUnits, sVarCoords, &
                             iTimeStep, sTimeSave)
        ! TIME
        sVarName = 'time'; sVarNameLong = 'time definition of output datasets'; 
        sVarDescription = 'synthesized time coordinate from Times(time)'; sVarUnits = 'secs since '//trim(sTimeRef); 
        sVarCoords = 't'; iTimeStep = 1
        call S3M_Tools_IO_PutTime_DBL_NC(iFileID, iID_Dim_Time_DBL, & 
                             sVarName, sVarNameLong, sVarDescription, &
                             sVarUnits, sVarCoords, &
                             iTimeStep, dTimeSeconds_Diff)
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Writing static variable(s) in netcdf output file
        ! LONGITUDE
        sVarName = 'Longitude'; sVarNameLong = 'Longitude Coordinate'; sVarDescription = 'longitude';
        sVarUnits = 'degree_east'; sVarGridMap = 'epsg:4326'; dVarMissingValue = -9E15;
        sVarCoords = '';
        call S3M_Tools_IO_Put2d_NC(iFileID, iID_Dim_Cols, iID_Dim_Rows, & 
                             sVarName, sVarNameLong, sVarDescription, &
                             sVarUnits, sVarCoords, sVarGridMap, dVarMissingValue, &
                             iCols, iRows, transpose(a2dVarLon))
        ! LATITUDE
        sVarName = 'Latitude'; sVarNameLong = 'Latitude Coordinate'; sVarDescription = 'latitude';
        sVarUnits = 'degree_north'; sVarGridMap = 'epsg:4326'; dVarMissingValue = -9E15;
        sVarCoords = '';
        call S3M_Tools_IO_Put2d_NC(iFileID, iID_Dim_Cols, iID_Dim_Rows, & 
                             sVarName, sVarNameLong, sVarDescription, &
                             sVarUnits, sVarCoords, sVarGridMap, dVarMissingValue, &
                             iCols, iRows, transpose(a2dVarLat))
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Writing dynamic variable(s) in netCDF output file

        if ( (iFlagOutputMode .eq. 1)) then                       
                             
            !SNOWFALL AND RAINFALL
            sVarName = 'SnowFall'; sVarNameLong = 'Snowfall Amount'; sVarDescription = 'snowfall amount';
            sVarUnits = 'mm'; sVarGridMap = 'epsg:4326'; dVarMissingValue = -9E15;
            sVarCoords = 'Longitude Latitude';
            call S3M_Tools_IO_Put2d_NC(iFileID, iID_Dim_Cols, iID_Dim_Rows, & 
                                 sVarName, sVarNameLong, sVarDescription, &
                                 sVarUnits, sVarCoords, sVarGridMap, dVarMissingValue, &
                                 iCols, iRows, transpose(a2dVarSnowFall))

            sVarName = 'Precip'; sVarNameLong = 'Total Precipitation Amount'; sVarDescription = 'total precipitation amount';
            sVarUnits = 'mm'; sVarGridMap = 'epsg:4326'; dVarMissingValue = -9E15;
            sVarCoords = 'Longitude Latitude';
            call S3M_Tools_IO_Put2d_NC(iFileID, iID_Dim_Cols, iID_Dim_Rows, & 
                                 sVarName, sVarNameLong, sVarDescription, &
                                 sVarUnits, sVarCoords, sVarGridMap, dVarMissingValue, &
                                 iCols, iRows, transpose(a2dVarPrecip))

            sVarName = 'RainFall'; sVarNameLong = 'Rainfall Amount'; sVarDescription = 'rainfall amount';
            sVarUnits = 'mm'; sVarGridMap = 'epsg:4326'; dVarMissingValue = -9E15;
            sVarCoords = 'Longitude Latitude';
            call S3M_Tools_IO_Put2d_NC(iFileID, iID_Dim_Cols, iID_Dim_Rows, & 
                                 sVarName, sVarNameLong, sVarDescription, &
                                 sVarUnits, sVarCoords, sVarGridMap, dVarMissingValue, &
                                 iCols, iRows, transpose(a2dVarRainFall))                             
                             
        endif
        
        !DRY WET SNOW
        sVarName = 'SWE_D'; sVarNameLong = 'Dry SWE'; sVarDescription = 'dry SWE';
        sVarUnits = 'mm'; sVarGridMap = 'epsg:4326'; dVarMissingValue = -9E15;
        sVarCoords = 'Longitude Latitude';
        call S3M_Tools_IO_Put2d_NC(iFileID, iID_Dim_Cols, iID_Dim_Rows, & 
                             sVarName, sVarNameLong, sVarDescription, &
                             sVarUnits, sVarCoords, sVarGridMap, dVarMissingValue, &
                             iCols, iRows, transpose(a2dVarSWE_D))                             

        sVarName = 'Rho_D'; sVarNameLong = 'Dry Snow Density'; sVarDescription = 'dry snow density';
        sVarUnits = 'kg/m^3'; sVarGridMap = 'epsg:4326'; dVarMissingValue = -9E15;
        sVarCoords = 'Longitude Latitude';
        call S3M_Tools_IO_Put2d_NC(iFileID, iID_Dim_Cols, iID_Dim_Rows, & 
                             sVarName, sVarNameLong, sVarDescription, &
                             sVarUnits, sVarCoords, sVarGridMap, dVarMissingValue, &
                             iCols, iRows, transpose(a2dVarRho_D))                               

        sVarName = 'SWE_W'; sVarNameLong = 'Wet SWE'; sVarDescription = 'wet SWE';
        sVarUnits = 'mm'; sVarGridMap = 'epsg:4326'; dVarMissingValue = -9E15;
        sVarCoords = 'Longitude Latitude';
        call S3M_Tools_IO_Put2d_NC(iFileID, iID_Dim_Cols, iID_Dim_Rows, & 
                             sVarName, sVarNameLong, sVarDescription, &
                             sVarUnits, sVarCoords, sVarGridMap, dVarMissingValue, &
                             iCols, iRows, transpose(a2dVarSWE_W))     
                                 
        if ( (iFlagOutputMode .eq. 1)) then       

            sVarName = 'Theta_W'; sVarNameLong = 'Bulk Vol. LWC'; sVarDescription = 'bulk vol. LWC';
            sVarUnits = '-'; sVarGridMap = 'epsg:4326'; dVarMissingValue = -9E15;
            sVarCoords = 'Longitude Latitude';
            call S3M_Tools_IO_Put2d_NC(iFileID, iID_Dim_Cols, iID_Dim_Rows, & 
                                 sVarName, sVarNameLong, sVarDescription, &
                                 sVarUnits, sVarCoords, sVarGridMap, dVarMissingValue, &
                                 iCols, iRows, transpose(a2dVarTheta_W))                              
                             
        endif
                                 
        !BULK SNOW 
        sVarName = 'SWE'; sVarNameLong = 'Snow Water Equivalent'; sVarDescription = 'total SWE';
        sVarUnits = 'mm'; sVarGridMap = 'epsg:4326'; dVarMissingValue = -9E15;
        sVarCoords = 'Longitude Latitude';
        call S3M_Tools_IO_Put2d_NC(iFileID, iID_Dim_Cols, iID_Dim_Rows, & 
                             sVarName, sVarNameLong, sVarDescription, &
                             sVarUnits, sVarCoords, sVarGridMap, dVarMissingValue, &
                             iCols, iRows, transpose(a2dVarSWE))  

        if ( (iFlagOutputMode .eq. 1)) then       

            sVarName = 'RhoS0'; sVarNameLong = 'Fresh-Snow Density'; sVarDescription = 'fresh-snow density';
            sVarUnits = 'kg/m^3'; sVarGridMap = 'epsg:4326'; dVarMissingValue = -9E15;
            sVarCoords = 'Longitude Latitude';
            call S3M_Tools_IO_Put2d_NC(iFileID, iID_Dim_Cols, iID_Dim_Rows, & 
                                sVarName, sVarNameLong, sVarDescription, &
                                sVarUnits, sVarCoords, sVarGridMap, dVarMissingValue, &
                                iCols, iRows, transpose(a2dVarRhoS0))                             
        
            sVarName = 'RhoS'; sVarNameLong = 'Bulk-Snow Density'; sVarDescription = 'bulk-snow density';
            sVarUnits = 'kg/m^3'; sVarGridMap = 'epsg:4326'; dVarMissingValue = -9E15;
            sVarCoords = 'Longitude Latitude';
            call S3M_Tools_IO_Put2d_NC(iFileID, iID_Dim_Cols, iID_Dim_Rows, & 
                                sVarName, sVarNameLong, sVarDescription, &
                                sVarUnits, sVarCoords, sVarGridMap, dVarMissingValue, &
                                iCols, iRows, transpose(a2dVarRhoS))

            sVarName = 'H_S'; 
            sVarNameLong = 'Bulk Snow Depth'; sVarDescription = 'bulk snow depth (control volume)';
            sVarUnits = 'm'; sVarGridMap = 'epsg:4326'; dVarMissingValue = -9E15;
            sVarCoords = 'Longitude Latitude';
            call S3M_Tools_IO_Put2d_NC(iFileID, iID_Dim_Cols, iID_Dim_Rows, & 
                                sVarName, sVarNameLong, sVarDescription, &
                                sVarUnits, sVarCoords, sVarGridMap, dVarMissingValue, &
                                iCols, iRows, transpose(a2dVarH_S))                            

        endif
                       
        sVarName = 'SnowMask'; 
        sVarNameLong = 'Snow Mask'; sVarDescription = 'snow mask';
        sVarUnits = '-'; sVarGridMap = 'epsg:4326'; dVarMissingValue = -9E15;
        sVarCoords = 'Longitude Latitude';
        call S3M_Tools_IO_Put2d_NC(iFileID, iID_Dim_Cols, iID_Dim_Rows, & 
                            sVarName, sVarNameLong, sVarDescription, &
                            sVarUnits, sVarCoords, sVarGridMap, dVarMissingValue, &
                            iCols, iRows, transpose(a2dVarSnowMask))  
          
        !SNOWMELT AND REFREEZING                     
        if ( (iFlagOutputMode .eq. 1)) then 
            
            sVarName = 'MeltingS'; sVarNameLong = 'Snow Melt'; sVarDescription = 'snow melt';
            sVarUnits = 'mm'; sVarGridMap = 'epsg:4326'; dVarMissingValue = -9E15;
            sVarCoords = 'Longitude Latitude';
            call S3M_Tools_IO_Put2d_NC(iFileID, iID_Dim_Cols, iID_Dim_Rows, & 
                                sVarName, sVarNameLong, sVarDescription, &
                                sVarUnits, sVarCoords, sVarGridMap, dVarMissingValue, &
                                iCols, iRows, transpose(a2dVarMeltingS))                             

            sVarName = 'RefreezingS'; sVarNameLong = 'Snow Refreezing'; sVarDescription = 'snow refreezing';
            sVarUnits = 'mm'; sVarGridMap = 'epsg:4326'; dVarMissingValue = -9E15;
            sVarCoords = 'Longitude Latitude';
            call S3M_Tools_IO_Put2d_NC(iFileID, iID_Dim_Cols, iID_Dim_Rows, & 
                                sVarName, sVarNameLong, sVarDescription, &
                                sVarUnits, sVarCoords, sVarGridMap, dVarMissingValue, &
                                iCols, iRows, transpose(a2dVarRefreezingS))                              
        
            sVarName = 'MeltingG'; sVarNameLong = 'Glacier Melt'; sVarDescription = 'glacier melt';
            sVarUnits = 'mm'; sVarGridMap = 'epsg:4326'; dVarMissingValue = -9E15;
            sVarCoords = 'Longitude Latitude';
            call S3M_Tools_IO_Put2d_NC(iFileID, iID_Dim_Cols, iID_Dim_Rows, & 
                                sVarName, sVarNameLong, sVarDescription, &
                                sVarUnits, sVarCoords, sVarGridMap, dVarMissingValue, &
                                iCols, iRows, transpose(a2dVarMeltingG))   

        endif 
                    
        ! Snow albedo
        sVarName = 'AlbedoS'; sVarNameLong = 'Snow Albedo'; sVarDescription = 'snow albedo';
        sVarUnits = '-'; sVarGridMap = 'epsg:4326'; dVarMissingValue = -9E15;
        sVarCoords = 'Longitude Latitude';
        call S3M_Tools_IO_Put2d_NC(iFileID, iID_Dim_Cols, iID_Dim_Rows, & 
                                 sVarName, sVarNameLong, sVarDescription, &
                                 sVarUnits, sVarCoords, sVarGridMap, dVarMissingValue, &
                                 iCols, iRows, transpose(a2dVarAlbedoS))                 
        ! Snow age
        sVarName = 'AgeS'; sVarNameLong = 'Snow Age'; sVarDescription = 'snow age';
        sVarUnits = 'timestep'; sVarGridMap = 'epsg:4326'; dVarMissingValue = -9E15;                
        sVarCoords = 'Longitude Latitude';
            call S3M_Tools_IO_Put2d_NC(iFileID, iID_Dim_Cols, iID_Dim_Rows, & 
                                 sVarName, sVarNameLong, sVarDescription, &
                                 sVarUnits, sVarCoords, sVarGridMap, dVarMissingValue, &
                                 iCols, iRows, transpose(real(a2iVarAgeS)))  

        ! T 10 days
        sVarName = 'T_10Days'; sVarNameLong = 'Average T 10 Days'; sVarDescription = 'avg T 10 days';
        sVarUnits = 'C'; sVarGridMap = 'epsg:4326'; dVarMissingValue = -9E15;                
        sVarCoords = 'Longitude Latitude';
        call S3M_Tools_IO_Put2d_NC(iFileID, iID_Dim_Cols, iID_Dim_Rows, & 
                                 sVarName, sVarNameLong, sVarDescription, &
                                 sVarUnits, sVarCoords, sVarGridMap, dVarMissingValue, &
                                 iCols, iRows, transpose(a2dVarTaC_MeanDays10))                                  
        
        if ( (iFlagOutputMode .eq. 1)) then   
            
            !SNOW HYDRAULICS                    
            sVarName = 'Outflow'; sVarNameLong = 'Snowpack Runoff'; sVarDescription = 'snowpack runoff';
            sVarUnits = 'mm'; sVarGridMap = 'epsg:4326'; dVarMissingValue = -9E15;
            sVarCoords = 'Longitude Latitude';
            call S3M_Tools_IO_Put2d_NC(iFileID, iID_Dim_Cols, iID_Dim_Rows, & 
                                 sVarName, sVarNameLong, sVarDescription, &
                                 sVarUnits, sVarCoords, sVarGridMap, dVarMissingValue, &
                                 iCols, iRows, transpose(a2dVarOutflow))  
                                 
        endif
        
        sVarName = 'REff'; sVarNameLong = 'Effective Rainfall'; sVarDescription = 'Reff';
        sVarUnits = 'mm'; sVarGridMap = 'epsg:4326'; dVarMissingValue = -9E15;
        sVarCoords = 'Longitude Latitude';
        call S3M_Tools_IO_Put2d_NC(iFileID, iID_Dim_Cols, iID_Dim_Rows, & 
                             sVarName, sVarNameLong, sVarDescription, &
                             sVarUnits, sVarCoords, sVarGridMap, dVarMissingValue, &
                             iCols, iRows, transpose(a2dVarReff))  
                             
        ! GLACIERS                     
        if ( (iFlagIceMassBalance.eq.1.0) .or. (iFlagIceMassBalance.eq.2.0) ) then
            
            sVarName = 'Ice_Thickness'; sVarNameLong = 'Ice Thickness'; sVarDescription = 'icethickness';
            sVarUnits = 'm'; sVarGridMap = 'epsg:4326'; dVarMissingValue = -9E15;
            sVarCoords = 'Longitude Latitude';
            call S3M_Tools_IO_Put2d_NC(iFileID, iID_Dim_Cols, iID_Dim_Rows, &
                                 sVarName, sVarNameLong, sVarDescription, &
                                 sVarUnits, sVarCoords, sVarGridMap, dVarMissingValue, &
                                 iCols, iRows, transpose(a2dVarIceThick))
                                 
            sVarName = 'Cum_WY_MeltingG'; sVarNameLong = 'Cumulative WY Glacier Melt'; sVarDescription = 'wyglaciermelt';
            sVarUnits = 'mm'; sVarGridMap = 'epsg:4326'; dVarMissingValue = -9E15;
            sVarCoords = 'Longitude Latitude';
            call S3M_Tools_IO_Put2d_NC(iFileID, iID_Dim_Cols, iID_Dim_Rows, &
                                 sVarName, sVarNameLong, sVarDescription, &
                                 sVarUnits, sVarCoords, sVarGridMap, dVarMissingValue, &
                                 iCols, iRows, transpose(a2dVarMeltingGCumWY))     
                                 
            sVarName = 'Ice_Thickness_Change'; sVarNameLong = 'Ice Thickness Change'; sVarDescription = 'icethicknesschange';
            sVarUnits = 'm'; sVarGridMap = 'epsg:4326'; dVarMissingValue = -9E15;
            sVarCoords = 'Longitude Latitude';
            call S3M_Tools_IO_Put2d_NC(iFileID, iID_Dim_Cols, iID_Dim_Rows, &
                                 sVarName, sVarNameLong, sVarDescription, &
                                 sVarUnits, sVarCoords, sVarGridMap, dVarMissingValue, &
                                 iCols, iRows, transpose(a2dVarChangeThickness))
                                 
        endif                             
                  
        ! Daily variable(s)     
        if( sTime(12:13) .eq. '23' ) then
                
            ! Snow melting daily cumulated           
            sVarName = 'MeltingSDayCum'; sVarNameLong = 'Daily Cumulative Snow Melt'; sVarDescription = 'daily cum snow melt';
            sVarUnits = 'mm'; sVarGridMap = 'epsg:4326'; dVarMissingValue = -9E15;
            sVarCoords = 'Longitude Latitude';
            call S3M_Tools_IO_Put2d_NC(iFileID, iID_Dim_Cols, iID_Dim_Rows, & 
                                 sVarName, sVarNameLong, sVarDescription, &
                                 sVarUnits, sVarCoords, sVarGridMap, dVarMissingValue, &
                                 iCols, iRows, transpose(a2dVarMeltingSDayCum))
                                 
            ! Snowfall daily cumulated                       
            sVarName = 'SnowfallCum'; sVarNameLong = 'Daily Cumulative Snowfall'; sVarDescription = 'daily cum snowfall';
            sVarUnits = 'mm'; sVarGridMap = 'epsg:4326'; dVarMissingValue = -9E15;
            sVarCoords = 'Longitude Latitude';
            call S3M_Tools_IO_Put2d_NC(iFileID, iID_Dim_Cols, iID_Dim_Rows, & 
                                 sVarName, sVarNameLong, sVarDescription, &
                                 sVarUnits, sVarCoords, sVarGridMap, dVarMissingValue, &
                                 iCols, iRows, transpose(a2dVarSnowFallDayCum))                                  
                                 
            ! Air temperature last 1 day(s)           
            sVarName = 'T_1Days'; sVarNameLong = 'Air Temperature Last 1 Day'; sVarDescription = 't_1days';
            sVarUnits = 'C'; sVarGridMap = 'epsg:4326'; dVarMissingValue = -9E15;
            sVarCoords = 'Longitude Latitude';
            call S3M_Tools_IO_Put3d_NC(iFileID, iID_Dim_Cols, iID_Dim_Rows, iID_Dim_Time_1Day, & 
                                 sVarName, sVarNameLong, sVarDescription, &
                                 sVarUnits, sVarCoords, sVarGridMap, dVarMissingValue, &
                                 iCols, iRows, iDaySteps, transpose3Dvar(a3dVarTaC_1Days))                                   
                                 
        endif
            
        !------------------------------------------------------------------------------------------
                
        !------------------------------------------------------------------------------------------
        ! Close
        call check( nf90_close(iFileID) )
        ! Info
        call mprintf(.true., iINFO_Verbose, ' Save filename (result gridded): '//trim(sFileNameData_Output)//' ... OK')
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Zip file
        call S3M_Tools_Generic_ZipFile(sCommandZipFile, &
                                       sFileNameData_Output//'.gz', sFileNameData_Output, .false.)
        ! Remove un-zipped file
        !call S3M_Tools_Generic_RemoveFile(oS3M_Namelist(iID)%sCommandRemoveFile, sFileNameData_Output, .false.)
                                       
        ! Info end      
        call mprintf(.true., iINFO_Extra, ' Data :: Output gridded :: NetCDF ... OK ' )
        !------------------------------------------------------------------------------------------
        
    end subroutine S3M_Data_Output_Gridded_NC
#endif
    !------------------------------------------------------------------------------------------
    
end module S3M_Module_Data_Output_Gridded
!------------------------------------------------------------------------------------------

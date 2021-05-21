!------------------------------------------------------------------------------------------    
! File:   S3M_Module_Data_Static_Gridded.f90
! Author(s): Fabio Delogu, Francesco Silvestro, Simone Gabellani, Francesco Avanzi.
!
! Created on April 30, 2015, 9:45 AM
! Last update on Dec 10, 2020 10:00 AM
!
! Module to read static data.
!------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
! Module Header
module S3M_Module_Data_Static_Gridded
    
    !------------------------------------------------------------------------------------
    ! External module(s) and implicit none
#ifdef LIB_NC
    use netcdf
#endif
    
    use S3M_Module_Namelist,        only:   oS3M_Namelist
    use S3M_Module_Vars_Loader,     only:   oS3M_Vars
    
    use S3M_Module_Tools_Debug
    use gnufor2
    
#ifdef LIB_NC
    use S3M_Module_Tools_IO,        only:   S3M_Tools_IO_Get1d_ASCII, &  
                                            S3M_Tools_IO_GetArcGrid_ASCII, &
                                            S3M_Tools_IO_Get1D_NC, &
                                            S3M_Tools_IO_Get2d_NC, &
                                            check
#else 
    use S3M_Module_Tools_IO,        only:   S3M_Tools_IO_Get1d_ASCII, &  
                                            S3M_Tools_IO_GetArcGrid_ASCII
#endif
                                                                             
    use S3M_Module_Tools_Generic,   only:   nullborder2DVar, check2Dvar, &
                                            S3M_Tools_Generic_UnzipFile
                                            
        
    implicit none
    !------------------------------------------------------------------------------------
    
contains

    !------------------------------------------------------------------------------------
    ! Subroutine to couple main.f90 and various subroutines related to static data
    subroutine S3M_Data_Static_Gridded_Cpl(iID, iRows, iCols, iRowsPivot, iColsPivot)

        !------------------------------------------------------------------------------------
        ! Variable(s) declaration
        integer(kind = 4)       :: iID, iRows, iCols, iRowsPivot, iColsPivot
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Info start
        call mprintf(.true., iINFO_Verbose, ' Data :: Static gridded ... ' )
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Calling subroutine to compute static data land
        call S3M_Data_Static_Gridded_Land(iID, iRows, iCols, iRowsPivot, iColsPivot)
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Calling subroutine to compute extra derived model parameter(s)
        call S3M_Data_Static_Gridded_Params(iID, iRows, iCols)
        !------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------
        ! Info end
        call mprintf(.true., iINFO_Verbose, ' Data :: Static gridded ... OK' )
        !------------------------------------------------------------------------------------
        
    
    end subroutine
    !------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------ 
    ! Subroutine to distribute model parameter(s) across the computational domain
    subroutine S3M_Data_Static_Gridded_Params(iID, iRows, iCols)
        
        !------------------------------------------------------------------------------------
        ! Variable(s)
        integer(kind = 4)       :: iID, iRows, iCols
        
        real(kind = 4)          :: dVarDEMMax
        
        integer(kind = 4), dimension (iRows, iCols)   :: a2iVarMask

        real(kind = 4), dimension (iRows, iCols)      :: a2dVarDEM 
        real(kind = 4), dimension (iRows, iCols)      :: a2dVarArctUp
        
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Variable(s) initialization
        dVarDEMMax = 0.0;
        a2dVarDEM = 0.0;
        a2dVarArctUp = 0.0; 
        !------------------------------------------------------------------------------------
      
        !------------------------------------------------------------------------------------ 
        ! Extracting parameters
        dVarDEMMax = oS3M_Vars(iID)%dDEMMax        
        a2dVarDEM = oS3M_Vars(iID)%a2dDem
        a2iVarMask = oS3M_Vars(iID)%a2iMask
        
        ! Info start
        call mprintf(.true., iINFO_Verbose, ' Data :: Static gridded :: Distribute parameter(s) ... ' )
        !------------------------------------------------------------------------------------ 

        !------------------------------------------------------------------------------------
        ! Distributing parameters by elevation bands
        where( a2dVarDEM .le. oS3M_Namelist(iID)%a1dAltRange(1) )
                a2dVarArctUp  = oS3M_Namelist(iID)%a1dArctUp(1)
        endwhere
        where( (a2dVarDEM .gt. oS3M_Namelist(iID)%a1dAltRange(1)) .and. &
                   (a2dVarDEM .le. oS3M_Namelist(iID)%a1dAltRange(2)) )
                a2dVarArctUp  = oS3M_Namelist(iID)%a1dArctUp(2)
        endwhere
        where( (a2dVarDEM .gt. oS3M_Namelist(iID)%a1dAltRange(2)) .and. &
                   (a2dVarDEM .le. oS3M_Namelist(iID)%a1dAltRange(3)) )
                a2dVarArctUp  = oS3M_Namelist(iID)%a1dArctUp(3)
        endwhere
        where( a2dVarDem .gt. oS3M_Namelist(iID)%a1dAltRange(3))
                a2dVarArctUp  = oS3M_Namelist(iID)%a1dArctUp(4)
        endwhere
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------        
        ! Passing parameter(s) to global declaration
        oS3M_Vars(iID)%a2dArctUp = a2dVarArctUp; 

        ! Info end
        call mprintf(.true., iINFO_Verbose, ' Data :: Static gridded :: Distribute parameter(s) ... OK' )
        !------------------------------------------------------------------------------------ 
        
    end subroutine S3M_Data_Static_Gridded_Params
    !------------------------------------------------------------------------------------ 
    
    !------------------------------------------------------------------------------------
    ! Subroutine for loading and initializing land variable(s)
    subroutine S3M_Data_Static_Gridded_Land(iID, iRows, iCols, iRowsPivot, iColsPivot)
        
        !------------------------------------------------------------------------------------------
        ! Variable(s)
        character(len = 256)    :: sDomainName, sPathData
        
        character(len = 700)    :: sFileName, sFileNameZip 
        character(len = 700)    :: sCommandUnzip
        
        integer(kind = 4)       :: iID, iI, iJ, iRows, iCols, iRowsPivot, iColsPivot
        integer(kind = 4)       :: iPixCount
        
        integer(kind = 4)       :: iFileID, iErr
        
        integer(kind = 4)       :: iNTime, iDomainPixels

        real(kind = 4)          :: dDEMMax, dDEMMin
        real(kind = 4)          :: dDomainArea
        
        character(len = 256)    :: sVarName
        character(len = 256)    :: sVarUnits
        
        real(kind = 4),     dimension (100)             :: a1dVar 
        real(kind = 4),     dimension (iCols, iRows)    :: a2dVar 
        
        integer(kind = 4),  dimension (iRows, iCols)    :: a2iVarMask
        integer(kind = 4),  dimension (iRows, iCols)    :: a2iVarGlacierMask, a2iVarGlaciers_ID
        
        real(kind = 4),     dimension (iRows, iCols)    :: a2dVarLon, a2dVarLat 
        real(kind = 4),     dimension (iRows, iCols)    :: a2dVarDEM, a2dVarGlacierDebris, a2dVarIceThick, a2dVarAreaCell
        
        real(kind = 4),     dimension (iRowsPivot, iColsPivot)    :: a2dVarGlacierPivotTable
        real(kind = 4),     dimension (iColsPivot, iRowsPivot)    :: a2dVar_PivotTable
        
        logical                                         :: bFileExist
        
        character(len = 256)                            :: sStrDomPix, sStrDomArea, sStrDemMax
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Variable(s) initialization
        bFileExist = .false.
        sVarName = ''; sVarUnits = ''; iNTime = 0;
        sDomainName = ''; sPathData = ''; sFileName = '';
        iI = 0; iJ = 0; iFileID = 0; iErr = 0; iDomainPixels = 0;
        
        dDEMMax = 0.0; dDEMMin = 0.0; dDomainArea = 0.0;
                
        a1dVar = 0.0; a2dVar = 0.0 
        
        a2dVarLon = 0.0; a2dVarLat = 0.0; 
        a2dVarDEM = 0.0; a2iVarMask = 0.0; a2iVarGlacierMask = 0; a2iVarGlaciers_ID = 0; a2dVarAreaCell = 0;
        a2dVarGlacierDebris = 0.0; a2dVarIceThick = 0.0; 
        a2dVarGlacierPivotTable = 0.0; a2dVar_PivotTable = 0.0;
        
        iPixCount = 0
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Get global info
        iNTime = oS3M_Namelist(iID)%iNTime
        sDomainName = oS3M_Namelist(iID)%sDomainName
        sPathData = oS3M_Namelist(iID)%sPathData_Static_Gridded
        sCommandUnzip = oS3M_Namelist(iID)%sCommandUnzipFile

        ! Info start
        call mprintf(.true., iINFO_Verbose, ' Data :: Static gridded :: Get land information ... ' )
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Land data type in netCDF format
#ifdef LIB_NC
            !------------------------------------------------------------------------------------------
            ! Info
            call mprintf(.true., iINFO_Extra, ' Data static gridded in netCDF format' )
            !------------------------------------------------------------------------------------------
            
            !------------------------------------------------------------------------------------------
            ! Filename nc
            sFileName = trim(sPathData)//'Terrain_Data.nc'

            ! Checking file input availability
            sFileNameZip = sFileName(1:len_trim(sFileName))//'.gz'
            inquire (file = trim(sFileNameZip), exist = bFileExist)
            if ( .not. bFileExist ) then
                ! Exit code file not found
                call mprintf(.true., iERROR, ' No compressed data static gridded file found: '//trim(sFileNameZip) )
            endif

            ! Unzip file
            call S3M_Tools_Generic_UnzipFile(sCommandUnzip, sFileNameZip, sFileName, .true.)

            ! Check file availability
            inquire (file = trim(sFileName), exist = bFileExist)
            if ( .not. bFileExist ) then
                ! Exit code file not found
                call mprintf(.true., iERROR, 'No data static gridded file found: '//trim(sFileName) )
            else         
                
                !------------------------------------------------------------------------------------------
                ! Open nc file
                call check(nf90_open(trim(sFileName), NF90_NOWRITE, iFileID) )
                !------------------------------------------------------------------------------------------
                
                !------------------------------------------------------------------------------------------
                ! Longitude
                sVarName = 'Longitude'
                call S3M_Tools_IO_Get2d_NC(sVarName, iFileID, a2dVar, sVarUnits, iCols, iRows, .true., iErr)
                a2dVarLon = transpose(a2dVar)
                !------------------------------------------------------------------------------------------
                
                !------------------------------------------------------------------------------------------
                ! Latitude
                sVarName = 'Latitude'
                call S3M_Tools_IO_Get2d_NC(sVarName, iFileID, a2dVar, sVarUnits, iCols, iRows, .true., iErr)
                a2dVarLat = transpose(a2dVar)
                !------------------------------------------------------------------------------------------
                
                !------------------------------------------------------------------------------------------
                ! DEM Terrain
                sVarName = 'Terrain'
                call S3M_Tools_IO_Get2d_NC(sVarName, iFileID, a2dVar, sVarUnits, iCols, iRows, .true., iErr)
                a2dVarDEM = transpose(a2dVar)           
                !------------------------------------------------------------------------------------------

                !------------------------------------------------------------------------------------------                
                ! Glacier mask
                sVarName = 'GlacierMask'
                call S3M_Tools_IO_Get2d_NC(sVarName, iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                if (iErr /= 0) then
                    call mprintf(.true., iWARN, ' Glacier mask not found. Initializing glacier mask with default values')
                    a2iVarGlacierMask = -9999
                else
                    a2iVarGlacierMask = int(transpose(a2dVar))
                endif
                !------------------------------------------------------------------------------------------
                
                !------------------------------------------------------------------------------------------                
                ! Glacier ID map
                sVarName = 'GlacierID'
                call S3M_Tools_IO_Get2d_NC(sVarName, iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                if (iErr /= 0) then
                    call mprintf(.true., iWARN, ' Glaciers ID map not found. Initializing Glaciers ID with undefined values')
                    a2iVarGlaciers_ID = -9999
                else
                    a2iVarGlaciers_ID = int(transpose(a2dVar))
                endif
                !------------------------------------------------------------------------------------------  
                
                !------------------------------------------------------------------------------------------                
                ! Glacier debris map
                sVarName = 'GlacierDebris'
                call S3M_Tools_IO_Get2d_NC(sVarName, iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                if (iErr /= 0) then
                    call mprintf(.true., iWARN, ' Glaciers debris map not found. Initializing debris with undefined values')
                    a2dVarGlacierDebris = -9999
                else
                    a2dVarGlacierDebris = transpose(a2dVar)
                endif
                !------------------------------------------------------------------------------------------ 

                !------------------------------------------------------------------------------------------                
                ! Glacier Thickness
                if (oS3M_Namelist(iID)%iFlagThickFromTerrData .eq. 1) then                
                    sVarName = 'Thickness'
                    call S3M_Tools_IO_Get2d_NC(sVarName, iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                    if (iErr /= 0) then
                        call mprintf(.true., iWARN, ' Glaciers debris map not found. Initializing debris with undefined values')
                        a2dVarIceThick = -9999
                    else
                        a2dVarIceThick = transpose(a2dVar)
                    endif
                endif
                !------------------------------------------------------------------------------------------ 
                
                !------------------------------------------------------------------------------------------                 
                ! Mask
                sVarName = 'Mask'
                call S3M_Tools_IO_Get2d_NC(sVarName, iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                
                if (iErr /= 0) then 
                    call mprintf(.true., iWARN, ' Mask data not found. Initializing Mask with default values')
                    where(a2dVarDEM.gt.0.0)
                        a2iVarMask = 1 
                    elsewhere
                        a2iVarMask = -9999
                    endwhere
                else
                    a2iVarMask = int(transpose(a2dVar))
                endif
                !------------------------------------------------------------------------------------------

                !------------------------------------------------------------------------------------------                
                ! Glacier deltaH pivot table
                sVarName = 'PivotTable'
                call S3M_Tools_IO_Get2d_NC(sVarName, iFileID, a2dVar_PivotTable, sVarUnits,&
                iColsPivot, iRowsPivot, .false., iErr)
                if (iErr /= 0) then
                    call mprintf(.true., iWARN, ' Pivot table not found. Initializing pivot table with undefined values')
                    a2dVarGlacierPivotTable = -9999
                else
                    a2dVarGlacierPivotTable = transpose(a2dVar_PivotTable)
                endif
                !------------------------------------------------------------------------------------------  

                !------------------------------------------------------------------------------------------                
                ! Area Cell
                sVarName = 'AreaCell'
                call S3M_Tools_IO_Get2d_NC(sVarName, iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                if (iErr /= 0) then
                    call mprintf(.true., iWARN, ' AreaCell not found. Initializing AreaCell with undefined values')
                    a2dVarAreaCell = -9999
                else
                    a2dVarAreaCell = int(transpose(a2dVar))
                endif
                !------------------------------------------------------------------------------------------                
                
                !------------------------------------------------------------------------------------------
                ! Closing netcdf file
                iErr = nf90_close(iFileID)
                !------------------------------------------------------------------------------------------
                
            endif
            !------------------------------------------------------------------------------------------
#endif
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Nullify map boundary
        a2dVarDEM = nullborder2DVar(a2dVarDEM, -9999.0)
        a2iVarMask = int(nullborder2DVar(float(a2iVarMask), -9999.0))

        where(a2dVarDEM.lt.0.0)
            a2iVarMask = -9999
            a2iVarGlacierMask = -9999
            a2iVarGlaciers_ID = -9999
            a2dVarGlacierDebris = -9999
            a2dVarIceThick = -9999
            a2dVarAreaCell = -9999
        endwhere
        !------------------------------------------------------------------------------------------
        
        ! Info end
        call mprintf(.true., iINFO_Verbose, ' Data :: Static gridded :: Get land information ... OK' )
        !------------------------------------------------------------------------------------ 
        
        !------------------------------------------------------------------------------------ 
        ! Land data derived fields
        call mprintf(.true., iINFO_Verbose, ' Data :: Static gridded :: Compute derived land information ... ' )
        
        ! DEM max and min values and step mean
        dDEMMax = maxval(maxval(a2dVarDEM,DIM = 1, MASK=a2dVarDEM.gt.0),DIM = 1)
        dDEMMin = minval(minval(a2dVarDEM,DIM = 1, MASK=a2dVarDEM.gt.0),DIM = 1)
       
        ! Computing total catchment pixels and area
        iDomainPixels = sum(sum(a2iVarMask,dim=1, mask=a2dVarDEM.gt.0.0))             
        
        ! Domain information
        write(sStrDomPix, *) iDomainPixels; write(sStrDomArea, *) dDomainArea; 
        write(sStrDemMax, *) dDEMMax;
        call mprintf(.true., iINFO_Main, ' DOMAIN INFO --- '// &
                    'NPixels : '//trim(sStrDomPix)//' [-] '//'TerrainHeigthMax: '//trim(sStrDemMax)//' [m] ')
        !------------------------------------------------------------------------------------ 
                
        !------------------------------------------------------------------------------------
        ! Debug
        if (iDEBUG.gt.0) then
            call mprintf(.true., iINFO_Extra, ' ========= STATIC GRIDDED START =========== ')   
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarDEM, a2iVarMask, 'DEM ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarLon, a2iVarMask, 'LON ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarLat, a2iVarMask, 'LAT ') )
            call mprintf(.true., iINFO_Extra, checkvar(float(a2iVarGlacierMask), a2iVarMask, 'GLACIER MASK ') )
            call mprintf(.true., iINFO_Extra, checkvar(float(a2iVarGlaciers_ID), a2iVarMask, 'GLACIER ID ') )  
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarGlacierDebris, a2iVarMask, 'GLACIER DEBRIS ') )
            if (oS3M_Namelist(iID)%iFlagThickFromTerrData .eq. 1) then
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarIceThick, a2iVarMask, 'GLACIER THICKNESS ') )
            endif
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarGlacierPivotTable, int(a2dVarGlacierPivotTable), 'PIVOT TABLE ') )  
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarAreaCell, a2iVarMask, 'AREA CELL ') )
            call mprintf(.true., iINFO_Extra, ' ========= STATIC GRIDDED END =========== ') 
        endif
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Pass local variable(s) to global workspace
        oS3M_Vars(iID)%iDomainPixels = iDomainPixels
        oS3M_Vars(iID)%dDEMMax = dDEMMax
        oS3M_Vars(iID)%dDEMMin = dDEMMin
        oS3M_Vars(iID)%dDomainArea = dDomainArea
        
        oS3M_Vars(iID)%a2iMask = a2iVarMask
        oS3M_Vars(iID)%a2iGlacierMask = a2iVarGlacierMask
        oS3M_Vars(iID)%a2iGlaciers_ID = a2iVarGlaciers_ID
        oS3M_Vars(iID)%a2dGlacierDebris = a2dVarGlacierDebris
        oS3M_Vars(iID)%a2dGlacierPivotTable = a2dVarGlacierPivotTable

        if (oS3M_Namelist(iID)%iFlagThickFromTerrData .eq. 1) then
            oS3M_Vars(iID)%a2dIceThick = a2dVarIceThick
        endif
        
        oS3M_Vars(iID)%a2dLon = a2dVarLon
        oS3M_Vars(iID)%a2dLat = a2dVarLat
        
        oS3M_Vars(iID)%a2dDem = a2dVarDEM
        oS3M_Vars(iID)%a2dAreaCell = a2dVarAreaCell
        !------------------------------------------------------------------------------------ 
                              
    end subroutine S3M_Data_Static_Gridded_Land
    !------------------------------------------------------------------------------------ 
   
end module S3M_Module_Data_Static_Gridded
!------------------------------------------------------------------------------------
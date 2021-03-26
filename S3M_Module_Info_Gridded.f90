!------------------------------------------------------------------------------------
! File:   S3M_Module_Info_Gridded.f90
! Author(s): Fabio Delogu, Francesco Silvestro, Simone Gabellani, Daniele Dolia, Francesco Avanzi.
!
! Created on May, 20 2014, 9:57 AM
! Last update on October 26, 2020 03:00 PM
!
! Module to get info gridded
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
! Module Header
module S3M_Module_Info_Gridded
    
    !------------------------------------------------------------------------------------
    ! External module(s) and implicit none
#ifdef LIB_NC
    use netcdf
#endif
   
    use S3M_Module_Namelist,        only:   oS3M_Namelist
    use S3M_Module_Vars_Loader,     only:   oS3M_Vars
   
    use S3M_Module_Tools_Debug
    use S3M_Module_Tools_Generic,   only:   S3M_Tools_Generic_ReplaceText, &
                                            S3M_Tools_Generic_CreateIndexGrid, &
                                            S3M_Tools_Generic_UnzipFile, &
                                            reshape2DVar, &
                                            filter_array_unique, &
                                            filter_array_condition
#ifdef LIB_NC
    use S3M_Module_Tools_IO,        only:   S3M_Tools_IO_Get2d_NC, &
                                            check
#endif
    
    implicit none
    !------------------------------------------------------------------------------------
    
contains 

    !------------------------------------------------------------------------------------
    ! Subroutine to get gridded-file dimension(s) -- Terrain Data
    subroutine S3M_Info_Gridded_GetDims_Ancillary(iID, iRows, iCols, iGlaciers_ID_DEF, a1iGlaciers_ID_DEF)

        !------------------------------------------------------------------------------------
        ! Variable(s) declaration
        integer(kind = 4)           :: iID
        integer(kind = 4)           :: iRows, iCols
    
        character(len = 256)        :: sDomainName, sPathData
        
        integer(kind = 4)           :: iFileID, iDimId
        
        character(len = 256)        :: sVarName
        character(len = 256)        :: sVarUnits
        character(len = 700)        :: sFileName, sFileNameZip
        character(len = 700)        :: sCommandUnzip
        
        logical                     :: bFileExist
        
        integer(kind = 4)           ::  iErr
        
        real(kind = 4),  dimension (iCols, iRows)       :: a2dVar_TMP
        real(kind = 4),  dimension (iRows, iCols)       :: a2dVar_DEF 
        real(kind = 4),  dimension (iRows * iCols)      :: a1dVar_DEF
        
        integer(kind = 4)                               :: iGlaciers_ID_DEF
        integer(kind = 4), dimension(:), allocatable    :: a1iGlaciers_ID_UNIQUE, a1iGlaciers_ID_DEF
        real(kind = 4), dimension(:), allocatable    :: a1dGlaciers_ID_DEF
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Initialize variable(s)
        sVarUnits = ''
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Get information
        sDomainName = oS3M_Namelist(iID)%sDomainName
        sPathData = oS3M_Namelist(iID)%sPathData_Static_Gridded
        sCommandUnzip = oS3M_Namelist(iID)%sCommandUnzipFile

        ! Get variable(s) dimension(s)
        call mprintf(.true., iINFO_Main, ' Define land geographical data  ... ')
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Get info from netCDF data
#ifdef LIB_NC
            !------------------------------------------------------------------------------------
            ! Info
            call mprintf(.true., iINFO_Extra, ' Land data in netCDF format ')
            !------------------------------------------------------------------------------------

            !------------------------------------------------------------------------------------
            ! Filename
            sFileName = trim(sPathData)//'Terrain_Data.nc'

            ! Checking file input availability
            sFileNameZip = sFileName(1:len_trim(sFileName))//'.gz'
            inquire (file = trim(sFileNameZip), exist = bFileExist)
            if ( .not. bFileExist ) then
                call mprintf(.true., iERROR, ' No compressed static file netCDF found '//trim(sFileNameZip) )
            endif

            ! Unzip file
            call S3M_Tools_Generic_UnzipFile(sCommandUnzip, sFileNameZip, sFileName, .true.)

            ! Check file availability
            inquire (file = trim(sFileName), exist = bFileExist)
            if ( .not. bFileExist ) then
                call mprintf(.true., iERROR, ' No static file netCDF found '//trim(sFileName) )
            else            
            
                ! Open netCDF file
                call check( nf90_open(trim(sFileName), NF90_NOWRITE, iFileID) )
               
                ! Glaciers ID (to estimate number of glaciers in GlacierID)
                sVarName = 'GlacierID'
                call S3M_Tools_IO_Get2d_NC(sVarName, iFileID, a2dVar_TMP, sVarUnits, iCols, iRows, .false., iErr)
                if (iErr /= 0) then
                    call mprintf(.true., iWARN, ' Glacier ID map not found. Initializing glacier ID map with undefined values')
                    a2dVar_DEF = -9999
                else
                    a2dVar_DEF = transpose(a2dVar_TMP)
                endif

                ! Close netCDF file
                call check( nf90_close(iFileID) )
                ! Info
                call mprintf(.true., iINFO_Main, ' Define land geographical data  ... OK')
            endif
            
            call filter_array_unique(int(reshape2DVar(a2dVar_DEF)), a1iGlaciers_ID_UNIQUE)
            call filter_array_condition(real(a1iGlaciers_ID_UNIQUE), a1dGlaciers_ID_DEF, 0.0, .false. , .false. , .true.)
            a1iGlaciers_ID_DEF = int(a1dGlaciers_ID_DEF)

            iGlaciers_ID_DEF = size(a1iGlaciers_ID_DEF)
            !------------------------------------------------------------------------------------
#endif
        !------------------------------------------------------------------------------------
        
    end subroutine S3M_Info_Gridded_GetDims_Ancillary
    !------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------
    ! Subroutine to get file gridded dimension(s)
    subroutine S3M_Info_Gridded_GetDims_Static(iID, iRows, iCols, iRows_Pivot, iCols_Pivot)
        
        !------------------------------------------------------------------------------------
        ! Variable(s) declaration
        integer(kind = 4)           :: iID
        integer(kind = 4)           :: iRows, iCols, iRows_Pivot, iCols_Pivot, iCheckPivot
        
        character(len = 256)        :: sDomainName, sPathData
        
        integer(kind = 4)           :: iNCid, iDimId
        character(len = 19)         :: sTime
        
        character(len = 256)        :: sText
        character(len = 700)        :: sFileName, sFileNameZip
        character(len = 700)        :: sCommandUnzip
        
        logical                     :: bFileExist
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Initialize variable(s)
        iRows = 0; iCols = 0; iRows_Pivot = 0; iCols_Pivot = 0;
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Get information
        sDomainName = oS3M_Namelist(iID)%sDomainName
        sPathData = oS3M_Namelist(iID)%sPathData_Static_Gridded
        sTime = oS3M_Namelist(iID)%sTimeStart
        sCommandUnzip = oS3M_Namelist(iID)%sCommandUnzipFile
        iCheckPivot = oS3M_Namelist(iID)%iFlagIceMassBalance
        
        ! Get variable(s) dimension(s)
        call mprintf(.true., iINFO_Main, ' Define land data dims ... ')
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Get info from netCDF file
#ifdef LIB_NC
            !------------------------------------------------------------------------------------
            ! Filename 
            sFileName = trim(sPathData)//'Terrain_Data.nc'
            
            ! Checking file input availability
            sFileNameZip = sFileName(1:len_trim(sFileName))//'.gz'
            inquire (file = trim(sFileNameZip), exist = bFileExist)
            if ( .not. bFileExist ) then
                call mprintf(.true., iERROR, ' No compressed static file netCDF found '//trim(sFileNameZip) )
            endif

            ! Unzip file
            call S3M_Tools_Generic_UnzipFile(sCommandUnzip, sFileNameZip, sFileName, .true.)

            ! Check file availability
            inquire (file = trim(sFileName), exist = bFileExist)
            if ( .not. bFileExist ) then
                call mprintf(.true., iERROR, ' No static file netCDF found '//trim(sFileName) )
            else
                
                ! Info
                call mprintf(.true., iINFO_Extra, ' Land data in netCDF format ')
                ! Open netCDF file
                call check( nf90_open(trim(sFileName), NF90_NOWRITE, iNCid) )
                ! Get dimension(s)
                call check( nf90_inq_dimid(iNCid, "X", iDimId) )
                call check( nf90_inquire_dimension(iNCid, iDimId, len = iCols) )
                call check( nf90_inq_dimid(iNCid, "Y", iDimId) )
                call check( nf90_inquire_dimension(iNCid, iDimId, len = iRows) )
                if (iCheckPivot .eq. 2) then !that is, if we are using a glacier-MB mode that requires the pivot table...
                    call check( nf90_inq_dimid(iNCid, "COL_PIVOT", iDimId) )
                    call check( nf90_inquire_dimension(iNCid, iDimId, len = iCols_Pivot) )
                    call check( nf90_inq_dimid(iNCid, "ROW_PIVOT", iDimId) )
                    call check( nf90_inquire_dimension(iNCid, iDimId, len = iRows_Pivot) )                    
                endif                
                ! Close netCDF file
                call check( nf90_close(iNCid) )
                ! Info
                call mprintf(.true., iINFO_Main, ' Define land data dims ... OK')
 
            endif
            !------------------------------------------------------------------------------------ 
#endif
        !------------------------------------------------------------------------------------
        ! Dims from local to global workspace
        oS3M_Namelist(iID)%iRowsL = iRows
        oS3M_Namelist(iID)%iColsL = iCols
        oS3M_Namelist(iID)%iColsPivot = iCols_Pivot  
        oS3M_Namelist(iID)%iRowsPivot = iRows_Pivot  
        !------------------------------------------------------------------------------------
        
    end subroutine S3M_Info_Gridded_GetDims_Static
    !------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------
    ! Subroutine to get file forcing dimension(s)
    subroutine S3M_Info_Gridded_GetDims_Forcing(iID, iRows, iCols)
        
        !------------------------------------------------------------------------------------
        ! Variable(s) declaration
        integer(kind = 4)           :: iID
        integer(kind = 4)           :: iRows, iCols
        
        character(len = 19)         :: sTime
        character(len = 256)        :: sPathData
        character(len = 700)        :: sFileName, sFileNameZip
        character(len = 700)        :: sCommandUnzip
        character(len = 256)        :: sDomainName
        
        integer(kind = 4)           :: iFileID, iDimId

        logical                     :: bFileExist
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Initialize variable(s)
        sFileName = ""; sFileNameZip = "";
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Get information
        sDomainName = oS3M_Namelist(iID)%sDomainName
        sPathData   = oS3M_Namelist(iID)%sPathData_Forcing_Gridded

        ! Get variable(s) dimension(s)
        call mprintf(.true., iINFO_Main, ' Define forcing data dims ... ')
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Get info from netCDF file
#ifdef LIB_NC
            !------------------------------------------------------------------------------------
            ! Info
            call mprintf(.true., iINFO_Extra, ' Forcing data in netCDF format ')
            !------------------------------------------------------------------------------------
            
            !------------------------------------------------------------------------------------
            ! Get information
            sTime = oS3M_Namelist(iID)%sTimeStart
            sCommandUnzip = oS3M_Namelist(iID)%sCommandUnzipFile
            
            ! Filename netCDF
            sFileName = trim(sPathData)//"MeteoData_"// &
                sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
                sTime(12:13)//sTime(15:16)// &
                ".nc"

            ! Replace general path with specific time feature(s)
            call S3M_Tools_Generic_ReplaceText(sFileName, '$yyyy', sTime(1:4))
            call S3M_Tools_Generic_ReplaceText(sFileName, '$mm', sTime(6:7))
            call S3M_Tools_Generic_ReplaceText(sFileName, '$dd', sTime(9:10))

            ! Checking file input availability
            sFileNameZip = sFileName(1:len_trim(sFileName))//'.gz'
            inquire (file = trim(sFileNameZip), exist = bFileExist)
            if ( .not. bFileExist ) then
                call mprintf(.true., iERROR, ' No compressed forcing file netCDF found '//trim(sFileNameZip) )
            endif

            ! Unzip file
            call S3M_Tools_Generic_UnzipFile(sCommandUnzip, sFileNameZip, sFileName, .true.)

            ! Check file availability
            inquire (file = trim(sFileName), exist = bFileExist)
            if ( .not. bFileExist ) then
                call mprintf(.true., iERROR, ' No forcing file netCDF found '//trim(sFileName) )
            else

                ! Open netCDF file
                call check( nf90_open(trim(sFileName), NF90_NOWRITE, iFileID) )
                ! Get global attribute(s)
                call check( nf90_get_att(iFileID, nf90_global, "ncols", iCols) )
                call check( nf90_get_att(iFileID, nf90_global, "nrows", iRows) )
                ! Close netCDF file
                call check( nf90_close(iFileID) )
                ! Info
                call mprintf(.true., iINFO_Main, ' Define forcing data dims ... OK')
                
            endif

#endif
            
        !------------------------------------------------------------------------------------
        ! Dims from local to global workspace
        oS3M_Namelist(iID)%iRowsF = iRows
        oS3M_Namelist(iID)%iColsF = iCols
        !------------------------------------------------------------------------------------
        
    end subroutine S3M_Info_Gridded_GetDims_Forcing
    !------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------
    ! Subroutine to get geographical information
    subroutine S3M_Info_Gridded_GetGeo_Static(iID)
    
        !------------------------------------------------------------------------------------
        ! Variable(s) declaration
        integer(kind = 4)           :: iID
        
        character(len = 256)        :: sDomainName, sPathData
        
        integer(kind = 4)           :: iNCid, iDimId
        
        character(len = 256)        :: sText
        character(len = 700)        :: sFileName, sFileNameZip
        character(len = 700)        :: sCommandUnzip
        
        logical                     :: bFileExist
        
        integer(kind = 4)           :: iVar
        real(kind = 4)              :: dVar
        real(kind = 4)              :: dVarXLLCorner, dVarYLLCorner, dVarCellSize, dVarNoData
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Initialize variable(s)
        iVar = 0; dVar = 0;
        dVarXLLCorner = -9999.0; dVarYLLCorner = -9999.0; dVarCellSize = -9999.0; 
        dVarNoData = -9999.0; 
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Get information
        sDomainName = oS3M_Namelist(iID)%sDomainName
        sPathData = oS3M_Namelist(iID)%sPathData_Static_Gridded
        sCommandUnzip = oS3M_Namelist(iID)%sCommandUnzipFile

        ! Get variable(s) dimension(s)
        call mprintf(.true., iINFO_Main, ' Define land geographical data  ... ')
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Get info from netCDF data
#ifdef LIB_NC
            !------------------------------------------------------------------------------------
            ! Info
            call mprintf(.true., iINFO_Extra, ' Land data in netCDF format ')
            !------------------------------------------------------------------------------------

            !------------------------------------------------------------------------------------
            ! Filename
            sFileName = trim(sPathData)//'Terrain_Data.nc'

            ! Checking file input availability
            sFileNameZip = sFileName(1:len_trim(sFileName))//'.gz'
            inquire (file = trim(sFileNameZip), exist = bFileExist)
            if ( .not. bFileExist ) then
                call mprintf(.true., iERROR, ' No compressed static file netCDF found '//trim(sFileNameZip) )
            endif

            ! Unzip file
            call S3M_Tools_Generic_UnzipFile(sCommandUnzip, sFileNameZip, sFileName, .true.)

            ! Check file availability
            inquire (file = trim(sFileName), exist = bFileExist)
            if ( .not. bFileExist ) then
                call mprintf(.true., iERROR, ' No static file netCDF found '//trim(sFileName) )
            else            
            
                ! Open netCDF file
                call check( nf90_open(trim(sFileName), NF90_NOWRITE, iNCid) )
                ! Get global attribute(s)
                call check( nf90_get_att(iNCid, nf90_global, "xllcorner",    dVarXLLCorner) )
                call check( nf90_get_att(iNCid, nf90_global, "yllcorner",    dVarYLLCorner) )
                call check( nf90_get_att(iNCid, nf90_global, "cellsize",     dVarCellSize) )
                call check( nf90_get_att(iNCid, nf90_global, "nodata_value", dVarNoData) )
                ! Close netCDF file
                call check( nf90_close(iNCid) )
                ! Info
                call mprintf(.true., iINFO_Main, ' Define land geographical data  ... OK')
            endif
            !------------------------------------------------------------------------------------

#endif
        !------------------------------------------------------------------------------------
        ! Pass local variable(s) to global workspace
        oS3M_Namelist(iID)%dXLLCornerL = dVarXLLCorner
        oS3M_Namelist(iID)%dYLLCornerL = dVarYLLCorner
        oS3M_Namelist(iID)%dXCellSizeL = dVarCellSize
        oS3M_Namelist(iID)%dYCellSizeL = dVarCellSize
        oS3M_Namelist(iID)%dNoDataL = dVarNoData
        !------------------------------------------------------------------------------------
        
    end subroutine S3M_Info_Gridded_GetGeo_Static
    !------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------
    ! Subroutine to create meteo grid (if static and forcing data have different grid)
    subroutine S3M_Info_Gridded_GetGeo_Forcing(iID)
        
        !------------------------------------------------------------------------------------
        ! Variable(s) declaration
        integer(kind = 4)           :: iID

        character(len = 19)         :: sTime
        character(len = 256)        :: sPathData
        character(len = 700)        :: sFileName, sFileNameZip
        character(len = 700)        :: sCommandUnzip
        
        integer(kind = 4)           :: iVarCols, iVarRows
        real(kind = 4)              :: dVarXLLCorner, dVarYLLCorner
        real(kind = 4)              :: dVarXCellSize, dVarYCellSize
        real(kind = 4)              :: dVarNoData
        
        integer(kind = 4)           :: iErr, iFileID
        
        logical                     :: bFileExist
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Initialize variable(s)
        iVarCols = -9999; iVarRows = -9999; 
        dVarXLLCorner = -9999.0; dVarYLLCorner = -9999.0;
        dVarXCellSize = -9999.0; dVarYCellSize = -9999.0;
        dVarNoData = -9999.0;
        
        sFileName = ""; sFileNameZip = "";
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Check if grid control is defined
        if (.not. oS3M_Namelist(iID)%bGridCheck) then
            
            !------------------------------------------------------------------------------------
            ! Get information
            sPathData = oS3M_Namelist(iID)%sPathData_Forcing_Gridded

            ! Data type
            call mprintf(.true., iINFO_Main, ' Define forcing geographical data  ... ')
            !------------------------------------------------------------------------------------
            
            !------------------------------------------------------------------------------------
            ! Get info from netCDF data
#ifdef LIB_NC
                !------------------------------------------------------------------------------------
                ! Info
                call mprintf(.true., iINFO_Extra, ' Forcing data in netCDF format ')
                !------------------------------------------------------------------------------------
                
                !------------------------------------------------------------------------------------ 
                ! Get information
                sTime = oS3M_Namelist(iID)%sTimeStart
                sCommandUnzip = oS3M_Namelist(iID)%sCommandUnzipFile
                !------------------------------------------------------------------------------------ 
                
                !------------------------------------------------------------------------------------   
                ! Filename netCDF
                sFileName = trim(sPathData)//"MeteoData_"// &
                    sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
                    sTime(12:13)//sTime(15:16)// &
                    ".nc"
                !------------------------------------------------------------------------------------
                
                !------------------------------------------------------------------------------------------
                ! Replace general path with specific time feature(s)
                call S3M_Tools_Generic_ReplaceText(sFileName, '$yyyy', sTime(1:4))
                call S3M_Tools_Generic_ReplaceText(sFileName, '$mm', sTime(6:7))
                call S3M_Tools_Generic_ReplaceText(sFileName, '$dd', sTime(9:10))
                !------------------------------------------------------------------------------------------
                    
                !------------------------------------------------------------------------------------------
                ! Checking file input availability
                sFileNameZip = trim(sFileName)//'.gz'
                inquire (file = trim(sFileNameZip), exist = bFileExist)
                if ( .not. bFileExist ) then
                    call mprintf(.true., iERROR, ' No compressed forcing file netCDF found '//trim(sFileNameZip) )
                endif
                !------------------------------------------------------------------------------------------

                !------------------------------------------------------------------------------------------
                ! Unzip file
                call S3M_Tools_Generic_UnzipFile(sCommandUnzip, sFileNameZip, sFileName, .true.)
                !------------------------------------------------------------------------------------------
                
                !------------------------------------------------------------------------------------------
                ! Check uncompressed file availability
                inquire (file = trim(sFileName), exist = bFileExist)
                if ( .not. bFileExist ) then
                    call mprintf(.true., iERROR, ' No forcing file netCDF found '//trim(sFileName) )
                else
                    
                    ! Open nc file
                    call check( nf90_open(trim(sFileName), NF90_NOWRITE, iFileID) )
                    
                    ! Get global attribute(s)
                    call check( nf90_get_att(iFileID, nf90_global, "xllcorner",    dVarXLLCorner) )
                    call check( nf90_get_att(iFileID, nf90_global, "yllcorner",    dVarYLLCorner) )
                    call check( nf90_get_att(iFileID, nf90_global, "cellsize",     dVarXCellSize) )
                    call check( nf90_get_att(iFileID, nf90_global, "nodata_value", dVarNoData) )
                    
                    dVarYCellSize = dVarXCellSize
                    
                    ! Close nc file
                    call check( nf90_close(iFileID) )
                    
                    ! Info
                    call mprintf(.true., iINFO_Main, ' Define forcing geographical data  ... OK ')

                endif
                !------------------------------------------------------------------------------------------
#endif
            !------------------------------------------------------------------------------------------
            !------------------------------------------------------------------------------------
            ! Pass local variable(s) to global workspace
            oS3M_Namelist(iID)%dXLLCornerF = dVarXLLCorner
            oS3M_Namelist(iID)%dYLLCornerF = dVarYLLCorner
            oS3M_Namelist(iID)%dXCellSizeF = dVarXCellSize
            oS3M_Namelist(iID)%dYCellSizeF = dVarYCellSize
            oS3M_Namelist(iID)%dNoDataF = dVarNoData
            !------------------------------------------------------------------------------------
            
            !------------------------------------------------------------------------------------
            ! Compute indexes grid(s) ---> if grid are different XIndex and YIndex are not equal to -9999
            call S3M_Tools_Generic_CreateIndexGrid(oS3M_Namelist(iID)%iRowsL, oS3M_Namelist(iID)%iColsL, &
                                                   oS3M_Namelist(iID)%iRowsF, oS3M_Namelist(iID)%iColsF, &
                                                   oS3M_Namelist(iID)%dYLLCornerL, oS3M_Namelist(iID)%dXLLCornerL, &
                                                   oS3M_Namelist(iID)%dYLLCornerF, oS3M_Namelist(iID)%dXLLCornerF, &
                                                   oS3M_Namelist(iID)%dYCellSizeL, oS3M_Namelist(iID)%dXCellSizeL, &
                                                   oS3M_Namelist(iID)%dYCellSizeF, oS3M_Namelist(iID)%dXCellSizeF, &
                                                   oS3M_Namelist(iID)%iFlagGrid, &
                                                   oS3M_Vars(iID)%a2iXIndex, oS3M_Vars(iID)%a2iYIndex)
            !------------------------------------------------------------------------------------
                               
            !------------------------------------------------------------------------------------
            ! Flag grid control
            oS3M_Namelist(iID)%bGridCheck = .true.
            !------------------------------------------------------------------------------------
            
        endif
        !------------------------------------------------------------------------------------
        
    end subroutine S3M_Info_Gridded_GetGeo_Forcing
    !--------------------------------------------------------------------------------
    
end module S3M_Module_Info_Gridded
!--------------------------------------------------------------------------------

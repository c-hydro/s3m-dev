!------------------------------------------------------------------------------------------    
! File:   S3M_Module_Data_Forcing_Gridded.f90
! Author(s): Fabio Delogu, Francesco Silvestro, Simone Gabellani, Francesco Avanzi
!
! Created on April 22, 2015, 5:19 PM
! Last update on Dec 15, 2023 09:30 AM
!
! Module to read forcing map.
!------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------
! Module Header
module S3M_Module_Data_Forcing_Gridded
    
    !------------------------------------------------------------------------------------------
    ! External module(s) for all subroutine in this module
#ifdef LIB_NC
    use netcdf
#endif

    use S3M_Module_Namelist,        only:   oS3M_Namelist
    use S3M_Module_Vars_Loader,     only:   oS3M_Vars
    
    use S3M_Module_Tools_Debug

#ifdef LIB_NC
    use S3M_Module_Tools_IO,        only:   S3M_Tools_IO_Get2d_Binary_INT, &
                                            S3M_Tools_IO_Get2d_Binary_DBL, &
                                            S3M_Tools_IO_Get2d_NC, &
                                            S3M_Tools_IO_CheckVar_NC, &
                                            check
#else
    use S3M_Module_Tools_IO,        only:   S3M_Tools_IO_Get2d_Binary_INT, &
                                            S3M_Tools_IO_Get2d_Binary_DBL
#endif                                    
                                                                                  
    use S3M_Module_Tools_Generic,   only:   S3M_Tools_Generic_ReplaceText, &
                                            S3M_Tools_Generic_SwitchGrid, &
                                            S3M_Tools_Generic_UnzipFile, &
                                            S3M_Tools_Generic_RemoveFile, &
                                            check2Dvar, getProcessID
                                            
    use S3M_Module_Tools_Time,      only:   S3M_Tools_Time_MonthVal
                             
    ! Implicit none for all subroutines in this module
    implicit none
    !------------------------------------------------------------------------------------------
    
contains
    
    !------------------------------------------------------------------------------------------
    ! Subroutine to read and manage forcing gridded data
    subroutine S3M_Data_Forcing_Gridded_Cpl( iID, sTime, &
                                     iRowsStartL, iRowsEndL, iColsStartL, iColsEndL, &
                                     iRowsStartF, iRowsEndF, iColsStartF, iColsEndF)
        
        !------------------------------------------------------------------------------------------
        ! Variable(s)
        integer(kind = 4)           :: iID
                                    
        integer(kind = 4)           :: iRowsStartL, iRowsEndL, iColsStartL, iColsEndL
        integer(kind = 4)           :: iRowsStartF, iRowsEndF, iColsStartF, iColsEndF
        integer(kind = 4)           :: iRowsL, iColsL, iRowsF, iColsF
        integer(kind = 4)           :: iFlagTypeData_Forcing
        integer(kind = 4)           :: iScaleFactor_Forcing
        
        character(len = 19)         :: sTime
        character(len = 12)         :: sTimeMonth
        
        character(len = 256)        :: sPathData_Forcing
        
        real(kind = 4), dimension(iRowsEndL - iRowsStartL + 1, &
                                  iColsEndL - iColsStartL + 1) ::   a2dVarDem
                                  
        real(kind = 4), dimension(iRowsEndL - iRowsStartL + 1, &
                                  iColsEndL - iColsStartL + 1) ::   a2dVarPrecipL, a2dVarTaL, &
                                                                    a2dVarIncRadL, & 
                                                                    a2dVarRelHumL
                                                                    
                                                                    
        real(kind = 4), dimension(iRowsEndF - iRowsStartF + 1, &
                                  iColsEndF - iColsStartF + 1) ::   a2dVarPrecipF, a2dVarTaF, &
                                                                    a2dVarIncRadF, & 
                                                                    a2dVarRelHumF
                                                                    
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Forcing data (mandatory):
        !   a2dVarPrecipF       : total precip [mm]
        !   a2dVarTaF           : air temperature [degC]
        !   a2dVarIncRadF       : incoming radiation [W/m^2]
        !   a2dVarRelHumF       : relative humidity [%]
        
        !------------------------------------------------------------------------------------------
                                                                                                                        
        !------------------------------------------------------------------------------------------
        ! Initialize variable(s)
        a2dVarPrecipF = -9999.0; a2dVarTaF = -9999.0; a2dVarIncRadF = -9999.0;  
        a2dVarRelHumF = -9999.0;
        
        a2dVarPrecipL = -9999.0; a2dVarTaL = -9999.0; a2dVarIncRadL = -9999.0;
        a2dVarRelHumL = -9999.0;
        !------------------------------------------------------------------------------------------
                                                                                                
        !------------------------------------------------------------------------------------------
        ! Defining iRows and iCols (Land data)
        iRowsL = iRowsEndL - iRowsStartL + 1
        iColsL = iColsEndL - iColsStartL + 1
        ! Defining iRows and iCols (forcing data)
        iRowsF = iRowsEndF - iRowsStartF + 1
        iColsF = iColsEndF - iColsStartF + 1
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Get global information
        sPathData_Forcing = oS3M_Namelist(iID)%sPathData_Forcing_Gridded
        iFlagTypeData_Forcing = oS3M_Namelist(iID)%iFlagTypeData_Forcing_Gridded
        iScaleFactor_Forcing = oS3M_Namelist(iID)%iScaleFactor_Forcing
        a2dVarDem = oS3M_Vars(iID)%a2dDem
        
        ! Info start
        call mprintf(.true., iINFO_Extra, ' Data :: Forcing gridded ... ' )
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Replace general path with specific time feature(s)
        call S3M_Tools_Generic_ReplaceText(sPathData_Forcing, '$yyyy', sTime(1:4))
        call S3M_Tools_Generic_ReplaceText(sPathData_Forcing, '$mm', sTime(6:7))
        call S3M_Tools_Generic_ReplaceText(sPathData_Forcing, '$dd', sTime(9:10))
        call S3M_Tools_Generic_ReplaceText(sPathData_Forcing, '$HH', sTime(12:13))
        call S3M_Tools_Generic_ReplaceText(sPathData_Forcing, '$MM', sTime(15:16))
        
        ! Checking date
        write(sTimeMonth,'(A,A,A)') sTime(1:4), sTime(6:7), sTime(9:10)
        !------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        ! Check time step (iT)
        if (oS3M_Vars(iID)%iTime .lt. oS3M_Namelist(iID)%iNTime) then
            
            !------------------------------------------------------------------------------------------
            ! Subroutine for reading sequential netCDF forcing data 
            if (iFlagTypeData_Forcing == 3) then

                !------------------------------------------------------------------------------------------
                ! Call subroutine to get forcing data in netCDF format
#ifdef LIB_NC
                call S3M_Data_Forcing_Gridded_NC(iID, &
                                        sPathData_Forcing, &
                                        iRowsF, iColsF, &
                                        sTime, &
                                        a2dVarPrecipF, a2dVarTaF, a2dVarIncRadF, &
                                        a2dVarRelHumF)
#else   
                ! Redefinition of forcing data flag (if netCDF library is not linked)
                iFlagTypeData_Forcing = 1 
                call mprintf(.true., iWARN, ' '// &
                                            'Forcing gridded data type selected was netCDF but library is not linked! '// &
                                            'Data in binary int format will be used!')
#endif
                !------------------------------------------------------------------------------------------
                           
            endif
            !------------------------------------------------------------------------------------------

            !------------------------------------------------------------------------------------------
            ! Subroutine for reading sequential binary forcing data
            if (iFlagTypeData_Forcing == 1 .or. iFlagTypeData_Forcing == 2) then
                
                !------------------------------------------------------------------------------------------
                ! Calling subroutine to read data in binary format
                call S3M_Data_Forcing_Gridded_Binary(iID, iFlagTypeData_Forcing, &
                                            sPathData_Forcing, &
                                            iRowsF, iColsF, &
                                            sTime, &
                                            a2dVarPrecipF, a2dVarTaF, a2dVarIncRadF, &
                                            a2dVarRelHumF, &
                                            iScaleFactor_Forcing)
                !------------------------------------------------------------------------------------------
                                        
            endif
            !------------------------------------------------------------------------------------------

            !------------------------------------------------------------------------------------------
            ! Debug
            if (iDEBUG.gt.0) then
                call mprintf(.true., iINFO_Extra, ' ========= FORCING GRIDDED START =========== ')
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarPrecipF, oS3M_Vars(iID)%a2iMask, 'PRECIP START') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarTaF, oS3M_Vars(iID)%a2iMask, 'TA START') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarIncRadF, oS3M_Vars(iID)%a2iMask, 'INCRAD START') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarRelHumF, oS3M_Vars(iID)%a2iMask, 'RELHUM START') )
		call mprintf(.true., iINFO_Extra, '')
            endif
            !------------------------------------------------------------------------------------------
            
            !------------------------------------------------------------------------------------------
            ! Grid switcher land-forcing
            call S3M_Tools_Generic_SwitchGrid(oS3M_Namelist(iID)%iFlagGrid, &
                                              iRowsL, iColsL, a2dVarPrecipL, &
                                              iRowsF, iColsF, a2dVarPrecipF, &
                                              oS3M_Vars(iID)%a2dDem, &
                                              oS3M_Vars(iID)%a2iXIndex, oS3M_Vars(iID)%a2iYIndex)
            call S3M_Tools_Generic_SwitchGrid(oS3M_Namelist(iID)%iFlagGrid, &
                                              iRowsL, iColsL, a2dVarTaL, &
                                              iRowsF, iColsF, a2dVarTaF, &
                                              oS3M_Vars(iID)%a2dDem, &
                                              oS3M_Vars(iID)%a2iXIndex, oS3M_Vars(iID)%a2iYIndex)
            call S3M_Tools_Generic_SwitchGrid(oS3M_Namelist(iID)%iFlagGrid, &
                                              iRowsL, iColsL, a2dVarIncRadL, &
                                              iRowsF, iColsF, a2dVarIncRadF, &
                                              oS3M_Vars(iID)%a2dDem, &
                                              oS3M_Vars(iID)%a2iXIndex, oS3M_Vars(iID)%a2iYIndex)
            call S3M_Tools_Generic_SwitchGrid(oS3M_Namelist(iID)%iFlagGrid, &
                                              iRowsL, iColsL, a2dVarRelHumL, &
                                              iRowsF, iColsF, a2dVarRelHumF, &
                                              oS3M_Vars(iID)%a2dDem, &
                                              oS3M_Vars(iID)%a2iXIndex, oS3M_Vars(iID)%a2iYIndex)                                               
            !------------------------------------------------------------------------------------------
                                            
            !------------------------------------------------------------------------------------------
            ! Check variable(s) limits and domain
            a2dVarPrecipL = check2Dvar(a2dVarPrecipL,           oS3M_Vars(iID)%a2iMask,     0.0,    850.0,  0.0)
            a2dVarTaL = check2Dvar(a2dVarTaL,                   oS3M_Vars(iID)%a2iMask,     -70.0,  60.0,   0.0)    
            a2dVarIncRadL = check2Dvar(a2dVarIncRadL,           oS3M_Vars(iID)%a2iMask,     0.0,    1412.0, 0.0) 
            a2dVarRelHumL = check2Dvar(a2dVarRelHumL,           oS3M_Vars(iID)%a2iMask,     0.0,    100.0,  0.0)
            !------------------------------------------------------------------------------------------
            
        else
            !------------------------------------------------------------------------------------------
            ! Extra steps condition
            a2dVarPrecipL = 0.0;
            a2dVarTaL = oS3M_Vars(iID)%a2dTa; 
            a2dVarIncRadL = oS3M_Vars(iID)%a2dK;  
            a2dVarRelHumL = oS3M_Vars(iID)%a2dRHum
            
            ! Info message for extra time step(s)
            call mprintf(.true., iINFO_Extra, ' Extra time step ---> Forcing data are set constant to last observed value')
            call mprintf(.true., iINFO_Extra, ' Extra time step ---> Precip data are set to 0.0')
            !------------------------------------------------------------------------------------------
        endif
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Nullify map boundary
        where(a2dVarDem.lt.0.0)
            a2dVarPrecipL = -9999
            a2dVarTaL = -9999
            a2dVarIncRadL = -9999
            a2dVarRelHumL = -9999
        endwhere
        !------------------------------------------------------------------------------------------        
        
        !------------------------------------------------------------------------------------------
        ! Check update and save forcing data to local variable(s) to global workspace
        ! Precip
        if ( .not. all(a2dVarPrecipL.eq.-9999.0) ) then
            oS3M_Vars(iID)%a2dPrecip = a2dVarPrecipL
        else
            call mprintf(.true., iWARN, ' All Precip values are undefined! Check forcing data!' )
        endif
        ! Air temperature
        if ( .not. all(a2dVarTaL.eq.-9999.0) ) then
            oS3M_Vars(iID)%a2dTa = a2dVarTaL
        else
            call mprintf(.true., iWARN, ' All air temperature values are undefined! Check forcing data!' )
        endif
        
        ! Incoming radiation
        if ( .not. all(a2dVarIncRadL.eq.-9999.0) ) then
            oS3M_Vars(iID)%a2dK = a2dVarIncRadL
        else
            call mprintf(.true., iWARN, ' All incoming radiation values are undefined! Check forcing data!' )
        endif
        
        ! Relative humidity
        if ( .not. all(a2dVarRelHumL.eq.-9999.0) ) then
            oS3M_Vars(iID)%a2dRHum = a2dVarRelHumL
        else
            call mprintf(.true., iWARN, ' All relative humidity values are undefined! Check forcing data!' )
        endif
     
        !------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        ! Debug
        if (iDEBUG.gt.0) then
            call mprintf(.true., iINFO_Extra, '')
            call mprintf(.true., iINFO_Extra, checkvar(oS3M_Vars(iID)%a2dPrecip, oS3M_Vars(iID)%a2iMask, 'PRECIP END') )
            call mprintf(.true., iINFO_Extra, checkvar(oS3M_Vars(iID)%a2dTa, oS3M_Vars(iID)%a2iMask, 'TA END') )
            call mprintf(.true., iINFO_Extra, checkvar(oS3M_Vars(iID)%a2dK, oS3M_Vars(iID)%a2iMask, 'INCRAD END') )
            call mprintf(.true., iINFO_Extra, checkvar(oS3M_Vars(iID)%a2dRHum, oS3M_Vars(iID)%a2iMask, 'RELHUM END') )
            call mprintf(.true., iINFO_Extra, ' ========= FORCING GRIDDED END =========== ')
        endif
        
        ! Info end
        call mprintf(.true., iINFO_Extra, ' Data :: Forcing gridded ... OK' )
        !------------------------------------------------------------------------------------------
        
    end subroutine S3M_Data_Forcing_Gridded_Cpl
    !------------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------------
    ! Subroutine to read NC data forcing
#ifdef LIB_NC
    subroutine S3M_Data_Forcing_Gridded_NC(iID,  &
                                  sPathData_Forcing, &
                                  iRows, iCols, sTime, &
                                  a2dVarPrecip, a2dVarTa, a2dVarIncRad, &
                                  a2dVarRelHum)
                                  
        !------------------------------------------------------------------------------------------
        ! Variable(s)
        integer(kind = 4)                       :: iID                  
        
        character(len = 256), intent(in)        :: sPathData_Forcing
        character(len = 700)                    :: sFileNameData_Forcing, sFileNameData_Forcing_Zip, sFileNameData_Temp
        character(len = 700)                    :: sCommandUnzipFile, sCommandRemoveFile
        character(len = 256)                    :: sVarName
        integer(kind = 4), intent(in)           :: iRows, iCols

        character(len = 19), intent(in)         :: sTime
        character(len = 12)                     :: sTimeMonth
        
        
        real(kind = 4), dimension(iCols, iRows)                 :: a2dVar
        
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarPrecip
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarTa
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarIncRad
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarRelHum       

        character(len = 256):: sVarUnits, sPID
        integer(kind = 4)   :: iErr
        integer(kind = 4)   :: iFileID
        
        logical             :: bFileExist
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Initialize variable(s)
        a2dVarPrecip = -9999.0; a2dVarTa = -9999.0; a2dVarIncRad = -9999.0; a2dVarRelHum = -9999.0; 

        sFileNameData_Forcing = ''; sFileNameData_Forcing_Zip = ''; sTimeMonth = ''
        sFileNameData_Temp = ''; sPid = ''
        
        ! Checking date
        write(sTimeMonth,'(A,A,A)') sTime(1:4), sTime(6:7), sTime(9:10)
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Get global information
        sCommandUnzipFile = oS3M_Namelist(iID)%sCommandUnzipFile
        sCommandRemoveFile = oS3M_Namelist(iID)%sCommandRemoveFile
        
        ! Info start
        call mprintf(.true., iINFO_Extra, ' Data :: Forcing gridded :: NetCDF ... ' )
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Get unique process ID
        sPID = adjustl(getProcessID())
        
        ! Filename forcing (example: MeteoData_201404300000.nc.gz)
        sFileNameData_Forcing = trim(sPathData_Forcing)//"MeteoData_"// &
        sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
        sTime(12:13)//sTime(15:16)// &
        ".nc"
        
        ! Create Filename with unique PID number to avoid simultaneously access to the same Forcing file       
        sFileNameData_Temp = trim(sPathData_Forcing)//"MeteoData_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)//'_'//trim(sPID)// &
            ".nc"  
        
        ! Info netCDF filename
        call mprintf(.true., iINFO_Verbose, ' Get filename (forcing gridded): '//trim(sFileNameData_Forcing)//' ... ' )
        !------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        ! Checking file input availability
        sFileNameData_Forcing_Zip = trim(sFileNameData_Forcing)//'.gz'
        inquire (file = trim(sFileNameData_Forcing)//'.gz', exist = bFileExist)
        if ( .not. bFileExist ) then
            !------------------------------------------------------------------------------------------
            ! Warning message
            call mprintf(.true., iWARN, ' No compressed forcing netCDF data found: '//trim(sFileNameData_Forcing_Zip) )
            ! Info netCDF filename
            call mprintf(.true., iINFO_Verbose, &
                         ' Get filename (forcing gridded): '//trim(sFileNameData_Forcing)//' ... FAILED' )
            ! Info end
            call mprintf(.true., iINFO_Extra, ' Data :: Forcing gridded :: NetCDF ... SKIPPED!' )
            !------------------------------------------------------------------------------------------
        else
            
            !------------------------------------------------------------------------------------------
            ! Unzip file
            call S3M_Tools_Generic_UnzipFile(oS3M_Namelist(iID)%sCommandUnzipFile, &
                                             sFileNameData_Forcing_Zip, &
                                             sFileNameData_Temp, .true.)
            !------------------------------------------------------------------------------------------

            !------------------------------------------------------------------------------------------
            ! Open netCDF file
            iErr = nf90_open(trim(sFileNameData_Temp), NF90_NOWRITE, iFileID)
            if (iErr /= 0) then
                call mprintf(.true., iWARN, ' Problem opening uncompressed netCDF file: '// &
                             trim(sFileNameData_Forcing)//' --> Undefined forcing data values' )
                call mprintf(.true., iINFO_Verbose, &
                             ' Get filename (forcing gridded): '//trim(sFileNameData_Forcing)//' ... FAILED' )
            else
                
                
                !------------------------------------------------------------------------------------------
                ! PRECIP
                call S3M_Tools_IO_CheckVar_NC('Precipitation;Rain', iFileID, sVarName) !We allow S3M to look for two possible 
                !sVarName for precipitation in the forcing file, mainly for compatibility issues. 
!                sVarName = 'Rain'
                call S3M_Tools_IO_Get2d_NC(sVarName, iFileID, a2dVar, sVarUnits, iCols, iRows, .true., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, ' Get forcing gridded data FAILED! Check forcing data for '//sVarName//'!')
                    a2dVarPrecip = -9999.0;
                else
                    a2dVarPrecip = transpose(a2dVar)
                endif
                !------------------------------------------------------------------------------------------
                
                !------------------------------------------------------------------------------------------
                ! AIR TEMPERATURE
                sVarName = 'AirTemperature'
                call S3M_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, ' Get forcing gridded data FAILED! Check forcing data for '//sVarName//'!')
                    a2dVarTa = -9999.0;
                else
                    a2dVarTa = transpose(a2dVar)
                endif
                !------------------------------------------------------------------------------------------
                
                !------------------------------------------------------------------------------------------
                ! INCOMING RADIATION
                sVarName = 'IncRadiation'
                call S3M_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, ' Get forcing gridded data FAILED! Check forcing data for '//sVarName//'!')
                    a2dVarIncRad = -9999.0;
                else
                    a2dVarIncRad = transpose(a2dVar)
                endif
                !------------------------------------------------------------------------------------------
                
                !------------------------------------------------------------------------------------------
                ! RELATIVE HUMIDITY
                sVarName = 'RelHumidity'
                call S3M_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, ' Get forcing gridded data FAILED! Check forcing data for '//sVarName//'!')
                    a2dVarRelHum = -9999.0;
                else
                    a2dVarRelHum = transpose(a2dVar)
                endif
                !------------------------------------------------------------------------------------------
                
                !------------------------------------------------------------------------------------------
                ! Closing netCDF file
                iErr = nf90_close(iFileID)
                ! Remove uncompressed file (to save space on disk)
                call S3M_Tools_Generic_RemoveFile(oS3M_Namelist(iID)%sCommandRemoveFile, &
                                                  sFileNameData_Temp, .false.)
                !------------------------------------------------------------------------------------------
                                        
                !------------------------------------------------------------------------------------------
                ! Info netCDF filename
                call mprintf(.true., iINFO_Verbose, ' Get filename (forcing gridded): '//trim(sFileNameData_Forcing)//' ... OK' )
                ! Info end
                call mprintf(.true., iINFO_Extra, ' Data :: Forcing gridded :: NetCDF ... OK' )
                !------------------------------------------------------------------------------------------
                
            endif
            !------------------------------------------------------------------------------------------
            
        endif
        !------------------------------------------------------------------------------------------
        
    end subroutine S3M_Data_Forcing_Gridded_NC
#endif
    !------------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------------
    ! Subroutine to read binary forcing data
    subroutine S3M_Data_Forcing_Gridded_Binary(iID, iFlagTypeData_Forcing, &
                                      sPathData_Forcing, &
                                      iRows, iCols, sTime, &
                                      a2dVarPrecip, a2dVarTa, a2dVarIncRad, &
                                      a2dVarRelHum, iScaleFactor_Forcing)
    
        !------------------------------------------------------------------------------------------
        ! Variable(s)
        integer(kind = 4)                   :: iID
        integer(kind = 4)                   :: iFlagTypeData_Forcing
                                      
        character(len = 256), intent(in)    :: sPathData_Forcing
        character(len = 700)                :: sFileNameData_Forcing, sFileNameData_Forcing_Zip, sFileNameData_Temp
        character(len = 700)                :: sCommandUnzipFile
        character(len = 256)                :: sVarName
        integer(kind = 4), intent(in)       :: iRows, iCols
        real(kind = 4)                      :: dVar
               
        character(len = 19), intent(in)     :: sTime
        character(len = 12)                 :: sTimeMonth
        
        real(kind = 4), dimension(iRows, iCols)                 :: a2dVar

        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarPrecip
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarTa
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarIncRad
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarRelHum       
       
        character(len = 256):: sVarUnits, sPID
        integer(kind = 4)   :: iErr
        integer(kind = 4)   :: iFileID, iScaleFactor_Forcing
        
        logical             :: bFileExist
        !------------------------------------------------------------------------------------------
	
        !------------------------------------------------------------------------------------------
        ! Initialize variable(s)
        a2dVarPrecip = -9999.0; a2dVarTa = -9999.0; a2dVarIncRad = -9999.0; a2dVarRelHum = -9999.0;

        sFileNameData_Forcing = ''; sFileNameData_Forcing_Zip = ''; sFileNameData_Temp = ''; sTimeMonth = ''; sPid = ''
        
        ! Checking date
        write(sTimeMonth,'(A,A,A)') sTime(1:4), sTime(6:7), sTime(9:10)
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Get global information
        sCommandUnzipFile = oS3M_Namelist(iID)%sCommandUnzipFile
        
        ! Info start
        call mprintf(.true., iINFO_Extra, ' Data :: Forcing gridded :: Binary ... ' )
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Info binary file(s) time step
        call mprintf(.true., iINFO_Verbose, ' Get (forcing gridded) at time '//trim(sTime)//' ... ')
        
        ! Get unique process ID
        sPID = adjustl(getProcessID())
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Precip  (example: Rain_201405010000.bin.gz)
        sFileNameData_Forcing = trim(sPathData_Forcing)//"Rain_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)// &
            ".bin"  
        sFileNameData_Temp = trim(sPathData_Forcing)//"Rain_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)//'_'//trim(sPID)// &
            ".bin"  
            
        call mprintf(.true., iINFO_Extra, ' Get filename (forcing gridded): '//trim(sFileNameData_Forcing) )

        ! Checking file input availability
        sFileNameData_Forcing_Zip = trim(sFileNameData_Forcing)//'.gz'
        inquire (file = sFileNameData_Forcing_Zip, exist = bFileExist)
        if ( .not. bFileExist ) then
            call mprintf(.true., iWARN, ' Problem opening uncompressed binary file: '// &
                         trim(sFileNameData_Forcing_Zip)//' --> Undefined forcing data values!' )
            a2dVar = -9999.0
        else
            ! Unzip file
            call S3M_Tools_Generic_UnzipFile(oS3M_Namelist(iID)%sCommandUnzipFile, &
                                             sFileNameData_Forcing_Zip, &
                                             sFileNameData_Temp, .true.)
                                             
            ! Read binary data
            if (iFlagTypeData_Forcing == 1) then
                call S3M_Tools_IO_Get2d_Binary_INT(sFileNameData_Temp, a2dVar, iRows, iCols, iScaleFactor_Forcing, .true., iErr) 
            elseif (iFlagTypeData_Forcing == 2) then
                call S3M_Tools_IO_Get2d_Binary_DBL(sFileNameData_Temp, a2dVar, iRows, iCols, iScaleFactor_Forcing, .true., iErr) 
            endif
   
            ! Remove uncompressed file (to save space on disk)
            call S3M_Tools_Generic_RemoveFile(oS3M_Namelist(iID)%sCommandRemoveFile, &
                                              sFileNameData_Temp, .false.)
            
        endif
        a2dVarPrecip = a2dVar
        !------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        ! Temperature  (example: Temperature_201405010000.bin.gz)
        sFileNameData_Forcing = trim(sPathData_Forcing)//"Temperature_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)// &
            ".bin"
        sFileNameData_Temp = trim(sPathData_Forcing)//"Temperature_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)//'_'//trim(sPID)// &
            ".bin"
        call mprintf(.true., iINFO_Extra, ' Get filename (forcing gridded): '//trim(sFileNameData_Forcing) )

        ! Checking file input availability
        sFileNameData_Forcing_Zip = trim(sFileNameData_Forcing)//'.gz'
        inquire (file = sFileNameData_Forcing_Zip, exist = bFileExist)
        if ( .not. bFileExist ) then
            call mprintf(.true., iWARN, ' Problem opening uncompressed binary file: '// &
                         trim(sFileNameData_Forcing_Zip)//' --> Undefined forcing data values!' )
            a2dVar = -9999.0
        else
            ! Unzip file
            call S3M_Tools_Generic_UnzipFile(oS3M_Namelist(iID)%sCommandUnzipFile, &
                                             sFileNameData_Forcing_Zip, &
                                             sFileNameData_Temp, .true.)
                                             
            ! Read binary data
            if (iFlagTypeData_Forcing == 1) then
                call S3M_Tools_IO_Get2d_Binary_INT(sFileNameData_Temp, a2dVar, iRows, iCols, iScaleFactor_Forcing, .true., iErr) 
            elseif (iFlagTypeData_Forcing == 2) then
                call S3M_Tools_IO_Get2d_Binary_DBL(sFileNameData_Temp, a2dVar, iRows, iCols, iScaleFactor_Forcing, .true., iErr) 
            endif
            
            ! Remove uncompressed file (to save space on disk)
            call S3M_Tools_Generic_RemoveFile(oS3M_Namelist(iID)%sCommandRemoveFile, &
                                              sFileNameData_Temp, .false.)
        endif
        a2dVarTa = a2dVar
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Radiation  (example: Radiation_201405010000.bin.gz)
        sFileNameData_Forcing = trim(sPathData_Forcing)//"Radiation_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)// &
            ".bin"
        sFileNameData_Temp = trim(sPathData_Forcing)//"Radiation_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)//'_'//trim(sPID)// &
            ".bin" 
        call mprintf(.true., iINFO_Extra, ' Get filename (forcing gridded): '//trim(sFileNameData_Forcing) )

        ! Checking file input availability
        sFileNameData_Forcing_Zip = trim(sFileNameData_Forcing)//'.gz'
        inquire (file = sFileNameData_Forcing_Zip, exist = bFileExist)
        if ( .not. bFileExist ) then
            call mprintf(.true., iWARN, ' Problem opening uncompressed binary file: '// &
                         trim(sFileNameData_Forcing_Zip)//' --> Undefined forcing data values!' )
            a2dVar = -9999.0
        else
            ! Unzip file
            call S3M_Tools_Generic_UnzipFile(oS3M_Namelist(iID)%sCommandUnzipFile, &
                                             sFileNameData_Forcing_Zip, &
                                             sFileNameData_Temp, .true.)
                                             
            ! Read binary data
            if (iFlagTypeData_Forcing == 1) then
                call S3M_Tools_IO_Get2d_Binary_INT(sFileNameData_Temp, a2dVar, iRows, iCols, iScaleFactor_Forcing, .true., iErr) 
            elseif (iFlagTypeData_Forcing == 2) then
                call S3M_Tools_IO_Get2d_Binary_DBL(sFileNameData_Temp, a2dVar, iRows, iCols, iScaleFactor_Forcing, .true., iErr) 
            endif
            
            ! Remove uncompressed file (to save space on disk)
            call S3M_Tools_Generic_RemoveFile(oS3M_Namelist(iID)%sCommandRemoveFile, &
                                              sFileNameData_Temp, .false.)
        endif
        a2dVarIncRad = a2dVar
        !------------------------------------------------------------------------------------------
 
        !------------------------------------------------------------------------------------------
        ! RelHum  (example: RelUmid_201405010000.bin.gz)
        sFileNameData_Forcing = trim(sPathData_Forcing)//"RelUmid_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)// &
            ".bin"
        sFileNameData_Temp = trim(sPathData_Forcing)//"RelUmid_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)//'_'//trim(sPID)// &
            ".bin"
        call mprintf(.true., iINFO_Extra, ' Get filename (forcing gridded): '//trim(sFileNameData_Forcing) )

        ! Checking file input availability
        sFileNameData_Forcing_Zip = trim(sFileNameData_Forcing)//'.gz'
        inquire (file = sFileNameData_Forcing_Zip, exist = bFileExist)
        if ( .not. bFileExist ) then
            call mprintf(.true., iWARN, ' Problem opening uncompressed binary file: '// &
                         trim(sFileNameData_Forcing_Zip)//' --> Undefined forcing data values!' )
            a2dVar = -9999.0
        else
            ! Unzip file
            call S3M_Tools_Generic_UnzipFile(oS3M_Namelist(iID)%sCommandUnzipFile, & 
                                             sFileNameData_Forcing_Zip, & 
                                             sFileNameData_Temp, .true.)
                                             
            ! Read binary data
            if (iFlagTypeData_Forcing == 1) then
                call S3M_Tools_IO_Get2d_Binary_INT(sFileNameData_Temp, a2dVar, iRows, iCols, iScaleFactor_Forcing, .true., iErr) 
            elseif (iFlagTypeData_Forcing == 2) then
                call S3M_Tools_IO_Get2d_Binary_DBL(sFileNameData_Temp, a2dVar, iRows, iCols, iScaleFactor_Forcing, .true., iErr) 
            endif
            
            ! Remove uncompressed file (to save space on disk)
            call S3M_Tools_Generic_RemoveFile(oS3M_Namelist(iID)%sCommandRemoveFile, &
                                              sFileNameData_Temp, .false.)
        endif
        a2dVarRelHum = a2dVar
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Info binary file(s) time step
        call mprintf(.true., iINFO_Verbose, ' Get (forcing gridded) at time '//trim(sTime)//' ... OK')
        ! Info end
        call mprintf(.true., iINFO_Extra, ' Data :: Forcing gridded :: Binary ... OK' )
        !------------------------------------------------------------------------------------------
        
    end subroutine S3M_Data_Forcing_Gridded_Binary
    !------------------------------------------------------------------------------------------
    
end module S3M_Module_Data_Forcing_Gridded
!-----------------------------------------------------------------------------------------

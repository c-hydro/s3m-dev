!------------------------------------------------------------------------------------------    
! File:   S3M_Module_Data_Restart_Gridded.f90
! Author(s): Fabio Delogu, Francesco Silvestro, Simone Gabellani, Francesco Avanzi.
!
! Created on May 7, 2015, 1:27 PM
! Last update on December 15, 2023 11:20 AM
!
! Module to read restart data.
!------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------
! Module Header
module S3M_Module_Data_Restart_Gridded
    
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
                                            S3M_Tools_IO_Get3d_Binary_INT, &
                                            S3M_Tools_IO_Get3d_Binary_DBL, &
                                            S3M_Tools_IO_Get2d_NC, &
                                            S3M_Tools_IO_Get3d_NC, &
                                            check
#else
    use S3M_Module_Tools_IO,        only:   S3M_Tools_IO_Get2d_Binary_INT, &
                                            S3M_Tools_IO_Get2d_Binary_DBL, &
                                            S3M_Tools_IO_Get3d_Binary_INT, &
                                            S3M_Tools_IO_Get3d_Binary_DBL                                     
#endif                                   
    
    use S3M_Module_Tools_Generic,   only:   S3M_Tools_Generic_ReplaceText, & 
                                            S3M_Tools_Generic_CreateFolder, &
                                            S3M_Tools_Generic_ZipFile, &
                                            S3M_Tools_Generic_UnzipFile, &
                                            S3M_Tools_Generic_RemoveFile, &
                                            transpose3Dvar, &
                                            checkdomainvar, getProcessID
    
    ! Implicit none for all subroutines in this module
    implicit none
    !------------------------------------------------------------------------------------------
    
contains
    !------------------------------------------------------------------------------------------
    ! Subroutine to manage restart gridded data
    subroutine S3M_Data_Restart_Gridded_Cpl( iID, sTime, &
                                             iRowsStart, iRowsEnd, &
                                             iColsStart, iColsEnd, &
                                             iDaySteps)

        !------------------------------------------------------------------------------------------
        ! Variable(s)                                    
        integer(kind = 4)           :: iID
        integer(kind = 4)           :: iFlagRestart
        integer(kind = 4)           :: iRows, iCols
        integer(kind = 4)           :: iRowsStart, iColsStart, iRowsEnd, iColsEnd
        integer(kind = 4)           :: iDaySteps

        integer(kind = 4)           :: iFlagIceMassBalance
        
        character(len = 19)         :: sTime
        character(len = 700)        :: sPathData_Restart
        character(len = 700)        :: sFileNameData_Restart, sFileNameData_Restart_Zip
        character(len = 700)        :: sCommandUnzipFile
        
        real(kind = 4), dimension(iRowsEnd - iRowsStart + 1, iColsEnd - iColsStart + 1) :: a2dVarDEM
        real(kind = 4), dimension(iRowsEnd - iRowsStart + 1, iColsEnd - iColsStart + 1) :: a2dVarLat, a2dVarLon
                                                  
        integer(kind = 4), dimension(iRowsEnd - iRowsStart + 1, iColsEnd - iColsStart + 1)  :: a2iVarAgeS        
        
        real(kind = 4), dimension(iRowsEnd - iRowsStart + 1, iColsEnd - iColsStart + 1)     ::  a2dVarSWE_D, &
                                                                                                a2dVarRho_D, &
                                                                                                a2dVarSWE_W, &
                                                                                                a2dVarSWE, &
                                                                                                a2dVarAlbedoS, &
                                                                                                a2dVarIceThick, &
                                                                                                a2dVarMeltingGCumWY, &
                                                                                                a2dVarTaC_MeanDaysSuppressMelt, &
                                                                                                a2dVarSnowFallDayCum, &
                                                                                                a2dVarMeltingSDayCum
           
        real(kind = 4), dimension(iRowsEnd - iRowsStart + 1, iColsEnd - iColsStart + 1, iDaySteps)      :: a3dVarTaC_1Days 
        
        logical                     :: bFileExist, bCheckRestart
        !------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        ! Initialize variable(s)
        a2dVarLat = 0.0; a2dVarLon = 0.0
        
        a2iVarAgeS = 0; a2dVarSWE_D = 0.0; a2dVarRho_D = 0.0; a2dVarSWE_W = 0.0; a2dVarSWE = 0.0
        a3dVarTaC_1Days = 0.0; a2dVarIceThick = 0.0; a2dVarMeltingGCumWY = 0.0;
        a2dVarAlbedoS = 0.0; a2dVarTaC_MeanDaysSuppressMelt = 0.0;
        a2dVarSnowFallDayCum = 0.0; a2dVarMeltingSDayCum = 0.0;
        
        bCheckRestart = .false.;
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Defining iRows and iCols (output data)
        iRows = iRowsEnd - iRowsStart + 1
        iCols = iColsEnd - iColsStart + 1
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Get global information
        iFlagRestart = oS3M_Namelist(iID)%iFlagRestart
        iFlagIceMassBalance = oS3M_Namelist(iID)%iFlagIceMassBalance
        sPathData_Restart = oS3M_Namelist(iID)%sPathData_Restart_Gridded
        sCommandUnzipFile = oS3M_Namelist(iID)%sCommandUnzipFile
        
        ! Get global variable(s)
        a2dVarDem = oS3M_Vars(iID)%a2dDem

        ! Info start
        call mprintf(.true., iINFO_Extra, ' Data :: Restart gridded ... ' )
        !------------------------------------------------------------------------------------------
          
        !------------------------------------------------------------------------------------------
        ! Check restart flag value
        if (iFlagRestart == 1) then
        
            !------------------------------------------------------------------------------------------
            ! Replace general path with specific time feature(s)
            call S3M_Tools_Generic_ReplaceText(sPathData_Restart, '$yyyy', sTime(1:4))
            call S3M_Tools_Generic_ReplaceText(sPathData_Restart, '$mm', sTime(6:7))
            call S3M_Tools_Generic_ReplaceText(sPathData_Restart, '$dd', sTime(9:10))
            call S3M_Tools_Generic_ReplaceText(sPathData_Restart, '$HH', sTime(12:13))
            call S3M_Tools_Generic_ReplaceText(sPathData_Restart, '$MM', sTime(15:16))
            !------------------------------------------------------------------------------------------

            !------------------------------------------------------------------------------------------
            ! Subroutine for reading netCDF restart data 
            
            call S3M_Data_Restart_Gridded_NC(iID, &
                                        sPathData_Restart, &
                                        iRows, iCols, &
                                        iDaySteps, &
                                        sTime, &
                                        a2iVarAgeS, a2dVarSWE_D, a2dVarRho_D, a2dVarSWE_W, a2dVarSWE, a2dVarAlbedoS, &
                                        a3dVarTaC_1Days, a2dVarTaC_MeanDaysSuppressMelt, &
                                        a2dVarLat, a2dVarLon, a2dVarIceThick, a2dVarMeltingGCumWY, iFlagIceMassBalance, &
                                        a2dVarMeltingSDayCum, a2dVarSnowFallDayCum, &
                                        bCheckRestart)
                                        
            !------------------------------------------------------------------------------------------
                                            
            !------------------------------------------------------------------------------------------
            ! Nullify map boundary
            where(a2dVarDEM.lt.0.0)
                a2iVarAgeS = -9999
                a2dVarSWE_D = -9999
                a2dVarRho_D = -9999
                a2dVarSWE_W = -9999
                a2dVarSWE = -9999                
                a2dVarIceThick = -9999
                a2dVarMeltingGCumWY = -9999
                a2dVarAlbedoS = -9999
                a2dVarTaC_MeanDaysSuppressMelt = -9999
                a2dVarMeltingSDayCum = -9999
                a2dVarSnowFallDayCum = -9999
            endwhere
            !------------------------------------------------------------------------------------------
            
            !------------------------------------------------------------------------------------------
            ! Check restart flag on data availability
            call mprintf(.true., iINFO_Extra, ' Data :: Restart gridded for snow physics ... ' )
            if (bCheckRestart .eqv. .true.) then
                
                !------------------------------------------------------------------------------------------
                ! Pass variable(s) to global workspace
                oS3M_Vars(iID)%a2dSWE_D = a2dVarSWE_D
                oS3M_Vars(iID)%a2dRho_D = a2dVarRho_D 
                oS3M_Vars(iID)%a2dSWE_W = a2dVarSWE_W
                oS3M_Vars(iID)%a2dSWE = a2dVarSWE
            
                if (oS3M_Namelist(iID)%iFlagThickFromTerrData .eq. 0) then
                    oS3M_Vars(iID)%a2dIceThick = a2dVarIceThick
                endif
            
                oS3M_Vars(iID)%a2dMeltingGCumWY = a2dVarMeltingGCumWY
            
                oS3M_Vars(iID)%a3dTaC_Days1 = a3dVarTaC_1Days
                oS3M_Vars(iID)%a2dAlbedo_Snow = a2dVarAlbedoS
                oS3M_Vars(iID)%a2dTaC_MeanDaysSuppressMelt = a2dVarTaC_MeanDaysSuppressMelt
                oS3M_Vars(iID)%a2dMeltingDayCum = a2dVarMeltingSDayCum
                oS3M_Vars(iID)%a2dSnowFallDayCum = a2dVarSnowFallDayCum           
                oS3M_Vars(iID)%a2iAge = a2iVarAgeS
                
                ! Info end
                call mprintf(.true., iINFO_Extra, ' Data :: Restart gridded for snow physics ... OK' )
                !------------------------------------------------------------------------------------------
                
            else
                !------------------------------------------------------------------------------------------
                ! Exit message for not using restart data
                call mprintf(.true., iINFO_Verbose, ' Restart flag selected but problems with reading snow data')
                ! Info end
                call mprintf(.true., iINFO_Extra, ' Data :: Restart gridded for snow physics ... SKIPPED' )
                !------------------------------------------------------------------------------------------
                
            endif
            !------------------------------------------------------------------------------------------  
            
        else
            !------------------------------------------------------------------------------------------
            ! Exit message for not using restart data
            call mprintf(.true., iINFO_Verbose, ' No restart run selected (gridded data)')
            ! Info end
            call mprintf(.true., iINFO_Extra, ' Data :: Restart gridded ... SKIPPED ' )
            !------------------------------------------------------------------------------------------
        endif
        !------------------------------------------------------------------------------------------            

        
        !------------------------------------------------------------------------------------------
        ! Debug
        if (iDEBUG.gt.0) then
            call mprintf(.true., iINFO_Extra, '')
            call mprintf(.true., iINFO_Extra, checkvar(real(oS3M_Vars(iID)%a2iAge), oS3M_Vars(iID)%a2iMask, 'AGES END') )
            call mprintf(.true., iINFO_Extra, checkvar(oS3M_Vars(iID)%a2dSWE_D, oS3M_Vars(iID)%a2iMask, 'SWED END') )
            call mprintf(.true., iINFO_Extra, checkvar(oS3M_Vars(iID)%a2dRho_D, oS3M_Vars(iID)%a2iMask, 'RHOD END') )
            call mprintf(.true., iINFO_Extra, checkvar(oS3M_Vars(iID)%a2dSWE_W, oS3M_Vars(iID)%a2iMask, 'SWEW END') )
            call mprintf(.true., iINFO_Extra, checkvar(oS3M_Vars(iID)%a2dSWE, oS3M_Vars(iID)%a2iMask, 'SWE END') )
            call mprintf(.true., iINFO_Extra, checkvar(oS3M_Vars(iID)%a2dAlbedo_Snow, oS3M_Vars(iID)%a2iMask, 'ALBEDO END') )
            
            if (oS3M_Namelist(iID)%iFlagThickFromTerrData .eq. 0) then            
              call mprintf(.true., iINFO_Extra, checkvar(oS3M_Vars(iID)%a2dIceThick, oS3M_Vars(iID)%a2iMask, 'ICE THICKNESS END') ) 
            endif
            
            call mprintf(.true., iINFO_Extra, checkvar(oS3M_Vars(iID)%a2dMeltingGCumWY, &
                                                                            oS3M_Vars(iID)%a2iMask, 'MELTING G CUM WY END') )
            call mprintf(.true., iINFO_Extra, checkvar(oS3M_Vars(iID)%a3dTaC_Days1(:,:,1), oS3M_Vars(iID)%a2iMask, 'TA 1DAYS END') )
            call mprintf(.true., iINFO_Extra, checkvar(oS3M_Vars(iID)%a2dTaC_MeanDaysSuppressMelt, & 
                                                                            oS3M_Vars(iID)%a2iMask, 'TA AVG SUPPRESSMELT END') )                            
            call mprintf(.true., iINFO_Extra, checkvar(oS3M_Vars(iID)%a2dMeltingDayCum, oS3M_Vars(iID)%a2iMask, 'MELTCUM END') )
            call mprintf(.true., iINFO_Extra, checkvar(oS3M_Vars(iID)%a2dSnowFallDayCum, oS3M_Vars(iID)%a2iMask, 'SNOWFCUM END') )            
            call mprintf(.true., iINFO_Extra, '')
            call mprintf(.true., iINFO_Extra, ' ========= RESTART GRIDDED END =========== ')
        endif
        !------------------------------------------------------------------------------------------
        
    end subroutine S3M_Data_Restart_Gridded_Cpl
    !------------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------------
    ! Subroutine to read netCDF data restart
#ifdef LIB_NC
    subroutine S3M_Data_Restart_Gridded_NC(iID, &
                                           sPathData_Restart, &
                                           iRows, iCols, &
                                           iDaySteps, &
                                           sTime, &
                                           a2iVarAgeS, a2dVarSWE_D, a2dVarRho_D, a2dVarSWE_W, a2dVarSWE, a2dVarAlbedoS, &
                                           a3dVarTaC_1Days, a2dVarTaC_MeanDaysSuppressMelt, &
                                           a2dVarLat, a2dVarLon, a2dVarIceThick, a2dVarMeltingGCumWY, iFlagIceMassBalance, &
                                           a2dVarMeltingSDayCum, a2dVarSnowFallDayCum, &
                                           bCheckRestart)
        !------------------------------------------------------------------------------------------
        ! Variable(s)
        integer(kind = 4)                       :: iID                  
                                  
        character(len = 256), intent(in)        :: sPathData_Restart
        character(len = 700)                    :: sFileNameData_Restart, sFileNameData_Restart_Zip, sFileNameData_Temp
        character(len = 700)                    :: sCommandUnzipFile
        character(len = 256)                    :: sVarName
        integer(kind = 4), intent(in)           :: iRows, iCols
        integer(kind = 4), intent(in)           :: iDaySteps

        character(len = 19), intent(in)         :: sTime

        real(kind = 4), dimension(iCols, iRows)                                :: a2dVar
        real(kind = 4), dimension(iCols, iRows, iDaySteps)                     :: a3dVar3
        real(kind = 4), dimension(iCols, iRows, iDaySteps*5)                   :: a3dVar4
        
        integer(kind = 4), dimension(iRows, iCols),             intent(out)    :: a2iVarAgeS
        real(kind = 4), dimension(iRows, iCols),                intent(out)    :: a2dVarSWE_D, a2dVarRho_D, a2dVarSWE_W, & 
                                                                                  a2dVarSWE, a2dVarAlbedoS
        real(kind = 4), dimension(iRows, iCols, iDaySteps),     intent(out)    :: a3dVarTaC_1Days
        real(kind = 4), dimension(iRows, iCols),                intent(out)    :: a2dVarIceThick, a2dVarMeltingGCumWY
        real(kind = 4), dimension(iRows, iCols),                intent(out)    :: a2dVarTaC_MeanDaysSuppressMelt, &
                                                                                  a2dVarMeltingSDayCum, &
                                                                                  a2dVarSnowFallDayCum                   
       
        real(kind = 4), dimension(iRows, iCols),                intent(out)    :: a2dVarLat
        real(kind = 4), dimension(iRows, iCols),                intent(out)    :: a2dVarLon

        character(len = 256)    :: sVarUnits, sPID
        integer(kind = 4)       :: iErr
        integer(kind = 4)       :: iFileID
        integer(kind = 4)       :: iFlagIceMassBalance

        logical                 :: bFileExist, bCheckRestart, bCheckVar
        
        
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Initialize variable(s)
        a2iVarAgeS = -9999; a2dVarSWE_D = -9999.0; a2dVarRho_D = -9999.0; a2dVarSWE_W = -9999.0; a2dVarSWE = -9999.0
        a3dVarTaC_1Days = -9999.0; a2dVarTaC_MeanDaysSuppressMelt = -9999.0; a2dVarAlbedoS = -9999
        a2dVarIceThick = -9999.0; a2dVarMeltingGCumWY = 0; a2dVarMeltingSDayCum = -9999; a2dVarSnowFallDayCum = -9999;
        a2dVarLat = -9999.0; a2dVarLon = -9999.0;
        bCheckRestart = .false.; 
        bCheckVar = .true.;
        
        sPID = ''; sFileNameData_Temp = ''
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Get global information
        sCommandUnzipFile = oS3M_Namelist(iID)%sCommandUnzipFile
        ! Info start
        call mprintf(.true., iINFO_Extra, ' Data :: Restart gridded :: NetCDF ... ' )
        
        ! Get unique process ID
        sPID = adjustl(getProcessID())
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Filename restart (example: S3M_201404300000.nc)
        sFileNameData_Restart = trim(sPathData_Restart)//"S3M_"// &
        sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
        sTime(12:13)//sTime(15:16)// &
        ".nc"
        ! Create Filename with unique PID number to avoid simultaneously access to the same Forcing file       
        sFileNameData_Temp = trim(sPathData_Restart)//"S3M_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)//'_'//trim(sPID)// &
            ".nc"  

        ! Info netCDF filename
        call mprintf(.true., iINFO_Basic, ' Get filename (restart gridded): '//trim(sFileNameData_Restart)//' ... ' )
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Checking file input availability
        sFileNameData_Restart_Zip = trim(sFileNameData_Restart)//'.gz'
        inquire (file = trim(sFileNameData_Restart)//'.gz', exist = bFileExist)
        if ( .not. bFileExist ) then
            !------------------------------------------------------------------------------------------
            ! Warning message
            call mprintf(.true., iWARN, ' No compressed restart netCDF data found: '//trim(sFileNameData_Restart_Zip) )
            call mprintf(.true., iINFO_Verbose, &
                         ' Get filename (restart gridded): '//trim(sFileNameData_Restart)//' ... FAILED' )
            bCheckVar = .false.
            !------------------------------------------------------------------------------------------
        else
            !------------------------------------------------------------------------------------------
            ! Unzip file
            call S3M_Tools_Generic_UnzipFile(sCommandUnzipFile, &
                                             sFileNameData_Restart_Zip, &
                                             sFileNameData_Temp, .true.)
            !------------------------------------------------------------------------------------------
        
            !------------------------------------------------------------------------------------------
            ! Opening netCDF file
            iErr = nf90_open(trim(sFileNameData_Temp), NF90_NOWRITE, iFileID)
            if (iErr /= 0) then
                
                !------------------------------------------------------------------------------------------
                ! Condition for no file restart found
                call mprintf(.true., iWARN, ' Problem opening uncompressed netCDF file: '// &
                             trim(sFileNameData_Restart)//' --> Undefined restart data values' ) 
                call mprintf(.true., iINFO_Verbose, &
                            ' Get filename (restart gridded): '//trim(sFileNameData_Restart)//' ... FAILED' )
                bCheckVar = .false.
                !------------------------------------------------------------------------------------------
                            
            else
                
                !------------------------------------------------------------------------------------------
                ! Condition for file restart found
                ! SWE_D
                sVarName = 'SWE_D'
                call S3M_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, ' Get restart gridded data FAILED! '// &
                        'Snow physics is activated! If needed check restart data for '//sVarName//'!')
                    a2dVarSWE_D = -9999.0;
                    bCheckVar = bCheckVar .and. .false. 
                else
                    a2dVarSWE_D = transpose(a2dVar)
                    bCheckVar = bCheckVar .and. .true. 
                endif                
                
                ! RHO_D
                sVarName = 'Rho_D'
                call S3M_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, ' Get restart gridded data FAILED! '// &
                        'Snow physics is activated! If needed check restart data for '//sVarName//'!')
                    a2dVarRho_D = -9999.0;
                    bCheckVar = bCheckVar .and. .false. 
                else
                    a2dVarRho_D = transpose(a2dVar)
                    bCheckVar = bCheckVar .and. .true. 
                endif                 
                
                ! SWE_W
                sVarName = 'SWE_W'
                call S3M_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, ' Get restart gridded data FAILED! '// &
                        'Snow physics is activated! If needed check restart data for '//sVarName//'!')
                    a2dVarSWE_W = -9999.0;
                    bCheckVar = bCheckVar .and. .false. 
                else
                    a2dVarSWE_W = transpose(a2dVar)
                    bCheckVar = bCheckVar .and. .true. 
                endif                  
                
                ! SWE
                sVarName = 'SWE'
                call S3M_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, ' Get restart gridded data FAILED! '// &
                        'Snow physics is activated! If needed check restart data for '//sVarName//'!')
                    a2dVarSWE = -9999.0;
                    bCheckVar = bCheckVar .and. .false. 
                else
                    a2dVarSWE = transpose(a2dVar)
                    bCheckVar = bCheckVar .and. .true. 
                endif
                  
                ! Albedo
                sVarName = 'AlbedoS'
                call S3M_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, ' Get restart gridded data FAILED! '// &
                        'Snow physics is activated! If needed check restart data for '//sVarName//'!')
                    a2dVarAlbedoS = -9999.0;
                    bCheckVar = bCheckVar .and. .false. 
                else
                    a2dVarAlbedoS = transpose(a2dVar)
                    bCheckVar = bCheckVar .and. .true. 
                endif            
                
                ! Snow age
                sVarName = 'AgeS';
                call S3M_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, ' Get restart gridded data FAILED! '// &
                        'Snow physics is activated! If needed check restart data for '//sVarName//'!')
                    a2iVarAgeS = -9999;
                    bCheckVar = bCheckVar .and. .false. 
                else
                    a2iVarAgeS = int(transpose(a2dVar))
                    bCheckVar = bCheckVar .and. .true. 
                endif
                
                ! Air T avg SuppressMelt Days
                sVarName = 'T_SuppressMeltDays'
                call S3M_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, ' Get restart gridded data FAILED! '// &
                        'Snow physics is activated! If needed check restart data for '//sVarName//'!')
                    a2dVarTaC_MeanDaysSuppressMelt = -9999.0;
                    bCheckVar = bCheckVar .and. .false. 
                else
                    a2dVarTaC_MeanDaysSuppressMelt = transpose(a2dVar)
                    bCheckVar = bCheckVar .and. .true. 
                endif              
                    
                ! Air temperature last 1 day(s) 
                sVarName = 'T_1Days';
                call S3M_Tools_IO_Get3d_NC((sVarName), iFileID, a3dVar3, sVarUnits, iDaySteps, iCols, iRows, .true., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, ' Get restart gridded data FAILED! '// &
                       'Snow physics is activated! If needed check restart data for '//sVarName//'!')
                    a3dVarTaC_1Days = -9999.0;
                    bCheckVar = bCheckVar .and. .false. 
                else
                    a3dVarTaC_1Days = transpose3Dvar(a3dVar3)
                    bCheckVar = bCheckVar .and. .true. 
                endif
                
                if ( (iFlagIceMassBalance.eq.1.0) .or. (iFlagIceMassBalance.eq.2.0) ) then
                    
                    ! Ice thickness
                    if (oS3M_Namelist(iID)%iFlagThickFromTerrData .eq. 0) then
                    !We load Thickness from Restart data only if we said S3M to do so in the infoFile. 
                        sVarName = 'Thickness'
                        call S3M_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                        if(iErr /= 0) then
                            call mprintf(.true., iWARN, ' Get restart gridded data FAILED! '// &
                                'Snow physics is activated! If needed check restart data for '//sVarName//'!')
                            a2dVarIceThick = -9999.0;
                            bCheckVar = bCheckVar .and. .false. 
                        else
                            a2dVarIceThick = transpose(a2dVar)
                            bCheckVar = bCheckVar .and. .true. 
                        endif

                    endif

                    ! Cum_WY_MeltingG
                    sVarName = 'Cum_WY_MeltingG'
                    call S3M_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                    if(iErr /= 0) then
                        call mprintf(.true., iWARN, ' Get restart gridded data FAILED! '// &
                            'Snow physics is activated! If needed check restart data for '//sVarName//'!')
                        a2dVarMeltingGCumWY = 0.0;
                        bCheckVar = bCheckVar .and. .false.                         
                    else
                        a2dVarMeltingGCumWY = transpose(a2dVar)
                        bCheckVar = bCheckVar .and. .true. 
                    endif 
                    
                endif
                
                ! MeltingDayCum
                sVarName = 'MeltingSDayCum'
                call S3M_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, ' Get restart gridded data FAILED! '// &
                        'Snow physics is activated! If needed check restart data for '//sVarName//'!')
                    a2dVarMeltingSDayCum = -9999.0;
                    bCheckVar = bCheckVar .and. .false. 
                else
                    a2dVarMeltingSDayCum = transpose(a2dVar)
                    bCheckVar = bCheckVar .and. .true. 
                endif 
                
                ! SnowFallDayCum
                sVarName = 'SnowfallCum'
                call S3M_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, ' Get restart gridded data FAILED! '// &
                        'Snow physics is activated! If needed check restart data for '//sVarName//'!')
                    a2dVarSnowFallDayCum = -9999.0;
                    bCheckVar = bCheckVar .and. .false. 
                else
                    a2dVarSnowFallDayCum = transpose(a2dVar)
                    bCheckVar = bCheckVar .and. .true. 
                endif                 
                
                
                ! Closing netcdf file (drops db)
                iErr = nf90_close(iFileID)
                
                ! Remove uncompressed file (to save space on disk)
                call S3M_Tools_Generic_RemoveFile(oS3M_Namelist(iID)%sCommandRemoveFile, &
                                                  sFileNameData_Temp, .false.)
                !------------------------------------------------------------------------------------------
                    
                !------------------------------------------------------------------------------------------
                ! Info filename
                call mprintf(.true., iINFO_Basic, ' Get filename (restart gridded): '//trim(sFileNameData_Restart)//' ... OK' )
                ! Info end
                call mprintf(.true., iINFO_Extra, ' Data :: Restart gridded :: NetCDF ... OK' )
                !------------------------------------------------------------------------------------------
                
            endif
            !------------------------------------------------------------------------------------------
            
        endif
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Check restart
        if (bCheckVar .eqv. .true.) then
            call mprintf(.true., iINFO_Basic, ' Data :: Restart gridded :: NetCDF :: All variable(s) are loaded! ' )
            bCheckRestart = .true.
        else
            call mprintf(.true., iINFO_Basic, ' Data :: Restart gridded :: NetCDF :: Some/All variable(s) are N/A! ' )
            call mprintf(.true., iWARN, ' Restart flag activated but some data restart are not available! ')
            call mprintf(.true., iWARN, ' Restart gridded conditions are null! ')
            bCheckRestart = .false.
        endif
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Debug
        if (iDEBUG.gt.0) then
            call mprintf(.true., iINFO_Extra, ' ========= CHECK FORCING GRIDDED NC =========== ')
            call mprintf(.true., iINFO_Extra, checkvar(real(a2iVarAgeS), oS3M_Vars(iID)%a2iMask, 'AGES NC') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarSWE_D, oS3M_Vars(iID)%a2iMask, 'SWE_D NC') )            
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarRho_D, oS3M_Vars(iID)%a2iMask, 'RHO_D NC') )             
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarSWE_W, oS3M_Vars(iID)%a2iMask, 'SWE_W NC') )              
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarSWE, oS3M_Vars(iID)%a2iMask, 'SWE NC') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarAlbedoS, oS3M_Vars(iID)%a2iMask, 'ALBEDO NC') )
            
            if (oS3M_Namelist(iID)%iFlagThickFromTerrData .eq. 0) then
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarIceThick, oS3M_Vars(iID)%a2iMask, 'ICE THICKNESS NC') )
            endif
            
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarMeltingGCumWY, oS3M_Vars(iID)%a2iMask, 'MELTING G CUM WY NC') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarTaC_MeanDaysSuppressMelt, & 
                                                        oS3M_Vars(iID)%a2iMask, 'TA AVG SUPPRESSMELTDAYS NC') )     
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarMeltingSDayCum, oS3M_Vars(iID)%a2iMask, 'MELTCUM NC') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarSnowFallDayCum, oS3M_Vars(iID)%a2iMask, 'SNOWFCUM NC') )             
            call mprintf(.true., iINFO_Extra, ' ========= CHECK FORCING GRIDDED NC =========== ')
        endif
        !------------------------------------------------------------------------------------------
        
    end subroutine S3M_Data_Restart_Gridded_NC
#endif
    !------------------------------------------------------------------------------------------
    
end module S3M_Module_Data_Restart_Gridded
!------------------------------------------------------------------------------------------

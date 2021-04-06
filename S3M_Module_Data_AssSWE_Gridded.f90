!------------------------------------------------------------------------------------------    
! File:   S3M_Module_Data_AssSWE_Gridded.f90
! Author(s): Francesco Avanzi
!
! Created on September 05, 2018, 11:00 AM
! Last update on October 26, 2020 08:55 AM
!
! Module to read assimilated SWE map.
!------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------
! Module Header
module S3M_Module_Data_AssSWE_Gridded
    
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
                                            check
#else
    use S3M_Module_Tools_IO,        only:   S3M_Tools_IO_Get2d_Binary_INT, &
                                            S3M_Tools_IO_Get2d_Binary_DBL                                      
#endif                                    
                                                                                  
    use S3M_Module_Tools_Generic,   only:   S3M_Tools_Generic_ReplaceText, &
                                            S3M_Tools_Generic_SwitchGrid, &
                                            S3M_Tools_Generic_UnzipFile, &
                                            S3M_Tools_Generic_RemoveFile, &
                                            check2Dvar
                                            
    use S3M_Module_Tools_Time,      only:   S3M_Tools_Time_MonthVal
                             
    ! Implicit none for all subroutines in this module
    implicit none
    !------------------------------------------------------------------------------------------
    
contains
    
    !------------------------------------------------------------------------------------------
    ! Subroutine to read SWE gridded data
    subroutine S3M_Data_AssSWE_Gridded_Cpl( iID, sTime, &
                                     iRowsStartL, iRowsEndL, iColsStartL, iColsEndL)
        
        !------------------------------------------------------------------------------------------
        ! Variable(s)
        integer(kind = 4)           :: iID
        
        integer(kind = 4)           :: iRowsStartL, iRowsEndL, iColsStartL, iColsEndL
        integer(kind = 4)           :: iRowsL, iColsL
        integer(kind = 4)           :: iFlagTypeData_Ass_SWE
        integer(kind = 4)           :: iVarDayNud
        integer(kind = 4)           :: iVarSWEassInfluence
                
        character(len = 19)         :: sTime
        character(len = 2)          :: sTimeVarSWEass
        character(len = 12)         :: sTimeMonth
        
        character(len = 256)        :: sPathData_SWE_Assimilation

        real(kind = 4), dimension(iRowsEndL - iRowsStartL + 1, &
                                  iColsEndL - iColsStartL + 1) :: a2dVarDem
                                  
        real(kind = 4), dimension(iRowsEndL - iRowsStartL + 1, &
                                  iColsEndL - iColsStartL + 1) :: a2dVarSWEass, a2dVarSWEassHist
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Initialize variable(s)
        iVarDayNud = 0;
        a2dVarSWEass = -9999.0; a2dVarSWEassHist = -9999.0;
        sTimeVarSWEass = '';
        !------------------------------------------------------------------------------------------
                                                                                                
        !------------------------------------------------------------------------------------------
        ! Defining iRows and iCols (Land data)
        iRowsL = iRowsEndL - iRowsStartL + 1
        iColsL = iColsEndL - iColsStartL + 1
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Get global information
        sPathData_SWE_Assimilation = oS3M_Namelist(iID)%sPathData_SWE_Assimilation_Gridded
        iFlagTypeData_Ass_SWE = oS3M_Namelist(iID)%iFlagTypeData_Ass_SWE_Gridded
        iVarSWEassInfluence = oS3M_Namelist(iID)%iSWEassInfluence
        a2dVarDem = oS3M_Vars(iID)%a2dDem

        ! Info start
        call mprintf(.true., iINFO_Extra, ' Data :: SWE gridded assimilation data... ' )
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Replace general path with specific time feature(s)
        call S3M_Tools_Generic_ReplaceText(sPathData_SWE_Assimilation, '$yyyy', sTime(1:4))
        call S3M_Tools_Generic_ReplaceText(sPathData_SWE_Assimilation, '$mm', sTime(6:7))
        call S3M_Tools_Generic_ReplaceText(sPathData_SWE_Assimilation, '$dd', sTime(9:10))
        call S3M_Tools_Generic_ReplaceText(sPathData_SWE_Assimilation, '$HH', sTime(12:13))
        call S3M_Tools_Generic_ReplaceText(sPathData_SWE_Assimilation, '$MM', sTime(15:16))
        
        ! Checking date
        write(sTimeMonth,'(A,A,A)') sTime(1:4), sTime(6:7), sTime(9:10)
        !------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        ! Check time step (iT)
        if (oS3M_Vars(iID)%iTime .lt. oS3M_Namelist(iID)%iNTime) then

            !------------------------------------------------------------------------------------------
            ! Subroutine for reading netCDF SWE data 
            if (iFlagTypeData_Ass_SWE == 3) then

                !------------------------------------------------------------------------------------------
                ! Call subroutine to get SWE data in netCDF format
#ifdef LIB_NC
                call S3M_Data_AssSWE_Gridded_NC(iID, &
                                        sPathData_SWE_Assimilation, &
                                        iRowsL, iColsL, &
                                        sTime, &
                                        a2dVarSWEass)
#else   
                ! Redefinition of SWE data flag (if netCDF library is not linked)
                iFlagTypeData_Ass_SWE = 1 
                call mprintf(.true., iWARN, ' '// &
                                            'SWE gridded data type selected was netCDF, but library is not linked! '// &
                                            'Data in binary int format will be used!')
#endif
                !------------------------------------------------------------------------------------------
                           
            endif
            !------------------------------------------------------------------------------------------

            !------------------------------------------------------------------------------------------
            ! Subroutine for reading binary SWE data
            if (iFlagTypeData_Ass_SWE == 1 .or. iFlagTypeData_Ass_SWE == 2) then

                !------------------------------------------------------------------------------------------
                ! Calling subroutine to read data in binary format
                call S3M_Data_AssSWE_Gridded_Binary(iID, iFlagTypeData_Ass_SWE, &
                                            sPathData_SWE_Assimilation, &
                                            iRowsL, iColsL, sTime, &
                                            a2dVarSWEass)
                !------------------------------------------------------------------------------------------

            endif
            !------------------------------------------------------------------------------------------

            !------------------------------------------------------------------------------------------
            ! Debug
            if (iDEBUG.gt.0) then
                call mprintf(.true., iINFO_Extra, ' ========= SWE ASSIMILATION GRIDDED START =========== ')
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarSWEass, oS3M_Vars(iID)%a2iMask, 'SWE ASSIMILATION START') )
                call mprintf(.true., iINFO_Extra, '')                
            endif
            !------------------------------------------------------------------------------------------
 
        else
            
            !------------------------------------------------------------------------------------------
            ! Extra steps condition
            a2dVarSWEass = -9999.0; 
            
            ! Info message for extra time step(s)
            call mprintf(.true., iINFO_Extra, ' Extra time step ---> SWE assimilation data are set to null')
            call mprintf(.true., iINFO_Extra, ' Extra time step ---> SWE assimilation data are set to -9999.0')
            !------------------------------------------------------------------------------------------
            
        endif
        !------------------------------------------------------------------------------------------
       
        !------------------------------------------------------------------------------------------
        ! Nullify map boundary
        where(a2dVarDem.lt.0.0)
            a2dVarSWEass = -9999
        endwhere
        !------------------------------------------------------------------------------------------        
        
        !------------------------------------------------------------------------------------------
        ! Save SWE assimilation data to local variable(s) and to global workspace
        ! Snow Water Equivalent assimilated (SWE assimilated)
        if ( .not. all(a2dVarSWEass.eq.-9999.0) ) then
            oS3M_Vars(iID)%a2dSWEass = a2dVarSWEass
            a2dVarSWEassHist = a2dVarSWEass
            oS3M_Vars(iID)%a2dSWEassHist = a2dVarSWEassHist
            
            sTimeVarSWEass = sTime(12:13)
            oS3M_Vars(iID)%sTimeSWEass = sTimeVarSWEass
            
            iVarDayNud = iVarSWEassInfluence
            oS3M_Vars(iID)%iDayNud = iVarDayNud
            
        else
            
            iVarDayNud = oS3M_Vars(iID)%iDayNud
            
            if(iVarDayNud.gt.0) then
                
                sTimeVarSWEass = oS3M_Vars(iID)%sTimeSWEass
                if(sTime(12:13).eq.sTimeVarSWEass) then
                    a2dVarSWEassHist = oS3M_Vars(iID)%a2dSWEassHist
                    oS3M_Vars(iID)%a2dSWEass = a2dVarSWEassHist
                    iVarDayNud = iVarDayNud -1
                    
                    if(iVarDayNud.lt.0) then
                        iVarDayNud = 0
                    endif
                    
                    oS3M_Vars(iID)%iDayNud = iVarDayNud
                    
                endif
                
                oS3M_Vars(iID)%a2dSWEass = -9999.0
                
            else
            
                oS3M_Vars(iID)%a2dSWEass = -9999.0
                oS3M_Vars(iID)%a2dSWEassHist = -9999.0
                
            endif
            
        endif
        !
        !------------------------------------------------------------------------------------------  

        !------------------------------------------------------------------------------------------
        ! Debug
        !
        if (iDEBUG.gt.0) then
            call mprintf(.true., iINFO_Extra, '')
            call mprintf(.true., iINFO_Extra, checkvar(oS3M_Vars(iID)%a2dSWEass, &
            oS3M_Vars(iID)%a2iMask, 'SWE ASSIMILATION END') )
            call mprintf(.true., iINFO_Extra, checkvar(oS3M_Vars(iID)%a2dSWEassHist, &
            oS3M_Vars(iID)%a2iMask, 'SWE hist ASSIMILATION END') )
            call mprintf(.true., iINFO_Extra, ' ========= SWE GRIDDED END =========== ')
        endif
        
        ! Info end
        call mprintf(.true., iINFO_Extra, ' Data :: SWE assimilation gridded ... OK' )
        !------------------------------------------------------------------------------------------
        
    end subroutine S3M_Data_AssSWE_Gridded_Cpl
    !------------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------------
    ! Subroutine to read NC data SWE assimilated
#ifdef LIB_NC
    subroutine S3M_Data_AssSWE_Gridded_NC(iID,  &
                                  sPathData_SWE_Assimilation, &
                                  iRows, iCols, sTime, &
                                  a2dVarSWEass)
                                  
        !------------------------------------------------------------------------------------------
        ! Variable(s)
        integer(kind = 4)                       :: iID                  
        
        character(len = 256), intent(in)        :: sPathData_SWE_Assimilation
        character(len = 700)                    :: sFileNameData_SWEass, sFileNameData_SWEass_Zip
        character(len = 700)                    :: sCommandUnzipFile, sCommandRemoveFile
        character(len = 256)                    :: sVarName
        integer(kind = 4), intent(in)           :: iRows, iCols

        character(len = 19), intent(in)         :: sTime
        character(len = 12)                     :: sTimeMonth
                
        real(kind = 4), dimension(iCols, iRows)                 :: a2dVar
        
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarSWEass
       
        character(len = 256):: sVarUnits
        integer(kind = 4)   :: iErr
        integer(kind = 4)   :: iFileID
        
        logical             :: bFileExist
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Initialize variable(s)
        a2dVarSWEass = -9999.0;
        sFileNameData_SWEass = ''; sFileNameData_SWEass_Zip = ''; sTimeMonth = ''
        
        ! Checking date
        write(sTimeMonth,'(A,A,A)') sTime(1:4), sTime(6:7), sTime(9:10)
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Get global information
        sCommandUnzipFile = oS3M_Namelist(iID)%sCommandUnzipFile
        sCommandRemoveFile = oS3M_Namelist(iID)%sCommandRemoveFile
        
        ! Info start
        call mprintf(.true., iINFO_Extra, ' Data :: SWE gridded :: NetCDF ... ' )
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Filename SWEass (example: S3M.SWEass-grid.201404300000.nc.gz)
        sFileNameData_SWEass = trim(sPathData_SWE_Assimilation)//"SWEass_"// &
        sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
        sTime(12:13)//sTime(15:16)// &
        ".nc"

        ! Info netCDF filename
        call mprintf(.true., iINFO_Verbose, ' Get filename (SWE gridded): '//trim(sFileNameData_SWEass)//' ... ' )
        !------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        ! Checking file input availability
        sFileNameData_SWEass_Zip = trim(sFileNameData_SWEass)//'.gz'
        inquire (file = trim(sFileNameData_SWEass)//'.gz', exist = bFileExist)
        if ( .not. bFileExist ) then
            !------------------------------------------------------------------------------------------
            ! Warning message
            call mprintf(.true., iWARN, ' No compressed SWE netCDF data found: '//trim(sFileNameData_SWEass_Zip) )
            ! Info netCDF filename
            call mprintf(.true., iINFO_Verbose, &
                         ' Get filename (SWE gridded): '//trim(sFileNameData_SWEass)//' ... FAILED' )
            ! Info end
            call mprintf(.true., iINFO_Extra, ' Data :: SWE assimilation gridded :: NetCDF ... SKIPPED!' )
            !------------------------------------------------------------------------------------------
        else
            
            !------------------------------------------------------------------------------------------
            ! Unzip file
            call S3M_Tools_Generic_UnzipFile(oS3M_Namelist(iID)%sCommandUnzipFile, &
                                             sFileNameData_SWEass_Zip, &
                                             sFileNameData_SWEass, .true.)
            !------------------------------------------------------------------------------------------
        
            !------------------------------------------------------------------------------------------
            ! Open netCDF file
            iErr = nf90_open(trim(sFileNameData_SWEass), NF90_NOWRITE, iFileID)
            if (iErr /= 0) then
                call mprintf(.true., iWARN, ' Problem opening uncompressed netCDF file: '// &
                             trim(sFileNameData_SWEass)//' --> Undefined SWE data values' )
                call mprintf(.true., iINFO_Verbose, &
                             ' Get filename (SWE gridded): '//trim(sFileNameData_SWEass)//' ... FAILED' )
            else
                               
                !------------------------------------------------------------------------------------------
                
                ! Snow variable(s)
                !------------------------------------------------------------------------------------------
                ! SWE assimilated 
                sVarName = 'SWEAssimilated'
                call S3M_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, &
                        ' Get SWE gridded data FAILED! Check SWE data for '//sVarName//'!'// &
                            ' Snow physics is activated! If needed check SWE data!')
                    a2dVarSWEass = -9999.0;
                else
                    a2dVarSWEass = transpose(a2dVar)
                endif
                !------------------------------------------------------------------------------------------
               
                !------------------------------------------------------------------------------------------
                ! Closing netCDF file
                iErr = nf90_close(iFileID)
                ! Remove uncompressed file (to save space on disk)
                call S3M_Tools_Generic_RemoveFile(oS3M_Namelist(iID)%sCommandRemoveFile, &
                                                  sFileNameData_SWEass, .false.)
                !------------------------------------------------------------------------------------------
                
                !------------------------------------------------------------------------------------------
                ! Info netCDF filename
                call mprintf(.true., iINFO_Verbose, ' Get filename (SWE gridded): '//trim(sFileNameData_SWEass)//' ... OK' )
                ! Info end
                call mprintf(.true., iINFO_Extra, ' Data :: SWE gridded :: NetCDF ... OK' )
                !------------------------------------------------------------------------------------------
                
            endif
            !------------------------------------------------------------------------------------------
            
        endif
        !------------------------------------------------------------------------------------------
        
    end subroutine S3M_Data_AssSWE_Gridded_NC
#endif
    !------------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------------
    ! Subroutine to read binary SWE assimilated data
    subroutine S3M_Data_AssSWE_Gridded_Binary(iID, iFlagTypeData_Ass_SWE, &
                                      sPathData_SWE_Assimilation, &
                                      iRows, iCols, sTime, &
                                      a2dVarSWEass)
    
        !------------------------------------------------------------------------------------------
        ! Variable(s)
        integer(kind = 4)                   :: iID
        integer(kind = 4)                   :: iFlagTypeData_Ass_SWE
                                      
        character(len = 256), intent(in)    :: sPathData_SWE_Assimilation
        character(len = 700)                :: sFileNameData_SWEass, sFileNameData_SWEass_Zip
        character(len = 700)                :: sCommandUnzipFile
        character(len = 256)                :: sVarName
        
        integer(kind = 4), intent(in)       :: iRows, iCols
        
        character(len = 19), intent(in)     :: sTime
        character(len = 12)                 :: sTimeMonth
        character(len = 256):: sVarUnits
       
        real(kind = 4), dimension(iRows, iCols)                 :: a2dVar
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarSWEass     
       
        integer(kind = 4)   :: iErr
        integer(kind = 4)   :: iFileID, iScaleFactor_SWEass
        
        logical             :: bFileExist
        !------------------------------------------------------------------------------------------
	
        !------------------------------------------------------------------------------------------
        ! Initialize variable(s)
        a2dVarSWEass = -9999.0;        
        
        sFileNameData_SWEass = ''; sFileNameData_SWEass_Zip = ''; sTimeMonth = ''
        
        ! Checking date
        write(sTimeMonth,'(A,A,A)') sTime(1:4), sTime(6:7), sTime(9:10)
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Get global information
        sCommandUnzipFile = oS3M_Namelist(iID)%sCommandUnzipFile
        iScaleFactor_SWEass = oS3M_Namelist(iID)%iScaleFactor_SWEass

        ! Info start
        call mprintf(.true., iINFO_Extra, ' Data :: SWE gridded :: Binary ... ' )
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Info binary file(s) time step
        call mprintf(.true., iINFO_Verbose, ' Get (SWE gridded) at time '//trim(sTime)//' ... ')
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Snow variable(s)
        !------------------------------------------------------------------------------------------
        ! SWE assimilated (example: SWEAssimilated_201405010000.bin.gz)
        sFileNameData_SWEass = trim(sPathData_SWE_Assimilation)//"SWEAssimilated_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)// &
            ".bin"
        call mprintf(.true., iINFO_Extra, ' Get filename (SWE assimilation gridded): '//trim(sPathData_SWE_Assimilation) )

        ! Checking file input availability
        sFileNameData_SWEass_Zip = trim(sFileNameData_SWEass)//'.gz'
        inquire (file = sFileNameData_SWEass_Zip, exist = bFileExist)
        
        if ( .not. bFileExist ) then

            call mprintf(.true., iWARN, ' Problem opening uncompressed binary file: '// &
                         trim(sFileNameData_SWEass_Zip)//' --> Undefined SWE data values.'// &
                         ' Snow physics is activated! If needed check updating data!')
            a2dVar = -9999.0; 

        else
            ! Unzip file
            call S3M_Tools_Generic_UnzipFile(oS3M_Namelist(iID)%sCommandUnzipFile, &
                                             sFileNameData_SWEass_Zip, &
                                             sFileNameData_SWEass, .true.)
            ! Read binary data
            if (iFlagTypeData_Ass_SWE == 1) then                                                                
                call S3M_Tools_IO_Get2d_Binary_INT(sFileNameData_SWEass, a2dVar, iRows, iCols, iScaleFactor_SWEass, .true., iErr)
            elseif (iFlagTypeData_Ass_SWE == 2) then
                call S3M_Tools_IO_Get2d_Binary_DBL(sFileNameData_SWEass, a2dVar, iRows, iCols, iScaleFactor_SWEass, .true., iErr)
            endif

            ! Remove uncompressed file (to save space on disk)
            call S3M_Tools_Generic_RemoveFile(oS3M_Namelist(iID)%sCommandRemoveFile, &
                                             sFileNameData_SWEass, .false.)               
        endif
        
        a2dVarSWEass = a2dVar
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Info binary file(s) time step
        call mprintf(.true., iINFO_Verbose, ' Get (SWE assimilation gridded) at time '//trim(sTime)//' ... OK')
        ! Info end
        call mprintf(.true., iINFO_Extra, ' Data :: SWE assimilation gridded :: Binary ... OK' )
        !------------------------------------------------------------------------------------------
        
    end subroutine S3M_Data_AssSWE_Gridded_Binary
    !------------------------------------------------------------------------------------------
    
end module S3M_Module_Data_AssSWE_Gridded
!-----------------------------------------------------------------------------------------

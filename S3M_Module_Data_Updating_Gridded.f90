!------------------------------------------------------------------------------------------    
! File:   S3M_Module_Data_Updating_Gridded.f90
! Author(s): Fabio Delogu, Valerio Basso, Francesco Avanzi.
!
! Created on December 19, 2017, 1:19 PM
! Last update on December 15, 2023 10:30 PM
!
! Module to read the snow-depth-assimilation data: SH (snow depth), SQA (quality), SCA, kernel.
!------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------
! Module Header
module S3M_Module_Data_Updating_Gridded
    
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
                                            check2Dvar, getProcessID
                                            
    use S3M_Module_Tools_Time,      only:   S3M_Tools_Time_MonthVal
                             
    ! Implicit none for all subroutines in this module
    implicit none
    !------------------------------------------------------------------------------------------
    
contains
    
    !------------------------------------------------------------------------------------------
    ! Subroutine to read and manage snow-depth-assimilation gridded data
    subroutine S3M_Data_Updating_Gridded_Cpl( iID, sTime, &
                                     iRowsStartL, iRowsEndL, iColsStartL, iColsEndL)
        
        !------------------------------------------------------------------------------------------
        ! Variable(s)
        integer(kind = 4)           :: iID
        
        integer(kind = 4)           :: iRowsStartL, iRowsEndL, iColsStartL, iColsEndL
        integer(kind = 4)           :: iRowsL, iColsL
        integer(kind = 4)           :: iFlagTypeData_Updating
        
        character(len = 19)         :: sTime
        character(len = 12)         :: sTimeMonth
        
        character(len = 256)        :: sPathData_Updating
        
        real(kind = 4), dimension(iRowsEndL - iRowsStartL + 1, &
                                  iColsEndL - iColsStartL + 1) ::   a2dVarSnowCAL, a2dVarSnowQAL, &
                                                                    a2dVarSnowHeightL, a2dVarSnowKernelL
        !------------------------------------------------------------------------------------------
                                                                                                                        
        !------------------------------------------------------------------------------------------
        ! Initialize variable(s)
        a2dVarSnowCAL = -9999.0; a2dVarSnowQAL = -9999.0;
        a2dVarSnowHeightL = -9999.9; a2dVarSnowKernelL = -9999.0; 
        !------------------------------------------------------------------------------------------
                                                                                                
        !------------------------------------------------------------------------------------------
        ! Defining iRows and iCols (Land data)
        iRowsL = iRowsEndL - iRowsStartL + 1
        iColsL = iColsEndL - iColsStartL + 1
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Get global information
        sPathData_Updating = oS3M_Namelist(iID)%sPathData_Updating_Gridded
        iFlagTypeData_Updating = oS3M_Namelist(iID)%iFlagTypeData_Updating_Gridded
                
        ! Info start
        call mprintf(.true., iINFO_Extra, ' Data :: Updating gridded assimilated data... ' )
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Replace general path with specific time feature(s)
        call S3M_Tools_Generic_ReplaceText(sPathData_Updating, '$yyyy', sTime(1:4))
        call S3M_Tools_Generic_ReplaceText(sPathData_Updating, '$mm', sTime(6:7))
        call S3M_Tools_Generic_ReplaceText(sPathData_Updating, '$dd', sTime(9:10))
        call S3M_Tools_Generic_ReplaceText(sPathData_Updating, '$HH', sTime(12:13))
        call S3M_Tools_Generic_ReplaceText(sPathData_Updating, '$MM', sTime(15:16))
        
        ! Checking date
        write(sTimeMonth,'(A,A,A)') sTime(1:4), sTime(6:7), sTime(9:10)
        !------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        ! Check time step (iT)
        if (oS3M_Vars(iID)%iTime .lt. oS3M_Namelist(iID)%iNTime) then
 
            !------------------------------------------------------------------------------------------
            ! Subroutine for reading sequential netCDF updating data 
            if (iFlagTypeData_Updating == 3) then

                !------------------------------------------------------------------------------------------
                ! Call subroutine to get updating data in netCDF format
#ifdef LIB_NC
                call S3M_Data_Updating_Gridded_NC(iID, &
                                        sPathData_Updating, &
                                        iRowsL, iColsL, &
                                        sTime, &
                                        a2dVarSnowCAL, a2dVarSnowQAL, &
                                        a2dVarSnowHeightL, a2dVarSnowKernelL)
#else   
                ! Redefinition of updating data flag (if netCDF library is not linked)
                iFlagTypeData_Updating = 1 
                call mprintf(.true., iWARN, ' '// &
                                            'Updating gridded data type selected was netCDF but library is not linked! '// &
                                            'Data in binary int format will be used!')
#endif
                !------------------------------------------------------------------------------------------
                           
            endif
            !------------------------------------------------------------------------------------------

            !------------------------------------------------------------------------------------------
            ! Subroutine for reading sequential binary updating data
            if (iFlagTypeData_Updating == 1 .or. iFlagTypeData_Updating == 2) then

                !------------------------------------------------------------------------------------------
                ! Calling subroutine to read data in binary format
                call S3M_Data_Updating_Gridded_Binary(iID, iFlagTypeData_Updating, &
                                            sPathData_Updating, &
                                            iRowsL, iColsL, sTime, &
                                            a2dVarSnowCAL, a2dVarSnowQAL, &
                                            a2dVarSnowHeightL, a2dVarSnowKernelL)
                !------------------------------------------------------------------------------------------

            endif
            !------------------------------------------------------------------------------------------
           
            !------------------------------------------------------------------------------------------
            ! Debug
            if (iDEBUG.gt.0) then
                call mprintf(.true., iINFO_Extra, ' ========= UPDATING GRIDDED START =========== ')
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarSnowCAL, oS3M_Vars(iID)%a2iMask, 'SNOWCA START') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarSnowQAL, oS3M_Vars(iID)%a2iMask, 'SNOWQA START') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarSnowHeightL, oS3M_Vars(iID)%a2iMask, 'SNOWHEIGHT START') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarSnowKernelL, oS3M_Vars(iID)%a2iMask, 'SNOWKERNEL START') )
                call mprintf(.true., iINFO_Extra, '')                
            endif
            !------------------------------------------------------------------------------------------
                     
            !------------------------------------------------------------------------------------------
            ! Check variable(s) limits and domain
            a2dVarSnowCAL = check2Dvar(a2dVarSnowCAL,           oS3M_Vars(iID)%a2iMask,     -1.0,   3.0,    -9999.0 )  
            a2dVarSnowQAL = check2Dvar(a2dVarSnowQAL,           oS3M_Vars(iID)%a2iMask,     0.0,    1.0 ,   -9999.0 ) 
            a2dVarSnowHeightL = check2Dvar(a2dVarSnowHeightL,   oS3M_Vars(iID)%a2iMask,     0.0,    10000.0,-9999.0 ) 
            a2dVarSnowKernelL = check2Dvar(a2dVarSnowKernelL,   oS3M_Vars(iID)%a2iMask,     0.0,    1.0,    -9999.0 )
            !------------------------------------------------------------------------------------------
            
        else
            
            !------------------------------------------------------------------------------------------
            ! Extra steps condition
            a2dVarSnowCAL = -9999.0; 
            a2dVarSnowQAL = -9999.0;
            a2dVarSnowHeightL = -9999.0; 
            a2dVarSnowKernelL = -9999.0; 
            
            ! Info message for extra time step(s)
            call mprintf(.true., iINFO_Extra, ' Extra time step ---> Snow-depth-assimilation data are set null')
            call mprintf(.true., iINFO_Extra, ' Extra time step ---> Snow-depth-assimilation data are set to -9999.0')
            !------------------------------------------------------------------------------------------
            
        endif
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Check and save snow-depth-assimilation data to local variable(s) to global workspace
        ! Snow Cover Area
        if ( .not. all(a2dVarSnowCAL.eq.-9999.0) ) then
            oS3M_Vars(iID)%a2dSCA = a2dVarSnowCAL
        else
            oS3M_Vars(iID)%a2dSCA = -9999.0
            
        endif
        !
        ! Snow Quality Assessment
        if ( .not. all(a2dVarSnowQAL.eq.-9999.0) ) then
            oS3M_Vars(iID)%a2dSQA = a2dVarSnowQAL
        else
            oS3M_Vars(iID)%a2dSQA = -9999.0
        endif
        !
        ! Snow Height
        if ( .not. all(a2dVarSnowHeightL.eq.-9999.0) ) then
            oS3M_Vars(iID)%a2dSHeight = a2dVarSnowHeightL
        else
            oS3M_Vars(iID)%a2dSHeight = -9999.0
        endif 
        !
        ! Snow Kernel
        if ( .not. all(a2dVarSnowKernelL.eq.-9999.0) ) then
            oS3M_Vars(iID)%a2dSKernel = a2dVarSnowKernelL
        else
            oS3M_Vars(iID)%a2dSKernel = -9999.0
        endif         
        !
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Debug
        if (iDEBUG.gt.0) then
            call mprintf(.true., iINFO_Extra, '')
            call mprintf(.true., iINFO_Extra, checkvar(oS3M_Vars(iID)%a2dSCA, oS3M_Vars(iID)%a2iMask, 'SNOWCA END') )
            call mprintf(.true., iINFO_Extra, checkvar(oS3M_Vars(iID)%a2dSQA, oS3M_Vars(iID)%a2iMask, 'SNOWQA END') )
            call mprintf(.true., iINFO_Extra, checkvar(oS3M_Vars(iID)%a2dSHeight, oS3M_Vars(iID)%a2iMask, 'SNOWHEIGHT END') )
            call mprintf(.true., iINFO_Extra, checkvar(oS3M_Vars(iID)%a2dSKernel, oS3M_Vars(iID)%a2iMask, 'SNOWKERNEL END') )
            call mprintf(.true., iINFO_Extra, ' ========= UPDATING GRIDDED END =========== ')
        endif
        
        ! Info end
        call mprintf(.true., iINFO_Extra, ' Data :: Updating gridded ... OK' )
        !------------------------------------------------------------------------------------------
        
    end subroutine S3M_Data_Updating_Gridded_Cpl
    !------------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------------
    ! Subroutine to read NC data updating
#ifdef LIB_NC
    subroutine S3M_Data_Updating_Gridded_NC(iID, &
                                  sPathData_Updating, &
                                  iRows, iCols, sTime, &
                                  a2dVarSnowCA, a2dVarSnowQA, &
                                  a2dVarSnowHeight, a2dVarSnowKernel)
                                  
        !------------------------------------------------------------------------------------------
        ! Variable(s)
        integer(kind = 4)                       :: iID                  
        
        character(len = 256), intent(in)        :: sPathData_Updating
        character(len = 700)                    :: sFileNameData_Updating, sFileNameData_Updating_Zip, sFileNameData_Temp
        character(len = 700)                    :: sCommandUnzipFile, sCommandRemoveFile
        character(len = 256)                    :: sVarName
        integer(kind = 4), intent(in)           :: iRows, iCols

        character(len = 19), intent(in)         :: sTime
        character(len = 12)                     :: sTimeMonth
                
        real(kind = 4), dimension(iCols, iRows)                 :: a2dVar
        
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarSnowCA
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarSnowQA
        
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarSnowHeight
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarSnowKernel
       
        character(len = 256):: sVarUnits, sPID
        integer(kind = 4)   :: iErr
        integer(kind = 4)   :: iFileID
        
        logical             :: bFileExist
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Initialize variable(s)
        a2dVarSnowCA = -9999.0; a2dVarSnowCA = -9999.0;
        a2dVarSnowHeight = -9999.0; a2dVarSnowKernel = -9999.0;
        
        sFileNameData_Updating = ''; sFileNameData_Updating_Zip = ''; sFileNameData_Temp = ''; sTimeMonth = ''
        sPID = ''
        
        ! Checking date
        write(sTimeMonth,'(A,A,A)') sTime(1:4), sTime(6:7), sTime(9:10)
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Get global information
        sCommandUnzipFile = oS3M_Namelist(iID)%sCommandUnzipFile
        sCommandRemoveFile = oS3M_Namelist(iID)%sCommandRemoveFile
        
        ! Info start
        call mprintf(.true., iINFO_Extra, ' Data :: Updating gridded :: NetCDF ... ' )
        
        ! Get unique process ID
        sPID = adjustl(getProcessID())
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Filename updating (example: Updating_201404300000.nc.gz)
        sFileNameData_Updating = trim(sPathData_Updating)//"Updating_"// &
        sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
        sTime(12:13)//sTime(15:16)// &
        ".nc"
        ! Create Filename with unique PID number to avoid simultaneously access to the same Forcing file       
        sFileNameData_Temp = trim(sPathData_Updating)//"Updating_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)//'_'//trim(sPID)// &
            ".nc" 

        ! Info netCDF filename
        call mprintf(.true., iINFO_Verbose, ' Get filename (updating gridded): '//trim(sFileNameData_Updating)//' ... ' )
        !------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        ! Checking file input availability
        sFileNameData_Updating_Zip = trim(sFileNameData_Updating)//'.gz'
        inquire (file = trim(sFileNameData_Updating)//'.gz', exist = bFileExist)
        if ( .not. bFileExist ) then
            !------------------------------------------------------------------------------------------
            ! Warning message
            call mprintf(.true., iWARN, ' No compressed updating netCDF data found: '//trim(sFileNameData_Updating_Zip) )
            ! Info netCDF filename
            call mprintf(.true., iINFO_Verbose, &
                         ' Get filename (updating gridded): '//trim(sFileNameData_Updating)//' ... FAILED' )
            ! Info end
            call mprintf(.true., iINFO_Extra, ' Data :: Updating gridded :: NetCDF ... SKIPPED!' )
            !------------------------------------------------------------------------------------------
        else
            
            !------------------------------------------------------------------------------------------
            ! Unzip file
            call S3M_Tools_Generic_UnzipFile(oS3M_Namelist(iID)%sCommandUnzipFile, &
                                             sFileNameData_Updating_Zip, &
                                             sFileNameData_Temp, .true.)
            !------------------------------------------------------------------------------------------
        
            !------------------------------------------------------------------------------------------
            ! Open netCDF file
            iErr = nf90_open(trim(sFileNameData_Temp), NF90_NOWRITE, iFileID)
            if (iErr /= 0) then
                call mprintf(.true., iWARN, ' Problem opening uncompressed netCDF file: '// &
                             trim(sFileNameData_Updating)//' --> Undefined updating data values' )
                call mprintf(.true., iINFO_Verbose, &
                             ' Get filename (updating gridded): '//trim(sFileNameData_Updating)//' ... FAILED' )
            else
                               
                !------------------------------------------------------------------------------------------
                ! Snow variable(s)

                !------------------------------------------------------------------------------------------
                ! SNOW HEIGHT
                sVarName = 'SnowHeight'
                call S3M_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, &
                        ' Get updating gridded data FAILED! Check updating data for '//sVarName//'!'// &
                            ' Snow physics is activated! If needed check updating data!')
                    a2dVarSnowHeight = -9999.0;
                else
                    a2dVarSnowHeight = transpose(a2dVar)
                endif
                !------------------------------------------------------------------------------------------
                
                !------------------------------------------------------------------------------------------
                ! SNOW KERNEL
                sVarName = 'Kernel'
                call S3M_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, &
                        ' Get updating gridded data FAILED! Check updating data for '//sVarName//'!'// &
                        ' Snow physics is activated! If needed check updating data!')
                    a2dVarSnowKernel = -9999.0;
                else
                    a2dVarSnowKernel = transpose(a2dVar)
                endif
                !------------------------------------------------------------------------------------------                    
                    
                !------------------------------------------------------------------------------------------
                ! SNOW COVER AREA
                sVarName = 'SCA'
                call S3M_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, &
                        ' Get updating gridded data FAILED! Check updating data for '//sVarName//'!'// &
                        ' Snow physics is activated! If needed check updating data!')
                    a2dVarSnowCA = -9999.0;
                else
                    a2dVarSnowCA = transpose(a2dVar)
                endif
                !------------------------------------------------------------------------------------------

                !------------------------------------------------------------------------------------------
                ! SNOW QUALITY ASSESSMENT
                sVarName = 'SQA'
                call S3M_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, &
                        ' Get updating gridded data FAILED! Check updating data for '//sVarName//'!'// &
                        ' Snow physics is activated! If needed check updating data!')
                    a2dVarSnowQA = -9999.0;
                else
                    a2dVarSnowQA = transpose(a2dVar)
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
                call mprintf(.true., iINFO_Verbose, ' Get filename (updating gridded): '//trim(sFileNameData_Updating)//' ... OK' )
                ! Info end
                call mprintf(.true., iINFO_Extra, ' Data :: Updating gridded :: NetCDF ... OK' )
                !------------------------------------------------------------------------------------------
                
            endif
            !------------------------------------------------------------------------------------------
            
        endif
        !------------------------------------------------------------------------------------------
        
    end subroutine S3M_Data_Updating_Gridded_NC
#endif
    !------------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------------
    ! Subroutine to read binary updating data
    subroutine S3M_Data_Updating_Gridded_Binary(iID, iFlagTypeData_Updating, &
                                      sPathData_Updating, &
                                      iRows, iCols, sTime, &
                                      a2dVarSnowCA, a2dVarSnowQA, &
                                      a2dVarSnowHeight, a2dVarSnowKernel)
    
        !------------------------------------------------------------------------------------------
        ! Variable(s)
        integer(kind = 4)                   :: iID
        integer(kind = 4)                   :: iFlagTypeData_Updating
                                      
        character(len = 256), intent(in)    :: sPathData_Updating
        character(len = 700)                :: sFileNameData_Updating, sFileNameData_Updating_Zip, sFileNameData_Temp
        character(len = 700)                :: sCommandUnzipFile
        character(len = 256)                :: sVarName
        integer(kind = 4), intent(in)       :: iRows, iCols
                
        character(len = 19), intent(in)     :: sTime
        character(len = 12)                 :: sTimeMonth
        
        real(kind = 4), dimension(iRows, iCols)                 :: a2dVar

        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarSnowCA
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarSnowQA
 
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarSnowHeight
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarSnowKernel        
       
        character(len = 256):: sVarUnits, sPID
        integer(kind = 4)   :: iErr
        integer(kind = 4)   :: iFileID, iScaleFactor_Update
        
        logical             :: bFileExist
        !------------------------------------------------------------------------------------------
	
        !------------------------------------------------------------------------------------------
        ! Initialize variable(s)
        a2dVarSnowCA = -9999.0; a2dVarSnowCA = -9999.0;
        a2dVarSnowKernel = -9999.0; a2dVarSnowHeight = -9999.0;        
                
        sFileNameData_Updating = ''; sFileNameData_Updating_Zip = ''; sFileNameData_Temp = ''; sTimeMonth = ''
        sPID = ''
        
        ! Checking date
        write(sTimeMonth,'(A,A,A)') sTime(1:4), sTime(6:7), sTime(9:10)
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Get global information
        sCommandUnzipFile = oS3M_Namelist(iID)%sCommandUnzipFile
        iScaleFactor_Update = oS3M_Namelist(iID)%iScaleFactor_Update

        ! Info start
        call mprintf(.true., iINFO_Extra, ' Data :: Updating gridded :: Binary ... ' )
        
        ! Get unique process ID
        sPID = adjustl(getProcessID())
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Info binary file(s) time step
        call mprintf(.true., iINFO_Verbose, ' Get (updating gridded) at time '//trim(sTime)//' ... ')
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Snow variable(s)

        !------------------------------------------------------------------------------------------
        ! SnowHeight (example: SnowHeight_201405010000.bin.gz)
        sFileNameData_Updating = trim(sPathData_Updating)//"SnowHeight_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)// &
            ".bin"
        sFileNameData_Temp = trim(sPathData_Updating)//"SnowHeight_"// &
                sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
                sTime(12:13)//sTime(15:16)//'_'//trim(sPID)// &
                ".bin"  
        call mprintf(.true., iINFO_Extra, ' Get filename (updating gridded): '//trim(sPathData_Updating) )

        ! Checking file input availability
        sFileNameData_Updating_Zip = trim(sFileNameData_Updating)//'.gz'
        inquire (file = sFileNameData_Updating_Zip, exist = bFileExist)
        if ( .not. bFileExist ) then
            call mprintf(.true., iWARN, ' Problem opening uncompressed binary file: '// &
                         trim(sFileNameData_Updating_Zip)//' --> Undefined updating data values.'// &
                         ' Snow physics is activated! If needed check updating data!')
            a2dVar = -9999.0; 
        else
            ! Unzip file
            call S3M_Tools_Generic_UnzipFile(oS3M_Namelist(iID)%sCommandUnzipFile, &
                                             sFileNameData_Updating_Zip, &
                                             sFileNameData_Temp, .true.)
                                             
            ! Read binary data
            if (iFlagTypeData_Updating == 1) then
                call S3M_Tools_IO_Get2d_Binary_INT(sFileNameData_Temp, a2dVar, iRows, iCols, iScaleFactor_Update, .true., iErr) 
            elseif (iFlagTypeData_Updating == 2) then
                call S3M_Tools_IO_Get2d_Binary_DBL(sFileNameData_Temp, a2dVar, iRows, iCols, iScaleFactor_Update, .true., iErr) 
            endif
            
            ! Remove uncompressed file (to save space on disk)
            call S3M_Tools_Generic_RemoveFile(oS3M_Namelist(iID)%sCommandRemoveFile, &
                                             sFileNameData_Temp, .false.)
        endif
        a2dVarSnowHeight = a2dVar
        !------------------------------------------------------------------------------------------
            
        !------------------------------------------------------------------------------------------
        ! Snow Kernel (example: Kernel_201405010000.bin.gz)
        sFileNameData_Updating = trim(sPathData_Updating)//"Kernel_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)// &
            ".bin"
        sFileNameData_Temp = trim(sPathData_Updating)//"Kernel_"// &
                sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
                sTime(12:13)//sTime(15:16)//'_'//trim(sPID)// &
                ".bin"  
        call mprintf(.true., iINFO_Extra, ' Get filename (updating gridded): '//trim(sFileNameData_Updating) )

        ! Checking file input availability
        sFileNameData_Updating_Zip = trim(sFileNameData_Updating)//'.gz'
        inquire (file = sFileNameData_Updating_Zip, exist = bFileExist)
        if ( .not. bFileExist ) then
            call mprintf(.true., iWARN, ' Problem opening uncompressed binary file: '// &
                         trim(sFileNameData_Updating_Zip)//' --> Undefined updating data values.'// &
                         ' Snow physics is activated! If needed check updating data!')
            a2dVar = -9999.0; 
        else
            ! Unzip file
            call S3M_Tools_Generic_UnzipFile(oS3M_Namelist(iID)%sCommandUnzipFile, &
                                             sFileNameData_Updating_Zip, &
                                             sFileNameData_Temp, .true.)
                                             
            ! Read binary data
            if (iFlagTypeData_Updating == 1) then
                call S3M_Tools_IO_Get2d_Binary_INT(sFileNameData_Temp, a2dVar, iRows, iCols, iScaleFactor_Update, .true., iErr)
            elseif (iFlagTypeData_Updating == 2) then
                call S3M_Tools_IO_Get2d_Binary_DBL(sFileNameData_Temp, a2dVar, iRows, iCols, iScaleFactor_Update, .true., iErr)
            endif
                
            ! Remove uncompressed file (to save space on disk)
            call S3M_Tools_Generic_RemoveFile(oS3M_Namelist(iID)%sCommandRemoveFile, &
                                              sFileNameData_Temp, .false.)
        endif
        a2dVarSnowKernel = a2dVar
        !------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        ! Snow Cover Area (example: SCA_201405010000.bin.gz)
        sFileNameData_Updating = trim(sPathData_Updating)//"SCA_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)// &
            ".bin"
        sFileNameData_Temp = trim(sPathData_Updating)//"SCA_"// &
                sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
                sTime(12:13)//sTime(15:16)//'_'//trim(sPID)// &
                ".bin"          
        call mprintf(.true., iINFO_Extra, ' Get filename (updating gridded): '//trim(sFileNameData_Updating) )

        ! Checking file input availability
        sFileNameData_Updating_Zip = trim(sFileNameData_Updating)//'.gz'
        inquire (file = sFileNameData_Updating_Zip, exist = bFileExist)
        if ( .not. bFileExist ) then
            call mprintf(.true., iWARN, ' Problem opening uncompressed binary file: '// &
                         trim(sFileNameData_Updating_Zip)//' --> Undefined updating data values.'// &
                         ' Snow physics is activated! If needed check updating data!')
            a2dVar = -9999.0; 
        else
            ! Unzip file
            call S3M_Tools_Generic_UnzipFile(oS3M_Namelist(iID)%sCommandUnzipFile, &
                                             sFileNameData_Updating_Zip, &
                                             sFileNameData_Temp, .true.)
                                             
            ! Read binary data
            if (iFlagTypeData_Updating == 1) then
                call S3M_Tools_IO_Get2d_Binary_INT(sFileNameData_Temp, a2dVar, iRows, iCols, iScaleFactor_Update, .true., iErr)
            elseif (iFlagTypeData_Updating == 2) then
                call S3M_Tools_IO_Get2d_Binary_DBL(sFileNameData_Temp, a2dVar, iRows, iCols, iScaleFactor_Update, .true., iErr)
            endif
                
            ! Remove uncompressed file (to save space on disk)
            call S3M_Tools_Generic_RemoveFile(oS3M_Namelist(iID)%sCommandRemoveFile, &
                                              sFileNameData_Temp, .false.)
        endif
        a2dVarSnowCA = a2dVar
        !------------------------------------------------------------------------------------------
            
        !------------------------------------------------------------------------------------------
        ! Snow Quality Assessment (example: SQA_201405010000.bin.gz)
        sFileNameData_Updating = trim(sPathData_Updating)//"SQA_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)// &
            ".bin"
        sFileNameData_Temp = trim(sPathData_Updating)//"SQA_"// &
                sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
                sTime(12:13)//sTime(15:16)//'_'//trim(sPID)// &
                ".bin"         
        call mprintf(.true., iINFO_Extra, ' Get filename (updating gridded): '//trim(sFileNameData_Updating) )

        ! Checking file input availability
        sFileNameData_Updating_Zip = trim(sFileNameData_Updating)//'.gz'
        inquire (file = sFileNameData_Updating_Zip, exist = bFileExist)
        if ( .not. bFileExist ) then
            call mprintf(.true., iWARN, ' Problem opening uncompressed binary file: '// &
                         trim(sFileNameData_Updating_Zip)//' --> Undefined updating data values.'// &
                         ' Snow physics is activated! If needed check updating data!')
            a2dVar = -9999.0; 
        else
            ! Unzip file
            call S3M_Tools_Generic_UnzipFile(oS3M_Namelist(iID)%sCommandUnzipFile, &
                                             sFileNameData_Updating_Zip, &
                                             sFileNameData_Temp, .true.)
            ! Read binary data
            if (iFlagTypeData_Updating == 1) then
                call S3M_Tools_IO_Get2d_Binary_INT(sFileNameData_Temp, a2dVar, iRows, iCols, iScaleFactor_Update, .true., iErr)
            elseif (iFlagTypeData_Updating == 2) then
                call S3M_Tools_IO_Get2d_Binary_DBL(sFileNameData_Temp, a2dVar, iRows, iCols, iScaleFactor_Update, .true., iErr)
            endif
                
            ! Remove uncompressed file (to save space on disk)
            call S3M_Tools_Generic_RemoveFile(oS3M_Namelist(iID)%sCommandRemoveFile, &
                                              sFileNameData_Temp, .false.)
        endif
        a2dVarSnowQA = a2dVar            
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Info binary file(s) time step
        call mprintf(.true., iINFO_Verbose, ' Get (updating gridded) at time '//trim(sTime)//' ... OK')
        ! Info end
        call mprintf(.true., iINFO_Extra, ' Data :: Updating gridded :: Binary ... OK' )
        !------------------------------------------------------------------------------------------
        
    end subroutine S3M_Data_Updating_Gridded_Binary
    !------------------------------------------------------------------------------------------
    
end module S3M_Module_Data_Updating_Gridded
!-----------------------------------------------------------------------------------------

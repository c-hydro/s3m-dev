!--------------------------------------------------------------------------------  
! File:   S3M_Module_Info_Time.f90
! Author(s): Fabio Delogu, Francesco Silvestro, Simone Gabellani.
!
! Created on November, 06 2015, 8:42 AM
!
! Module to get info time
!--------------------------------------------------------------------------------

!--------------------------------------------------------------------------------
! Module Header
module S3M_Module_Info_Time
    
    !--------------------------------------------------------------------------------
    ! External module(s) for all subroutine in this module
#ifdef LIB_NC
    use netcdf
#endif

    use S3M_Module_Namelist,        only:   oS3M_Namelist
    use S3M_Module_Vars_Loader,     only:   oS3M_Vars
    
    use S3M_Module_Tools_Debug
    
#ifdef LIB_NC
    use S3M_Module_Tools_IO,        only:   S3M_Tools_IO_Get2d_NC, &
                                            S3M_Tools_IO_GetArcGrid_ASCII, &
                                            check
#else
    use S3M_Module_Tools_IO,        only:   S3M_Tools_IO_GetArcGrid_ASCII                                    
#endif
                                            
    ! Implicit none for all subroutines in this module
    implicit none
    !--------------------------------------------------------------------------------

contains 

    !------------------------------------------------------------------------------------------
    ! Subroutine to get time dims
    subroutine S3M_Info_Time_GetDims(iID, iDaySteps)
        
        !------------------------------------------------------------------------------------
        ! Variable(s) declaration
        integer(kind = 4)               :: iID
        integer(kind = 4)               :: iFileID, iErr
         
        integer(kind = 4)               :: iDtModel
        integer(kind = 4)               :: iNData, iNTime
        integer(kind = 4)               :: iDaySteps
        
        logical                         :: bFileExist
        
        character(len = 256)            :: sStrDaySteps
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Variable(s) initialization
        bFileExist = .false.    
        iDaySteps = 0;
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Get information
        iDtModel = oS3M_Namelist(iID)%iDtModel
        
        iNData = oS3M_Namelist(iID)%iNData
        iNTime = oS3M_Namelist(iID)%iNTime
        
        ! Info
        call mprintf(.true., iINFO_Main, ' Define time dims ... ')
        !------------------------------------------------------------------------------------
 
        !------------------------------------------------------------------------------------
        ! Day steps
        iDaySteps = 24*3600/oS3M_Namelist(iID)%iDtData_Forcing
        !------------------------------------------------------------------------------------
              
        !------------------------------------------------------------------------------------------
        ! Pass local variable(s) to global workspace
        oS3M_Namelist(iID)%iDaySteps = iDaySteps
        
        ! Time info
        write(sStrDaySteps, *) iDaySteps;
        call mprintf(.true., iINFO_Main, ' TIME INFO --- '// &
                'DaySteps: '//sStrDaySteps)
        
        ! Info
        call mprintf(.true., iINFO_Main, ' Define time dims ... OK')
        !------------------------------------------------------------------------------------------
        
    end subroutine S3M_Info_Time_GetDims
    !------------------------------------------------------------------------------------------
    
end module S3M_Module_Info_Time
!------------------------------------------------------------------------------------------
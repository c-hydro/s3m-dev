!------------------------------------------------------------------------------------------    
! File:   S3M_Module_Phys.f90
! Author: Fabio Delogu.
!
! Created on April 2, 2014, 5:19 PM
! Last update on Sep 06, 2018 11:00 AM
!
! Snow physics module
!------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------
! Module Header
module S3M_Module_Phys
    
    !------------------------------------------------------------------------------------------
    ! External module(s) for all subroutine in this module
    use S3M_Module_Namelist,            only:   oS3M_Namelist
    use S3M_Module_Vars_Loader,         only:   oS3M_Vars
    
    use S3M_Module_Tools_Debug
    
    use S3M_Module_Phys_Snow,           only:   S3M_Phys_Snow_Cpl

    ! Implicit none for all subroutines in this module
    implicit none
    !------------------------------------------------------------------------------------------
    
contains

    !------------------------------------------------------------------------------------------
    ! Subroutine to run model physics
    subroutine S3M_Phys_Cpl(iID, &
                            iRowsStart, iRowsEnd, iColsStart, iColsEnd, &
                            iTime, iNTime, sTime, &
                            iNData, &
                            iDaySteps)
        !------------------------------------------------------------------------------------------
        ! Variable(s)
        integer(kind = 4)           :: iID
        integer(kind = 4)           :: iRows, iRowsStart, iRowsEnd, iCols, iColsStart, iColsEnd
        integer(kind = 4)           :: iNData
        integer(kind = 4)           :: iTime, iNTime
        integer(kind = 4)           :: iDaySteps
        
        character(len = 19)         :: sTime
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Defining iRows and iCols
        iRows = iRowsEnd - iRowsStart + 1
        iCols = iColsEnd - iColsStart + 1
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Subroutine to compute SNOW
        call S3M_Phys_Snow_Cpl(iID, iRows, iCols)
        !------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        ! Info message for extra time step(s)
        call mprintf(.true., iINFO_Extra, ' Extra time step ---> Swow routine(s) are skipped!')
        !------------------------------------------------------------------------------------------
        
    end subroutine S3M_Phys_Cpl
    !------------------------------------------------------------------------------------------
    
end module S3M_Module_Phys
!------------------------------------------------------------------------------------------ù


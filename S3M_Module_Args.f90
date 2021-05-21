!--------------------------------------------------------------------------------  
! File:   S3M_Module_Args.f90
! Author(s): Fabio Delogu, Francesco Silvestro, Simone Gabellani, Francesco Avanzi
!
! Created on February 15, 2017, 4:40 PM
! Last update on Oct 27, 2020 11:45 AM
!
! Module to read argument(s) defined by command-line
!--------------------------------------------------------------------------------

!--------------------------------------------------------------------------------
! Module Args
module S3M_Module_Args
    
    !--------------------------------------------------------------------------------
    ! External module(s) and implicit none
    use S3M_Module_Tools_Debug          ! to import global variable(s) declaration

    implicit none
    !--------------------------------------------------------------------------------
    
contains 
    
    !--------------------------------------------------------------------------------
    ! Subroutine to read argument(s)
    subroutine S3M_Args_Read(sFileInfo, iArgsType) 
        
        !------------------------------------------------------------------------------------------
        ! Variable(s)
        integer(kind = 4)               :: iArgsN, iArgsType
        character(len = 700)            :: sLineBuffer
        character(len = 700)            :: sFileInfo
        !------------------------------------------------------------------------------------------                    
        
        !------------------------------------------------------------------------------------------
        ! Variable(s) initialization
        iArgsN = -9999; 
        iArgsType = -9999
        sFileInfo = "";
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Get number of argument(s)
        iArgsN = iargc()   
        !------------------------------------------------------------------------------------------
        
        if (iArgsN == 1) then
            
            !------------------------------------------------------------------------------------------
            ! Get information file name 
            call getarg(1, sLineBuffer); sFileInfo = sLineBuffer

            ! Define argument(s) type
            iArgsType = 2
            !------------------------------------------------------------------------------------------
        
        else
            
            call mprintf(.true., iERROR, 'Too many arguments in the command line, see S3M_Args_Read subroutine')
            
        endif
        !------------------------------------------------------------------------------------------
        
    end subroutine S3M_Args_Read
    !--------------------------------------------------------------------------------
    
end module S3M_Module_Args
!--------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
! File:      S3M_Module_Vars_Loader.f90
! Author(s): Fabio Delogu, Francesco Silvestro, Simone Gabellani.
!
! Created on March, 4 2015, 9:57 AM
!
! Module to import global variable(s)
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
! Module Header
module S3M_Module_Vars_Loader
    
    !------------------------------------------------------------------------------------
    ! External module(s) and implicit none
    implicit none
    integer, parameter :: iMaxDomain = 1
    !------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------
    ! Defining global variables type
    include "S3M_Type_Vars.inc"
    type (S3M_Type_Vars), dimension (iMaxDomain) :: oS3M_Vars
    save oS3M_Vars
    !------------------------------------------------------------------------------------
    
end module S3M_Module_Vars_Loader
!------------------------------------------------------------------------------------
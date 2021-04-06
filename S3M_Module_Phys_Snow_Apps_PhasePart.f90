!------------------------------------------------------------------------------------
! File:   S3M_Module_Phys_Snow_Apps_PhasePart.f90
! Author:   Simone Gabellani, Fabio Delogu, Francesco Avanzi.
!
! Created on June 10, 2020 5:00 PM
! Last update on October 27, 2020 11:06 AM
!
! Precipitation phase-partitioning application module
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
! Module Header
module S3M_Module_Phys_Snow_Apps_PhasePart

    !------------------------------------------------------------------------------------
    ! External module(s) 
    use S3M_Module_Namelist,            only:   oS3M_Namelist
    use S3M_Module_Vars_Loader,         only:   oS3M_Vars
   
    use S3M_Module_Tools_Debug
    
    
    ! Implicit none for all subroutines in this module
    implicit none
    !------------------------------------------------------------------------------------------

contains

    !------------------------------------------------------------------------------------------
    ! Subroutine to compute precipitation phase-partitioning based on Froidurot et al. (2014)
    ! DOI: https://doi.org/10.1175/JHM-D-13-073.1
    subroutine S3M_Phys_Snow_Apps_PhasePart_Froidurot(iID, iRows, iCols, &
                                            a2dVarDEM, &
                                            a2dVarTa, a2dVarRelHum, a2dVarPrecip, a2dVarSepCoeff, a2dVarSnowFall, &
                                            a2dVarRainfall)

        !------------------------------------------------------------------------------------------
        ! Variable(s) declaration 
        integer(kind = 4)   :: iID, iRows, iCols
               
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarDem
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarTa, a2dVarRelHum, a2dVarPrecip
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarSepCoeff, a2dVarSnowFall, a2dVarRainfall    
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------ 
        ! Info start
        call mprintf(.true., iINFO_Extra, ' Phys :: Snow :: Phase Part Froidurot ... ' )
        !------------------------------------------------------------------------------------------ 
        
        !------------------------------------------------------------------------------------------        
        ! Compute phase-separation coefficient
        where(a2dVarDEM.ge.0.0)
            a2dVarSepCoeff = 1/(1 + exp(22 - 2.7*a2dVarTa - 0.2*a2dVarRelHum))
        endwhere   
          
        ! Compute snowfall and rainfall from total precipitation
        where( (a2dVarDEM.ge.0.0) .and. (a2dVarPrecip.gt.0.0) )
            a2dVarSnowFall = (1 - a2dVarSepCoeff)*a2dVarPrecip 
            a2dVarRainfall = a2dVarSepCoeff*a2dVarPrecip  
        elsewhere( (a2dVarDEM.ge.0.0) )
            a2dVarSnowFall = 0.0 
            a2dVarRainfall = 0.0
        endwhere
            
        where(a2dVarSnowFall.lt.0.01) a2dVarSnowFall = 0.0        
        where(a2dVarRainfall.lt.0.01) a2dVarRainfall = 0.0  
        !------------------------------------------------------------------------------------------
            
        !------------------------------------------------------------------------------------------
        ! Info end
        call mprintf(.true., iINFO_Extra, ' Phys :: Snow :: Phase Part Froidurot... OK' )
        !------------------------------------------------------------------------------------------
        
    end subroutine S3M_Phys_Snow_Apps_PhasePart_Froidurot
    !------------------------------------------------------------------------------------------   


end module S3M_Module_Phys_Snow_Apps_PhasePart
!------------------------------------------------------------------------------------

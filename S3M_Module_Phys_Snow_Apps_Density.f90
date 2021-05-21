!------------------------------------------------------------------------------------
! File:   S3M_Module_Phys_Snow_Apps_Density.f90
! Author:   Simone Gabellani, Fabio Delogu, Francesco Avanzi.
!
! Created on June 10, 2020 5:00 PM
! Last update on October 26, 2020 5:04 PM
!
! Snow-density applications module
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
! Module Header
module S3M_Module_Phys_Snow_Apps_Density

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
    ! Subroutine to compute snow density according to De Michele et al. 2013: https://tc.copernicus.org/articles/7/433/2013/
    subroutine S3M_Phys_Snow_Apps_Rho(iID, iRows, iCols, &
                                       sTime, iTime, iDtForcing, &
                                       a2dVarDem, &
                                       a2dVarTa, a2dVarSnowFall, a2dVarSWE_D, &
                                       a2dVarH_D, a2dVarRho_D, a2dVarRhoS0)

        !------------------------------------------------------------------------------------------
        ! Variable(s) declaration 
        integer(kind = 4)   :: iID, iRows, iCols
        integer(kind = 4)   :: iTime, iDtForcing
        character(len = 19) :: sTime
       
        real(kind = 4)      :: dVarRhoW, dVarRhoSMax, dVarRhoSMin, dVarRhoS0Max 
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarDem
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarTa, a2dVarSnowFall, a2dVarSWE_D
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarRhoS0, a2dVarRho_D, a2dVarH_D
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarSnowTemp
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Initialize variable(s)
        dVarRhoSMax = -9999.0; dVarRhoSMin = -9999.0; dVarRhoS0Max = -9999.0
        dVarRhoW = -99990; a2dVarSnowTemp = -9999.0
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Get information
        dVarRhoSMax = oS3M_Namelist(iID)%dRhoSnowMax
        dVarRhoS0Max = oS3M_Namelist(iID)%dRhoSnowFresh
        dVarRhoSMin = oS3M_Namelist(iID)%dRhoSnowMin
        dVarRhoW = oS3M_Namelist(iID)%dRhoW  
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Info start
        call mprintf(.true., iINFO_Extra, ' Phys :: Snow :: Rho ... ' )
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Compute fresh-snow density
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSnowFall.gt.0.0) )
            a2dVarRhoS0 = 67.9 + 51.3*exp(a2dVarTa/2.6)
        endwhere
        
        ! Check limits for fresh-snow density
        where( (a2dVarDem.ge.0.0) .and. (a2dVarRhoS0.gt.dVarRhoS0Max) )
            a2dVarRhoS0 = dVarRhoS0Max
        endwhere
        
        where( (a2dVarDem.ge.0.0) .and. (a2dVarRhoS0.lt.dVarRhoSMin) )
            a2dVarRhoS0 = dVarRhoSMin
        endwhere
        
        ! Set to 0 where there's no snowfall
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSnowFall.eq.0.0) )
            a2dVarRhoS0 = 0.0
        endwhere   
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Update bulk-snow density to account for fresh-snow density                  
        where( (a2dVarDem.ge.0.0) .and. &
            ((a2dVarSWE_D - a2dVarSnowFall).gt.1.0) .and. (a2dVarSnowFall.gt.1.0) .and. (a2dVarRho_D.gt.dVarRhoSMin) )

               a2dVarRho_D = (a2dVarSWE_D)/((a2dVarSnowFall/a2dVarRhoS0) + ((a2dVarSWE_D - a2dVarSnowFall)/a2dVarRho_D))
               ! This equation assumes that a2dVarSWE_D already includes a2dVarSnowFall, per an update that is supposed to 
               ! take place in Module_Phys_Snow right after the computation of snowfall and rainfall. Accordingly, we must
               ! subtract a2dVarSnowFall from a2dVarSWE_D here to correctly weight pre-existing a2dVarRho_D
        endwhere
        ! This first "where" has updated dry-snow density where DEM > 0, where there were at least some snow on the ground, 
        ! and some substantial snowfall. In all other cases, a2dVarRho_D has not changed yet. Where there is little snowfall
        ! over a large snowpack, this is OK, we assume fresh snow density plays a minimal role. But we need to include some 
        ! provisions when there is very little snow on the ground, for example at the beginning of the snow season... 
        
        where( (a2dVarDem.ge.0.0) .and. ((a2dVarSWE_D - a2dVarSnowFall).le.1.0) .and. (a2dVarSnowFall.gt.0.0) )
            a2dVarRho_D = a2dVarRhoS0
        endwhere
        ! Now we made sure that a2dVarRho_D has been updated wherever there was sone snowfall and no pre-existing SWE. 
        
        ! Check rhos limit(s)
        where( (a2dVarDem.ge.0.0) .and. (a2dVarRho_D.gt.dVarRhoSMax) )
            a2dVarRho_D = dVarRhoSMax
        endwhere
        where( (a2dVarDem.ge.0.0) .and. (a2dVarRho_D.lt.dVarRhoSMin) )
            a2dVarRho_D = dVarRhoSMin
        endwhere
        !------------------------------------------------------------------------------------------ 
        
        !------------------------------------------------------------------------------------------        
        ! Compute updated dry-snow height    
        where( (a2dVarDem.ge.0.0) .and. (a2dVarRho_D.gt.0.0) )
            a2dVarH_D = ((a2dVarSWE_D/1000)*dVarRhoW)/a2dVarRho_D
        elsewhere(a2dVarDem.ge.0.0) 
            a2dVarH_D = 0.0
        endwhere
        !------------------------------------------------------------------------------------------ 
        
        !------------------------------------------------------------------------------------------        
        ! Compute linearly approximated snow temperature         
        where( (a2dVarDem.ge.0.0) .and. (a2dVarTa.ge.0.0) )
            a2dVarSnowTemp = 0.0
        elsewhere( (a2dVarDem.ge.0.0) ) 
            a2dVarSnowTemp = 0.5*a2dVarTa
        endwhere
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------        
        ! Compute dry-snow-density compaction now          
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE_D.gt.0.0) )  
            a2dVarRho_D = a2dVarRho_D + 0.66*(iDtForcing/3600)*0.001*a2dVarH_D*(a2dVarRho_D**2)*EXP(0.08*a2dVarSnowTemp - &
                                                    0.021*a2dVarRho_D)
        endwhere 
        
        ! Check rhos limit(s)
        where( (a2dVarDem.ge.0.0) .and. (a2dVarRho_D.gt.dVarRhoSMax) )
            a2dVarRho_D = dVarRhoSMax
        endwhere
        where( (a2dVarDem.ge.0.0) .and. (a2dVarRho_D.lt.dVarRhoSMin) )
            a2dVarRho_D = dVarRhoSMin
        endwhere
        !------------------------------------------------------------------------------------------    
        
        !------------------------------------------------------------------------------------------        
        ! Re-compute updated dry-snow height    
        where( (a2dVarDem.ge.0.0) .and. (a2dVarRho_D.gt.0.0) )
            a2dVarH_D = ((a2dVarSWE_D/1000)*dVarRhoW)/a2dVarRho_D
        elsewhere(a2dVarDem.ge.0.0) 
            a2dVarH_D = 0.0
        endwhere       
        !------------------------------------------------------------------------------------------
              
        !------------------------------------------------------------------------------------------
        ! Info end
        call mprintf(.true., iINFO_Extra, ' Phys :: Snow :: Rho ... OK' )
        !------------------------------------------------------------------------------------------
        
    end subroutine S3M_Phys_Snow_Apps_Rho
    !------------------------------------------------------------------------------------------    

end module S3M_Module_Phys_Snow_Apps_Density
!------------------------------------------------------------------------------------
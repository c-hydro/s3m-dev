!------------------------------------------------------------------------------------
! File:   S3M_Module_Phys_Snow_Apps_Hydraulics.f90
! Author:   Francesco Avanzi
!
! Created on June 10, 2020 5:00 PM
! Last update on October 27, 2020 10:45 AM
!
! Snow-hydraulics application module
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
! Module Header
module S3M_Module_Phys_Snow_Apps_Hydraulics

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
    ! Subroutine to compute snowpack outflow based on Darcy's 
    !(Avanzi et al. 2015, https://doi.org/10.1016/j.advwatres.2015.09.021)
    subroutine S3M_Phys_Snow_Apps_Outflow(iID, iRows, iCols, iDtForcing, &
        dVarRhoW, &
        a2dVarDem, &
        a2dVarSWE_D, a2dVarRho_D, a2dVarSWE_W, & 
        a2dVarRainfall, a2dVarMeltingS, a2dVarRefreezingS, &
        a2dVarOutflow_K, a2dVarH_D, a2dVarH_S, a2dVarSSA, a2dVarPerm, a2dVarTheta_W)    

        !------------------------------------------------------------------------------------------
        ! Variable(s) declaration 
        integer(kind = 4)   :: iID, iRows, iCols, iDtForcing
        real(kind = 4)      :: dVarRhoW
        
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarDem
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarTheta_W, a2dVarSWE_D, a2dVarRho_D, a2dVarSWE_W
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarOutflow_K, a2dVarRainfall, a2dVarMeltingS, a2dVarRefreezingS
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarPorosity, a2dVarH_D, a2dVarH_S 
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarSr, a2dVarSr_irr, a2dVarSr_star 
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarSSA, a2dVarr_e, a2dVarPerm, a2dVarCond
        !------------------------------------------------------------------------------------------
                
        !------------------------------------------------------------------------------------------
        ! Initialize variable(s)
        a2dVarPorosity = -9999; a2dVarSr = -9999; a2dVarSr_irr = -9999; a2dVarSr_star = -9999; a2dVarr_e = -9999;
        a2dVarCond = -9999; 

        ! Info start
        call mprintf(.true., iINFO_Extra, ' Phys :: Snow :: Outflow ... ' )
        !------------------------------------------------------------------------------------------   
        
        !------------------------------------------------------------------------------------------        
        ! Compute porosity
        where( (a2dVarDem.ge.0.0) .and. (a2dVarH_D.gt.0.0))
            a2dVarPorosity = 1 - a2dVarRho_D/917        
        elsewhere(a2dVarDem.ge.0.0)    
            a2dVarPorosity = 0.0
        endwhere
        !------------------------------------------------------------------------------------------ 
        
        !------------------------------------------------------------------------------------------        
        ! Compute dry-snow depth
        where( (a2dVarDem.ge.0.0) .and. (a2dVarRho_D.gt.0.0) )    
            a2dVarH_D = ((a2dVarSWE_D/1000)*dVarRhoW)/a2dVarRho_D
        elsewhere     
            a2dVarH_D = 0.0
        endwhere
        !------------------------------------------------------------------------------------------  

        !------------------------------------------------------------------------------------------        
        ! Compute control volume and saturation degree        
        where( (a2dVarDem.ge.0.0) .and. (a2dVarH_D.gt.0.0) & 
            .and. ((a2dVarSWE_W/1000 - (a2dVarPorosity*a2dVarH_D)) .ge. 0.0) ) !we divide by 1000 to convert mm to m 
                                                            !and so compare a2dVarSWE_W with a2dVarH_D
            
            a2dVarH_S = a2dVarH_D + ((a2dVarSWE_W/1000) - (a2dVarPorosity*a2dVarH_D))
            a2dVarSr = 1.0
        
        elsewhere( (a2dVarDem.ge.0.0) .and. (a2dVarH_D.gt.0.0) )
            
            a2dVarH_S = a2dVarH_D
            a2dVarSr = (a2dVarSWE_W/1000)/(a2dVarPorosity*a2dVarH_D)
        
        elsewhere( (a2dVarDem.ge.0.0) )
                
            a2dVarH_S = 0.0
            a2dVarSr = 0.0
        
        endwhere
        
        where( (a2dVarDem.ge.0.0) .and. (a2dVarPorosity.gt.0.0))
            
            a2dVarSr_irr = 0.02*((a2dVarRho_D/dVarRhoW)/a2dVarPorosity)
            a2dVarSr_star = (a2dVarSr - a2dVarSr_irr)/(1 - a2dVarSr_irr)
            
        elsewhere 
            
            a2dVarSr_irr = 0.0
            a2dVarSr_star = 0.0
    
        endwhere
        !------------------------------------------------------------------------------------------ 
        
        !------------------------------------------------------------------------------------------        
        ! Compute SSA, r_e, permeability, and conductivity
        where( (a2dVarDem.ge.0.0) .and. (a2dVarRho_D.gt.0.0))

            a2dVarSSA = -308.2*LOG(a2dVarRho_D/1000) - 206 ! SSA here will be in cm2/g as in Domine et al 07. 
                                                           ! a2dVarRho_D must be in g/cm3
            a2dVarSSA = a2dVarSSA/10 !now a2dVarSSA is in m2/kg
            a2dVarr_e = 3/(a2dVarSSA*917)
            a2dVarPerm = 3*(a2dVarr_e**2)*EXP((-0.013*a2dVarRho_D))
            a2dVarCond = a2dVarPerm*(a2dVarSr_star**3)
        endwhere    
            
        where(a2dVarSr .lt. a2dVarSr_irr) a2dVarCond = 0.0
        ! Here we force a2dVarCond to 0 where a2dVarSr lt a2dVarSr_irr; the implication is the a2dVarOutflow will be 0            
        !------------------------------------------------------------------------------------------ 
            
        !------------------------------------------------------------------------------------------        
        ! Compute a2dVarOutflow_K, but some conditions apply....
        ! Note that a2dVarRainfall has been already added to SWE_W at this point.
            
        where( (a2dVarSr.ge.0.5) .or. (a2dVarSWE_D .lt. 10) ) 

                a2dVarOutflow_K = a2dVarSWE_W
                ! This is a pragmatic solution that is particularly helpful given that S3M is employed for flood forecasting: 
                ! wherever saturation is extremely high OR the snowpack is very shallow, we assume snow has no retention capability 
                ! and all wet SWE (here possibly including a2dVarRainfall already, although in principle it should 
                ! have been allocated to a2dVarOutflow_ExcessRain) is immediately discharged. 

        elsewhere(a2dVarSWE_W .le. (5.47*(10**5)*a2dVarCond*(iDtForcing)*1000))     
                ! 5.47*(10**5) includes viscosity and density. 

                a2dVarOutflow_K = a2dVarSWE_W
                
        elsewhere    
                
                a2dVarOutflow_K = 5.47*(10**5)*a2dVarCond*(iDtForcing)*1000
                ! 5.47*(10**5) includes viscosity and density. 

        endwhere
            
        where( a2dVarOutflow_K.le.0.0 )  a2dVarOutflow_K = 0.0 !just in case....
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Info end
        call mprintf(.true., iINFO_Extra, ' Phys :: Snow :: Phase Part OUtflow... OK' )
        !------------------------------------------------------------------------------------------
        
    end subroutine S3M_Phys_Snow_Apps_Outflow
    !------------------------------------------------------------------------------------------    


end module S3M_Module_Phys_Snow_Apps_Hydraulics
!------------------------------------------------------------------------------------
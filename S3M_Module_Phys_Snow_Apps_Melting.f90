!------------------------------------------------------------------------------------
! File:   S3M_Module_Phys_Snow_Apps_Melting.f90
! Author:   Simone Gabellani, Fabio Delogu, Daniele Dolia, Francesco Avanzi.
!
! Created on June 10, 2020 5:00 PM
! Last update on October 27, 2020 10:56 AM
!
! Snow-melt application module
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
! Module Header
module S3M_Module_Phys_Snow_Apps_Melting
    
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
    ! Subroutine to compute snow albedo
    subroutine S3M_Phys_Snow_Apps_Albedo(iID, iRows, iCols, &
                                         sTime, iTime, iGlacierValue, &
                                         a2dVarDem, a2iVarGlacierMask, &
                                         a2dVarTaC_MeanDays1, a2iVarAgeS, &
                                         a2dVarAlbedoS, iFlagIceMassBalance, a2dVarSWE_D, a2dVarIceThick_WE)

        !------------------------------------------------------------------------------------------
        ! Variable(s) declaration 
        integer(kind = 4)   :: iID, iRows, iCols, iK, iI, iJ

        integer(kind = 4)   :: iTime, iGlacierValue, iFlagIceMassBalance
        character(len = 19) :: sTime
        
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarDem
        
        integer(kind = 4), dimension(iRows, iCols)      :: a2iVarAgeS, a2iVarGlacierMask, a2iIndex_for_pivot_table
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarTaC_MeanDays1, a2dVarAlbedoS, a2dVarSWE_D, a2dVarIceThick_WE, &
                                                           a2dVarAlbedoS_old
        real(kind = 4), dimension(500)                  :: a1dAlbedoPivotDry, a1dAlbedoPivotWet
        integer(kind = 4), dimension(500)               :: a1iAlbedoPivotDays
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        !Initialization
        iI = 0; iJ = 0; iK = 0;
        a1iAlbedoPivotDays = (/ (iK, iK = 1, 500) /) !We initialize this as an array of days, assuming 500 days are enough
        a1dAlbedoPivotDry = -9999.0; a1dAlbedoPivotWet = -9999.0; a2dVarAlbedoS_old = -9999.0; a2iIndex_for_pivot_table = -9999;
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Info start
        call mprintf(.true., iINFO_Extra, ' Phys :: Snow :: Albedo ... ' )
        !------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        ! Update snow albedo once each day
        ! Laramie and Schaake 1972.
        if ( (sTime(12:13).eq.'23') ) then
            
            ! Compute pivot albedo arrays 
            a1dAlbedoPivotWet = 0.5 + 0.45*exp(dble(-a1iAlbedoPivotDays)*0.12)
            a1dAlbedoPivotDry = 0.5 + 0.45*exp(dble(-a1iAlbedoPivotDays)*0.05)
            
            !we store current albedo values BEFORE update
            a2dVarAlbedoS_old = a2dVarAlbedoS
            
            !Now we update albedo
            do iI = 1, iRows
                do iJ = 1, iCols
                    
                    !WET-CONDITION UPDATE
                    if( (a2dVarDem(iI,iJ).ge.0.0) .and. (a2dVarTaC_MeanDays1(iI,iJ).gt.0.0) ) then
                    !We do this in two steps: we first compute where we stand in the pivot table; then update albedo accordingly.
                        
                        a2iIndex_for_pivot_table(iI,iJ) = MINLOC(ABS(a1dAlbedoPivotWet - a2dVarAlbedoS_old(iI,iJ) ), DIM = 1)
                        
                        if( a2iIndex_for_pivot_table(iI,iJ) .lt. SIZE(a1dAlbedoPivotWet, DIM = 1) ) then
                            a2dVarAlbedoS(iI,iJ) = a1dAlbedoPivotWet(a2iIndex_for_pivot_table(iI,iJ) + 1)
                        endif
                        
                    endif
                    
                    !DRY-CONDITION UPDATE
                    if( (a2dVarDem(iI,iJ).ge.0.0) .and. (a2dVarTaC_MeanDays1(iI,iJ).le.0.0) ) then
                    !We do this in two steps: we first compute where we stand in the pivot table; then update albedo accordingly.
                        
                        a2iIndex_for_pivot_table(iI,iJ) = MINLOC(ABS(a1dAlbedoPivotDry - a2dVarAlbedoS_old(iI,iJ) ), DIM = 1)
                        
                        if( a2iIndex_for_pivot_table(iI,iJ) .lt. SIZE(a1dAlbedoPivotDry, DIM = 1) ) then
                            a2dVarAlbedoS(iI,iJ) = a1dAlbedoPivotDry(a2iIndex_for_pivot_table(iI,iJ) + 1)
                        endif
                        
                    endif                    
                    
                enddo
            enddo
            
            ! Check snow albedo boundaries conditions
            where( (a2dVarDem.ge.0.0) .and. (a2dVarAlbedoS.le.0.5) )
                a2dVarAlbedoS = 0.5
            endwhere
            
            where( (a2dVarDem.ge.0.0) .and. (a2dVarAlbedoS.gt.0.95) )
                a2dVarAlbedoS = 0.95  
            endwhere
            
            ! Update albedo where Age is 0!
            where( (a2dVarDem.ge.0.0) .and. (a2iVarAgeS.eq.0) )
                a2dVarAlbedoS = 0.95 
            endwhere            
                        
            ! Set snow albedo glaciers condition
            if ( (iFlagIceMassBalance.eq.1.0) .or. (iFlagIceMassBalance.eq.2.0) ) then
                where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE_D.eq.0.0) .and. (a2dVarIceThick_WE.gt.0.0) ) 
                    a2dVarAlbedoS = 0.4
                endwhere

            else
                where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE_D.eq.0.0) .and. (a2iVarGlacierMask.eq.iGlacierValue) ) 
                    a2dVarAlbedoS = 0.4
                endwhere
            endif     
            
        endif
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Info end
        call mprintf(.true., iINFO_Extra, ' Phys :: Snow :: Albedo ... OK' )
        !------------------------------------------------------------------------------------------

    end subroutine S3M_Phys_Snow_Apps_Albedo
    !------------------------------------------------------------------------------------------    
        
    !------------------------------------------------------------------------------------------
    ! Subroutine to compute snow age
    subroutine S3M_Phys_Snow_Apps_Age(iID, iRows, iCols, &
                                       sTime, iTime, &
                                       a2dVarDem, &
                                       a2dVarSnowFallDayCum, a2dVarSWE, &
                                       a2iVarAgeS)

        !------------------------------------------------------------------------------------------
        ! Variable(s) declaration 
        integer(kind = 4)   :: iID, iRows, iCols

        integer(kind = 4)   :: iTime
        character(len = 19) :: sTime
        
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarDem
        
        integer(kind = 4), dimension(iRows, iCols)      :: a2iVarAgeS
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarSnowFallDayCum, a2dVarSWE
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Info start
        call mprintf(.true., iINFO_Extra, ' Phys :: Snow-Surface :: Age ... ' )
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Initialize variable(s) every day at 23.00
        if ( (sTime(12:13).eq.'23') .and. (iTime.gt.1) ) then
            ! Snow age updating
            where( (a2dVarDem.ge.0.0) .and. (a2dVarSnowFallDayCum.le.3.0) )
                a2iVarAgeS = a2iVarAgeS + 1
            endwhere
            
            where( (a2dVarDem.ge.0.0) .and. (a2dVarSnowFallDayCum.gt.3.0) )
                a2iVarAgeS = 0
            endwhere
            
            ! Snow age re-initialization with SWE equal to 0.0 (without snow --> no data value)
            where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.eq.0.0) )
                a2iVarAgeS = 0
            endwhere

        endif
        !-----------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        ! Info end
        call mprintf(.true., iINFO_Extra, ' Phys :: Snow-Surface :: Age ... OK' )
        !------------------------------------------------------------------------------------------
        
    end subroutine S3M_Phys_Snow_Apps_Age
    !------------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------------
    ! Subroutine to compute average temperature over n days
    subroutine S3M_Phys_Snow_Apps_TMean(iID, iRows, iCols, &
                                        iDaySteps, &    
                                        a2iVarMask, &
                                        a3dVarT, &
                                        a2dVarT, a2dVarTMean)
                                        
        !------------------------------------------------------------------------------------------
        ! Variable(s) declaration 
        integer(kind = 4)   :: iID, iRows, iCols
        
        integer(kind = 4)   :: iStep, iDaySteps
        
        integer(kind = 4), dimension(iRows, iCols)          :: a2iVarMask
        real(kind = 4), dimension(iRows, iCols)             :: a2dVarT, a2dVarTMean
        real(kind = 4), dimension(iRows, iCols, iDaySteps)  :: a3dVarT
        !------------------------------------------------------------------------------------------           
        
        !------------------------------------------------------------------------------------------
        ! Info start
        call mprintf(.true., iINFO_Extra, ' Phys :: Snow :: TMean ... ' )
        !------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        ! Compute and update temperature 3d mean field(s)
        if (all(a3dVarT.eq.0.0))then

            ! Debug
            call mprintf(.true., iINFO_Extra, ' Phys :: Snow :: TMean over days :: '// &
                                              ' First mean temperature 3d field storing step ... ')

            ! Update with a starting field all temporal steps
            do iStep = 1, int(iDaySteps)
                where(a2iVarMask.gt.0.0)
                    a3dVarT(:,:,int(iStep)) = a2dVarT
                elsewhere
                    a3dVarT(:,:,int(iStep)) = 0.0
                endwhere
            enddo

            call mprintf(.true., iINFO_Extra, ' Phys :: Snow :: TMean over days :: '// &
                                              ' First mean temperature 3d field storing step ... OK')
        else
            ! Re-initialize temperature 3d array
            do iStep=2, int(iDaySteps)
                a3dVarT(:,:,int(iStep-1)) = a3dVarT(:,:,int(iStep))
            enddo

            ! Update with new field
            where(a2iVarMask.gt.0.0)
                a3dVarT(:,:,int(iDaySteps)) =  a2dVarT
            elsewhere
                a3dVarT(:,:,int(iDaySteps)) = 0.0
            endwhere

        endif

        ! Calculate mean temperature over n days
        where(a2iVarMask.gt.0.0)
            a2dVarTMean = sum(a3dVarT(:,:,1:int(iDaySteps)),3)/(iDaySteps)
        elsewhere
            a2dVarTMean = 0.0
        endwhere
        !------------------------------------------------------------------------------------------ 

        !------------------------------------------------------------------------------------------
        ! Info end
        call mprintf(.true., iINFO_Extra, ' Phys :: Snow :: TMean ... OK' )
        !------------------------------------------------------------------------------------------
                                        
    end subroutine S3M_Phys_Snow_Apps_TMean
    !------------------------------------------------------------------------------------------
        
   !------------------------------------------------------------------------------------------
    ! Subroutine to compute open-loop melting-refreezing
    subroutine S3M_Phys_Snow_Apps_MeltingRefr(iID, iRows, iCols, iDtForcing, iDaySteps1Days, &
                                                sTime, iTime, &
                                                iGlacierValue, dVarRhoW, dVarMeltingTRef, dVarIceMeltingCoeff, &
                                                a2dVarDem, a2iVarGlacierMask, a2dVarArctUp, &
                                                a2dVarTa, a2dVarIncRad, &
                                                a2dVarAlbedoS, a2dVarSWE_D, a2dVarMeltingS, a2dVarMeltingSc, &
                                                a2dVarSWE_W, a2dVarRefreezingS, dVarRefreezingSc, & 
                                                iFlagIceMassBalance, a2dVarIceThick_WE, a2dVarMeltingG, a2dVarMeltSIncRad, &
                                                a2dVarMeltSTemp, a2dVarMeltGIncRad, a2dVarMeltGTemp, &
                                                a2dVarTaC_MeanDays10, dVarModFactorRadS, a2dVarMeltingSRadc, iFlagGlacierDebris)
        
        !------------------------------------------------------------------------------------------
        ! Variable(s) declaration 
        integer(kind = 4)   :: iID, iRows, iCols
        
        integer(kind = 4)   :: iTime, iDaySteps1Days, iDtForcing
        integer(kind = 4)   :: iGlacierValue, iFlagIceMassBalance, iFlagGlacierDebris
        real(kind = 4)      :: dVarRhoW, dVarMeltingTRef, dVarLamba, dVarIceMeltingCoeff, dVarRefreezingSc, dVarModFactorRadS, &
                               dDebrisThreshold
        character(len = 19) :: sTime
        
        integer(kind = 4), dimension(iRows, iCols)      :: a2iVarGlacierMask, a2iVarMask
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarDem
        
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarTa, a2dVarIncRad
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarMeltingSc, a2dVarMeltingScG, a2dVarMeltingSRadc, a2dVarMeltingGRadc
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarAlbedoS, a2dVarMeltingS, a2dVarRefreezingS, a2dVarMeltingG, &
                                                           a2dVarMeltGIncRad, a2dVarMeltGTemp, &
                                                           a2dVarMeltSIncRad, a2dVarMeltSTemp, a2dVarGlacierDebris, &
                                                           a2dVarTaC_MeanDays10
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarSWE_D, a2dVarSWE_W     
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarIceThick_WE   
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarArctUp
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Initialize variable(s)
        a2dVarMeltingSc = 0.0; a2dVarMeltingScG = 0.0; a2dVarMeltingSRadc = 0.0; a2dVarMeltingGRadc = 0.0;
        a2dVarGlacierDebris = oS3M_Vars(iID)%a2dGlacierDebris
        dDebrisThreshold = oS3M_Namelist(iID)%dDebrisThreshold
        a2iVarMask = oS3M_Vars(iID)%a2iMask

        ! Fusion latent heat [MJ/kg] 
        dVarLamba = 0.334	
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Info start
        call mprintf(.true., iINFO_Extra, ' Phys :: Snow :: MeltingOL ... ' )
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Compute snow melting first coefficient
        where ( (a2dVarDem.ge.0.0) )
            a2dVarMeltingSc = 0.598862*atan(0.27439*a2dVarTaC_MeanDays10 - 0.5988)-0.598862*3.14/2 + a2dVarArctUp 
            a2dVarMeltingSRadc = 0.49338*atan(0.27439*a2dVarTaC_MeanDays10 - 0.5988) -0.49338*3.14/2 + dVarModFactorRadS
        endwhere
        
        ! Melting coefficient lower limit
        where(a2dVarMeltingSc.lt.0.0) a2dVarMeltingSc = 0.0 
        where(a2dVarMeltingSRadc.lt.0.0) a2dVarMeltingSRadc = 0.0 
        
        ! Compute snow melting second coefficient (for glacier(s))
        a2dVarMeltingScG = dVarIceMeltingCoeff*a2dVarMeltingSc !It is recommended 
        ! NOT to increase the DDF of ice. This means setting dVarIceMeltingCoeff = 1 in the namelist.
        ! The reason is that increasing the DDF of snow on ice surfaces should take into account a decrease in albedo, 
        ! which we are already accounting for elsewhere. So a2dVarMeltingScG = a2dVarMeltingScd
        a2dVarMeltingGRadc = a2dVarMeltingSRadc
        ! We also assume that the melt coefficient for radiation is the same for both snow and ice, given again that the main 
        ! difference between snow and ice is in albedo, which we already capture in a2dVarAlbedoS. 
        !------------------------------------------------------------------------------------------
         
        !------------------------------------------------------------------------------------------
        ! Compute snow melting
        where ( (a2dVarDem.ge.0.0) .and. (a2dVarTa.ge.dVarMeltingTRef) &
                .and. (a2dVarSWE_D.gt.0.0) .and. (a2dVarTaC_MeanDays10.ge.dVarMeltingTRef) )
                
            a2dVarMeltSIncRad = (1000.0/(dVarRhoW*dVarLamba))*((a2dVarIncRad*(1.0 - a2dVarAlbedoS))*(iDtForcing)/1000000)
            ! Assuming a2dVarIncRad is in W/m2 and iDtForcing is in seconds, then a2dVarIncRad is first converted 
            ! in MJ/m2 and integrated across the timestep, then it is converted in mm w.e. using density and latent heat.
            a2dVarMeltSIncRad = a2dVarMeltingSRadc*a2dVarMeltSIncRad
            
            a2dVarMeltSTemp = a2dVarMeltingSc*(a2dVarTa - dVarMeltingTRef)/iDaySteps1Days
            ! a2dVarMeltingSc*(a2dVarTa - dVarMeltingTRef)/iDaySteps1Days is 
            ! originally in mm w.e./d (a2dVarMeltingSc*(a2dVarTa - dVarMeltingTRef)) and we then convert to mm w.e. 
            ! over the time step by
            ! dividing by iDaySteps1Days 
            
            a2dVarMeltingS = a2dVarMeltSIncRad + a2dVarMeltSTemp
            
        elsewhere
            
            a2dVarMeltingS = 0.0
            a2dVarMeltSIncRad = 0.0
            a2dVarMeltSTemp = 0.0
            
        endwhere
        
        ! Check snow melting limit
        where (a2dVarMeltingS.lt.0.01) a2dVarMeltingS = 0.0
        where (a2dVarMeltingS.lt.0.01) a2dVarMeltSIncRad = 0.0
        where (a2dVarMeltingS.lt.0.01) a2dVarMeltSTemp = 0.0
        !------------------------------------------------------------------------------------------
            
        !------------------------------------------------------------------------------------------
        ! Compute glacier melting
        if ( (iFlagIceMassBalance.eq.1.0) .or. (iFlagIceMassBalance.eq.2.0) ) then
            
            where ( (a2dVarDem.ge.0.0) .and. (a2dVarTa.ge.dVarMeltingTRef) & 
                    .and. (a2dVarSWE_D.eq.0.0) .and. (a2dVarIceThick_WE.gt.0.0) )
                    
                a2dVarMeltGIncRad = (1000.0/(dVarRhoW*dVarLamba))*(a2dVarIncRad*(1.0 - a2dVarAlbedoS))*(iDtForcing)/1000000
                ! Assuming a2dVarIncRad is in W/m2 and iDtForcing is in seconds, then a2dVarIncRad is first converted 
                ! in MJ/m2 and integrated across the timestep, then it is converted in mm w.e. using density and latent heat.
                a2dVarMeltGIncRad = a2dVarMeltingGRadc*a2dVarMeltGIncRad
            
                a2dVarMeltGTemp = a2dVarMeltingScG*(a2dVarTa - dVarMeltingTRef)/iDaySteps1Days
                ! a2dVarMeltingScG*(a2dVarTa - dVarMeltingTRef)/iDaySteps1Days is 
                ! originally in mm w.e./d (a2dVarMeltingScG*(a2dVarTa - dVarMeltingTRef)) and we then convert to mm w.e. 
                ! over the time step by
                ! dividing by iDaySteps1Days                
                a2dVarMeltingG = a2dVarMeltGIncRad + a2dVarMeltGTemp
            
            elsewhere
                
                a2dVarMeltingG = 0.0
                a2dVarMeltGTemp = 0.0
                a2dVarMeltGIncRad = 0.0
                
            endwhere
            
        else
            
            where ( (a2dVarDem.ge.0.0) .and. (a2dVarTa.ge.dVarMeltingTRef) & 
                    .and. (a2dVarSWE_D.eq.0.0) .and. (a2iVarGlacierMask.eq.iGlacierValue) )
                    
                a2dVarMeltGIncRad = (1000.0/(dVarRhoW*dVarLamba))*(a2dVarIncRad*(1.0 - a2dVarAlbedoS))*(iDtForcing)/1000000
                ! Assuming a2dVarIncRad is in W/m2 and iDtForcing is in seconds, then a2dVarIncRad is first converted 
                ! in MJ/m2 and integrated across the timestep, then it is converted in mm w.e. using density and latent heat.
                a2dVarMeltGIncRad = a2dVarMeltingGRadc*a2dVarMeltGIncRad
            
                a2dVarMeltGTemp = a2dVarMeltingScG*(a2dVarTa - dVarMeltingTRef)/iDaySteps1Days
                ! a2dVarMeltingScG*(a2dVarTa - dVarMeltingTRef)/iDaySteps1Days is 
                ! originally in mm w.e./d (a2dVarMeltingScG*(a2dVarTa - dVarMeltingTRef)) and we then convert to mm w.e. 
                ! over the time step by
                ! dividing by iDaySteps1Days                
                a2dVarMeltingG = a2dVarMeltGIncRad + a2dVarMeltGTemp
            
            elsewhere
                
                a2dVarMeltingG = 0.0  
                a2dVarMeltGTemp = 0.0
                a2dVarMeltGIncRad = 0.0                
            
            endwhere
            
        endif      
        
        !Apply glacier-debris mask if needed
        if (iFlagGlacierDebris .eq. 1) then
            
            where( (a2dVarDem.ge.0.0) .and. (a2dVarGlacierDebris .gt. dDebrisThreshold))
            
            a2dVarMeltingG = (1 - a2dVarGlacierDebris)*a2dVarMeltingG
            a2dVarMeltGTemp = (1 - a2dVarGlacierDebris)*a2dVarMeltGTemp
            a2dVarMeltGIncRad = (1 - a2dVarGlacierDebris)*a2dVarMeltGIncRad
            
            endwhere

        endif
        
        where (a2dVarMeltingG.lt.0.01) a2dVarMeltingG = 0.0        
        where (a2dVarMeltingG.lt.0.01) a2dVarMeltGTemp = 0.0        
        where (a2dVarMeltingG.lt.0.01) a2dVarMeltGIncRad = 0.0  
        !------------------------------------------------------------------------------------------
            
        !------------------------------------------------------------------------------------------
        ! Compute snow refreezing
        where ( (a2dVarDem.ge.0.0) .and. (a2dVarTa.lt.dVarMeltingTRef) &
                .and. (a2dVarSWE_W.gt.0.0) )
                
            a2dVarRefreezingS = -dVarRefreezingSc*a2dVarMeltingSc*(a2dVarTa - dVarMeltingTRef)/iDaySteps1Days
            ! here, we compute refreezing using a simplified approach following Avanzi et al. 2015.
            ! https://doi.org/10.1016/j.advwatres.2015.09.021
            ! We assume a2dVarRefreezingS to be in mm w.e.
            ! Given that refreezing is computed with negative temperatures, we need extra care to make sure sign is correct.
            ! a2dVarTa - dVarMeltingTRef is negative if 2dVarTa.lt.dVarMeltingTRef, which is why we have included a 
            ! minus sign at the beginning of this equation. 
            
        elsewhere
            
            a2dVarRefreezingS = 0.0
            
        endwhere
        
        ! Check snow & ice melting limit
        where (a2dVarRefreezingS.lt.0.01) a2dVarRefreezingS = 0.0            
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Info end
        call mprintf(.true., iINFO_Extra, ' Phys :: Snow :: MeltingOL ... OK' )
        !------------------------------------------------------------------------------------------
              
    end subroutine S3M_Phys_Snow_Apps_MeltingRefr
    !------------------------------------------------------------------------------------------        

end module S3M_Module_Phys_Snow_Apps_Melting
!------------------------------------------------------------------------------------
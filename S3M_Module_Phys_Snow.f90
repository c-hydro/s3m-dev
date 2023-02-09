!------------------------------------------------------------------------------------
! File:   S3M_Module_Phys_Snow.f90
! Author:   Francesco Avanzi, Fabio Delogu, Simone Gabellani.
!
! Created on Jul 15, 2015 11:00 AM
! Last update on February 09, 2023 10:30 AM
!
! Snow- and glacier-physics subroutine for S3M
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
! Module Header
module S3M_Module_Phys_Snow
    
    !------------------------------------------------------------------------------------
    ! External module(s)
    use S3M_Module_Namelist,                    only:  oS3M_Namelist
    use S3M_Module_Vars_Loader,                 only:  oS3M_Vars
    
    use S3M_Module_Tools_Debug
    
    use S3M_Module_Phys_Snow_Apps_Assimilation, only:  S3M_Phys_Snow_Apps_AssimSH, &
                                                       S3M_Phys_Snow_Apps_AssimSWE
                                                
    use S3M_Module_Phys_Snow_Apps_Density,      only:  S3M_Phys_Snow_Apps_Rho  
                                                   
    use S3M_Module_Phys_Snow_Apps_Hydraulics,   only:  S3M_Phys_Snow_Apps_Outflow  
    
    use S3M_Module_Phys_Snow_Apps_Melting,      only:  S3M_Phys_Snow_Apps_TMean, &
                                                       S3M_Phys_Snow_Apps_Age, &
                                                       S3M_Phys_Snow_Apps_Albedo, &
                                                       S3M_Phys_Snow_Apps_MeltingRefr

    use S3M_Module_Phys_Snow_Apps_PhasePart,    only:  S3M_Phys_Snow_Apps_PhasePart_Froidurot
    
    use S3M_Module_Phys_Snow_Apps_Glaciers,     only:  S3M_Phys_Snow_Apps_GlacierDeltaH

    use gnufor2 
    ! Implicit none for all subroutines in this module
    implicit none
    !------------------------------------------------------------------------------------------

contains

    !------------------------------------------------------------------------------------------
    ! Subroutine to calculate snow- and glacier-physics
    subroutine S3M_Phys_Snow_Cpl(iID, iRows, iCols)
        
        !-------------------------------------------------------------------------------------
        ! Variable(s) declaration
        integer(kind = 4)   :: iID, iRows, iCols
        integer(kind = 4)   :: iRows_Pivot, iCols_Pivot 
        integer(kind = 4)   :: iDaySteps, iDaySteps1Days, iTime, iDtDataForcing, iDt, iHour 
        integer(kind = 4)   :: iGlacierValue 
        integer(kind = 4)   :: iFlagSnowAssim, iFlagSnowAssim_SWE, iFlagIceMassBalance, iFlagGlacierDebris, iFlagAssOnlyPos
        integer(kind = 4)   :: iDaysAvgTSuppressMelt 
        
        real(kind = 4)      :: dVarRhoW 
        real(kind = 4)      :: dVarMeltingTRef, dVarIceMeltingCoeff, dVarRefreezingSc, dVarModFactorRadS
        
        character(len = 19) :: sTime
        character(len = 2) :: sWYstart
        
        integer(kind = 4), dimension(iRows, iCols)      :: a2iVarMask, a2iVarGlacierMask 
        integer(kind = 4), dimension(iRows, iCols)      :: a2iVarAgeS 

        !topography
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarDEM
        
        !weather inputs
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarTa, a2dVarRelHum, a2dVarIncRad, a2dVarPrecip
        
        !assimilation variables: snow-depth and SCA information
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarSnowHeight, a2dVarSnowKernel, a2dVarSnowCA, a2dVarSnowQA
        
        !assimilation variables: SWE information & data-assimilation-related variables
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarSWEass, a2dVarSWE_preass, a2dVarUass   
        
        !parameters
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarArctUp 
        
        !variables: snowfall and rainfall
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarSnowFall, a2dVarSnowFallDayCum, a2dVarRainfall
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarSepCoeff
        
        !variables: dry and wet snow
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarSWE_D, a2dVarRho_D, a2dVarH_D, a2dVarSWE_W, a2dVarTheta_W         
        
        !variables: bulk snow        
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarSWE, a2dVarRhoS, a2dVarRhoS0, a2dVarH_S, a2dVarSnowMask
        
        !variables: snow melt and refreezing
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarMeltingS, a2dVarMeltingSDayCum, a2dVarMeltingSc, a2dVarMeltingSRadc 
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarMeltSIncRad, a2dVarMeltSTemp, a2dVarMeltGIncRad, a2dVarMeltGTemp        
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarRefreezingS, a2dVarAlbedoS        
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarTaC_MeanDays1, &
                                                           a2dVarTaC_MeanDaysSuppressMelt
        !variables: snow hydraulics    
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarOutflow, a2dVarSSA, a2dVarPerm, a2dVarReff, &
                                                           a2dVarOutflow_ExcessRain, a2dVarOutflow_ExcessMelt, a2dVarOutflow_K

        !variables: glaciers
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarIceThick, a2dVarIceThick_WE, a2dVarIceFlag, a2dVarMeltingG
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarMeltingGCumWY , a2dVarChangeThickness      
        !-------------------------------------------------------------------------------------
        
        !-------------------------------------------------------------------------------------
        ! Local variable(s) initialization  
        ! Topography        
        a2iVarGlacierMask = -9999; a2dVarDEM = -9999.0; 
        
        ! Parameters
        a2dVarArctUp = -9999.0; 
        
        ! Weather inputs
        a2dVarPrecip = -9999.0; a2dVarTa = -9999.0; a2dVarRelHum = -9999.0; a2dVarIncRad = -9999.0; 
        
        ! Assimilation variables: snow-depth and SCA information
        a2dVarSnowHeight = -9999.0; a2dVarSnowKernel = -9999.0; a2dVarSnowCA = -9999.0; a2dVarSnowQA = -9999.0;
        ! Assimilation variables: SWE information & data-assimilation-related variables
        a2dVarSWEass = -9999.0; a2dVarSWE_preass = -9999.0; a2dVarUass = -9999.0;
        
        ! Variables: snowfall and rainfall
        a2dVarSnowFall = 0.0; a2dVarSnowFallDayCum = -9999; a2dVarRainfall = 0.0;
        a2dVarSepCoeff = -9999.0; 
        ! Variables: dry and wet snow
        a2dVarSWE_D = -9999.0; a2dVarRho_D = -9999.0; a2dVarH_D = -9999.0; a2dVarSWE_W = -9999.0; a2dVarTheta_W  = -9999.0; 
        ! Variables: bulk snow  
        a2dVarSWE = -9999.0; a2dVarRhoS = -9999.0; a2dVarRhoS0 = -9999; a2dVarH_S = -9999.0;      
        a2dVarSnowMask = -9999.0;
        ! Variables: snow melt and refreezing
        a2dVarMeltingS = 0; a2dVarMeltingSDayCum = -9999.0; a2dVarMeltingSc = -9999.0; a2dVarMeltingSRadc = -9999.0;
        a2dVarAlbedoS = -9999.0; a2iVarAgeS = 0; 
        dVarRefreezingSc = -9999.0; dVarModFactorRadS = -9999.0;
        a2dVarRefreezingS = 0;
        a2dVarMeltingG = 0;
        a2dVarTaC_MeanDays1 = -9999.0; a2dVarTaC_MeanDaysSuppressMelt = -9999.0;
        a2dVarMeltSIncRad = 0.0; a2dVarMeltSTemp = 0.0; a2dVarMeltGIncRad = 0.0; a2dVarMeltGTemp = 0.0;
        ! Variables: snow hydraulics   
        a2dVarOutflow = 0; a2dVarSSA = -9999.0; a2dVarPerm = -9999.0; a2dVarReff = 0.0;
        a2dVarOutflow_ExcessRain = 0.0; a2dVarOutflow_ExcessMelt = 0.0; a2dVarOutflow_K = 0.0;
        ! Variables: glaciers
        a2dVarIceThick = -9999.0; a2dVarIceThick_WE = -9999.0; a2dVarIceFlag = 1.0; a2dVarMeltingGCumWY = 0.0;
        a2dVarChangeThickness = 0.0;
        !-------------------------------------------------------------------------------------                         
                                     
        !-------------------------------------------------------------------------------------   
        ! Read values from oS3M_Namelist
        ! Flags
        iFlagSnowAssim = oS3M_Namelist(iID)%iFlagSnowAssim 
        iFlagSnowAssim_SWE = oS3M_Namelist(iID)%iFlagSnowAssim_SWE 
        iFlagIceMassBalance = oS3M_Namelist(iID)%iFlagIceMassBalance 
        iFlagGlacierDebris = oS3M_Namelist(iID)%iFlagGlacierDebris
        iFlagAssOnlyPos = oS3M_Namelist(iID)%iFlagAssOnlyPos 
        
        ! Constants
        dVarRhoW = oS3M_Namelist(iID)%dRhoW 
        iGlacierValue = oS3M_Namelist(iID)%iGlacierValue 
        sWYstart = oS3M_Namelist(iID)%sWYstart
        iRows_Pivot = oS3M_Namelist(iID)%iRowsPivot  
        iCols_Pivot = oS3M_Namelist(iID)%iColsPivot     
        
        ! Snow Parameters
        dVarMeltingTRef = oS3M_Namelist(iID)%dMeltingTRef 
        dVarIceMeltingCoeff = oS3M_Namelist(iID)%dIceMeltingCoeff 
        dVarRefreezingSc = oS3M_Namelist(iID)%dRefreezingSc 
        dVarModFactorRadS = oS3M_Namelist(iID)%dModFactorRadS
        iDaysAvgTSuppressMelt = oS3M_Namelist(iID)%iDaysAvgTSuppressMelt        
        
        ! Time
        iDaySteps1Days = oS3M_Namelist(iID)%iDaySteps ! Define days steps [-]
        iDtDataForcing = oS3M_Namelist(iID)%iDtData_Forcing ! Dt data forcing
        !------------------------------------------------------------------------------------- 
        
        !-------------------------------------------------------------------------------------  
        ! Static variable(s)
        a2dVarDEM = oS3M_Vars(iID)%a2dDem
        a2iVarMask = oS3M_Vars(iID)%a2iMask
        a2iVarGlacierMask = oS3M_Vars(iID)%a2iGlacierMask
        a2dVarArctUp = oS3M_Vars(iID)%a2dArctUp
        sTime = oS3M_Vars(iID)%sTimeStep ! Time information
        iTime = oS3M_Vars(iID)%iTime ! Time information
        read (sTime(12:13),*) iHour ! Time information
        !------------------------------------------------------------------------------------- 
        
        !-------------------------------------------------------------------------------------         
        ! Extracting dynamic forcing variable(s)
        a2dVarPrecip = oS3M_Vars(iID)%a2dPrecip
        a2dVarTa = oS3M_Vars(iID)%a2dTa
        a2dVarRelHum = oS3M_Vars(iID)%a2dRHum        
        if (iHour .ge. int(7) .and. iHour.le.int(19)) then 
            a2dVarIncRad = oS3M_Vars(iID)%a2dK
        else
            a2dVarIncRad = 0.0
        endif
        !------------------------------------------------------------------------------------- 

        !------------------------------------------------------------------------------------- 
        ! Extracting snow-depth and SWE assimilation variable(s)    
        a2dVarSnowHeight = oS3M_Vars(iID)%a2dSHeight
        a2dVarSnowKernel = oS3M_Vars(iID)%a2dSKernel
        a2dVarSnowCA = oS3M_Vars(iID)%a2dSCA
        a2dVarSnowQA = oS3M_Vars(iID)%a2dSQA
        a2dVarSWEass = oS3M_Vars(iID)%a2dSWEass
        !------------------------------------------------------------------------------------- 

        !-------------------------------------------------------------------------------------         
        !Re-initialize snowfall and melting cumulated variable(s)
        if ( (sTime(12:13).eq.'00') ) then
            ! Update cumulative daily variables
            oS3M_Vars(iID)%a2dMeltingDayCum = 0.0
            oS3M_Vars(iID)%a2dSnowFallDayCum = 0.0
        endif
        !-------------------------------------------------------------------------------------         

        !-------------------------------------------------------------------------------------         
        ! Extracting state variable(s)
        ! Here we only extract those variables that depend on previous time-step values. Other fluxes are only defined and 
        ! initialized above. a2dVarSWE is an exception: it could be computed just as a2dVarSWE_D + a2dVarSWE_W, but we need it 
        ! right away to perform a first set of sanity checks. So we simply take it from previous time step.
        a2dVarSWE_D = oS3M_Vars(iID)%a2dSWE_D
        a2dVarRho_D = oS3M_Vars(iID)%a2dRho_D
        a2dVarSWE_W = oS3M_Vars(iID)%a2dSWE_W
        a2dVarSWE = oS3M_Vars(iID)%a2dSWE
        a2iVarAgeS = oS3M_Vars(iID)%a2iAge
        a2dVarMeltingSDayCum = oS3M_Vars(iID)%a2dMeltingDayCum 
        a2dVarSnowFallDayCum = oS3M_Vars(iID)%a2dSnowFallDayCum 
        a2dVarAlbedoS = oS3M_Vars(iID)%a2dAlbedo_Snow  
        a2dVarTaC_MeanDaysSuppressMelt = oS3M_Vars(iID)%a2dTaC_MeanDaysSuppressMelt
        a2dVarMeltingGCumWY = oS3M_Vars(iID)%a2dMeltingGCumWY
        
        ! Check ice melting flag and variable(s)
        if ( (iFlagIceMassBalance.eq.1.0) .or. (iFlagIceMassBalance.eq.2.0) ) then
            
            a2dVarIceThick = oS3M_Vars(iID)%a2dIceThick
            
            ! Conversion in mm of water equivalent for calculation
            where(a2dVarDEM.gt.0.0)
                a2dVarIceThick_WE = a2dVarIceThick*917
            endwhere
        endif
        !------------------------------------------------------------------------------------- 
        
        !---------------------------------------------------------------------------------------- 
        ! Info start
        call mprintf(.true., iINFO_Verbose, ' Phys :: Snow model ... ' )
        !----------------------------------------------------------------------------------------       
        
        !-----------------------------------------------------------------------------------------
        ! Debug
        if (iDEBUG.gt.0) then
            call mprintf(.true., iINFO_Extra, ' ========= SNOW START =========== ')  
            ! METEO
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarPrecip, a2iVarMask, 'PRECIP START ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarTa, a2iVarMask, 'TA START ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarRelHum, a2iVarMask, 'RELHUM START ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarIncRad, a2iVarMask, 'INCRAD START ') )
            
            !ASSIMILATION
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarSnowHeight, a2iVarMask, 'SNOWH START ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarSnowKernel, a2iVarMask, 'SNOWK START ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarSnowCA, a2iVarMask, 'SNOWCA START ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarSnowQA, a2iVarMask, 'SNOWQA START ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarSWEass, a2iVarMask, 'SWEASS START ') )
                
            !SNOWFALL AND RAINFALL
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarSepCoeff, a2iVarMask, 'SEPCOEFF START ') )            
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarSnowFall, a2iVarMask, 'SNOWFALL START ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarSnowFallDayCum, a2iVarMask, 'SNOWFALL DAYCUM START ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarRainFall, a2iVarMask, 'RAINFALL START ') )

            !DRY WET SNOW
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarSWE_D, a2iVarMask, 'SWE_D START ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarRho_D, a2iVarMask, 'RHO_D START ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarSWE_W, a2iVarMask, 'SWE_W START ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarTheta_W, a2iVarMask, 'THETA_W START ') )              
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarH_D, a2iVarMask, 'HD START ') )             

            !BULK SNOW
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarSWE, a2iVarMask, 'SWE START ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarRhoS, a2iVarMask, 'RHOS START ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarH_S, a2iVarMask, 'HS START ') )             
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarSnowMask, a2iVarMask, 'SNOWMASK START ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarRhoS0, a2iVarMask, 'RHOS0 START ') )

            !SNOWMELT AND REFREEZING
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarMeltingS, a2iVarMask, 'MELTINGS START ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarMeltingSDayCum, a2iVarMask, 'MELTINGS DAYCUM START ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarMeltingSc, a2iVarMask, 'MELTINGS COEFF START ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarAlbedoS, a2iVarMask, 'ALBEDOS START ') )
            call mprintf(.true., iINFO_Extra, checkvar(real(a2iVarAgeS), a2iVarMask, 'AGES START ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarRefreezingS, a2iVarMask, 'REFREEZING START ') ) 
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarMeltingG, a2iVarMask, 'MELTINGG START ') ) 
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarMeltSIncRad, a2iVarMask, 'MELTSINCRAD START ') ) 
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarMeltSTemp, a2iVarMask, 'MELTSTEMP START ') ) 
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarMeltGIncRad, a2iVarMask, 'MELTGINCRAD START ') ) 
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarMeltGTemp, a2iVarMask, 'MELTGTEMP START ') )   
            
            !SNOW HYDRAULICS
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarOutflow, a2iVarMask, 'OUTFLOW START ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarSSA, a2iVarMask, 'SSA START ') ) 
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarPerm, a2iVarMask, 'PERMEABILITY START ') )   
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarReff, a2iVarMask, 'REFF START ') )              
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarOutflow_ExcessRain, a2iVarMask, 'OUTFLOW EXCESS RAIN START ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarOutflow_ExcessMelt, a2iVarMask, 'OUTFLOW EXCESS MELT START ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarOutflow_K, a2iVarMask, 'OUTFLOW PERMEABILITY START ') )
            
            !GLACIERS
            if ( (iFlagIceMassBalance.eq.1.0) .or. (iFlagIceMassBalance.eq.2.0) ) then
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarIceThick_WE, a2iVarMask, 'ICE W.E. START ') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarIceThick, a2iVarMask, 'ICE START ') )                
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarMeltingGCumWY, a2iVarMask, 'ICE MELT CUM START ') )  
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarChangeThickness, a2iVarMask, 'ICE CHANGE THICKNESS ') ) 
            endif
            
            call mprintf(.true., iINFO_Extra, ' ') 
        endif
        !-----------------------------------------------------------------------------------------
        
        !-----------------------------------------------------------------------------------------
        ! Check variables
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarSWE = 0.0 ! just in case
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarSWE_D = 0.0 ! just in case  
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarRho_D = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarH_D = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarSWE_W = 0.0 ! just in case              
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarTheta_W = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarRhoS = 0.0 ! just in case 
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarH_S = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2iVarAgeS = 0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarSSA = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarPerm = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarIceThick_WE.lt.0.01) ) a2dVarIceThick_WE = 0.0 ! just in case                 
        !-----------------------------------------------------------------------------------------
            
        !-----------------------------------------------------------------------------------------  
        ! Call subroutine to compute precipitation phase partitioning  
        call S3M_Phys_Snow_Apps_PhasePart_Froidurot(iID, iRows, iCols, &
                                            a2dVarDEM, &
                                            a2dVarTa, a2dVarRelHum, a2dVarPrecip, a2dVarSepCoeff, a2dVarSnowFall, &
                                            a2dVarRainfall)
        ! Compute cumulative daily snowfall             
        a2dVarSnowFallDayCum = a2dVarSnowFallDayCum + a2dVarSnowFall
        !-----------------------------------------------------------------------------------------  
        
        !-----------------------------------------------------------------------------------------  
        ! Update SWE (1): inputs
        where( (a2dVarDEM.ge.0.0) .and. (a2dVarSnowFall .gt. 0.0))
            a2dVarSWE_D = a2dVarSWE_D + a2dVarSnowFall
        endwhere
        
        where( (a2dVarDEM.ge.0.0) .and. (a2dVarRainfall .gt. 0.0) .and. (a2dVarSWE_D.ge.10.0))
            a2dVarSWE_W = a2dVarSWE_W + a2dVarRainfall
        endwhere 
        
        where( (a2dVarDEM.ge.0.0) .and. (a2dVarRainfall .gt. 0.0) .and. (a2dVarSWE_D.lt.10.0))
            a2dVarOutflow_ExcessRain = a2dVarRainfall + a2dVarSWE_W
            a2dVarSWE_W = 0.0
        endwhere         
        a2dVarSWE = a2dVarSWE_D + a2dVarSWE_W
     
        ! Check variables
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarSWE = 0.0 ! just in case
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarSWE_D = 0.0 ! just in case  
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarRho_D = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarH_D = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarSWE_W = 0.0 ! just in case              
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarTheta_W = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarRhoS = 0.0 ! just in case 
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarH_S = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2iVarAgeS = 0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarSSA = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarPerm = 0.0 ! just in case   
        where( (a2dVarDem.ge.0.0) .and. (a2dVarIceThick_WE.lt.0.01) ) a2dVarIceThick_WE = 0.0 ! just in case                  
        !-----------------------------------------------------------------------------------------   
            
        !-----------------------------------------------------------------------------------------  
        ! Call subroutine to compute snow density
        call S3M_Phys_Snow_Apps_Rho(iID, iRows, iCols, &
                                       sTime, iTime, iDtDataForcing, &
                                       a2dVarDem, &
                                       a2dVarTa, a2dVarSnowFall, a2dVarSWE_D, &
                                       a2dVarH_D, a2dVarRho_D, a2dVarRhoS0)
   
        ! Check variables
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarSWE = 0.0 ! just in case
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarSWE_D = 0.0 ! just in case  
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarRho_D = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarH_D = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarSWE_W = 0.0 ! just in case              
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarTheta_W = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarRhoS = 0.0 ! just in case 
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarH_S = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2iVarAgeS = 0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarSSA = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarPerm = 0.0 ! just in case   
        where( (a2dVarDem.ge.0.0) .and. (a2dVarIceThick_WE.lt.0.01) ) a2dVarIceThick_WE = 0.0 ! just in case                                
        !-----------------------------------------------------------------------------------------                                        
            
        !-----------------------------------------------------------------------------------------  
        ! Call subroutine to compute average temperature over 1 days
        call S3M_Phys_Snow_Apps_TMean(iID, iRows, iCols, & 
                                        iDaySteps1Days, &
                                        a2iVarMask, &
                                        oS3M_Vars(iID)%a3dTaC_Days1, &
                                        a2dVarTa, a2dVarTaC_MeanDays1)
        !-----------------------------------------------------------------------------------------  
                                        
        !-----------------------------------------------------------------------------------------  
        ! Call subroutine to compute surface-snow age
        call S3M_Phys_Snow_Apps_Age(iID, iRows, iCols, &
                                    sTime, iTime, &
                                    a2dVarDem, &
                                    a2dVarSnowFallDayCum, a2dVarSWE, &
                                    a2iVarAgeS) 
        !-----------------------------------------------------------------------------------------     
                                    
        !-----------------------------------------------------------------------------------------     
        ! We now compute average temperature over the last iDaysAvgTSuppressMelt days...
        where( (a2dVarDem.ge.0.0) )
            a2dVarTaC_MeanDaysSuppressMelt = & 
                    (a2dVarTaC_MeanDaysSuppressMelt*(iDaySteps1Days*iDaysAvgTSuppressMelt - 1) &
                    + a2dVarTa)/(iDaySteps1Days*iDaysAvgTSuppressMelt)
        !We subtract - 1 from iDaySteps1Days*iDaysAvgTSuppressMelt here to exclude the current timestep, which has T = a2dVarTa             
        endwhere   
        !------------------------------------------------------------------------------------- 
        
        !-------------------------------------------------------------------------------------
        ! Call subroutine to compute snow albedo
        call S3M_Phys_Snow_Apps_Albedo(iID, iRows, iCols, &
                                       sTime, iTime, iGlacierValue, &
                                       a2dVarDem, a2iVarGlacierMask, &
                                       a2dVarTaC_MeanDays1, a2iVarAgeS, &
                                       a2dVarAlbedoS, iFlagIceMassBalance, a2dVarSWE_D, a2dVarIceThick_WE)                                       

        !-------------------------------------------------------------------------------------
                                       
        !-------------------------------------------------------------------------------------
        ! Call subroutine to compute open-loop melting
        call S3M_Phys_Snow_Apps_MeltingRefr(iID, iRows, iCols, iDtDataForcing, iDaySteps1Days, &
                                                sTime, iTime, &
                                                iGlacierValue, dVarRhoW, dVarMeltingTRef, dVarIceMeltingCoeff, &
                                                a2dVarDem, a2iVarGlacierMask, a2dVarArctUp, &
                                                a2dVarTa, a2dVarIncRad, &
                                                a2dVarAlbedoS, a2dVarSWE_D, a2dVarMeltingS, a2dVarMeltingSc, &
                                                a2dVarSWE_W, a2dVarRefreezingS, dVarRefreezingSc, & 
                                                iFlagIceMassBalance, a2dVarIceThick_WE, a2dVarMeltingG,  a2dVarMeltSIncRad, &
                                                a2dVarMeltSTemp, a2dVarMeltGIncRad, a2dVarMeltGTemp, &
                                                a2dVarTaC_MeanDaysSuppressMelt, dVarModFactorRadS, a2dVarMeltingSRadc, &
                                                iFlagGlacierDebris)                            
        !-------------------------------------------------------------------------------------
        
        !-------------------------------------------------------------------------------------
        ! Update SWE (2): melt. Here we must consider glaciers too (if any)        
        if (iFlagIceMassBalance.eq.1) then
        ! This first case regards a simulation with mass balance but no movement, where we simply subtract glacier melt from 
        ! a2dVarIceThick_WE as simulation time passes (of course, if no snow is on the ground).     
            !-------------------------------------------------------------------------------------
            ! PIXELS W/O GLACIERS AND WITH SNOW
            ! where Melt > SWE_D: deplete all snowpack and send to Outflow for mass conservation
            where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE_D.gt.0.0) .and. &
                    (a2dVarSWE_D.le.a2dVarMeltingS) .and. (a2dVarIceThick_WE.le.0.0) )
                    
                a2dVarMeltingS = a2dVarSWE_D
                a2dVarOutflow_ExcessMelt = a2dVarSWE_D + a2dVarSWE_W                
                a2dVarSWE_D = 0.0
                a2dVarSWE_W = 0.0
                a2dVarSWE = 0.0
                
            ! elsewhere Melt < SWE_D: reduce SWE_D and increase SWE_W, update SWE, Outflow not involved here
            elsewhere( (a2dVarDem.ge.0.0) .and. (a2dVarSWE_D.gt.0.0) .and. (a2dVarIceThick_WE.le.0.0))

                a2dVarSWE_D = a2dVarSWE_D - a2dVarMeltingS
                a2dVarSWE_W = a2dVarSWE_W + a2dVarMeltingS
                a2dVarSWE = a2dVarSWE_D + a2dVarSWE_W
                
            endwhere
            !-------------------------------------------------------------------------------------

            !-------------------------------------------------------------------------------------
            ! PIXELS W GLACIERS AND WITH SNOW
            ! where Melt > SWE_D: deplete all snowpack, send to Outflow for mass conservation, erode glacier
            where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE_D.gt.0.0) .and. &
                    (a2dVarSWE_D.le.a2dVarMeltingS) .and. (a2dVarIceThick_WE.gt.0.0) )
                    
                a2dVarOutflow_ExcessMelt = a2dVarSWE_D + a2dVarSWE_W 
                a2dVarIceThick_WE = a2dVarIceThick_WE - (a2dVarMeltingS - a2dVarSWE_D)  
                a2dVarMeltingG = a2dVarMeltingG + (a2dVarMeltingS - a2dVarSWE_D) !original a2dVarMeltingG should be 0.0
                a2dVarSWE_D = 0.0
                a2dVarSWE_W = 0.0
                a2dVarSWE = 0.0                    
                a2dVarIceFlag = 0.0 !this a2dVarIceFlag is used to ''turn off'' ice melting and so avoid that a2dVarMeltingG 
                                    !is subtracted twice below in the ''PIXELS W GLACIERS BUT NO SNOW'' case

            ! elsewhere Melt < SWE_D: reduce SWE_D and increase SWE_W, update SWE, Outflow not involved here               
            elsewhere( (a2dVarDem.ge.0.0) .and. (a2dVarSWE_D.gt.0.0) .and. (a2dVarIceThick_WE.gt.0.0))
    
                a2dVarSWE_D = a2dVarSWE_D - a2dVarMeltingS
                a2dVarSWE_W = a2dVarSWE_W + a2dVarMeltingS
                a2dVarSWE = a2dVarSWE_D + a2dVarSWE_W

            endwhere
            !-------------------------------------------------------------------------------------

            !-------------------------------------------------------------------------------------
            ! PIXELS W GLACIERS BUT NO SNOW          
            where ( (a2dVarDem.ge.0.0) .and. (a2dVarSWE_D.eq.0.0) .and. (a2dVarIceThick_WE.gt.0.0) &
                 .and. (a2dVarMeltingG.gt.0.0) .and. (a2dVarIceFlag.gt.0.0))
            
                    a2dVarIceThick_WE = a2dVarIceThick_WE - a2dVarMeltingG
                 
            endwhere
            !-------------------------------------------------------------------------------------
        elseif (iFlagIceMassBalance.eq.2) then
        ! This second case regards a simulation with mass balance AND movement according to the deltaH parametrization, 
        ! so here WE DO NOT subtract glacier melt from a2dVarIceThick_WE as simulation time passes. We store melt into  
        ! a2dVarMeltingGCumWY    
            !-------------------------------------------------------------------------------------
            ! PIXELS W/O GLACIERS AND WITH SNOW
            ! where Melt > SWE_D: deplete all snowpack and send to Outflow for mass conservation
            where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE_D.gt.0.0) .and. &
                    (a2dVarSWE_D.le.a2dVarMeltingS) .and. (a2dVarIceThick_WE.le.0.0) )
                    
                a2dVarMeltingS = a2dVarSWE_D
                a2dVarOutflow_ExcessMelt = a2dVarSWE_D + a2dVarSWE_W                
                a2dVarSWE_D = 0.0
                a2dVarSWE_W = 0.0
                a2dVarSWE = 0.0
                
            ! elsewhere Melt < SWE_D: reduce SWE_D and increase SWE_W, update SWE, Outflow not involved here
            elsewhere( (a2dVarDem.ge.0.0) .and. (a2dVarSWE_D.gt.0.0) .and. (a2dVarIceThick_WE.le.0.0))

                a2dVarSWE_D = a2dVarSWE_D - a2dVarMeltingS
                a2dVarSWE_W = a2dVarSWE_W + a2dVarMeltingS
                a2dVarSWE = a2dVarSWE_D + a2dVarSWE_W
                
            endwhere
            !-------------------------------------------------------------------------------------

            !-------------------------------------------------------------------------------------
            ! PIXELS W GLACIERS AND WITH SNOW
            ! where Melt > SWE_D: deplete all snowpack, send to Outflow for mass conservation, erode glacier
            where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE_D.gt.0.0) .and. &
                    (a2dVarSWE_D.le.a2dVarMeltingS) .and. (a2dVarIceThick_WE.gt.0.0) )
                    
                a2dVarOutflow_ExcessMelt = a2dVarSWE_D + a2dVarSWE_W 
                a2dVarMeltingG = a2dVarMeltingG + (a2dVarMeltingS - a2dVarSWE_D) !original a2dVarMeltingG should be 0.0                
                a2dVarMeltingGCumWY = a2dVarMeltingG + a2dVarMeltingGCumWY
                a2dVarSWE_D = 0.0
                a2dVarSWE_W = 0.0
                a2dVarSWE = 0.0                    
                a2dVarIceFlag = 0.0 !this a2dVarIceFlag is used to ''turn off'' ice melting and so avoid that a2dVarMeltingG 
                                    !is subtracted twice below in the ''PIXELS W GLACIERS BUT NO SNOW'' case

            ! elsewhere Melt < SWE_D: reduce SWE_D and increase SWE_W, update SWE, Outflow not involved here               
            elsewhere( (a2dVarDem.ge.0.0) .and. (a2dVarSWE_D.gt.0.0) .and. (a2dVarIceThick_WE.gt.0.0))
    
                a2dVarSWE_D = a2dVarSWE_D - a2dVarMeltingS
                a2dVarSWE_W = a2dVarSWE_W + a2dVarMeltingS
                a2dVarSWE = a2dVarSWE_D + a2dVarSWE_W

            endwhere
            !-------------------------------------------------------------------------------------

            !-------------------------------------------------------------------------------------
            ! PIXELS W GLACIERS BUT NO SNOW          
            where ( (a2dVarDem.ge.0.0) .and. (a2dVarSWE_D.eq.0.0) .and. (a2dVarIceThick_WE.gt.0.0) &
                 .and. (a2dVarMeltingG.gt.0.0) .and. (a2dVarIceFlag.gt.0.0))
            
                    a2dVarMeltingGCumWY = a2dVarMeltingG + a2dVarMeltingGCumWY
                 
            endwhere
            !-------------------------------------------------------------------------------------            
        else            
            !-------------------------------------------------------------------------------------
            ! PIXELS where Melt > SWE_D: deplete all snowpack and send to Outflow for mass conservation            
            where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE_D.gt.0.0) .and. &
                    (a2dVarSWE_D.le.a2dVarMeltingS) ) 
                
                a2dVarMeltingS = a2dVarSWE_D
                a2dVarOutflow_ExcessMelt = a2dVarSWE_D + a2dVarSWE_W                
                a2dVarSWE_D = 0.0
                a2dVarSWE_W = 0.0
                a2dVarSWE = 0.0   
                
            ! elsewhere Melt < SWE_D: reduce SWE_D and increase SWE_W, update SWE, Outflow not involved here  
            elsewhere( (a2dVarDem.ge.0.0) .and. (a2dVarSWE_D.gt.0.0) )

                a2dVarSWE_D = a2dVarSWE_D - a2dVarMeltingS
                a2dVarSWE_W = a2dVarSWE_W + a2dVarMeltingS
                a2dVarSWE = a2dVarSWE_D + a2dVarSWE_W
                
            endwhere
            !-------------------------------------------------------------------------------------
        
        endif
        !-------------------------------------------------------------------------------------
        
        !-------------------------------------------------------------------------------------
        ! Update SWE 
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE_D.lt.0.0) ) a2dVarSWE_D = 0.0 !just in case....
        a2dVarSWE = a2dVarSWE_D + a2dVarSWE_W !an update just in case...

        ! Compute daily cumulated melting
        where( (a2dVarDem.ge.0.0) .and. (a2dVarMeltingS.gt.0.0) )
            a2dVarMeltingSDayCum = a2dVarMeltingSDayCum + a2dVarMeltingS
        endwhere
    
        !-------------------------------------------------------------------------------------
        
        !-----------------------------------------------------------------------------------------        
        ! Check variables
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarSWE = 0.0 ! just in case
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarSWE_D = 0.0 ! just in case  
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarRho_D = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarH_D = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarSWE_W = 0.0 ! just in case              
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarTheta_W = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarRhoS = 0.0 ! just in case 
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarH_S = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2iVarAgeS = 0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarSSA = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarPerm = 0.0 ! just in case               
        where( (a2dVarDem.ge.0.0) .and. (a2dVarIceThick_WE.lt.0.01) ) a2dVarIceThick_WE = 0.0 ! just in case                  
        !-----------------------------------------------------------------------------------------        
        
        !----------------------------------------------------------------------------------------- 
        ! Update SWE (3): Refreezing (better to decouple this from melt for readability). 
        ! Here we also update dry-snow density
            
        ! PIXELS with a2dVarSWE_W and where Refreezing > SWE_W: set SWE_W to 0 and send it to SWE_D, correct Refreezing          
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE_W.gt.0.0) .and. (a2dVarRefreezingS.gt.0.0) .and. &
                (a2dVarSWE_W.le.a2dVarRefreezingS) ) 

            a2dVarRefreezingS = a2dVarSWE_W
            a2dVarRho_D = (a2dVarSWE_D + a2dVarRefreezingS)/((a2dVarRefreezingS/917) + (a2dVarSWE_D/a2dVarRho_D))                
            a2dVarSWE_D = a2dVarSWE_D + a2dVarRefreezingS
            a2dVarSWE_W = 0.0
            a2dVarSWE = a2dVarSWE_D  

        ! elsewhere with a2dVarSWE_W where Refreezing < SWE_W: reduce SWE_W and increase SWE_D, update SWE  
        elsewhere( (a2dVarDem.ge.0.0) .and. (a2dVarSWE_W.gt.0.0) .and. (a2dVarRefreezingS.gt.0.0))

            a2dVarRho_D = (a2dVarSWE_D + a2dVarRefreezingS)/((a2dVarRefreezingS/917) + (a2dVarSWE_D/a2dVarRho_D))
            a2dVarSWE_D = a2dVarSWE_D + a2dVarRefreezingS
            a2dVarSWE_W = a2dVarSWE_W - a2dVarRefreezingS
            a2dVarSWE = a2dVarSWE_D + a2dVarSWE_W

        endwhere
        !----------------------------------------------------------------------------------------- 
        
        !----------------------------------------------------------------------------------------- 
        ! Update SWE 
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE_W.lt.0.0) ) a2dVarSWE_W = 0.0
        a2dVarSWE = a2dVarSWE_D + a2dVarSWE_W
        !----------------------------------------------------------------------------------------- 
                
        !-----------------------------------------------------------------------------------------        
        ! Check variables
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) ) a2dVarSWE = 0.0 ! just in case
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarSWE_D = 0.0 ! just in case  
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarRho_D = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarH_D = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarSWE_W = 0.0 ! just in case              
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarTheta_W = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarRhoS = 0.0 ! just in case 
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarH_S = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2iVarAgeS = 0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarSSA = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarPerm = 0.0 ! just in case                        
        where( (a2dVarDem.ge.0.0) .and. (a2dVarIceThick_WE.lt.0.01) ) a2dVarIceThick_WE = 0.0 ! just in case                  
        !-----------------------------------------------------------------------------------------          
            
        !-----------------------------------------------------------------------------------------  
        ! Call subroutine to compute outflow      
        call S3M_Phys_Snow_Apps_Outflow(iID, iRows, iCols, iDtDataForcing, &
        dVarRhoW, &
        a2dVarDem, &
        a2dVarSWE_D, a2dVarRho_D, a2dVarSWE_W, & 
        a2dVarRainfall, a2dVarMeltingS, a2dVarRefreezingS, &
        a2dVarOutflow_K, a2dVarH_D, a2dVarH_S, a2dVarSSA, a2dVarPerm, a2dVarTheta_W) 
        !-----------------------------------------------------------------------------------------  
        
        !-----------------------------------------------------------------------------------------  
        ! Update SWE (4): outflow
        ! Here we must bear in mind that a2dVarOutflow is the result of two computations: 
        ! (1) complete snow depletion due to melt or rain on shallow snowpack; in that case, SWE_W has already been updated
        ! (2) Outflow due to drainage of liquid water from the snowpack. In this case, SWE_W has not been updated yet. 
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE_W.gt.0.0) )

            a2dVarSWE_W = a2dVarSWE_W - a2dVarOutflow_K !the only flux we have not taken into account yet
            a2dVarSWE = a2dVarSWE_D + a2dVarSWE_W
            !Note that the S3M_Phys_Snow_Apps_Outflow routine already includes a number of sanity checks to make sure that 
            !a2dVarOutflow_K should be .le. a2dVarSWE_W

        endwhere

        !we now compute total outflow from all sources...
        a2dVarOutflow = a2dVarOutflow_K + a2dVarOutflow_ExcessRain + a2dVarOutflow_ExcessMelt

        !sanity check
        where( (a2dVarSWE_W.lt.0.0) ) 

            a2dVarOutflow = a2dVarOutflow + a2dVarSWE_W
            a2dVarSWE_W = 0.0 
            a2dVarSWE = a2dVarSWE_D
            !this should never happen, but if this was the case then remove the negative a2dVarSWE_W from a2dVarOutflow
            !and correct a2dVarSWE_W.

        endwhere 
        !--------------------------------------------------------------------        
        
        !-----------------------------------------------------------------------------------------        
        ! Check variables
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) ) a2dVarSWE = 0.0 ! just in case
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarSWE_D = 0.0 ! just in case  
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarRho_D = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarH_D = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarSWE_W = 0.0 ! just in case              
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarTheta_W = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarRhoS = 0.0 ! just in case 
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarH_S = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2iVarAgeS = 0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarSSA = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarPerm = 0.0 ! just in case                        
        where( (a2dVarDem.ge.0.0) .and. (a2dVarIceThick_WE.lt.0.01) ) a2dVarIceThick_WE = 0.0 ! just in case                 
        !-----------------------------------------------------------------------------------------           
        
        !-----------------------------------------------------------------------------------------        
        ! Compute Reff
        ! Reff is the sum of any snowpack runoff and glacier melt resulting at this point. 
        where( (a2dVarDem.ge.0.0) )    
            a2dVarReff = a2dVarOutflow + a2dVarMeltingG
            ! note that here we add a2dVarMeltingG to outflow regardless of iFlagIceMassBalance. If iFlagIceMassBalance = 0 or 
            ! 1, this choice has no implication; if iFlagIceMassBalance = 2, this choice implies that we disconnect the update of
            ! glacier thickness (which happens once per yr) from the sub-annual patterns of a2dVarMeltingG. 
        endwhere
        !-----------------------------------------------------------------------------------------  
        
        !-----------------------------------------------------------------------------------------        
        ! In preparation of the assimilation part, we update control volumes, bulk-snow density, 
        ! and some ancillary state variables       
        where( (a2dVarDem.ge.0.0) .and. (a2dVarRho_D .gt. 0.0) )    
            
            a2dVarH_D = ((a2dVarSWE_D/1000)*dVarRhoW)/a2dVarRho_D
        
            where( (a2dVarSWE_W/1000 .ge. ((1 - a2dVarRho_D/917)*a2dVarH_D)) ) 
            ! We divide by 1000 to convert mm to m and so compare a2dVarSWE_W with a2dVarH_D
            ! (1 -  - a2dVarRho_D/917) is porosity
                a2dVarH_S = a2dVarH_D + ((a2dVarSWE_W/1000) - (1 - a2dVarRho_D/917)*a2dVarH_D)
           
            elsewhere 
                
                a2dVarH_S = a2dVarH_D
                
            endwhere 
            
            where( (a2dVarH_S .gt. 0.0) )
                
                a2dVarTheta_W = a2dVarSWE_W/1000/a2dVarH_S 
                a2dVarRhoS = (a2dVarRho_D*a2dVarH_S + dVarRhoW*(a2dVarSWE_W/1000))/(a2dVarH_S)            
            
            elsewhere 
            
                a2dVarTheta_W = 0.0
                a2dVarRhoS = 0.0
            
            endwhere
            
        endwhere        
        !-----------------------------------------------------------------------------------------    
        
        !-----------------------------------------------------------------------------------------  
        ! Call subroutine to compute SWE assimilation from snow height interpolated observations
        call mprintf(.true., iINFO_Verbose, ' Phys :: Snow :: Assimilation of SH and MODIS... ' )
        if (iFlagSnowAssim.eq.1) then 
        
            !-----------------------------------------------------------------------------------------  
            ! Check forcing(s) to use assimilation method
            if ( any(a2dVarSnowHeight.ne.-9999.0) .and. any(a2dVarSnowKernel.ne.-9999.0) .and. &
                any(a2dVarSnowCA.ne.-9999.0) .and. any(a2dVarSnowQA.ne.-9999.0) ) then
                
               !-----------------------------------------------------------------------------------------  
                ! We store the modeled map of a2dVarSWE into a matrix, that we will then use to update SWE_D and SWE_W
                a2dVarSWE_preass = a2dVarSWE

                ! Call subroutine to compute SWE assimilation from snow height interpolated observations                
                call S3M_Phys_Snow_Apps_AssimSH(iID, iRows, iCols, &
                                                    dVarRhoW, &
                                                    a2iVarMask, &
                                                    a2dVarSnowHeight, a2dVarSnowKernel, &
                                                    a2dVarSnowCA, a2dVarSnowQA, &
                                                    a2dVarSWE, a2iVarAgeS, a2dVarAlbedoS, a2dVarRhoS, iFlagAssOnlyPos)
                                                    
                ! Now we update SWE_D and SWE_W as well
                where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE_preass.gt.10) .and. (a2dVarTheta_W .lt. 0.1) )                        
                                                    
                    a2dVarUass = a2dVarSWE/a2dVarSWE_preass
                    a2dVarSWE_D = a2dVarUass*a2dVarSWE_D
                    a2dVarSWE_W = a2dVarUass*a2dVarSWE_W
                    
                elsewhere( (a2dVarDem.ge.0.0) .and. (a2dVarSWE_preass.gt.10) ) 
                    
                    a2dVarUass = a2dVarSWE - a2dVarSWE_preass
                    a2dVarSWE_D = a2dVarSWE + a2dVarUass
                    a2dVarSWE_W = a2dVarSWE_W
                    
                elsewhere( (a2dVarDem.ge.0.0) ) 
                    
                    a2dVarSWE_D = a2dVarSWE 
                    a2dVarSWE_W = 0.0                    
                    
                endwhere
                
                ! Info end assimilation
                call mprintf(.true., iINFO_Verbose, ' Phys :: Snow :: Assimilation of SH and MODIS... OK' )
                !-----------------------------------------------------------------------------------------  
 
            else
              
                !-----------------------------------------------------------------------------------------  
                ! Info assimilation no data available
                call mprintf(.true., iINFO_Verbose, ' Phys :: Snow :: Assimilation of SH and MODIS... DATA NO AVAILABLE ' )
                !-----------------------------------------------------------------------------------------  
                    
            endif
            !-----------------------------------------------------------------------------------------  

        else
                
            !-----------------------------------------------------------------------------------------  
            ! Info assimilation no data available
            call mprintf(.true., iINFO_Verbose, ' Phys :: Snow :: Assimilation of SH and MODIS... NOT ACTIVATED ' )
            !-----------------------------------------------------------------------------------------  
                
        endif
        !-----------------------------------------------------------------------------------------  
        
        !-----------------------------------------------------------------------------------------        
        ! Check variables
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) ) a2dVarSWE = 0.0 ! just in case
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarSWE_D = 0.0 ! just in case  
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarRho_D = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarH_D = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarSWE_W = 0.0 ! just in case              
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarTheta_W = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarRhoS = 0.0 ! just in case 
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarH_S = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2iVarAgeS = 0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarSSA = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarPerm = 0.0 ! just in case                        
        where( (a2dVarDem.ge.0.0) .and. (a2dVarIceThick_WE.lt.0.01) ) a2dVarIceThick_WE = 0.0 ! just in case                
        !-----------------------------------------------------------------------------------------      
            
        !-----------------------------------------------------------------------------------------        
        ! In preparation of the assimilation part, we update control volumes, bulk-snow density, and some ancillary state variables       
        where( (a2dVarDem.ge.0.0) .and. (a2dVarRho_D .gt. 0.0) )    
            
            a2dVarH_D = ((a2dVarSWE_D/1000)*dVarRhoW)/a2dVarRho_D
        
            where( (a2dVarSWE_W/1000 .ge. ((1 - a2dVarRho_D/917)*a2dVarH_D)) ) 
            !we divide by 1000 to convert mm to m and so compare a2dVarSWE_W with a2dVarH_D
            !(1 -  - a2dVarRho_D/917) is porosity
                a2dVarH_S = a2dVarH_D + ((a2dVarSWE_W/1000) - (1 - a2dVarRho_D/917)*a2dVarH_D)
           
            elsewhere 
                
                a2dVarH_S = a2dVarH_D
                
            endwhere 
            
            where( (a2dVarH_S .gt. 0.0) )
                
                a2dVarTheta_W = a2dVarSWE_W/1000/a2dVarH_S 
                a2dVarRhoS = (a2dVarRho_D*a2dVarH_S + dVarRhoW*(a2dVarSWE_W/1000))/(a2dVarH_S)            
            
            elsewhere 
            
                a2dVarTheta_W = 0.0
                a2dVarRhoS = 0.0
            
            endwhere
            
            
        endwhere           
        !-----------------------------------------------------------------------------------------    

        !-------------------------------------------------------------------------------------
        ! Call subroutine to compute SWE assimilation map, ex. from ARPA VdA
        call mprintf(.true., iINFO_Verbose, ' Phys :: Snow :: Assimilation of SWE map... ' )
        if (iFlagSnowAssim_SWE.eq.1) then 

            ! Check forcing(s) to use assimilation method
            if ( any(a2dVarSWEass.ne.-9999.0) ) then
                
                ! We store the modeled map of a2dVarSWE into a matrix, that we will then use to update SWE_D and SWE_W
                a2dVarSWE_preass = a2dVarSWE
                
                ! Call subroutine to compute SWE assimilation from external map
                call S3M_Phys_Snow_Apps_AssimSWE(iID, iRows, iCols, &
                                                    a2iVarMask, &
                                                    a2dVarSWE, a2dVarSWEass, iFlagAssOnlyPos)

                ! Now we update SWE_D and SWE_W as well
                where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE_preass.gt.10) .and. (a2dVarTheta_W .lt. 0.1) )                        
                                                    
                    a2dVarUass = a2dVarSWE/a2dVarSWE_preass
                    a2dVarSWE_D = a2dVarUass*a2dVarSWE_D
                    a2dVarSWE_W = a2dVarUass*a2dVarSWE_W
                    
                elsewhere( (a2dVarDem.ge.0.0) .and. (a2dVarSWE_preass.gt.10) ) 
                    
                    a2dVarUass = a2dVarSWE - a2dVarSWE_preass
                    a2dVarSWE_D = a2dVarSWE + a2dVarUass
                    a2dVarSWE_W = a2dVarSWE_W
                    
                elsewhere( (a2dVarDem.ge.0.0) ) 
                    
                    a2dVarSWE_D = a2dVarSWE 
                    a2dVarSWE_W = 0.0                    
                    
                endwhere                            

                ! Info start assimilation
                call mprintf(.true., iINFO_Verbose, ' Phys :: Snow :: Assimilation of SWE map ... OK' )
                !-------------------------------------------------------------------------------------
 
            else
              
                !-------------------------------------------------------------------------------------
                ! Info assimilation no data available
                call mprintf(.true., iINFO_Verbose, ' Phys :: Snow :: Assimilation of SWE map ... DATA NO AVAILABLE ' )
                !-------------------------------------------------------------------------------------
                    
            endif
            !-------------------------------------------------------------------------------------

        else
                
            !-------------------------------------------------------------------------------------
            ! Info assimilation no data available
            call mprintf(.true., iINFO_Verbose, ' Phys :: Snow :: Assimilation of SWE map... NOT ACTIVATED ' )
            !-------------------------------------------------------------------------------------
                
        endif
        !-------------------------------------------------------------------------------------
        
        !-----------------------------------------------------------------------------------------        
        ! Check variables
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) ) a2dVarSWE = 0.0 ! just in case
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarSWE_D = 0.0 ! just in case  
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarRho_D = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarH_D = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarSWE_W = 0.0 ! just in case              
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarTheta_W = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarRhoS = 0.0 ! just in case 
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarH_S = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2iVarAgeS = 0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarSSA = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarPerm = 0.0 ! just in case                        
        where( (a2dVarDem.ge.0.0) .and. (a2dVarIceThick_WE.lt.0.01) ) a2dVarIceThick_WE = 0.0 ! just in case                
        !-----------------------------------------------------------------------------------------      
            
        !-----------------------------------------------------------------------------------------        
        ! We update again control volumes, bulk-snow density, and some ancillary state variables       
        where( (a2dVarDem.ge.0.0) .and. (a2dVarRho_D .gt. 0.0) )    
            
            a2dVarH_D = ((a2dVarSWE_D/1000)*dVarRhoW)/a2dVarRho_D
        
            where( (a2dVarSWE_W/1000 .ge. ((1 - a2dVarRho_D/917)*a2dVarH_D)) ) 
            !we divide by 1000 to convert mm to m and so compare a2dVarSWE_W with a2dVarH_D
            !(1 -  - a2dVarRho_D/917) is porosity
                a2dVarH_S = a2dVarH_D + ((a2dVarSWE_W/1000) - (1 - a2dVarRho_D/917)*a2dVarH_D)
           
            elsewhere 
                
                a2dVarH_S = a2dVarH_D
                
            endwhere 
            
            where( (a2dVarH_S .gt. 0.0) )
                
                a2dVarTheta_W = a2dVarSWE_W/1000/a2dVarH_S 
                a2dVarRhoS = (a2dVarRho_D*a2dVarH_S + dVarRhoW*(a2dVarSWE_W/1000))/(a2dVarH_S)            
            
            elsewhere 
            
                a2dVarTheta_W = 0.0
                a2dVarRhoS = 0.0
            
            endwhere
            
        endwhere        
        !-----------------------------------------------------------------------------------------      
        
        !-----------------------------------------------------------------------------------------        
        ! We apply the deltaH parametrization for glacier-thickness update, as well as snow-to-ice conversion        
        !-----------------------------------------------------------------------------------------  
        if ( (iFlagIceMassBalance.eq.2) .and. (sTime(9:10).eq.'01') .and. (sTime(6:7).eq.sWYstart) & 
            .and. (sTime(12:13).eq.'00')) then
                
            where((a2dVarDEM .ge. 0.0) .and. (a2dVarSWE.gt.0.0))
                
                a2dVarIceThick_WE = a2dVarIceThick_WE + a2dVarSWE
                a2dVarSWE = 0.0
                a2dVarSWE_D = 0.0
                a2dVarSWE_W = 0.0
                
            endwhere                
                
            call   S3M_Phys_Snow_Apps_GlacierDeltaH(iID, iRows, iCols, iRows_Pivot, iCols_Pivot, a2dVarDEM, &
                                                    a2dVarIceThick_WE, a2dVarMeltingGCumWY, a2iVarMask)
                                                    
            a2dVarMeltingGCumWY = 0.0                                        
                                                    
        endif
        
        if ( (iFlagIceMassBalance.eq.1.0) .or. (iFlagIceMassBalance.eq.2.0) ) then
            
            ! Conversion from mm of water equivalent to m of ice
            where(a2dVarDEM.gt.0.0)
                a2dVarChangeThickness = a2dVarIceThick_WE/917 - a2dVarIceThick
                a2dVarIceThick = a2dVarIceThick_WE/917
            endwhere
        endif
        
        !-----------------------------------------------------------------------------------------        
        ! Check variables
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) ) a2dVarSWE = 0.0 ! just in case
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarSWE_D = 0.0 ! just in case  
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarRho_D = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarH_D = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarSWE_W = 0.0 ! just in case              
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarTheta_W = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarRhoS = 0.0 ! just in case 
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarH_S = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2iVarAgeS = 0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarSSA = 0.0 ! just in case            
        where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.01) )  a2dVarPerm = 0.0 ! just in case                        
        where( (a2dVarDem.ge.0.0) .and. (a2dVarIceThick_WE.lt.0.01) ) a2dVarIceThick_WE = 0.0 ! just in case                
        !-----------------------------------------------------------------------------------------         
        
        !-----------------------------------------------------------------------------------------  
        ! Nullify variables where DEM < 0
         where(a2dVarDem.lt.0.0)
            
            !variables: snowfall and rainfall
            a2dVarRainFall = -9999.0; a2dVarSnowFall = -9999.0; a2dVarSnowFallDayCum = -9999.0; 
            
            !variables: dry and wet snow
            a2dVarSWE_D = -9999; a2dVarRho_D = -9999; a2dVarH_D = -9999.0;
            a2dVarSWE_W = -9999; a2dVarTheta_W = -9999;
            
            !variables: bulk snow 
            a2dVarSWE = -9999.0; a2dVarRhos = -9999.0; a2dVarH_S = -9999.0; a2dVarRhoS0 = -9999.0;
            
            !variables: snow melt and refreezing
            a2dVarMeltingS = -9999.0; a2dVarMeltingSDayCum = -9999.0; 
            a2dVarMeltingG = -9999; 
            a2dVarRefreezingS = -9999;
            a2dVarAlbedoS = -9999; a2iVarAgeS = -9999; 
            a2dVarTaC_MeanDays1 = -9999.0; a2dVarTaC_MeanDaysSuppressMelt = -9999.0;
            a2dVarMeltSIncRad = -9999.0; a2dVarMeltSTemp = -9999.0; a2dVarMeltGIncRad = -9999.0; a2dVarMeltGTemp = -9999.0;
            a2dVarMeltingSc = -9999.0; a2dVarMeltingSRadc = -9999.0;
            
            !variables: snow hydraulics 
            a2dVarOutflow = -9999.0;  a2dVarReff = -9999; 
            a2dVarSSA = -9999; a2dVarPerm = -9999;
            a2dVarOutflow_ExcessRain = -9999.0; a2dVarOutflow_ExcessMelt = -9999.0; a2dVarOutflow_K  = -9999.0; 

            !variables: glaciers
            a2dVarIceThick = -9999.0; a2dVarIceThick_WE = -9999.0; a2dVarMeltingGCumWY = -9999.0;
            a2dVarChangeThickness = -9999.0
        
        endwhere
        !-----------------------------------------------------------------------------------------  
         
        !-----------------------------------------------------------------------------------------  
        ! Snow-mask computation           
        where( (a2dVarDem.gt.0.0) .and. (a2dVarSWE.gt.0.1) )
            a2dVarSnowMask = 1 
        elsewhere
            a2dVarSnowMask = 0
        endwhere
        !-----------------------------------------------------------------------------------------  
        
        !-----------------------------------------------------------------------------------------
        ! Debug
        if (iDEBUG.gt.0) then

            call mprintf(.true., iINFO_Extra, ' ')

            ! METEO
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarPrecip, a2iVarMask, 'PRECIP END ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarTa, a2iVarMask, 'TA END ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarRelHum, a2iVarMask, 'RELHUM END ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarIncRad, a2iVarMask, 'INCRAD END ') )
            
            ! ASSIMILATION
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarSnowHeight, a2iVarMask, 'SNOWH END ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarSnowKernel, a2iVarMask, 'SNOWK END ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarSnowCA, a2iVarMask, 'SNOWCA END ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarSnowQA, a2iVarMask, 'SNOWQA END ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarSWEass, a2iVarMask, 'SWEASS END ') )
                
            ! SNOWFALL AND RAINFALL
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarSepCoeff, a2iVarMask, 'SEPCOEFF END ') )            
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarSnowFall, a2iVarMask, 'SNOWFALL END ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarSnowFallDayCum, a2iVarMask, 'SNOWFALL DAYCUM END ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarRainFall, a2iVarMask, 'RAINFALL END ') )

            ! DRY WET SNOW
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarSWE_D, a2iVarMask, 'SWE_D END ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarRho_D, a2iVarMask, 'RHO_D END ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarSWE_W, a2iVarMask, 'SWE_W END ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarTheta_W, a2iVarMask, 'THETA_W END ') )              
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarH_D, a2iVarMask, 'HD END ') )             

            ! BULK SNOW
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarSWE, a2iVarMask, 'SWE END ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarRhoS, a2iVarMask, 'RHOS END ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarH_S, a2iVarMask, 'HS END ') )             
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarSnowMask, a2iVarMask, 'SNOWMASK END ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarRhoS0, a2iVarMask, 'RHOS0 END ') )

            ! SNOWMELT AND REFREEZING
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarMeltingS, a2iVarMask, 'MELTINGS END ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarMeltingSDayCum, a2iVarMask, 'MELTINGS DAYCUM END ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarMeltingSc, a2iVarMask, 'MELTINGS COEFF END ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarAlbedoS, a2iVarMask, 'ALBEDOS END ') )
            call mprintf(.true., iINFO_Extra, checkvar(real(a2iVarAgeS), a2iVarMask, 'AGES END ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarRefreezingS, a2iVarMask, 'REFREEZING END ') ) 
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarMeltingG, a2iVarMask, 'MELTINGG END ') ) 
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarMeltSIncRad, a2iVarMask, 'MELTSINCRAD EMD ') ) 
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarMeltSTemp, a2iVarMask, 'MELTSTEMP END ') ) 
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarMeltGIncRad, a2iVarMask, 'MELTGINCRAD END ') ) 
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarMeltGTemp, a2iVarMask, 'MELTGTEMP END ') )      
                !OTHER MELT RELATED VARIABLES ARE NOT 2D AND SO WE CANNOT DEBUG THEM HERE
            
            ! SNOW HYDRAULICS
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarOutflow, a2iVarMask, 'OUTFLOW END ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarSSA, a2iVarMask, 'SSA END ') ) 
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarPerm, a2iVarMask, 'PERMEABILITY END ') )   
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarReff, a2iVarMask, 'REFF END ') )     
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarOutflow_ExcessRain, a2iVarMask, 'OUTFLOW EXCESS RAIN END ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarOutflow_ExcessMelt, a2iVarMask, 'OUTFLOW EXCESS MELT END ') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarOutflow_K, a2iVarMask, 'OUTFLOW PERMEABILITY END ') )            
            
            ! GLACIERS
            if ( (iFlagIceMassBalance.eq.1.0) .or. (iFlagIceMassBalance.eq.2.0)) then
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarIceThick_WE, a2iVarMask, 'ICE W.E. END ') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarIceThick, a2iVarMask, 'ICE END ') ) 
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarMeltingGCumWY, a2iVarMask, 'MELTING G CUM WY END ') ) 
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarChangeThickness, a2iVarMask, 'ICE CHANGE THICKNESS END ') ) 
            endif    
                
            call mprintf(.true., iINFO_Extra, ' ========= SNOW END =========== ')  
        endif
        !-----------------------------------------------------------------------------------------

        !-----------------------------------------------------------------------------------------
        ! Update state variables
        ! SNOWFALL AND RAINFALL
        oS3M_Vars(iID)%a2dRainFall = a2dVarRainFall
        oS3M_Vars(iID)%a2dSnowFall = a2dVarSnowFall        
        oS3M_Vars(iID)%a2dSnowFallDayCum = a2dVarSnowFallDayCum        
        
        ! DRY WET SNOW
        oS3M_Vars(iID)%a2dSWE_D = a2dVarSWE_D
        oS3M_Vars(iID)%a2dRho_D = a2dVarRho_D        
        oS3M_Vars(iID)%a2dH_D = a2dVarH_D        
        oS3M_Vars(iID)%a2dSWE_W = a2dVarSWE_W
        oS3M_Vars(iID)%a2dTheta_W = a2dVarTheta_W        
        
        ! BULK SNOW
        oS3M_Vars(iID)%a2dSWE = a2dVarSWE
        oS3M_Vars(iID)%a2dRhoS = a2dVarRhoS
        oS3M_Vars(iID)%a2dRhoS0 = a2dVarRhoS0
        oS3M_Vars(iID)%a2dH_S = a2dVarH_S        
        oS3M_Vars(iID)%a2dMaskS = a2dVarSnowMask
        
        ! SNOWMELT AND REFREEZING
        oS3M_Vars(iID)%a2dMelting = a2dVarMeltingS
        oS3M_Vars(iID)%a2dMeltingSc = a2dVarMeltingSc
        oS3M_Vars(iID)%a2dMeltingDayCum = a2dVarMeltingSDayCum
        oS3M_Vars(iID)%a2dRefreezingS = a2dVarRefreezingS
        oS3M_Vars(iID)%a2dMeltingG = a2dVarMeltingG
        oS3M_Vars(iID)%a2dAlbedo_Snow = a2dVarAlbedoS
        oS3M_Vars(iID)%a2iAge = a2iVarAgeS  
        oS3M_Vars(iID)%a2dMeltSIncRad = a2dVarMeltSIncRad 
        oS3M_Vars(iID)%a2dMeltSTemp = a2dVarMeltSTemp 
        oS3M_Vars(iID)%a2dMeltGIncRad = a2dVarMeltGIncRad 
        oS3M_Vars(iID)%a2dMeltGTemp = a2dVarMeltGTemp         
        oS3M_Vars(iID)%a2dTaC_MeanDaysSuppressMelt = a2dVarTaC_MeanDaysSuppressMelt
        oS3M_Vars(iID)%a2dMeltingSRadc = a2dVarMeltingSRadc
        
        ! SNOW HYDRAULICS
        oS3M_Vars(iID)%a2dOutflow = a2dVarOutflow
        oS3M_Vars(iID)%a2dReff = a2dVarReff        
        oS3M_Vars(iID)%a2dSSA = a2dVarSSA
        oS3M_Vars(iID)%a2dPerm = a2dVarPerm
        
        ! GLACIER
        oS3M_Vars(iID)%a2dIceThick = a2dVarIceThick
        oS3M_Vars(iID)%a2dMeltingGCumWY = a2dVarMeltingGCumWY 
        oS3M_Vars(iID)%a2dChangeThickness = a2dVarChangeThickness

        ! MISCELLANEOUS
        oS3M_Vars(iID)%a2dTaC_MeanDays1 = a2dVarTaC_MeanDays1

        ! Info end
        call mprintf(.true., iINFO_Verbose, ' Phys :: Snow model ... OK' )
        !-----------------------------------------------------------------------------------------
       
    !-----------------------------------------------------------------------------------------

    end subroutine S3M_Phys_Snow_Cpl
    !------------------------------------------------------------------------------------------

end module S3M_Module_Phys_Snow
!------------------------------------------------------------------------------------
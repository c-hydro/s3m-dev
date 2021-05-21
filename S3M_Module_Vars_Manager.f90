!------------------------------------------------------------------------------------
! File:      S3M_Module_Vars_Manager.f90
! Author(s): Fabio Delogu, Francesco Silvestro, Simone Gabellani, Francesco Avanzi
!
! Created on February 11 2015, 9:57 AM
! Last update on May 20, 2021 11:00 AM
!
! Module to allocate and initialize global variable(s)
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
! Module Header
module S3M_Module_Vars_Manager
    
    !------------------------------------------------------------------------------------
    ! External module(s) and implicit none
    use S3M_Module_Namelist,    only: oS3M_Namelist
    use S3M_Module_Vars_Loader,  only: oS3M_Vars
    
    use S3M_Module_Tools_Debug
    
    implicit none
    !------------------------------------------------------------------------------------
    
contains
    
    !------------------------------------------------------------------------------------
    ! Subroutine for allocating S3M static variable(s)
    subroutine S3M_Vars_Allocate(iID, &
                                 iRows, iCols, & 
                                 iRowsPivot, iColsPivot, &
                                 iDaySteps, &
                                 iNData, iNGlaciers)
        
        !------------------------------------------------------------------------------------ 
        ! Variable(s)
        integer(kind = 4)    :: iID
        integer(kind = 4)    :: iRows, iCols, iRowsPivot, iColsPivot
        integer(kind = 4)    :: iDaySteps 
        integer(kind = 4)    :: iNData
        integer(kind = 4)    :: iNGlaciers
        !------------------------------------------------------------------------------------ 
        
        !------------------------------------------------------------------------------------ 
        ! Start to allocate static variable(s)
        call mprintf(.true., iINFO_Main, ' Allocating static variable(s) ... ')

        ! Static Land variable(s)
        allocate( oS3M_Vars(iID)%a2dLon             (iRows, iCols) )
        allocate( oS3M_Vars(iID)%a2dLat             (iRows, iCols) )
        allocate( oS3M_Vars(iID)%a2dDEM             (iRows, iCols) )
        allocate( oS3M_Vars(iID)%a2dAreaCell        (iRows, iCols) )        
        allocate( oS3M_Vars(iID)%a2iMask            (iRows, iCols) )
        allocate( oS3M_Vars(iID)%a2iGlacierMask     (iRows, iCols) )
        allocate( oS3M_Vars(iID)%a2iGlaciers_ID     (iRows, iCols) )
        allocate( oS3M_Vars(iID)%a1iGlaciers_ID     (iNGlaciers)   )
        allocate( oS3M_Vars(iID)%a2dGlacierDebris   (iRows, iCols) )
        allocate( oS3M_Vars(iID)%a2dGlacierPivotTable  ( iRowsPivot, iColsPivot ) )
        
        allocate( oS3M_Vars(iID)%a2dArctUp          (iRows, iCols) )
      
        allocate( oS3M_Vars(iID)%a2iXIndex          (iRows, iCols) )
        allocate( oS3M_Vars(iID)%a2iYIndex          (iRows, iCols) )
 
        ! Finish to allocate static variable(s)
        call mprintf(.true., iINFO_Main, ' Allocating static variable(s) ... OK')
        !------------------------------------------------------------------------------------ 

        !------------------------------------------------------------------------------------ 
        ! Start to allocate dynamic variable(s)
        call mprintf(.true., iINFO_Main, ' Allocating dynamic variable(s) ... ')

        ! Dynamic forcing variable(s)
        allocate( oS3M_Vars(iID)%a2dPrecip          (iRows, iCols) )
        allocate( oS3M_Vars(iID)%a2dTa              (iRows, iCols) )
        allocate( oS3M_Vars(iID)%a2dK               (iRows, iCols) )
        allocate( oS3M_Vars(iID)%a2dRHum            (iRows, iCols) )
        allocate( oS3M_Vars(iID)%a2dSHeight         (iRows, iCols) )
        allocate( oS3M_Vars(iID)%a2dSKernel         (iRows, iCols) )
        
        ! Dynamic updating variable(s)
        allocate( oS3M_Vars(iID)%a2dSCA             (iRows, iCols) )
        allocate( oS3M_Vars(iID)%a2dSQA             (iRows, iCols) )

        ! Dynamic (monthly) Vegetation variable(s)
        allocate( oS3M_Vars(iID)%a2dAlbedo          (iRows, iCols) )
        
        ! Dynamic snow variable(s)
        allocate( oS3M_Vars(iID)%a2iAge             (iRows, iCols) )
        allocate( oS3M_Vars(iID)%a2dSWE             (iRows, iCols) )
        allocate( oS3M_Vars(iID)%a2dRhoS            (iRows, iCols) )
        allocate( oS3M_Vars(iID)%a2dRhoS0           (iRows, iCols) )
        allocate( oS3M_Vars(iID)%a2dAlbedo_Snow     (iRows, iCols) )
        allocate( oS3M_Vars(iID)%a3dTaC_Days1       (iRows, iCols, iDaySteps) )
        allocate( oS3M_Vars(iID)%a2dMelting         (iRows, iCols) )  
        allocate( oS3M_Vars(iID)%a2dMeltingSc       (iRows, iCols) ) 
        allocate( oS3M_Vars(iID)%a2dMeltingDayCum   (iRows, iCols) )
        allocate( oS3M_Vars(iID)%a2dSnowFall        (iRows, iCols) )  
        allocate( oS3M_Vars(iID)%a2dSnowFallDayCum  (iRows, iCols) )
        allocate( oS3M_Vars(iID)%a2dMaskS           (iRows, iCols) )
        allocate( oS3M_Vars(iID)%a2dIceThick        (iRows, iCols) ) 
        allocate( oS3M_Vars(iID)%a2dSWEass          (iRows, iCols) ) 
        allocate( oS3M_Vars(iID)%a2dSWEassHist      (iRows, iCols) ) 
        
        allocate( oS3M_Vars(iID)%a2dOutflow         (iRows, iCols) ) 
        allocate( oS3M_Vars(iID)%a2dRainFall        (iRows, iCols) ) 
        allocate( oS3M_Vars(iID)%a2dSWE_D           (iRows, iCols) ) 
        allocate( oS3M_Vars(iID)%a2dRho_D           (iRows, iCols) ) 
        allocate( oS3M_Vars(iID)%a2dSWE_W           (iRows, iCols) ) 
        allocate( oS3M_Vars(iID)%a2dTheta_W         (iRows, iCols) ) 
        allocate( oS3M_Vars(iID)%a2dRefreezingS     (iRows, iCols) ) 
        allocate( oS3M_Vars(iID)%a2dMeltingG        (iRows, iCols) )         
        allocate( oS3M_Vars(iID)%a2dH_D             (iRows, iCols) )       
        allocate( oS3M_Vars(iID)%a2dH_S             (iRows, iCols) )       
        allocate( oS3M_Vars(iID)%a2dSSA             (iRows, iCols) )       
        allocate( oS3M_Vars(iID)%a2dPerm            (iRows, iCols) )       
        allocate( oS3M_Vars(iID)%a2dReff            (iRows, iCols) )       
        allocate( oS3M_Vars(iID)%a2dMeltSIncRad     (iRows, iCols) )   
        allocate( oS3M_Vars(iID)%a2dMeltSTemp       (iRows, iCols) )   
        allocate( oS3M_Vars(iID)%a2dMeltGIncRad     (iRows, iCols) )   
        allocate( oS3M_Vars(iID)%a2dMeltGTemp       (iRows, iCols) )   
        allocate( oS3M_Vars(iID)%a2dMeltingSRadc    (iRows, iCols) )   
        allocate( oS3M_Vars(iID)%a2dMeltingGCumWY   (iRows, iCols) )   
        allocate( oS3M_Vars(iID)%a2dChangeThickness (iRows, iCols) ) 
        
        allocate( oS3M_Vars(iID)%a2dTaC_MeanDays1   (iRows, iCols) )   
        allocate( oS3M_Vars(iID)%a2dTaC_MeanDaysSuppressMelt  (iRows, iCols) ) 
        
        ! Finish to allocate dynamic variable(s)
        call mprintf(.true., iINFO_Main, ' Allocating dynamic variable(s) ... OK')
        !------------------------------------------------------------------------------------ 
        
    end subroutine S3M_Vars_Allocate
    !------------------------------------------------------------------------------------
        
    !------------------------------------------------------------------------------------
    ! Subroutine for initializing var(s) with default value
    subroutine S3M_Vars_InitDefault(iID, a1iGlaciers_ID, iNGlaciers)

        !------------------------------------------------------------------------------------ 
        ! Variable(s)
        integer(kind = 4)               :: iID
        integer(kind = 4), dimension(:) :: a1iGlaciers_ID
        integer(kind = 4)    :: iNGlaciers
        !------------------------------------------------------------------------------------
	
        !------------------------------------------------------------------------------------
        ! Start to initialize scalar variable(s)
        call mprintf(.true., iINFO_Main, ' Initialize scalar variable(s) ... ')
        oS3M_Vars(iID)%bFileForcingTimeSeries = .false.
        
        oS3M_Vars(iID)%dVErr = 0.0
        
        ! Finish to initialize scalar variable(s)
        call mprintf(.true., iINFO_Main, ' Initialize scalar variable(s) ... OK ')
        !------------------------------------------------------------------------------------ 

        !------------------------------------------------------------------------------------
        ! Start to initialize static variable(s)
        call mprintf(.true., iINFO_Main, ' Initialize static variable(s) ... ')

        ! Static land-data variable(s)
        oS3M_Vars(iID)%a2dLon = 0.0
        oS3M_Vars(iID)%a2dLat = 0.0
        oS3M_Vars(iID)%a2iMask = 0.0
        oS3M_Vars(iID)%a2dDEM = 0.0
        oS3M_Vars(iID)%a2dAreaCell = 0.0
        oS3M_Vars(iID)%a2iGlacierMask = 0
        oS3M_Vars(iID)%a2iGlaciers_ID = 0
        oS3M_Vars(iID)%a1iGlaciers_ID = a1iGlaciers_ID
        oS3M_Vars(iID)%a2dGlacierDebris = 0
        oS3M_Vars(iID)%a2dGlacierPivotTable = 0
        oS3M_Vars(iID)%iNumberGlaciers = iNGlaciers
        
        ! Static snow variable(s)
        oS3M_Vars(iID)%a2dArctUp = 0.0
        
        ! Static indexes grids (land-forcing)
        oS3M_Vars(iID)%a2iXIndex = 0
        oS3M_Vars(iID)%a2iYIndex = 0
        
        ! Finish to initialize static variable(s)
        call mprintf(.true., iINFO_Main, ' Initialize static variable(s) ... OK ')
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Start to initialize dynamic model variable(s)
        call mprintf(.true., iINFO_Main, ' Initialize dynamic variable(s) ... ')
        
        ! Time information
        oS3M_Vars(iID)%sTimeStep = ''
        oS3M_Vars(iID)%iTime = -1
        
        ! Dynamic forcing variable(s)
        oS3M_Vars(iID)%a2dPrecip = 0.0          ! Total precipitation   [mm]        --> Forcing
        oS3M_Vars(iID)%a2dTa = 0.0              ! Air Temperature       [C]         --> Forcing
        oS3M_Vars(iID)%a2dK = 0.0               ! Incoming radiation    [W/m^2]     --> Forcing
        oS3M_Vars(iID)%a2dRHum = 0.0            ! Relative humidity     [%]         --> Forcing
        oS3M_Vars(iID)%a2dSHeight = -9999.0     ! Snow Height           [cm]        --> Forcing for snow physics
        oS3M_Vars(iID)%a2dSKernel = -9999.0     ! Snow Kernel           [0, 1]      --> Forcing for snow physics
        
        ! Dynamic updating variable(s)
        oS3M_Vars(iID)%a2dSCA = -9999.0         ! Snow Cover Area       [-1, -3]    --> Updating for snow physics
        oS3M_Vars(iID)%a2dSQA = -9999.0         ! Snow Quality          [0, -1]     --> Updating for snow physics

        ! Dynamic snow variable(s)
        oS3M_Vars(iID)%a2iAge = 0           
        oS3M_Vars(iID)%a2dSWE = 0.0   
        oS3M_Vars(iID)%a2dRhoS = 0.0
        oS3M_Vars(iID)%a2dRhoS0 = 0.0
        oS3M_Vars(iID)%a2dAlbedo_Snow = 0.0
        oS3M_Vars(iID)%a2dMelting = 0.0
        oS3M_Vars(iID)%a2dMeltingSc = 0.0
        oS3M_Vars(iID)%a2dMeltingDayCum = 0.0
        oS3M_Vars(iID)%a3dTaC_Days1 = 0.0
        oS3M_Vars(iID)%a2dSnowFall = 0.0
        oS3M_Vars(iID)%a2dSnowFallDayCum = 0.0
        oS3M_Vars(iID)%a2dMaskS = 0.0
        oS3M_Vars(iID)%a2dIceThick = 0.0
        oS3M_Vars(iID)%a2dSWEass = 0.0
        oS3M_Vars(iID)%a2dSWEassHist = 0.0
        
        oS3M_Vars(iID)%a2dOutflow = 0.0
        oS3M_Vars(iID)%a2dRainFall = 0.0
        oS3M_Vars(iID)%a2dSWE_D = 0.0
        oS3M_Vars(iID)%a2dRho_D = 0.0        
        oS3M_Vars(iID)%a2dSWE_W = 0.0
        oS3M_Vars(iID)%a2dTheta_W = 0.0
        oS3M_Vars(iID)%a2dRefreezingS = 0.0
        oS3M_Vars(iID)%a2dMeltingG = 0.0
        oS3M_Vars(iID)%a2dH_D = 0.0
        oS3M_Vars(iID)%a2dH_S = 0.0
        oS3M_Vars(iID)%a2dSSA = 0.0
        oS3M_Vars(iID)%a2dPerm = 0.0
        oS3M_Vars(iID)%a2dReff = 0.0 
        oS3M_Vars(iID)%a2dMeltSIncRad = 0.0 
        oS3M_Vars(iID)%a2dMeltSTemp = 0.0 
        oS3M_Vars(iID)%a2dMeltGIncRad = 0.0 
        oS3M_Vars(iID)%a2dMeltGTemp = 0.0 
        oS3M_Vars(iID)%a2dMeltingSRadc = 0.0 
        oS3M_Vars(iID)%a2dMeltingGCumWY = 0.0
        oS3M_Vars(iID)%a2dChangeThickness = 0.0
        
        oS3M_Vars(iID)%a2dTaC_MeanDays1 = 0.0 
        oS3M_Vars(iID)%a2dTaC_MeanDaysSuppressMelt = 0.0          
        
        ! Finish to initialize dynamic variable(s)
        call mprintf(.true., iINFO_Main, ' Initialize dynamic variable(s) ... OK ')
        !------------------------------------------------------------------------------------ 
        
    end subroutine S3M_Vars_InitDefault
    !------------------------------------------------------------------------------------ 
    
end module S3M_Module_Vars_Manager
!------------------------------------------------------------------------------------

! S3M VARIABLE(S) DEFINITION

TYPE S3M_Type_Vars
    
    ! S3M STATIC VARIABLE(S) - SCALAR
    ! Land variable(s)
    integer(kind = 4)                   :: iDomainPixels
    real(kind = 4)                      :: dDxM, dDyM
    real(kind = 4)                      :: dDEMMax, dDEMMin, dDEMStepMean
    real(kind = 4)                      :: dDomainArea
    
    ! Integration variable(s)
    integer(kind = 4)                   :: iDtIntegr, iDtIntegrPStep
    
    ! S3M STATIC VARIABLE(S) - SCALAR
    integer(kind = 4)                   :: iNumberGlaciers    
    
    ! S3M STATIC VARIABLE(S) - 1D, 2D
    ! Land variable(s)
    integer(kind = 4),   allocatable, dimension (:, :)   :: a2iXIndex
    integer(kind = 4),   allocatable, dimension (:, :)   :: a2iYIndex
    real(kind = 4),      allocatable, dimension (:)      :: a1dFCN
    real(kind = 4),      allocatable, dimension (:, :)   :: a2dArctUp
    real(kind = 4),      allocatable, dimension (:, :)   :: a2dLon, a2dLat
    real(kind = 4),      allocatable, dimension (:, :)   :: a2dDem    
    real(kind = 4),      allocatable, dimension (:, :)   :: a2dAreaCell      
    integer(kind = 4),   allocatable, dimension (:, :)   :: a2iGlaciers_ID
    integer(kind = 4),   allocatable, dimension (:)      :: a1iGlaciers_ID
    integer(kind = 4),   allocatable, dimension (:, :)   :: a2iMask
    integer(kind = 4),   allocatable, dimension (:, :)   :: a2iGlacierMask
    real(kind = 4),      allocatable, dimension (:, :)   :: a2dGlacierDebris
    real(kind = 4),      allocatable, dimension (:, :)   :: a2dGlacierPivotTable

    ! S3M DYNAMIC VARIABLE(S) - SCALAR 
    ! Time variable(s)
    integer(kind = 4)                   :: iTime
    integer(kind = 4)                   :: iFlagSWEassExistingMap, iDayNud
    character(len = 19)                 :: sTimeStep
    character(len = 2)                  :: sTimeSWEass

    ! Volume variable(s)
    real(kind = 4)                      :: dVErr, dVarETTot
    
    ! S3M DYNAMIC VARIABLE(S) - 1D, 2D, 3D
    ! Forcing variable(s)
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dPrecip, a2dTa
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dK, a2dRHum
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dAlbedo
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dSHeight
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dSKernel

    ! Updating variable(s)
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dSCA, a2dSQA
    
    ! S3M variable(s)
    integer(kind = 4), allocatable,  dimension(:, :)     :: a2iAge
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dSWE, a2dRhoS
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dRhoS0
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dAlbedo_Snow
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dSnowFall
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dSnowFallDayCum
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dMeltingDayCum
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dMelting
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dMeltingSc
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dMaskS
    real(kind = 4), allocatable,     dimension(:, :, :)  :: a3dTaC_Days1
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dIceThick
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dSWEass
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dSWEassHist
    
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dOutflow
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dRainFall
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dSWE_D
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dRho_D
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dSWE_W
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dTheta_W
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dRefreezingS
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dMeltingG    
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dH_D  
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dH_S  
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dSSA   
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dPerm  
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dReff
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dMeltSIncRad
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dMeltSTemp
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dMeltGIncRad
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dMeltGTemp
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dMeltingSRadc
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dMeltingGCumWY
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dChangeThickness
    
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dTaC_MeanDays1
    real(kind = 4), allocatable,     dimension(:, :)     :: a2dTaC_MeanDays10
    
    ! Time-Series variable(s)
    logical                             :: bFileForcingTimeSeries
    integer(kind = 4)                   :: iFileUnitTSQ, iFileUnitTSVDam
    logical                             :: bFileUnitTSQ, bFileUnitTSVDam
    
    
END TYPE S3M_Type_Vars
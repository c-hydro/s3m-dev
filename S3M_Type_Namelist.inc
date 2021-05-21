
TYPE S3M_Type_Namelist
  
    ! Logical
    logical                 :: bGridCheck

    
    ! Init
    integer(kind = 4)       :: iFlagDebugSet, iFlagDebugLevel
    integer(kind = 4)       :: iFlagTypeData_IceThickness_Gridded
    integer(kind = 4)       :: iFlagTypeData_Forcing_Gridded
    integer(kind = 4)       :: iFlagTypeData_Updating_Gridded
    integer(kind = 4)       :: iFlagTypeData_Ass_SWE_Gridded
    integer(kind = 4)       :: iFlagRestart
    integer(kind = 4)       :: iFlagSnowAssim, iFlagSnowAssim_SWE
    integer(kind = 4)       :: iFlagIceMassBalance, iFlagThickFromTerrData
    integer(kind = 4)       :: iFlagGlacierDebris
    integer(kind = 4)       :: iFlagGrid
    integer(kind = 4)       :: iSimLength, iDtModel
    integer(kind = 4)       :: iNData, iDaySteps, iDtData_Forcing
    integer(kind = 4)       :: iDtData_Updating, iDtData_AssSWE
    integer(kind = 4)       :: iDtData_Output
    integer(kind = 4)       :: iScaleFactor_Forcing
    integer(kind = 4)       :: iScaleFactor_Update
    integer(kind = 4)       :: iScaleFactor_SWEass
    integer(kind = 4)       :: iRowsL, iColsL
    integer(kind = 4)       :: iRowsF, iColsF 
    integer(kind = 4)       :: iRowsPivot, iColsPivot     
    integer(kind = 4)       :: iGlacierValue
    integer(kind = 4)       :: iSWEassInfluence
    integer(kind = 4)       :: iNTime
    integer(kind = 4), allocatable,         dimension(:) :: a1iDimsForcing
    integer(kind = 4)       :: iFlagOutputMode
    integer(kind = 4)       :: iFlagAssOnlyPos
    integer(kind = 4)       :: iDaysAvgTSuppressMelt
    
    ! Real
    real(kind = 4)          :: dXLLCornerL, dYLLCornerL
    real(kind = 4)          :: dXCellSizeL, dYCellSizeL, dNoDataL
    real(kind = 4)          :: dXLLCornerF, dYLLCornerF
    real(kind = 4)          :: dXCellSizeF, dYCellSizeF, dNoDataF
    real(kind = 4)          :: dRhoSnowMax, dRhoSnowMin
    real(kind = 4)          :: dRhoSnowFresh, dSnowQualityThr
    real(kind = 4)          :: dMeltingTRef, dIceMeltingCoeff, dWeightSWEass
    real(kind = 4)          :: dRhoW
    real(kind = 4)          :: dRefreezingSc
    real(kind = 4)          :: dModFactorRadS
    real(kind = 4)          :: dDebrisThreshold
    
    real(kind = 4), allocatable,            dimension(:) :: a1dGeoForcing
    real(kind = 4), allocatable,            dimension(:) :: a1dResForcing
    real(kind = 4), allocatable,            dimension(:) :: a1dArctUp
    real(kind = 4), allocatable,            dimension(:) :: a1dAltRange
    
    ! Character
    character(len = 19)     :: sTimeStart
    character(len = 19)     :: sTimeRestart
    character(len = 256)    :: sDomainName
    character(len = 256)    :: sPathData_Static_Gridded 
    character(len = 256)    :: sPathData_Forcing_Gridded
    character(len = 256)    :: sPathData_Output_Gridded
    character(len = 256)    :: sPathData_Restart_Gridded
    character(len = 256)    :: sPathData_Updating_Gridded
    character(len = 256)    :: sPathData_IceThickness_Gridded
    character(len = 256)    :: sPathData_SWE_Assimilation_Gridded
    character(len = 1)      :: sPathBar
    character(len = 700)    :: sCommandUnzipFile
    character(len = 700)    :: sCommandZipFile
    character(len = 700)    :: sCommandRemoveFile
    character(len = 700)    :: sCommandCreateFolder
    character(len = 10)     :: sReleaseDate
    character(len = 700)    :: sAuthorNames
    character(len = 5)      :: sReleaseVersion
    character(len = 2)       :: sWYstart
    
END TYPE S3M_Type_Namelist
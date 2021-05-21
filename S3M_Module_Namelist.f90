!--------------------------------------------------------------------------------  
! File:   S3M_Module_Namelist.f90
! Author(s): Fabio Delogu, Francesco Silvestro, Simone Gabellani, Francesco Avanzi.
!
! Created on May, 20 2014, 9:57 AM
! Last update on October 26, 2020 03:14 PM
!
! Module to allocate and read namelist
!--------------------------------------------------------------------------------

!--------------------------------------------------------------------------------
! Module Header
module S3M_Module_Namelist
    
    !--------------------------------------------------------------------------------
    ! External module(s) and implicit none
    use S3M_Module_Tools_Debug          ! to import global variable(s) declaration

    implicit none
    
    include "S3M_Type_Namelist.inc"
    integer, parameter :: iMaxDomain = 1
    type(S3M_Type_Namelist), dimension(iMaxDomain) :: oS3M_Namelist
    save oS3M_Namelist
    !--------------------------------------------------------------------------------
    
contains 
    
    !--------------------------------------------------------------------------------
    ! Subroutine to read namelist
    subroutine S3M_Namelist_Read(sFileInfo, iArgsType, oS3M_Namelist_Init) 
        
        !--------------------------------------------------------------------------------
        ! External module(s) and implicit none
        implicit none

        type(S3M_Type_Namelist) oS3M_Namelist_Init
        !--------------------------------------------------------------------------------
        
        !--------------------------------------------------------------------------------
        ! Variable(s) declaration    
        integer(kind = 4)       :: iErr
        
        integer(kind = 4)       :: iArgsType
        character(len = 256)    :: sFileInfo
        logical                 :: bFileExist
        
        integer(kind = 4)       :: iFlagDebugSet, iFlagDebugLevel, iDebugUnitInit
        integer(kind = 4)       :: iFlagTypeData_Forcing_Gridded
        integer(kind = 4)       :: iFlagTypeData_Updating_Gridded
        integer(kind = 4)       :: iFlagTypeData_Ass_SWE_Gridded
        integer(kind = 4)       :: iFlagRestart
        integer(kind = 4)       :: iFlagSnowAssim
        integer(kind = 4)       :: iFlagSnowAssim_SWE
        integer(kind = 4)       :: iFlagIceMassBalance, iFlagThickFromTerrData, iFlagGlacierDebris
        integer(kind = 4)       :: iFlagGrid
        integer(kind = 4)       :: iFlagOutputMode
        integer(kind = 4)       :: iFlagAssOnlyPos
        integer(kind = 4)       :: iDaysAvgTSuppressMelt

        logical                 :: bGridCheck

        integer(kind = 4)       :: iSimLength, iDtModel
        integer(kind = 4)       :: iDtData_Forcing, iDtData_Output
        integer(kind = 4)       :: iDtData_Updating, iDtData_AssSWE
        integer(kind = 4)       :: iScaleFactor_Forcing, iScaleFactor_Update, iScaleFactor_SWEass

        integer(kind = 4)       :: iRowsL, iColsL
        
        real(kind = 4)          :: dXLLCornerL, dYLLCornerL, dXCellSizeL, dYCellSizeL, dNoDataL
        
        integer(kind = 4)       :: iRowsF, iColsF
        real(kind = 4)          :: dXLLCornerF, dYLLCornerF, dXCellSizeF, dYCellSizeF, dNoDataF
        integer(kind = 4)       :: iNTime
        real(kind = 4)          :: dRhoW
        integer(kind = 4)       :: iNData, iDaySteps
                
        integer(kind = 4)       :: iGlacierValue
        real(kind = 4)          :: dRhoSnowMax, dRhoSnowFresh, dSnowQualityThr, dRhoSnowMin
        real(kind = 4)          :: dMeltingTRef, dIceMeltingCoeff, dWeightSWEass
        integer(kind = 4)       :: iSWEassInfluence
        real(kind = 4)          :: dRefreezingSc, dModFactorRadS, dDebrisThreshold
        
        integer(kind = 4), target, dimension(2)     :: a1iDimsForcing
        real(kind = 4), target, dimension(2)        :: a1dGeoForcing
        real(kind = 4), target, dimension(2)        :: a1dResForcing
        
        real(kind = 4), target, dimension(4)        :: a1dArctUp
        real(kind = 4), target, dimension(3)        :: a1dAltRange
        
        character(len = 19)     :: sTimeStart
        character(len = 19)     :: sTimeRestart
        
        character(len = 19)     :: sTimeStartLong
        character(len = 19)     :: sTimeRestartLong
        
        character(len = 256)    :: sDomainName

        character(len = 256)    :: sPathData_Static_Gridded
        character(len = 256)    :: sPathData_Forcing_Gridded
        character(len = 256)    :: sPathData_Output_Gridded
        character(len = 256)    :: sPathData_Restart_Gridded
        character(len = 256)    :: sPathData_Updating_Gridded
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
        
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Read namelist(s)
        namelist /S3M_Namelist/         sDomainName, &
                                        iFlagDebugSet, iFlagDebugLevel, &
                                        iFlagTypeData_Forcing_Gridded, &
                                        iFlagTypeData_Updating_Gridded, &
                                        iFlagTypeData_Ass_SWE_Gridded, &
                                        iFlagRestart, &
                                        iFlagSnowAssim, iFlagSnowAssim_SWE, iFlagIceMassBalance, iFlagThickFromTerrData, &
                                        iFlagGlacierDebris, iFlagOutputMode, iFlagAssOnlyPos, &
                                        a1dGeoForcing, a1dResForcing, a1iDimsForcing, &
                                        iSimLength, iDtModel, iDtData_Forcing, iDtData_Output, &
                                        iDtData_Updating, iDtData_AssSWE, &
                                        iScaleFactor_Forcing, iScaleFactor_Update, iScaleFactor_SWEass, &
                                        sTimeStart, sTimeRestart, &
                                        sPathData_Static_Gridded, &
                                        sPathData_Forcing_Gridded, &
                                        sPathData_Output_Gridded, &
                                        sPathData_Restart_Gridded, &
                                        sPathData_Updating_Gridded, &
                                        sPathData_SWE_Assimilation_Gridded
                                        
        namelist /S3M_Snow/             a1dArctUp, a1dAltRange, &
                                        iGlacierValue, dRhoSnowFresh, dRhoSnowMax, dRhoSnowMin, dSnowQualityThr, &
                                        dMeltingTRef, dIceMeltingCoeff, iSWEassInfluence, dWeightSWEass, dRefreezingSc, &
                                        dModFactorRadS, sWYstart, dDebrisThreshold, iDaysAvgTSuppressMelt
                                       
        namelist /S3M_Constants/        dRhoW
                                        
        namelist /S3M_Command/          sCommandZipFile, &
                                        sCommandUnzipFile, &
                                        sCommandRemoveFile, &
                                        sCommandCreateFolder
        
        namelist /S3M_Info/             sReleaseDate, &
                                        sAuthorNames, &
                                        sReleaseVersion
        !--------------------------------------------------------------------------------
        
        !--------------------------------------------------------------------------------
        ! Set defaults for namelist variables
        iFlagTypeData_Forcing_Gridded = -9999;
        iFlagTypeData_Updating_Gridded = -9999;
        iFlagTypeData_Ass_SWE_Gridded = -9999;
        iFlagDebugSet = -9999; iFlagDebugLevel = -9999;
        iFlagRestart = -9999;
        iScaleFactor_Forcing = -9999; 
        iScaleFactor_Update = -9999;
        iScaleFactor_SWEass = -9999; 
        iFlagIceMassBalance = -9999; iFlagThickFromTerrData = -9999; iFlagGlacierDebris = -9999;
        iFlagOutputMode = -9999;
        iFlagAssOnlyPos = -9999;
        iFlagSnowAssim = -9999; iFlagSnowAssim_SWE = -9999;
        a1dGeoForcing = -9999.0; a1dResForcing = -9999.0; a1iDimsForcing = -9999; 
        iSimLength = -9999; iDtModel = -9999; iDtData_Forcing = -9999; iDtData_Output = -9999;
        iDtData_Updating = -9999; iDtData_AssSWE = -9999;
        sTimeStart = ""; sTimeRestart = "";
        sDomainName = "";
        sPathData_Static_Gridded = "";
        sPathData_Forcing_Gridded = "";
        sPathData_Output_Gridded = ""; 
        sPathData_Restart_Gridded = "";
        sPathData_Updating_Gridded = "";
        sPathData_SWE_Assimilation_Gridded = "";
        sWYstart = "";
        
        iNTime = 0; iNData = 0;
        
        a1dArctUp = -9999.0; a1dAltRange = -9999.0;
        iGlacierValue = -9999; dRhoSnowFresh = -9999.0; dRhoSnowMax = -9999.0; dRhoSnowMin = -9999.0; 
        dSnowQualityThr = -9999.0; dMeltingTRef = -9999.0; dIceMeltingCoeff = -9999.0;
        iSWEassInfluence = 0; dWeightSWEass = -9999.0; dRefreezingSc = -9999.0; dModFactorRadS = -9999.0;
        dDebrisThreshold = -9999; iDaysAvgTSuppressMelt = -9999;
        
        dRhoW = -9999.0;
        
        sCommandZipFile = ""; sCommandUnzipFile = ""; sCommandRemoveFile = ""; sCommandCreateFolder = ""
        
        sReleaseDate = ""; sAuthorNames = ""; sReleaseVersion = "";
        !--------------------------------------------------------------------------------  
                                            
        !--------------------------------------------------------------------------------
        ! Checking information file availability
        inquire (file = trim(sFileInfo), exist = bFileExist)
        if ( .not. bFileExist ) then
            write(6, *) "No file info found ", trim(sFileInfo)
            stop "Stopped"
        else
            ! Read namelist section(s)
            open(20,file = trim(sFileInfo))
            read(20, S3M_Namelist, iostat=iErr); rewind(20)
            read(20, S3M_Snow, iostat=iErr); rewind(20)
            read(20, S3M_Command, iostat=iErr); rewind(20)
            read(20, S3M_Constants, iostat=iErr); rewind(20)
            read(20, S3M_Info, iostat=iErr); rewind(20)
 
            close(20) 
        endif
        !-----------------------------------------------------------------------------------
                
        !--------------------------------------------------------------------------------
        ! Set debug level
        call S3M_Tools_Debug_SetLevel(iFlagDebugSet, iFlagDebugLevel)
        ! Set file unit debug
        call S3M_Tools_Debug_SetUnit(80, 100, iDebugUnitInit)
        !--------------------------------------------------------------------------------

        !-----------------------------------------------------------------------------------
        ! S3M Version
        call mprintf(.true., iINFO_Basic, '************************************************************************')
        call mprintf(.true., iINFO_Basic, ' S3M Snow Multidata Mapping and Modeling (Version: '//trim(sReleaseVersion)//')')
        call mprintf(.true., iINFO_Basic, ' Author(s): '//trim(sAuthorNames))
        call mprintf(.true., iINFO_Basic, ' Release Date: '//trim(sReleaseDate))
        call mprintf(.true., iINFO_Basic, '************************************************************************')
        !-----------------------------------------------------------------------------------        
        
        !-----------------------------------------------------------------------------------
        ! Time long definition
        sTimeStartLong = sTimeStart(1:4)//"-"//sTimeStart(5:6)//"-"//sTimeStart(7:8)//"_" &
                         //sTimeStart(9:10)//":"//sTimeStart(11:12)//":"//"00"
        
        sTimeRestartLong = sTimeRestart(1:4)//"-"//sTimeRestart(5:6)//"-"//sTimeRestart(7:8)//"_" &
                           //sTimeRestart(9:10)//":"//sTimeRestart(11:12)//":"//"00"
        
        !----------------------------------------------------------------------------------
                                   
        !------------------------------------------------------------------------------------------
        ! Compute simulation time length
        iNTime = (iSimLength)*3600./nint(real(iDtModel))
        ! Compute data time steps
        iNData = int((iNTime*3600))/int(iDtData_Forcing)
        !------------------------------------------------------------------------------------------
 
        !--------------------------------------------------------------------------------
        ! Namelist definition
        call mprintf(.true., iINFO_Main, ' Read Namelist ... ')
        
        ! Model dims (initialization)
        oS3M_Namelist_Init%iRowsL = 0
        oS3M_Namelist_Init%iColsL = 0
        oS3M_Namelist_Init%iRowsF = 0
        oS3M_Namelist_Init%iColsF = 0
        oS3M_Namelist_Init%iRowsPivot = 0
        oS3M_Namelist_Init%iColsPivot = 0
        
        ! LSM constant(s)
        oS3M_Namelist_Init%dRhoW = dRhoW

        oS3M_Namelist_Init%iDaySteps = 0

        ! Flag(s) info
        oS3M_Namelist_Init%iFlagTypeData_Forcing_Gridded = iFlagTypeData_Forcing_Gridded
        oS3M_Namelist_Init%iFlagTypeData_Updating_Gridded = iFlagTypeData_Updating_Gridded
        oS3M_Namelist_Init%iFlagTypeData_Ass_SWE_Gridded = iFlagTypeData_Ass_SWE_Gridded
        oS3M_Namelist_Init%iFlagRestart = iFlagRestart
        oS3M_Namelist_Init%iFlagSnowAssim = iFlagSnowAssim
        oS3M_Namelist_Init%iFlagSnowAssim_SWE = iFlagSnowAssim_SWE
        oS3M_Namelist_Init%iFlagIceMassBalance = iFlagIceMassBalance
        oS3M_Namelist_Init%iFlagThickFromTerrData = iFlagThickFromTerrData
        oS3M_Namelist_Init%iFlagGlacierDebris = iFlagGlacierDebris
        oS3M_Namelist_Init%iFlagOutputMode = iFlagOutputMode
        oS3M_Namelist_Init%iFlagAssOnlyPos = iFlagAssOnlyPos
        oS3M_Namelist_Init%iFlagGrid = -9999
        
        ! Geographical land and forcing info
        oS3M_Namelist_Init%bGridCheck = .false.
        
        oS3M_Namelist_Init%dXLLCornerL = 0.0
        oS3M_Namelist_Init%dYLLCornerL = 0.0
        oS3M_Namelist_Init%dXCellSizeL = 0.0
        oS3M_Namelist_Init%dYCellSizeL = 0.0
        oS3M_Namelist_Init%dNoDataL    = 0.0
        oS3M_Namelist_Init%dXLLCornerF = 0.0
        oS3M_Namelist_Init%dYLLCornerF = 0.0
        oS3M_Namelist_Init%dXCellSizeF = 0.0
        oS3M_Namelist_Init%dYCellSizeF = 0.0
        oS3M_Namelist_Init%dNoDataF    = 0.0
        
        oS3M_Namelist_Init%a1dGeoForcing  = a1dGeoForcing
        oS3M_Namelist_Init%a1dResForcing  = a1dResForcing
        oS3M_Namelist_Init%a1iDimsForcing = a1iDimsForcing
        
        ! S3M parameter(s) and constant(s)
        oS3M_Namelist_Init%a1dArctUp     = a1dArctUp
        oS3M_Namelist_Init%a1dAltRange   = a1dAltRange
        oS3M_Namelist_Init%iGlacierValue = iGlacierValue
        oS3M_Namelist_Init%sWYstart = sWYstart
        
        ! backward compatibility with older version of info file (without fresh snow density reference)
        if (dRhoSnowFresh .eq. -9999) then                  
            oS3M_Namelist_Init%dRhoSnowFresh  = 100
        else
            oS3M_Namelist_Init%dRhoSnowFresh  = dRhoSnowFresh
        endif   
        oS3M_Namelist_Init%dRhoSnowMax     = dRhoSnowMax
        oS3M_Namelist_Init%dRhoSnowMin     = dRhoSnowMin
        oS3M_Namelist_Init%dSnowQualityThr = dSnowQualityThr
        
        ! backward compatibility with older version of info file (without melting t reference)
        if (dMeltingTRef .eq. -9999) then                   
            oS3M_Namelist_Init%dMeltingTRef = 1
        else
            oS3M_Namelist_Init%dMeltingTRef = dMeltingTRef
        endif
        
        ! backward compatibility with older version of info file (without melting t reference)
        if (dIceMeltingCoeff .eq. -9999) then                   
            oS3M_Namelist_Init%dIceMeltingCoeff = 3
        else
            oS3M_Namelist_Init%dIceMeltingCoeff = dIceMeltingCoeff
        endif
        
        oS3M_Namelist_Init%iSWEassInfluence = iSWEassInfluence
        oS3M_Namelist_Init%dWeightSWEass    = dWeightSWEass
        oS3M_Namelist_Init%dRefreezingSc    = dRefreezingSc
        oS3M_Namelist_Init%dModFactorRadS   = dModFactorRadS
        oS3M_Namelist_Init%dDebrisThreshold   = dDebrisThreshold
        oS3M_Namelist_Init%iDaysAvgTSuppressMelt = iDaysAvgTSuppressMelt
        
        ! Time, dt and step(s) info
        oS3M_Namelist_Init%iSimLength       = iSimLength
        oS3M_Namelist_Init%iDtModel         = iDtModel
        oS3M_Namelist_Init%iDtData_Forcing  = iDtData_Forcing
        oS3M_Namelist_Init%iDtData_Updating = iDtData_Updating
        oS3M_Namelist_Init%iDtData_AssSWE   = iDtData_AssSWE
        
        ! Time reference info
        oS3M_Namelist_Init%sTimeStart     = sTimeStartLong
        oS3M_Namelist_Init%sTimeRestart   = sTimeRestartLong
        oS3M_Namelist_Init%iDtData_Output = iDtData_Output 
        
        ! Data info
        oS3M_Namelist_Init%iNData               = iNData
        oS3M_Namelist_Init%iNTime               = iNTime
        oS3M_Namelist_Init%iScaleFactor_Forcing = iScaleFactor_Forcing
        oS3M_Namelist_Init%iScaleFactor_Update  = iScaleFactor_Update
        oS3M_Namelist_Init%iScaleFactor_SWEass  = iScaleFactor_SWEass
        
        ! Domain name
        oS3M_Namelist_Init%sDomainName = sDomainName

        ! Path(s) info
        oS3M_Namelist_Init%sPathData_Static_Gridded           = sPathData_Static_Gridded
        oS3M_Namelist_Init%sPathData_Forcing_Gridded          = sPathData_Forcing_Gridded
        oS3M_Namelist_Init%sPathData_Output_Gridded           = sPathData_Output_Gridded
        oS3M_Namelist_Init%sPathData_Restart_Gridded          = sPathData_Restart_Gridded
        oS3M_Namelist_Init%sPathData_SWE_Assimilation_Gridded = sPathData_SWE_Assimilation_Gridded
        
        ! backward compatibility with older version of info file (without updating gridded)
        if (sPathData_Updating_Gridded == "") then       
            oS3M_Namelist_Init%sPathData_Updating_Gridded = sPathData_Forcing_Gridded
        else
            oS3M_Namelist_Init%sPathData_Updating_Gridded = sPathData_Updating_Gridded
        endif
        
        ! Command line
        oS3M_Namelist_Init%sPathBar             = sPathBar
        oS3M_Namelist_Init%sCommandZipFile      = sCommandZipFile
        oS3M_Namelist_Init%sCommandUnzipFile    = sCommandUnzipFile
        oS3M_Namelist_Init%sCommandRemoveFile   = sCommandRemoveFile
        oS3M_Namelist_Init%sCommandCreateFolder = sCommandCreateFolder
    
        ! Model info
        oS3M_Namelist_Init%sReleaseDate    = sReleaseDate
        oS3M_Namelist_Init%sAuthorNames    = sAuthorNames
        oS3M_Namelist_Init%sReleaseVersion = sReleaseVersion
         
        ! Info
        call mprintf(.true., iINFO_Main, ' Read Namelist ... OK')
        !--------------------------------------------------------------------------------
        
    end subroutine S3M_Namelist_Read
    !--------------------------------------------------------------------------------

end module S3M_Module_Namelist
!--------------------------------------------------------------------------------
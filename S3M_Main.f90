!***********************************************************************************
! S3M
!***********************************************************************************
!
! brief		S3M Snow Multidata Mapping and Modeling (Boni et al. 2010, Avanzi et al. 2022)
!
! history	FRANCESCO AVANZI (CIMAFOUNDATION) 
!+ 		09/02/2023
!+		v5p2r0
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!-----------------------------------------------------------------------------------
! DESCRIPTION
!-----------------------------------------------------------------------------------
! CIMA Foundation's cryospheric model.
! S3M is a spatially distributed model combining several sources of data to estimate snowpack and glacier state.
! The snow component follows a one-layer, modified-degree-day approach and explicitly separates the dry and wet constituents of 
! snow (De Michele et al. 2013, Avanzi et al. 2015). Thus, it provides spatially explicitly estimates of all snowpack bulk 
! properties (SWE, density, snow depth), bulk liquid water content, melt and albedo, as well as snowpack runoff. 
! The glacier component optionally considers mass balance through the deltah parametrization (Huss et al. 2010). 

! Input data are: 
! 1) Land data, (mandatory)
! 2) Meteorological observations, (mandatory)
! 3) SCA satellite images, (optional)
! 4) SWE independent estimates for assimilation, (optional)
! 5) Ice thickness and a number of other ancillary glacier data, (optional).
! Output results are (among others):
! 1) Snowpack runoff,
! 2) SWE (dry and wet), snow density (dry and wet), snow depth, bulk liquid water content
! 3) Snowfall, Rainfall, and Precipitation rates, as well as fresh-snow density
! 4) Snow age, albedo, melt, and snowpack runoff
! 5) Ice thickness.

! S3M 5.2.0 is distributed through the CIMA Research Foundation - Department of Hydrology and Hydraulics GitHub repo at
! https://github.com/c-hydro 
! A complete manual regarding model installation, run, as well as pre- and post-processing can be found at: 
! https://gmd.copernicus.org/articles/15/4853/2022/
! 
!-----------------------------------------------------------------------------------
! COMMAND LINE (example)
!-----------------------------------------------------------------------------------
! ./S3M_v5p2.x --> parameter(s): S3M_info_domain_DataDa_DataA.txt 
!
! All compiling and debugging options for S3M are described at https://github.com/c-hydro/hmc-dev
!
! 3) DEBUG USING GNUPLOT (gnufor2):
! call surf(dble(VAR),pm3d='pm3d implicit map', palette='rgbformulae 31, -11, 32')
!
! 4) DEBUG USING FORTRAN+MATLAB: [ var*8 = dble(var*4) ]
! call debug_2dVar(dble(a2dVarTa), iRows, iCols, 1) + debug_2dVar.m
!
!******************************************************************************************
!
!------------------------------------------------------------------------------------------
! S3M main 
program S3M_Main
    
    !------------------------------------------------------------------------------------------
    ! Use and implicit none
    use S3M_Module_Args,                    only:   S3M_Args_Read

    use S3M_Module_Namelist,                only:   oS3M_Namelist, S3M_Namelist_Read

    use S3M_Module_Info_Gridded,            only:   S3M_Info_Gridded_GetDims_Static, &
                                                    S3M_Info_Gridded_GetDims_Forcing, &
                                                    S3M_Info_Gridded_GetDims_Ancillary, &
                                                    S3M_Info_Gridded_GetGeo_Static, &
                                                    S3M_Info_Gridded_GetGeo_Forcing
                                                    
    use S3M_Module_Info_Time

    use S3M_Module_Vars_Loader,             only:   S3M_Type_Vars, oS3M_Vars
    use S3M_Module_Vars_Manager,            only:   S3M_Vars_InitDefault, S3M_Vars_Allocate

    use S3M_Module_Tools_Time,              only:   S3M_Tools_Time_GetNewDate, &
                                                    S3M_Tools_Time_Printer
    use S3M_Module_Tools_Debug
    
    use S3M_Module_Data_Static_Gridded,     only:   S3M_Data_Static_Gridded_Cpl
    use S3M_Module_Data_Forcing_Gridded,    only:   S3M_Data_Forcing_Gridded_Cpl
    use S3M_Module_Data_Output_Gridded,     only:   S3M_Data_Output_Gridded_Cpl
    use S3M_Module_Data_Restart_Gridded,    only:   S3M_Data_Restart_Gridded_Cpl
    use S3M_Module_Data_Updating_Gridded,   only:   S3M_Data_Updating_Gridded_Cpl
    use S3M_Module_Data_AssSWE_Gridded,     only:   S3M_Data_AssSWE_Gridded_Cpl

    use S3M_Module_Phys,                    only:   S3M_Phys_Cpl

    implicit none
    !------------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------------
    ! Variable(s) declaration
    integer(kind = 4)               :: iArgsType
    character(len = 700)            :: sFileInfo

    integer(kind = 4)               :: iID
    integer(kind = 4)               :: iColsL, iRowsL, iColsF, iRowsF, iColsPivot, iRowsPivot, iTime, iT
    integer(kind = 4)               :: iDaySteps

    integer(kind = 4)               :: iDtModel
    integer(kind = 4)               :: iDtData_Forcing
    integer(kind = 4)               :: iDtData_Output
    integer(kind = 4)               :: iDtData_Updating, iDtData_AssSWE
    integer(kind = 4)               :: iSimLength, iNTime, iNData, iNGlaciers
    integer(kind = 4)               :: iFlagSnowAssim, iFlagSnowAssim_SWE
    
    character(len = 256)            :: sDomainName
    character(len = 19)             :: sTimeOld, sTimeForcing, sTimeRestart, sTimeNew
    character(len = 19)             :: sTimeOutput
    character(len = 19)             :: sTimeUpdating, sTimeAssSWE
    character(len = 256)            :: sTime, sNTime, sDtModel

    character(len=10)               :: sDateRunStart, sDateRunStep, sDateRunEnd
    character(len=12)               :: sTimeRunStart, sTimeRunStep, sTimeRunEnd
    
    integer(kind = 4), dimension(:), allocatable    :: a1iGlaciers_ID
    !------------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------------
    ! Initialize variable(s)
    sFileInfo = ""; iArgsType = -9999; iID = 1;
    iColsL = -9999; iRowsL = -9999; iColsF = -9999; iRowsF = -9999; iColsPivot = -9999; iRowsPivot = -9999;
    iTime = -9999; iDaySteps = -9999; iSimLength = 0;
    iNTime = -9999; iT = -9999; iNGlaciers = -9999;
    sTimeOld = ""; sTimeForcing = ""; sTimeRestart = ""; sTimeNew = ""; sTimeUpdating = "";
    sTimeAssSWE = ""; sTimeOutput = ""; iFlagSnowAssim = 0; iFlagSnowAssim_SWE = 0;
    !------------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------------
    ! Call time start
    call S3M_Tools_Time_Printer(sDateRunStart, sTimeRunStart)
    !------------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------------
    ! Read command line argument(s)
    call S3M_Args_Read(sFileInfo, iArgsType)
    !------------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------------
    ! Read namelist file
    call S3M_Namelist_Read(sFileInfo, iArgsType, oS3M_Namelist(iID))

    ! Start simulation message
    call mprintf(.true., iINFO_Basic, ' SIMULATION START :: Time: ' //trim(sDateRunStart)//' '//trim(sTimeRunStart) )
    
    ! Domain name
    sDomainName = oS3M_Namelist(iID)%sDomainName
              
    ! Flag
    iFlagSnowAssim = oS3M_Namelist(iID)%iFlagSnowAssim
    iFlagSnowAssim_SWE = oS3M_Namelist(iID)%iFlagSnowAssim_SWE
    
    ! Time definition(s)                
    sTimeOld = oS3M_Namelist(iID)%sTimeStart
    sTimeForcing = oS3M_Namelist(iID)%sTimeStart
    !sTimeOutput = oS3M_Namelist(iID)%sTimeStart
    sTimeRestart = oS3M_Namelist(iID)%sTimeRestart
    sTimeUpdating = oS3M_Namelist(iID)%sTimeStart
    sTimeAssSWE = oS3M_Namelist(iID)%sTimeStart
    
    call S3M_Tools_Time_GetNewDate(sTimeOutput, oS3M_Namelist(iID)%sTimeStart, &
        nint(real(oS3M_Namelist(iID)%iDtData_Output - real(oS3M_Namelist(iID)%iDtData_Forcing))))

    ! Model dt(s)
    iDtModel = oS3M_Namelist(iID)%iDtModel
    iDtData_Forcing = oS3M_Namelist(iID)%iDtData_Forcing
    iDtData_Output = oS3M_Namelist(iID)%iDtData_Output
    iDtData_Updating = oS3M_Namelist(iID)%iDtData_Updating
    iDtData_AssSWE = oS3M_Namelist(iID)%iDtData_AssSWE
    iSimLength = oS3M_Namelist(iID)%iSimLength
    
    ! Simulation step(s)
    iNTime = oS3M_Namelist(iID)%iNTime
    ! Data step(s)
    iNData = oS3M_Namelist(iID)%iNData
    !------------------------------------------------------------------------------------------
  
    !------------------------------------------------------------------------------------------
    ! Get static data dimension(s)
    call S3M_Info_Gridded_GetDims_Static(iID, iRowsL, iColsL, iRowsPivot, iColsPivot)
    ! Get forcing data dimension(s)
    call S3M_Info_Gridded_GetDims_Forcing(iID, iRowsF, iColsF)
    ! Get ancillary data dimension(s)
    call S3M_Info_Gridded_GetDims_Ancillary(iID, iRowsL, iColsL, iNGlaciers, a1iGlaciers_ID)

    ! Get time dimension(s)
    call S3M_Info_Time_GetDims(iID, iDaySteps)
    !------------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------------
    ! Allocate static and dynamic global variable(s)
    call S3M_Vars_Allocate( iID, & 
                            iRowsL, iColsL, & 
                            iRowsPivot, iColsPivot, &
                            iDaySteps, &
                            iNData, iNGlaciers)
                            
    ! Initialize static and dynamic global variable(s) using default value(s)
    call S3M_Vars_InitDefault(iID, a1iGlaciers_ID, iNGlaciers)
    !------------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------------
    ! Get land data static geographical information
    call S3M_Info_Gridded_GetGeo_Static(iID)
    ! Get land data forcing geographical information
    call S3M_Info_Gridded_GetGeo_Forcing(iID)
    !------------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------------
    ! Initialize static gridded variable(s)
    call S3M_Data_Static_Gridded_Cpl(iID, iRowsL, iColsL, iRowsPivot, iColsPivot)
    !------------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------------
    ! Initialize restart data
    call S3M_Data_Restart_Gridded_Cpl(iID, sTimeRestart, &
                                      1, iRowsL, 1, iColsL, &
                                      iDaySteps)
    !------------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------------
    ! Cycling on time period (and extra time steps)
    TimeLoop : do iTime = 1, iNTime
        
        !------------------------------------------------------------------------------------------
        ! Time step start
        write(sTime, *) iTime; write(sNTime, *) iNTime;  write(sDtModel, *) iDtModel
        call mprintf(.true., iINFO_Basic, ' ===== TIME STEP START :: Date, iT, iNTime, iDt, Clock :: ' &
                                          //sTimeOld//' '//trim(sTime)//' '// &
                                          trim(sNTime)//' '//trim(sDtModel)// ' ===== ')
        ! Clock step start                         
        call S3M_Tools_Time_Printer(sDateRunStep, sTimeRunStep)
        call mprintf(.true., iINFO_Basic, ' STEP CLOCK START :: Time: ' //trim(sDateRunStep)//' '//trim(sTimeRunStep) )
        
        ! Update time step
        oS3M_Vars(iID)%iTime = iTime; oS3M_Vars(iID)%sTimeStep = sTimeOld
        !------------------------------------------------------------------------------------------
       
        !------------------------------------------------------------------------------------------
        ! Forcing data step
        if(sTimeOld == sTimeForcing) then

            ! Get forcing gridded data
            call S3M_Data_Forcing_Gridded_Cpl( iID, sTimeForcing, &
                                       1, iRowsL, 1, iColsL, &
                                       1, iRowsF, 1, iColsF)
                                
            ! Update data forcing time
            call S3M_Tools_Time_GetNewDate(sTimeNew, sTimeForcing, nint(real(iDtData_Forcing)))
            sTimeForcing = sTimeNew
        endif
        !------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        ! Updating data step: used for snow depth assimilation + MODIS
        if(iFlagSnowAssim .eq. 1) then
            if(sTimeOld == sTimeUpdating) then
                
                ! Get updating gridded data
                call S3M_Data_Updating_Gridded_Cpl( iID, sTimeUpdating, &
                                       1, iRowsL, 1, iColsL)

                ! Update data forcing time
                call S3M_Tools_Time_GetNewDate(sTimeNew, sTimeUpdating, nint(real(iDtData_Updating)))
                sTimeUpdating = sTimeNew
            endif
        endif
        !------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        ! SWE assimilation from other sources
        if(iFlagSnowAssim_SWE .eq. 1) then
            if(sTimeOld == sTimeAssSWE) then
                
                ! Get updating gridded data
                call S3M_Data_AssSWE_Gridded_Cpl( iID, sTimeAssSWE, &
                                       1, iRowsL, 1, iColsL)
                
                ! SWE data forcing time
                call S3M_Tools_Time_GetNewDate(sTimeNew, sTimeAssSWE, nint(real(iDtData_AssSWE)))
                sTimeAssSWE = sTimeNew
                
            endif
        endif
        !------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        ! Physics step
        call S3M_Phys_Cpl(iID, &
                          1, iRowsL, 1, iColsL, &
                          iTime, iNTime, sTimeOld, &
                          iNData, &
                          iDaySteps)
        !------------------------------------------------------------------------------------------
                
        !------------------------------------------------------------------------------------------
        ! Output data step
        ! Save data output gridded
        if( (sTimeOld .eq. sTimeOutput) .and. (iDtData_Output .gt. 0) ) then
          
            ! Save output gridded data
            call S3M_Data_Output_Gridded_Cpl( iID, sTimeOutput, &
                                              1, iRowsL, 1, iColsL, &
                                              iTime, iDaySteps)
                                                          
            ! Update data output time
            call S3M_Tools_Time_GetNewDate(sTimeNew, sTimeOutput, nint(real(iDtData_Output)))
            sTimeOutput = sTimeNew
        endif
        

        !------------------------------------------------------------------------------------------
        ! Clock step end                          
        call S3M_Tools_Time_Printer(sDateRunStep, sTimeRunStep)
        call mprintf(.true., iINFO_Basic, ' STEP CLOCK END :: Time: ' //trim(sDateRunStep)//' '//trim(sTimeRunStep) )
        
        ! Global time updating
        call mprintf(.true., iINFO_Basic, ' ===== TIME STEP END :: Date, iT, iNTime, iDt :: ' &
                                          //sTimeOld//' '//trim(sTime)//' '// &
                                          trim(sNTime)//' '//trim(sDtModel)// ' ===== ')
                  
        call S3M_Tools_Time_GetNewDate(sTimeNew, sTimeOld, nint(real(iDtModel)))
        sTimeOld = sTimeNew
        !------------------------------------------------------------------------------------------

    enddo TimeLoop
    !------------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------------
    ! Call time end
    call S3M_Tools_Time_Printer(sDateRunEnd, sTimeRunEnd)
    ! Ending simulation message
    call mprintf(.true., iINFO_Basic, ' SIMULATION END :: Time: ' //trim(sDateRunEnd)//' '//trim(sTimeRunEnd) )
    call mprintf(.true., iINFO_Basic, '************************************************************************')
    !------------------------------------------------------------------------------------------

end program S3M_Main
!-----------------------------------------------------------------------------------

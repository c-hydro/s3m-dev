!------------------------------------------------------------------------------------
! File:   S3M_Module_Phys_Snow_Apps_Glaciers.f90
! Author:   Francesco Avanzi, Fabio Delogu.
!
! Created on October 07, 2020 6:00 PM
! Last update on October 27, 2020 10:16 AM
!
! Glacier-physics application module
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
! Module Header
module S3M_Module_Phys_Snow_Apps_Glaciers

    !------------------------------------------------------------------------------------
    ! External module(s) 
    use S3M_Module_Namelist,            only:   oS3M_Namelist
    use S3M_Module_Vars_Loader,         only:   oS3M_Vars
    use gnufor2
   
    use S3M_Module_Tools_Debug
    
    
    ! Implicit none for all subroutines in this module
    implicit none
    !------------------------------------------------------------------------------------------

contains

    !------------------------------------------------------------------------------------------
    ! Delta-H application routine
    ! DOI: https://hess.copernicus.org/articles/14/815/2010/
    subroutine S3M_Phys_Snow_Apps_GlacierDeltaH(iID, iRows, iCols, iRows_Pivot, iCols_Pivot, &
                                            a2dVarDEM, a2dVarIceThick_WE, a2dVarMeltingGCumWY, a2iVarMask)

        !------------------------------------------------------------------------------------------
        ! Variable(s) declaration 
        integer(kind = 4)   :: iID, iRows, iCols, iRows_Pivot, iCols_Pivot, iI, iI_Rows, iI_Cols
        
        real(kind = 4)      :: dMinElevGlacier, dMaxElevGlacier, dMeanElevGlacier, dTemp_MaxdeltaH, dTemp_GlacierArea, dTemp_f_s, &
                               dMaxInitialThicknessGlacier, dTempMatric_MeltingGCumWY_this_glacier, &
                               dTempMatrix_deltaH_by_Area, a2dVarIceThick_WE_change_dueDeltaH
        
        real(kind = 4), dimension(iRows_Pivot - 1)         :: a1dVarDiscStepsDeltaH
        real(kind = 4), dimension(iCols_Pivot - 1)         :: a1dVarGlacierIDinPivotTable
        real(kind = 4), dimension(iRows_Pivot - 1, iCols_Pivot - 1)  :: a2dVarPivotTableValues
        real(kind = 4), dimension(iRows_Pivot, iCols_Pivot)  :: a2dVarPivotTable
        
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarDem, a2dVarAreaCell
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarIceThick_WE, a2dVarMeltingGCumWY, a2dVarGlaciers_ID
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarTempMatrix_DEM, a2dVarTempMatrix_DEMnorm, &
                                                           a2dVarTempMatrix_deltaH
                                                           
                                                           
        integer(kind = 4), dimension(iRows, iCols)         :: a2iVarMask, a2iIndex_for_deltaH
        
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        !Initialization        
        iI = 0; iI_Rows = 0; iI_Cols = 0; dMinElevGlacier = -9999.0; dMaxElevGlacier = -9999.0; dTemp_MaxdeltaH = -9999.0;      
        dTemp_GlacierArea = -9999.0; dMaxInitialThicknessGlacier = -9999.0; dTempMatric_MeltingGCumWY_this_glacier = -9999.0; &
        dTempMatrix_deltaH_by_Area = -9999.0; a2dVarIceThick_WE_change_dueDeltaH = -9999.0; 
        dMeanElevGlacier = -9999.0;
        
        a1dVarDiscStepsDeltaH = -9999.0; a1dVarGlacierIDinPivotTable = -9999.0; a2dVarPivotTableValues = -9999.0;
        a2dVarPivotTable = -9999.0; 
        a2dVarTempMatrix_DEM = -9999.0; a2dVarTempMatrix_DEMnorm = -9999.0; a2dVarTempMatrix_deltaH = -9999;
        a2iIndex_for_deltaH = -9999;

        a2dVarPivotTable = oS3M_Vars(iID)%a2dGlacierPivotTable
        a2dVarGlaciers_ID = oS3M_Vars(iID)%a2iGlaciers_ID
        a2dVarAreaCell = oS3M_Vars(iID)%a2dAreaCell
        !------------------------------------------------------------------------------------------ 
        ! Info start
        call mprintf(.true., iINFO_Extra, ' Phys :: Snow :: Glacier DeltaH ... ' )
        !------------------------------------------------------------------------------------------ 
        
        !------------------------------------------------------------------------------------------        
        !We start by rearranging the pivot table
        a1dVarDiscStepsDeltaH = a2dVarPivotTable(2:iRows_Pivot, 1)
        a1dVarGlacierIDinPivotTable = a2dVarPivotTable(1, 2:iCols_Pivot)
        a2dVarPivotTableValues = a2dVarPivotTable(2:iRows_Pivot, 2:iCols_Pivot)
        
        !Now we perform a loop on a1dVarGlacierIDinPivotTable and for each glacier in that pivot table we: 
        !    (1) compute max, mean, and min elevation and determine normalized elevation for each pixel of that glacier
        !    (2) leveraging the pivot table, we convert normalized elevation to deltaH
        !    (3) we then compute the correction factor f_s and apply that to the deltaH above to get the actual change in thickness
        
        do iI = 1, (iCols_Pivot - 1) !we loop over the number of glaciers in the pivot table
        
            !We move the DEM to a temporary matrix that we will use for computations
            a2dVarTempMatrix_DEM = a2dVarDEM
        
            dMaxInitialThicknessGlacier = MAXVAL(a2dVarIceThick_WE, &
            MASK =( (a2dVarGlaciers_ID .eq. a1dVarGlacierIDinPivotTable(iI)) .and. (a2dVarDEM .ge. 0.0) ) )
            
            if ( (dMaxInitialThicknessGlacier .gt. 0.0) ) then
                !we include this condition to protect this code in case one glacier disappeared during the simulation. 
                !In that case we skip the entire routine since there is no point to update melting...
                !Note that if a glacier changed area during a simulation, a2dVarGlaciers_ID would not be updated; this layer
                !should thus be interpreted as INITIAL extent. 
                
                !We compute min, max and mean elevations
                dMinElevGlacier = MINVAL(a2dVarTempMatrix_DEM, &
                MASK =( (a2dVarGlaciers_ID .eq. a1dVarGlacierIDinPivotTable(iI)) & 
                .and. (a2dVarDEM .ge. 0.0) .and. (a2dVarIceThick_WE .gt. 0.0) ) )
                
                dMaxElevGlacier = MAXVAL(a2dVarTempMatrix_DEM, &
                MASK =( (a2dVarGlaciers_ID .eq. a1dVarGlacierIDinPivotTable(iI)) & 
                .and. (a2dVarDEM .ge. 0.0) .and. (a2dVarIceThick_WE .gt. 0.0) ) )
                
                dMeanElevGlacier = SUM(a2dVarTempMatrix_DEM, &
                MASK =( (a2dVarGlaciers_ID .eq. a1dVarGlacierIDinPivotTable(iI)) & 
                .and. (a2dVarDEM .ge. 0.0) .and. (a2dVarIceThick_WE .gt. 0.0) ) )/&
                COUNT((a2dVarGlaciers_ID .eq. a1dVarGlacierIDinPivotTable(iI)) & 
                .and. (a2dVarDEM .ge. 0.0) .and. (a2dVarIceThick_WE .gt. 0.0) )
                
                if (dMaxElevGlacier .gt. dMinElevGlacier) then    
                    !that is, if the glacier occupies more than one pixel, do deltaH; otherwise, it's time for non-dynamic MB
                    
                    !We compute normalized glacier elevation for pixels belonging to that glacier ID -- others to -9999
                    where( (a2dVarGlaciers_ID .eq. a1dVarGlacierIDinPivotTable(iI)) .and. (a2dVarDEM .ge. 0.0) &
                        .and. (a2dVarIceThick_WE .gt. 0.0) )

                        a2dVarTempMatrix_DEMnorm = (dMaxElevGlacier - a2dVarTempMatrix_DEM)/(dMaxElevGlacier - dMinElevGlacier)

                    elsewhere

                        a2dVarTempMatrix_DEMnorm = -9999.0

                    endwhere

                    !Now we compute deltaH for each non negative pixel
                    do iI_Rows = 1, iRows

                        do iI_Cols = 1, iCols

                            if ( (a2dVarTempMatrix_DEMnorm(iI_Rows, iI_Cols) .ge. 0.0) ) then

                                a2iIndex_for_deltaH(iI_Rows, iI_Cols) = &
                                            MINLOC(ABS(a1dVarDiscStepsDeltaH - a2dVarTempMatrix_DEMnorm(iI_Rows, iI_Cols)), DIM = 1)
                                 a2dVarTempMatrix_deltaH(iI_Rows, iI_Cols) = & 
                                                a2dVarPivotTableValues(a2iIndex_for_deltaH(iI_Rows, iI_Cols), iI)   

                            endif

                        enddo

                    enddo

                    !Now we want to compute f_s, which is a scalar correction factor. 
                    !Note that some glaciers might have -1 in the pivot
                    !table, meaning we DO NOT WANT to apply the deltaH approach for those glaciers. 
                    !So here we handle this possibility...
                    dTemp_MaxdeltaH = MAXVAL(a2dVarTempMatrix_deltaH, &
                    MASK =( (a2dVarGlaciers_ID .eq. a1dVarGlacierIDinPivotTable(iI)) & 
                    .and. (a2dVarDEM .ge. 0.0) .and. (a2dVarIceThick_WE .gt. 0.0) ) )

                    if (dTemp_MaxdeltaH .lt. 0.0) then

                        !then this is a glacier for which we want non-dynamic mass balance and we simply apply that
                        where( (a2dVarGlaciers_ID .eq. a1dVarGlacierIDinPivotTable(iI)) .and. (a2dVarDEM .ge. 0.0) &
                            .and. (a2dVarIceThick_WE .gt. 0.0) )
                            
                            a2dVarIceThick_WE = a2dVarIceThick_WE - a2dVarMeltingGCumWY

                        endwhere

                    else     

                        !then this is a glacier for which we want to apply the deltaH. So we compute f_s and apply it 
                        !to get the new thickness
                        !NOTE: the implementation here followed Seibert et al. 2018
                        !https://www.zora.uzh.ch/id/eprint/161601/1/2018_hess-22-2211-2018.pdf
                        !because this paper expresses all equations in mm w.e. Apart from this, the essence of this implementation
                        !is the same as the original Huss paper. 
                        dTemp_GlacierArea = SUM(a2dVarAreaCell, &
                                MASK =( (a2dVarGlaciers_ID .eq. a1dVarGlacierIDinPivotTable(iI)) &
                                .and. (a2dVarDEM .ge. 0.0) .and. (a2dVarIceThick_WE .gt. 0.0) ) )

                        dTempMatric_MeltingGCumWY_this_glacier = &
                                SUM(a2dVarMeltingGCumWY*(a2dVarAreaCell/dTemp_GlacierArea), &
                                MASK =( (a2dVarGlaciers_ID .eq. a1dVarGlacierIDinPivotTable(iI)) &
                                .and. (a2dVarDEM .ge. 0.0) .and. (a2dVarIceThick_WE .gt. 0.0) ) )

                        dTempMatrix_deltaH_by_Area  = &
                                SUM((a2dVarAreaCell/dTemp_GlacierArea)*a2dVarTempMatrix_deltaH, &
                                MASK =( (a2dVarGlaciers_ID .eq. a1dVarGlacierIDinPivotTable(iI)) &
                                .and. (a2dVarDEM .ge. 0.0) .and. (a2dVarIceThick_WE .gt. 0.0) ) )

                        dTemp_f_s = dTempMatric_MeltingGCumWY_this_glacier/(dTempMatrix_deltaH_by_Area)
                        
                        !Now: according to Huss et al. (2010) and our own experiments, the deltaH approach might return 
                        !"unrealistically high surface lowering on the glacier tongue", especially for 
                        !"years with strongly negative mass balance". There is likely no definitive approach to date re: how to 
                        !consistently handle this issue, so we consider a pragmatic strategy here and force 
                        !dTemp_f_s*a2dVarTempMatrix_deltaH (expected surface lowering) to be .eq. to a2dVarMeltingGCumWY
                        !(local melt) wherever the former is larger than the latter (in absolute terms) AND elevation is 
                        !lower than dMeanElevGlacier (so we are somewhere around the tongue)
                        
                        where( (a2dVarGlaciers_ID .eq. a1dVarGlacierIDinPivotTable(iI)) .and. (a2dVarDEM .ge. 0.0) &
                            .and. (a2dVarIceThick_WE .gt. 0.0) .and. (a2dVarDEM .lt. dMeanElevGlacier) .and. &
                            (dTemp_f_s*a2dVarTempMatrix_deltaH .gt. a2dVarMeltingGCumWY))                        
                        
                        a2dVarIceThick_WE = a2dVarIceThick_WE - a2dVarMeltingGCumWY
                        
                        elsewhere( (a2dVarGlaciers_ID .eq. a1dVarGlacierIDinPivotTable(iI)) .and. (a2dVarDEM .ge. 0.0) &
                            .and. (a2dVarIceThick_WE .gt. 0.0) )
                            
                            a2dVarIceThick_WE = a2dVarIceThick_WE - dTemp_f_s*a2dVarTempMatrix_deltaH
                            
                        endwhere
                        
                    endif
                    
                else !if the glacier occupies only one pixel, then non-dynamic MB
                    
                    where( (a2dVarGlaciers_ID .eq. a1dVarGlacierIDinPivotTable(iI)) &
                        .and. (a2dVarDEM .ge. 0.0) .and. (a2dVarIceThick_WE .gt. 0.0) )

                        a2dVarIceThick_WE = a2dVarIceThick_WE - a2dVarMeltingGCumWY

                    endwhere                    
                
                endif !this is the end of the if clause on dMaxElevGlacier .gt. dMinElevGlacier  
                
            endif !this is the end of the if clause on (dMaxInitialThicknessGlacier .gt. 0.0) 
            
        enddo 
            
            
        !There might be some pixels with positive thickness, but that are not part of any glacier for which we provided a deltaH
        !parametrization. For this pixels, we simply apply the local MB.
        where( (a2dVarGlaciers_ID .le. 0.0) .and. (a2dVarIceThick_WE .gt. 0.0) .and. (a2dVarDEM .ge. 0.0) )

            a2dVarIceThick_WE = a2dVarIceThick_WE - a2dVarMeltingGCumWY

        endwhere                            
    
        !Finally, we set to 0 any negative thickness: the glacier at that location has gone
        where( (a2dVarIceThick_WE .lt. 0.0) .and. (a2dVarDEM .ge. 0.0) )        
            
            a2dVarIceThick_WE = 0.0
            
        endwhere
        !------------------------------------------------------------------------------------------
            
        !------------------------------------------------------------------------------------------
        ! Info end
        call mprintf(.true., iINFO_Extra, ' Phys :: Snow :: Glacier DeltaH... OK' )
        !------------------------------------------------------------------------------------------
        
    end subroutine S3M_Phys_Snow_Apps_GlacierDeltaH
    !------------------------------------------------------------------------------------------   


end module S3M_Module_Phys_Snow_Apps_Glaciers
!------------------------------------------------------------------------------------
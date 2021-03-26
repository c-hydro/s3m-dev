!------------------------------------------------------------------------------------
! File:   S3M_Module_Phys_Snow_Apps_Assimilation.f90
! Author:   Simone Gabellani, Fabio Delogu, Daniele Dolia, Francesco Avanzi.
!
! Created on June 10, 2020 5:00 PM
! Last update on October 26, 2020 5:00 PM
!
! Assimilation application module
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
! Module Header
module S3M_Module_Phys_Snow_Apps_Assimilation

    !------------------------------------------------------------------------------------
    ! External module(s) 
    use S3M_Module_Namelist,            only:   oS3M_Namelist
    use S3M_Module_Vars_Loader,         only:   oS3M_Vars
   
    use S3M_Module_Tools_Debug
    
    use S3M_Module_Tools_Generic,       only:   assimNudging
    
    ! Implicit none for all subroutines in this module
    implicit none
    !------------------------------------------------------------------------------------------

contains

    !------------------------------------------------------------------------------------------
    ! Subroutine to compute SWE assimilation from snow-height interpolated observations
    subroutine S3M_Phys_Snow_Apps_AssimSH(iID, iRows, iCols, &
                                            dVarRhoW, &
                                            a2iVarMask, &
                                            a2dVarSnowHeight, a2dVarSnowKernel, &
                                            a2dVarSnowCA, a2dVarSnowQA, &
                                            a2dVarSWE, a2iVarAge, a2dVarAlbedoS, a2dVarRhoS, iFlagAssOnlyPos)
        
        !------------------------------------------------------------------------------------------
        ! Variable(s) declaration 
        integer(kind = 4)   :: iID, iRows, iCols, iFlagAssOnlyPos
        
        real(kind = 4)      :: dVarRhoW, dVarSnowQThr
        
        integer(kind = 4), dimension(iRows, iCols)      :: a2iVarMask 
        
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarSnowHeight, a2dVarSnowKernel
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarSnowCA, a2dVarSnowQA
        
        integer(kind = 4), dimension(iRows, iCols)      :: a2iVarAge 
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarSWE, a2dVarSWE_Assim
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarAlbedoS, a2dVarRhoS
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarCorr
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Snow Cover Area (SCA) value(s)
        ! NoSnow = 0; Snow = 1; Cloud = 2; NoData = -1, NoDecision = 3
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Get information
        dVarSnowQThr = oS3M_Namelist(iID)%dSnowQualityThr 
        
        ! Info start
        call mprintf(.true., iINFO_Extra, ' Phys :: Snow :: SWE Assim from snow height and satellite ... ' )
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Snow Height condition(s)
        
        ! Check snow quality value(s) using a threshold
        where( (a2dVarSnowQA.le.dVarSnowQThr ) ) a2dVarSnowHeight = -9999.0
                              
        ! Terrain condition --> use SCA value instead of snow height 
        where( (a2dVarSnowCA.eq.0.0) .and. (a2dVarSnowQA.gt.dVarSnowQThr) .and. (a2dVarSnowQA.lt.0.99) )  a2dVarSnowHeight = 0.0                    
              
        ! Cloud condition
        where( (a2dVarSWE.LT.20.0) .and. (a2dVarSnowCA.eq.2.0) .and. (a2dVarSnowQA.gt.dVarSnowQThr) ) a2dVarSnowHeight = -9999.0    
                                       
        ! NoDecision condition
        where( (a2dVarSWE.LT.20.0) .and. (a2dVarSnowCA.eq.3.0) .and. (a2dVarSnowQA.gt.dVarSnowQThr) ) a2dVarSnowHeight = -9999.0 
                             
        ! NoData condition
        where( (a2dVarSWE.LT.20.0) .and. (a2dVarSnowCA.eq.-1.0) .and. (a2dVarSnowQA.gt.dVarSnowQThr) ) a2dVarSnowHeight = -9999.0  
        !-----------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Snow Kernel condition(s)
            
        ! Terrain condition
        where( (a2dVarSnowCA.eq.0.0) .and. (a2dVarSnowQA.gt.dVarSnowQThr*2.5) .and. (a2dVarSnowQA.lt.0.99) ) 
            a2dVarSnowKernel = a2dVarSnowKernel*2.0
        endwhere
        ! Check upper boundary
        where( a2dVarSnowKernel.gt.1) a2dVarSnowKernel = 0.9;
        !------------------------------------------------------------------------------------------   

        !------------------------------------------------------------------------------------------
        ! Compute observed SWE starting from snow height interpolated observations
        where(a2dVarSnowHeight.gt.0.0) a2dVarSnowHeight = a2dVarSnowHeight*a2dVarRhoS/dVarRhoW
        !------------------------------------------------------------------------------------------
                      
        !------------------------------------------------------------------------------------------
        ! Nudging assimilation method
        call assimNudging(a2iVarMask, a2dVarSWE, a2dVarSnowHeight, a2dVarSnowKernel, &
                          a2dVarSWE_Assim, a2dVarCorr, iFlagAssOnlyPos) 
        !------------------------------------------------------------------------------------------
                            
        !------------------------------------------------------------------------------------------  
        ! Update SWE with nudging assimilation
        a2dVarSWE = a2dVarSWE_Assim              
        !------------------------------------------------------------------------------------------  
                          
        !------------------------------------------------------------------------------------------
        ! Info end
        call mprintf(.true., iINFO_Extra, ' Phys :: Snow :: SWE Assim from snow height and satellite ... OK' )
        !------------------------------------------------------------------------------------------
        
    end subroutine S3M_Phys_Snow_Apps_AssimSH
    !------------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------------
    ! Subroutine to compute SWE assimilation map
    subroutine S3M_Phys_Snow_Apps_AssimSWE(iID, iRows, iCols, &
                                            a2iVarMask, &
                                            a2dVarSWE, a2dVarSWEass, iFlagAssOnlyPos)

        !------------------------------------------------------------------------------------------
        ! Variable(s) declaration 
        integer(kind = 4)   :: iID, iRows, iCols
        integer(kind = 4)   :: iVarDayNud, iVarSWEassInfluence, iFlagAssOnlyPos

        real(kind = 4)      :: dSigma, dVarWeightSWEass

        integer(kind = 4), dimension(iRows, iCols)      :: a2iVarMask 

        real(kind = 4), dimension(iRows, iCols)         :: a2dVarSWEass
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarSWE
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarW
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarCorr, a2dVarSWE_Assim
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Get information
        iVarSWEassInfluence = oS3M_Namelist(iID)%iSWEassInfluence
        dVarWeightSWEass = oS3M_Namelist(iID)%dWeightSWEass
        iVarDayNud = oS3M_Vars(iID)%iDayNud
                
        ! Gaussian variance to determine dWt
        dSigma = iVarSWEassInfluence/1.5
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Initialize variable(s)
        a2dVarW = -9999.0;

        ! Info start
        call mprintf(.true., iINFO_Extra, ' Phys :: Snow :: SWE Assim map... ' )
        
        ! Debug start
        call mprintf(.true., iINFO_Extra, checkvar(a2dVarSWE, oS3M_Vars(iID)%a2iMask, 'SWE START') ) 
        call mprintf(.true., iINFO_Extra, checkvar(a2dVarSWEass, oS3M_Vars(iID)%a2iMask, 'SWE ASSIM START') ) 
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Compute SWE weight matrix
        a2dVarW = EXP(-(iVarDayNud*1.0-iVarSWEassInfluence)**2/0.5/dSigma**2)*dVarWeightSWEass
        
        ! Nudging assimilation method
        call assimNudging(a2iVarMask, a2dVarSWE, a2dVarSWEass, a2dVarW, a2dVarSWE_Assim, a2dVarCorr, iFlagAssOnlyPos)
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------  
        ! Update SWE with SWE with nudging assimilation
        a2dVarSWE = a2dVarSWE_Assim              
        !------------------------------------------------------------------------------------------  
                          
        !------------------------------------------------------------------------------------------
        ! Info end
        call mprintf(.true., iINFO_Extra, ' Phys :: Snow :: SWE Assim map... OK' )
        !------------------------------------------------------------------------------------------
        
    end subroutine S3M_Phys_Snow_Apps_AssimSWE
    !------------------------------------------------------------------------------------------

end module S3M_Module_Phys_Snow_Apps_Assimilation
!------------------------------------------------------------------------------------

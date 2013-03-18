!> Computes CMB and Deflection field power spectra from CAMB.
!! Also computes instrumental noise, 
!! lensing reconstruction noise spectra using LensingNoise
!! and derivatives used in Fisher calculation
!!@author Sudeep Das, Princeton University 
!!@version 1, 2008

module GetSpectra
  use Fisher
  use CAMB 
  use LambdaGeneral
  use miscUtilsForCAMB
  use LensingNoise
  implicit none
  
 
  

contains

  

  !> Gets CMB and deflection power spectra from CAMB 
  subroutine getCosmoPS(params,lmax,doLensing,allcls)
    real(dl)                 :: params(NPARAMSMAX)
    integer                  :: lmax 
    real(dl)                 :: clfac(2:lmax),clfacPhi(2:lmax),clfacPhiTemp(2:lmax)
    type(CAMBAndNoiseCLS)    :: allcls
    type(CAMBParams)         :: P
    real(dl)                 :: output_factor, AGE
    integer                  :: i 
    logical                  :: doLensing
    
    allcls%cl(2:lmax,:) = 0._dl
    allcls%clLensed(2:lmax,:) = 0._dl
    call CAMB_SetDefParams(P)

    P%WantScalars = .true.
    P%WantVectors = .false.
    P%WantTensors = .false. 

    P%OutputNormalization=outNone
    P%WantCls= P%WantScalars .or. P%WantTensors .or. P%WantVectors
    P%WantTransfer=.false.
    P%DoLensing = doLensing
    P%Max_l = lmax + 1000
    P%Max_eta_k = 2._dl*P%Max_l
    
    output_factor = 1.0_dl
    
    w_lam = params(I_w)
    cs2_lam = 1.0
    AccuracyBoost=1
    LAccuracyBoost=2
    P%AccurateBB = .true.
    P%NonLinear=NonLinear_Lens
    call solveThetaAandSetParams(params,P)
    
    P%tcmb = 2.726_dl
    P%yhe = params(I_HELIUM_FRACTION)
    
    P%Scalar_initial_condition = initial_adiabatic
    
    P%reion%Reionization = .true.
    
    P%Reion%use_optical_depth = .true.
    P%Reion%optical_depth = params(I_TAU)
    P%Reion%delta_redshift =  0.5_dl!for Oliver
    
    P%InitPower%nn  = 1
    P%InitPower%an(1) = params(I_N_S)
    P%InitPower%n_run(1) = params(I_N_RUN)
    P%InitPower%ScalarPowerAmp(1) = params(I_A_S)/1.d9
    !P%InitPower%k_0_scalar= 0.002_dl
    P%InitPower%k_0_scalar= 0.05_dl !Changing to compare with Oliver Zahn
    !print*, P%InitPower%ScalarPowerAmp(1), params(I_A_S)
    P%DoLensing = doLensing
    
    if (.not. CAMB_ValidateParams(P)) stop 'Stopped due to parameter error'
    
    FeedbackLevel = verbosity
    
    if (FeedbackLevel > 1) then
       Age = CAMB_GetAge(P)
       write (*,'("Age of universe/GYr  = ",f7.3)') Age
    end if
    
    P%AccuratePolarization = .true.
    P%AccurateReionization = .true.    
   
    call CAMB_GetResults(P)
    if (debugging) print*, "thetaA Final" , cosmomctheta(),params(I_THETA_A)
    if (debugging) print*, "Amplitude" ,CP%InitPower%ScalarPowerAmp(1) 
    if (debugging) print*, "Amplitude", CP%Reion%optical_depth
    if (debugging) print*, "ombh2,omch2,omnuh2", CP%Omegab*(P%H0/100._dl)**2,CP%Omegac*(P%H0/100._dl)**2,CP%Omegan*(P%H0/100._dl)**2
    if (debugging) print*, "omLambda, H0 ",CP%Omegav,CP%H0
    if (debugging) write(117,*) CP
    do i = 2,lmax
       clfac(i)=(real(i,dl)*(real(i,dl)+1._dl))/(2._dl*Pi)
       clfacPhi(i) = (real(i,dl)**4)/(real(i,dl)*(real(i,dl)+1._dl))
       clfacPhiTemp(i) = (real(i,dl)**3)/dsqrt(real(i,dl)*(real(i,dl)+1.))
    end do
    
    allcls%cl(2:lmax, INDEX_TT) = Cl_scalar(2:lmax,1,1)/clfac(2:lmax)
    allcls%cl(2:lmax, INDEX_EE) = Cl_scalar(2:lmax,1,2)/clfac(2:lmax)
    allcls%cl(2:lmax, INDEX_TE) = Cl_scalar(2:lmax,1,3)/clfac(2:lmax)
    !allcls%cl(2:lmax, INDEX_BB) = Cl_tensor(2:lmax,1,1)/clfac(2:lmax)
    !this BB is useless for now, but with tensor modes may become useful later 
    If (doLensing) then
       allcls%cl(2:lmax, INDEX_TD) = real(Cl_scalar(2:lmax,1,C_PhiTemp),dl)/clfacPhiTemp(2:lmax)
       allcls%cl(2:lmax, INDEX_DD) = real(Cl_scalar(2:lmax,1,C_Phi),dl)/clfacPhi(2:lmax)
       allcls%cllensed(2:lmax, INDEX_TT) = Cl_lensed(2:lmax,1,CT_Temp)/clfac(2:lmax)
       allcls%cllensed(2:lmax, INDEX_EE) = Cl_lensed(2:lmax,1,CT_E)/clfac(2:lmax)
       allcls%cllensed(2:lmax, INDEX_TE) = Cl_lensed(2:lmax,1,CT_Cross)/clfac(2:lmax)
       allcls%cllensed(2:lmax, INDEX_BB) = Cl_lensed(2:lmax,1,CT_B)/clfac(2:lmax)
    end If
        
  end subroutine getCosmoPS
  
  !> Solves for \f$\Omega_\Lambda\f$ given the size of the sound horizon
  !! and other parameters and sets up the calculations for CAMB
  subroutine solveThetaAandSetParams(params,P)
    real(dl)                    :: params(NPARAMSMAX) 
    type(CAMBParams)            :: P
    real(dl)                    :: omvleft, omvright, omv
    real(dl)                    :: thetaA, diff, diff0, nu_frac
    
    real(dl), parameter         :: tolerance = 1.d-7
    integer                     :: i 
    integer, parameter          :: MaxIter = 50
    P%Num_Nu_massless = params(I_NU_EFF)
    P%Num_Nu_massive  = params(I_NU_MASSIVE)
    if (P%Num_nu_massive .gt. 0) then 
       P%MassiveNuMethod = 3
       P%Nu_mass_eigenstates = 1
       P%Nu_mass_degeneracies(1) = P%num_nu_massive
       P%Nu_mass_splittings = .true.
       P%Nu_mass_fractions(1)=1  
    end if
    
    omvleft = 0.5
    omvright = 0.85
    do i = 1, MaxIter
       
       if (i .eq. 1) then
          omv = omvright 
       else 
          omv = (omvleft+omvright)/2._dl
       end if
       
       P%omegav = omv
       P%H0 = 100._dl*sqrt((params(I_OMEGA_B_H2) + params(I_OMEGA_C_H2))&
            &/(1-P%omegav -params(I_OMEGA_K)))
       P%omegab = params(I_OMEGA_B_H2)/(P%H0/100._dl)**2
       nu_frac = params(I_NU_FRAC)
       P%omegac = (1._dl-Nu_frac)*params(I_OMEGA_C_H2)/(P%H0/100._dl)**2
       P%omegan = NU_frac*params(I_OMEGA_C_H2)/(P%H0/100._dl)**2 
       FeedBackLevel=0
       call CAMBParams_Set(P) 
       if (.not. CAMB_ValidateParams(P)) stop 'Stopped due to parameter error'
       thetaA = cosmomctheta()
!       if (debugging) print*, "thetaA" , thetaA,params(I_THETA_A)
!       if(debugging) print*, "CP%r, omegak", CP%r,P%omegak, params(I_Omega_k),&
!            &(1-P%omegab-P%omegac-P%omegan-P%omegav),CP%omegak
       if (i .eq. 1) then 
          diff0 =  (thetaA - params(I_THETA_A))/params(I_THETA_A)
       else
          diff =  (thetaA - params(I_THETA_A))/params(I_THETA_A)
          !print*, "thetaA" , thetaA,params(I_THETA_A)
          !print*, "diff :", diff,i
          if (sign(1._dl,diff) .ne. sign(1._dl,diff0)) then 
             omvleft=omv   
          else
             omvright=omv 
          end if
          if(abs(diff) .lt. tolerance) exit
       end if
    end do
    if((i-1) .eq. MaxIter ) STOP 'Maximum Iterations Exceeded in: solveThetaAandSetParams'
       
  end subroutine solveThetaAandSetParams

  !> Gets the instrumental and Lensing noise power spetra
  subroutine GetNoisePS(fp,allcls)
    type(FisherParams)          :: fp
    type(CAMBAndNoiseCls)       :: allcls
    real(dl) :: wBlsq_T(2:fp%lmax,1:fp%Pexp%nchannels), wBlsq_T_comb(2:fp%lmax)
    real(dl) :: wBlsq_P(2:fp%lmax,1:fp%Pexp%nchannels), wblsq_P_comb(2:fp%lmax)
    integer  :: ichannel,il, polcombi
    real(dl) :: arcmin2radian,thetaFWHM(fp%Pexp%nchannels)
    ! w_c = 1/((\Delta T/T)_c * Theta_FWHM_c)^2 !c refers to channel
    
    ! B^2_{l,c} = exp(-l(l+1)theta_FWHM_c/(8ln2))
    
    ! w B^2_l = \sum w_c B^2_{l,c}
    arcmin2radian = 1._dl/60._dl*PI/180._dl

    thetaFWHM = fp%Pexp%thetaFWHM*arcmin2radian

    do ichannel = 1,fp%Pexp%nchannels
       do il = 2,fp%lmax
          wBlsq_T(il,ichannel) = 1._dl/(fp%Pexp%DeltaTemp(ichannel)*thetaFWHM(ichannel))**2*&
               & exp(-(real(il,dl)*(real(il,dl)+1)*thetaFWHM(ichannel)**2/(8._dl*log(2._dl))))
          wBlsq_P(il,ichannel) = 1._dl/(fp%Pexp%DeltaPol(ichannel)*thetaFWHM(ichannel))**2*&
               & exp(-(real(il,dl)*(real(il,dl)+1)*thetaFWHM(ichannel)**2/(8._dl*log(2._dl))))
       end do
    end do

    wBlsq_T_comb(2:fp%lmax) = 0.0_dl
    wBlsq_P_comb(2:fp%lmax) = 0.0_dl
    
    do ichannel = 1,fp%Pexp%nchannels
       wBlsq_T_comb(2:fp%lmax) = wBlsq_T_comb(2:fp%lmax) +  WBlsq_T(2:fp%lmax,ichannel)
       wBlsq_P_comb(2:fp%lmax) = wBlsq_P_comb(2:fp%lmax) + WBlsq_P(2:fp%lmax,ichannel)
    end do
    
    allcls%nl = 0.0_dl
    allcls%nl(2:fp%lmax,INDEX_TT) = 1./(wBlsq_T_comb(2:fp%lmax))
    allcls%nl(2:fp%lmax,INDEX_EE) = 1./(wBlsq_P_comb(2:fp%lmax))
    allcls%nl(2:fp%lmax,INDEX_BB) = allcls%nl(2:fp%lmax,INDEX_EE)
    
    !Now get the lensing noise
    if (fp%IncludeLensing) then 
       polcombi = 6
       if(fp%nspectra .eq. 3) polcombi = 1
       call GetLensingNoise(fp%lmax,polcombi,allcls)
    end if
    
  end subroutine GetNoisePS

  !> Allocates memeory for PowerSpectraAndDerivs
  subroutine AllocPowerSpectraAndDerivs(fp,PS)
    type(FisherParams) :: fp 
    type(PowerSpectraAndDerivs) :: PS 
    
    call assert(fp%isSet,'FisherParams has not be set yet') 
    
    allocate(PS%clFiducial(2:fp%lmax,1:fp%nspectra))
    allocate(PS%clDeriv(2:fp%lmax,1:fp%nspectra,1:fp%nvary))
    allocate(PS%clNoise(2:fp%lmax,1:fp%nspectra))
    PS%IsAllocated = .true. 
  end subroutine AllocPowerSpectraAndDerivs
  
  !> Gets the Derivatives and Noise power spectra 
  subroutine GetDerivsAndNoise(fp, PS)
    type(FisherParams) :: fp
    type(PowerSpectraAndDerivs) :: PS
    real(dl)                    :: params(NPARAMSMAX)
    integer                     :: i,kvary, verb, iL
    type(CAMBAndNoiseCls)       :: clfid,cltemp0,cltemp1,clCurvPlus,clCurvMinus
    real(dl), parameter         :: curvatureFloor = 1.d-6
    
    call AllocPowerSpectraAndDerivs(fp,PS)
    !Get the fiducial power spectra
    If (verbosity > 0) write(*,'("Getting the fiducial power spectra")')
    
    call GetCosmoPS(fp%fiducial,fp%lmax,fp%includeLensing,clfid)
    PS%clFiducial(2:fp%lmax,1:fp%nspectra) = clfid%cl(2:fp%lmax,fp%SpectrumIndex(1:fp%nspectra))
    
    verb = verbosity
    If (verbosity > 0) write(*,'("Getting the power spectra derivatives")')

    do i = 1, fp%nvary
       verbosity = 0
       !Get the left side 
      
       params = fp%fiducial
       
       kvary = fp%VaryingParamIndex(i)
       params(kvary) = params(kvary) + fp%stepsize(kvary)
       write(*,'("Getting the power spectra derivative w.r.t ",a)')paramNames(kVary)
       !print*, params(kVary) 
       debugging = .false.
       
       !if(ParamNames(kVary) .eq. 'w') debugging = .true. 
       
       call getCosmoPS(params,fp%lmax,fp%includeLensing,cltemp1) 
       
       !Get the right side 
       
       params = fp%fiducial
       params(kvary) = params(kvary) - fp%stepsize(kvary)
       !print *, params(kVary) 
       
       call getCosmoPS(params,fp%lmax,fp%includeLensing,cltemp0)
       
       !Get the derivative but treat Omega_k specially

       if((kVary .eq. I_OMEGA_K).and.&
            &(sign(1._dl,fp%fiducial(kVary)+fp%stepsize(kvary)) &
            &.ne. sign(1._dl, fp%fiducial(kVary)-fp%stepsize(kvary))))&
            & then 
          if (verbosity>0) print*,'Doing the fancy curvature derivative'
          params = fp%fiducial
          params(kVary) = -curvatureFloor
          call getCosmoPS(params,fp%lmax,fp%includeLensing,clCurvMinus)
          params = fp%fiducial
          params(kVary) = curvatureFloor
          call getCosmoPS(params,fp%lmax,fp%includeLensing,clCurvPlus)
          PS%ClDeriv(2:fp%lmax,1:fp%nspectra,i) = &
               & (cltemp1%cl(2:fp%lmax,fp%SpectrumIndex(1:fp%nspectra))&
               & -clCurvPlus%cl(2:fp%lmax,fp%SpectrumIndex(1:fp%nspectra))&
               & - cltemp0%cl(2:fp%lmax,fp%SpectrumIndex(1:fp%nspectra))&
               & + clCurvMinus%cl(2:fp%lmax,fp%SpectrumIndex(1:fp%nspectra)))&
               &/(2._dl*(fp%stepsize(kvary)-curvatureFloor))
          
       else
          !print*, 'Normal Deriv'
          PS%ClDeriv(2:fp%lmax,1:fp%nspectra,i) = &
               & (cltemp1%cl(2:fp%lmax,fp%SpectrumIndex(1:fp%nspectra))&
               & - cltemp0%cl(2:fp%lmax,fp%SpectrumIndex(1:fp%nspectra)))&
               &/(2._dl*fp%stepsize(kvary))
          if (kvary .eq. I_A_S) PS%ClDeriv(:,:,i) = PS%ClDeriv(:,:,i)*1.d9
          
          if(debugging) then 
             open(unit=101,file=trim(outdir)//'/ClTTleftRight_'//trim(paramNames(kVary))//'.'//TAG)
             do iL = 2,fp%lmax
                write(101,'(I10,x,e20.8,x,e20.8)') iL,cltemp0%cl(iL,INDEX_TT),cltemp1%cl(iL,INDEX_TT)
             end do
             close(101)
          end if
          
       end if
    end do
    verbosity = verb
    
    call GetNoisePS(fp,clfid)
    
    PS%clNoise(2:fp%lmax,1:fp%nspectra) = clfid%nl(2:fp%lmax,fp%SpectrumIndex(1:fp%nspectra))
    
    
  end subroutine GetDerivsAndNoise

 

end module GetSpectra


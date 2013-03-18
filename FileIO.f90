!> Module for input and output.
!! @author Sudeep Das, Princeton University
!! @version 1, 2008
module FileIO
  use Fisher
  use IniFile 
  use miscUtilsForCAMB 
  implicit none
  contains
    !> Parses the parameter file and sets up the Fisher calculation
    subroutine SetFisherParamsFromFile(InputFile,fp)
      character(len=ini_max_string_len)     :: InputFile  !< Input file name
      type(FisherParams)                    :: fp !< the FisherParams type that gets populated
      logical                               :: bad
      character                             :: junk
      integer                               :: i, j 
      character(len = ini_max_string_len)   :: freqs,thetas,deltaTs,deltaPs
      
      call Ini_Open(InputFile, 1, bad, .false.)
      if (bad) stop 'Error opening parameter file'
      Ini_fail_on_not_found = .false.
      OUTDIR = Ini_Read_String('out_dir')
      fp%lMax = ini_read_int('max_l')
      fp%fiducial(I_OMEGA_B_H2) = Ini_Read_Double('Ombh2')
      fp%fiducial(I_OMEGA_C_H2) = Ini_Read_Double('Omdh2')
      fp%fiducial(I_NU_FRAC) = Ini_Read_Double('Nu_frac')
      fp%fiducial(I_THETA_A) = Ini_Read_Double('ThetaA')
      fp%fiducial(I_OMEGA_K) = Ini_Read_Double('Omk')
      fp%fiducial(I_TAU) = Ini_Read_Double('tau')
      fp%fiducial(I_N_S) = Ini_Read_Double('ns')
      fp%fiducial(I_N_RUN) = Ini_Read_Double('nrun')
      fp%fiducial(I_A_S) = Ini_Read_Double('scalar_amplitude')
      fp%fiducial(I_W) = Ini_Read_Double('wlam')
      fp%fiducial(I_HELIUM_FRACTION) = Ini_Read_Double('helium_fraction')
      fp%fiducial(I_NU_EFF) = Ini_Read_Double('num_massless_neutrinos',3.04_dl)
      fp%fiducial(I_NU_MASSIVE) = Ini_Read_Double('num_massive_neutrinos',0._dl)
      
      
      fp%stepsize(I_OMEGA_B_H2) = Ini_Read_Double('Ombh2_step',0._dl)
      fp%stepsize(I_OMEGA_C_H2) = Ini_Read_Double('Omdh2_step',0._dl)
      fp%stepsize(I_NU_FRAC) = Ini_Read_Double('Nu_frac_step',0._dl)
      fp%stepsize(I_THETA_A) = Ini_Read_Double('ThetaA_step',0._dl)
      fp%stepsize(I_OMEGA_K) = Ini_Read_Double('Omk_step',0._dl)
      fp%stepsize(I_TAU) = Ini_Read_Double('tau_step',0._dl)
      fp%stepsize(I_N_S) = Ini_Read_Double('ns_step',0._dl)
      fp%stepsize(I_N_RUN) = Ini_Read_Double('nrun_step',0._dl)
      fp%stepsize(I_A_S) = Ini_Read_Double('scalar_amplitude_step',0._dl)
      fp%stepsize(I_W) = Ini_Read_Double('wlam_step',0._dl)
      fp%stepsize(I_HELIUM_FRACTION) = Ini_Read_Double('helium_fraction_step',0._dl)      
      fp%stepsize(I_NU_EFF) = Ini_Read_Double('num_massless_neutrinos_step',0._dl)
      fp%stepsize(I_NU_MASSIVE) = Ini_Read_Double('num_massive_neutrinos_step',0._dl)
      
      fp%includeLensing  = Ini_Read_Logical('include_lensing',.false.)
      

      fp%nspectra = Ini_Read_Int('n_spectra')

      allocate(fp%SpectrumIndex(fp%nspectra))

      If (not( fp%includeLensing )) then 
         call assert(fp%nspectra .eq. 1 .or. fp%nspectra .eq. 3 , 'Number of spectra has &
              & to be  1 or 3')
         If (fp%nspectra == 1) fp%SpectrumIndex(1) = INDEX_TT
         If (fp%nspectra == 3) fp%SpectrumIndex(:) = (/INDEX_TT,INDEX_TE,INDEX_EE/)
      else
         
         call assert(fp%nspectra .eq. 3 .or. fp%nspectra .eq. 5 , 'Number of spectra has &
              & to be  3 or 5 as lensing is Turned On')
         If (fp%nspectra == 3) fp%SpectrumIndex = (/INDEX_TT,INDEX_Td,INDEX_dd/)
         If (fp%nspectra == 5) fp%SpectrumIndex = (/INDEX_TT,INDEX_TE,INDEX_EE,INDEX_Td,INDEX_dd/)
      end If
      
      
      

      fp%pExp%fsky = Ini_Read_Double('fsky')
      fp%pExp%nchannels = Ini_Read_Int('nchannels')
      
      if (verbosity .gt. 1) then 
         print*, "Fsky for the experiment: ", fp%pExp%fsky 
         print*, "Number of channels:", fp%pExp%nchannels
      end if

      allocate(fp%pExp%frequency(fp%pExp%nchannels))
      allocate(fp%pExp%thetaFWHM(fp%pExp%nchannels))
      allocate(fp%pExp%DeltaTemp(fp%pExp%nchannels))
      allocate(fp%pExp%DeltaPol(fp%pExp%nchannels))
      
      freqs = Ini_Read_String('frequencies')
      read(freqs,*) fp%pExp%frequency(1:fp%PExp%nchannels)
      !print*, fp%pExp%frequency
      
      
      thetas = Ini_Read_String('beamFWHM') 
      read(thetas,*) fp%pExp%thetaFwhm(1:fp%pExp%nchannels)
      
      deltaTs = Ini_Read_String('deltaT')
      read(deltaTs,*) fp%pExp%DeltaTemp(1:fp%pExp%nchannels)

      
      deltaPs = Ini_Read_String('deltaP')
      read(deltaPs,*) fp%pExp%DeltaPol(1:fp%pExp%nchannels)
      
      fp%pExp%DeltaTemp = fp%pExp%DeltaTemp*1.d-6
      fp%pExp%DeltaPol = fp%pExp%DeltaPol*1.d-6
      
      
      call ini_close
      
      fp%nvary = 0
      do i = 1, NPARAMSMAX
         if (fp%stepsize(i) .ne. 0.0_dl) fp%nvary = fp%nvary +1 
      end do
      if (verbosity .gt. 1) write(*,'(" Number of parameters to be varied:",x,I3)') fp%nvary
      allocate(fp%varyingParamIndex(fp%nvary))
      j=1
      do i = 1, NPARAMSMAX
         if (fp%stepsize(i) .ne. 0.0_dl) then 
            fp%VaryingParamIndex(j) = i 
            j=j+1
         end if
      end do
     
      if (verbosity .gt. 1) then 
         write(*,'("The parameters being varied are")') 
         do i = 1,fp%nvary 
            print*, trim(paramNames(fp%VaryingParamIndex(i)))
         end do
      end if
      fp%isSet = .true.
    end subroutine SetFisherParamsFromFile
    
    !> Writes fiducial spectra and derivatives
    subroutine WritePS(fp,PS)
      type(FisherParams)          :: fp
      type(PowerSpectraAndDerivs) :: PS
      !      character(len=*)            :: TAG
      integer                     :: il
      character(len = 3)          :: cnvary
      
      call assert(PS%IsAllocated)
      
      open(unit=22, file = trim(OUTDIR)//'/PSAndDerivs.'//TAG)
      if ((fp%IncludeLensing .and. fp%nspectra > 3).or.(not(fp%IncludeLensing) .and. fp%nspectra>1)) &
           &open(unit=23, file = trim(OUTDIR)//'/PSEEAndDerivs.'//TAG) 
      if(fp%IncludeLensing) open(unit=25, file = trim(OUTDIR)//'/DeflectionPSAndDerivs.'//TAG)
      write(cnvary,'(I3)') fp%nvary+2
      
      write(22,'('//trim(cnvary)//'(a10,x))') "l", "Cl_fid", &
           & paramNames(fp%VaryingParamIndex(1:fp%nVary))
      
      if (fp%IncludeLensing) write(25,'('//trim(cnvary)//'(a10,x))') "l", "Cl_fid", &
           & paramNames(fp%VaryingParamIndex(1:fp%nVary))
      if ((fp%IncludeLensing .and. fp%nspectra > 3).or.(not(fp%IncludeLensing) .and. fp%nspectra>1)) &
           &write(23,'('//trim(cnvary)//'(a10,x))') "l", "Cl_fid", &
           & paramNames(fp%VaryingParamIndex(1:fp%nVary))
      write(cnvary,'(I3)') fp%nvary

      do il =  2, fp%lmax
         write(22,'(I10,x,e10.4,x'//trim(cnvary)//'(e10.4,x))') il, PS%ClFiducial(il,index_TT), PS%ClDeriv(il,index_TT,1:fp%nvary)
         if(fp%IncludeLensing .and. (fp%nspectra .eq. 3)) then 
            write(25,'(I10,x,e10.4,x'//trim(cnvary)//'(e10.4,x))') il, PS%ClFiducial(il,3), PS%ClDeriv(il,3,1:fp%nvary)
         else if (fp%IncludeLensing .and. (fp%nspectra .eq. 5)) then
            write(25,'(I10,x,e10.4,x'//trim(cnvary)//'(e10.4,x))') il, PS%ClFiducial(il,5), PS%ClDeriv(il,5,1:fp%nvary)
         end if
         if ((fp%IncludeLensing .and. fp%nspectra > 3).or.(not(fp%IncludeLensing) .and. fp%nspectra>1)) &
              & write(23,'(I10,x,e10.4,x'//trim(cnvary)//'(e10.4,x))') &
              &il, PS%ClFiducial(il,index_EE), PS%ClDeriv(il,index_EE,1:fp%nvary)
      end do
      
      close(22)
      if(fp%IncludeLensing) CLOSE(25)
      if ((fp%IncludeLensing .and. fp%nspectra > 3).or.(not(fp%IncludeLensing) .and. fp%nspectra>1))  close(23)
    end subroutine WritePS


    !> Writes fiducial spectra and noise
     subroutine WritePSAndNoise(fp,PS)
      type(FisherParams)          :: fp
      type(PowerSpectraAndDerivs) :: PS
     ! character(len=*)            :: TAG
      integer                     :: il
      character(len = 3)          :: cnvary
      call assert(PS%IsAllocated)
      !print*, "1"
      open(unit=125, file = trim(OUTDIR)//'/CMBPSAndNoise.'//TAG, status='unknown')
      
      !print*, "2"
      if (fp%IncludeLensing) then 
         open(unit=23, file = trim(OUTDIR)//'/DeflectionPSAndNoise.'//TAG)
         if(fp%nspectra .eq. 3) then
            do il = 2,fp%lmax
               write(125,'(I10,x,2(e10.4,x))') il, PS%ClFiducial(il,index_TT), PS%ClNoise(il,index_TT)
               write(23,'(I10,x,3(e10.4,x))') il, PS%ClFiducial(il,2),  PS%ClFiducial(il,3), PS%ClNoise(il,3)
            end do
         end if
         if(fp%nspectra .eq. 5) then
            do il = 2,fp%lmax
               write(125,'(I10,x,5(e10.4,x))') il, PS%ClFiducial(il,index_TT:index_EE),&
                    & PS%ClNoise(il,index_TT), PS%ClNoise(il,index_EE)
               write(23,'(I10,x,3(e10.4,x))') il, PS%ClFiducial(il,4),  PS%ClFiducial(il,5), PS%ClNoise(il,5)
            end do
         end if
         close(23)
      end if
      
      if (not(fp%IncludeLensing)) then 
         
         if(fp%nspectra .eq. 1) then
            do il = 2,fp%lmax
               write(125,'(I10,x,2(e10.4,x))') il, PS%ClFiducial(il,index_TT), PS%ClNoise(il,index_TT)
            end do
         end if
      
         if(fp%nspectra .eq. 3) then
            do il = 2,fp%lmax
               write(125,'(I10,x,5(e10.4,x))') il, PS%ClFiducial(il,index_TT:index_EE),&
                    & PS%ClNoise(il,index_TT), PS%ClNoise(il,index_EE)
             
            end do
         end if
         
      end if
      
      close(125)
    end subroutine WritePSAndNoise


    !> Writes out the Fisher Matrix results
    subroutine WriteFisherResult(fp,fullTag)
      type(FisherParams)          :: fp
      character (len = *)         :: fullTag
      integer                     :: i 
      open(unit=125, file = trim(OUTDIR)//'/FisherResults.'//fullTag, status='unknown')
      write(125,'(a16,2x,f9.6)') 'fsky',fp%pExp%fsky
            
      do i = 1, NPARAMSMAX
         write(125,'(a16,2x,3f9.6)') trim(paramNames(i)), fp%fiducial(i), fp%error(i)
      end do
      
      close(125)
      
    end subroutine WriteFisherResult

    !> Writes out the Fisher Matrix
    subroutine WriteFisherMatrix(fp,FM,fullTag)
      type(FisherParams)          :: fp
      type(FisherMatrix)          :: FM
      character (len = *)         :: fullTag
      integer                     :: i,j 
      open(unit=125, file = trim(OUTDIR)//'/FisherMatrix.'//fullTag, status='unknown')
      write(125,*) fp%nvary
      
      
      
            
      do i = 1,fp%nvary
         write(125,'(a16,2x,3f9.6)') trim(paramNames(fp%varyingParamIndex(i))),&
              & fp%fiducial(fp%varyingParamIndex(i)), fp%stepsize(fp%varyingParamIndex(i))  
      end do
      
      do i = 1, fp%nvary
         do j = 1,fp%nvary 
            write(125,'(e20.8)') FM%theFisherMatrix(i,j)
         end do
      end do

      close(125)
         
      
    end subroutine WriteFisherMatrix
end module FileIO

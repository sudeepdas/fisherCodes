!> \mainpage Documentation for fisherCodes
!! @author Sudeep Das, Princeton University
!! \section intro Introduction
!! This codes computes the CMB Fisher Matrix errors on a set of cosmological parameters,
!! optionally including CMB lensing. Currently, these parameters are supported:
!! \f[(\Omega_b h^2,\Omega_d h^2, \theta_A, f_\nu,N_\nu^{\mathrm{massive}},N_\nu^{\mathrm{masssless}},Y_{\mathrm{He}},n_s,\alpha,A_s,\tau, w,\Omega_k)\f]
!! Here \f$f_\nu\f$ is the fraction of dark matter in massive neutrinos, so that 
!! \f$\Omega_\nu h^2 = f_\nu \Omega_d h^2\f$ and the cold dark matter density is 
!! \f$ \Omega_c h^2 = (1-f_\nu) \Omega_d h^2\f$. Also, note that we use the angular size of the sound horizon
!! \f$\theta_A\f$ as a variable rather than \f$\Omega_\Lambda\f$. This lets us take the other derivatives
!! along the angular diameter distance degeneracy and makes these derivative more stable.  
!!
!! Report bugs to sudeep _at_ astro.princeton.edu
!!
!! \section compilation Compilation 
!! The code requires an exising installation of CAMB library libcamb.a . Currently  
!! tested only on March 2008 CAMB version.
!! Make necessary changes in the Makefile provided and it should compile without any problem.
!! \section basic_use Basic Usage
!! The executable generated by the code in \c FisherCombine.x which runs off of one (or more, see below)
!! parameter file, 
!! which is a text file to be named as \c params.your_experiment.
!! For example if I were doing a Fisher Matrix forecast for Planck, I could name 
!! the parameter
!! file, \c params.PLANCK . Note that \c PLANCK will be used as a tag for all result files produced 
!! by the code.
!! A typical parameter file will look like this (this my params.PLANCK): 
!! \verbinclude params.test
!! Note that each field is quite self explanatory and there is a sentence 
!! with instructions against most entries. Note also that you can combine channels. 
!! \subsection execution Execution
!! To execute the code, run \n
!!
!!  <tt> ./FisherCombine.x  params.PLANCK </tt>\n
!!
!! More than one experiment can be combined, but their parameter files must be consistent (i.e. same 
!! variables must be varied with same step sizes) 
!!
!! <tt> ./FisherCombine.x params.WMAP5 params.ACT </tt>\n
!!
!! 
!! \subsection output Output
!! The results will be stored in the directory defined by \c out_dir in the params file. Depending on whether
!! you did temperature only forecasts and included lensing, various files will be written out. Only the 
!! \c FisherResults.* file will be useful for most users, but for cross checking and lensing afficionados
!! the other files will be interesting as well.
!! 
!! <b> A. Temperature only no lensing </b> <tt> (n_spectra = 1; include_lensing=F) </tt>
!!
!! Suppose my param file was called \c params.PLANCK_TT_NOLENS \n
!! Then, upon execution four files will be produced:
!! \arg \c FisherResults.PLANCK_TT_NOLENS: Results (errorbars) from the Fisher Analysis (self explanatory).
!! \arg \c CMBPSAndNoise.PLANCK_TT_NOLENS: This file contains the fiducial CMB TT power spectrum 
!! \f$C_\ell^{TT}\f$ and the 
!! instrumental noise power spectrum \f$N_\ell^{TT}\f$. The columns are \f$\ell\f$, \f$C_\ell^{TT}\f$,
!! \f$N_\ell^{TT}\f$.
!! \arg \c FisherMatrix.PLANCK_TT_NOLENS: The actual Fisher Matrix (can be used to combine with other exps).
!! \arg \c PSAndDerivs.PLANCK_TT_NOLENS: A file with derivatives of the temperature power spectra:
!! \f$\ell,C_\ell^{TT},\partial{C_\ell^{TT}}/\partial{p_1},...,\partial{C_\ell^{TT}}/\partial{p_N}\f$
!! where p are the parameters.
!!
!! <b> B. Temperature only with lensing </b> <tt> (n_spectra = 3; include_lensing=T) </tt> 
!! 
!! Suppose my param file was called \c params.PLANCK_TT_LENS \n
!! In this case, 7 files are written out:
!! \arg \c FisherResults.PLANCK_TT_LENS
!! \arg \c CMBPSAndNoise.PLANCK_TT_LENS          
!! \arg \c FisherMatrix.PLANCK_TT_LENS
!! \arg \c PSAndDerivs.PLANCK_TT_LENS
!!
!! The above 4 are the same as the no lensing case. The next  are related to the lensing:
!! \arg \c IndivLensNoiseSpectra.PLANCK_TT_LENS: Noise bias in lensing reconstruction 
!! \f$ \ell, N_\ell^{dd}(TT)\f$  
!! \arg \c DeflectionPSAndDerivs.PLANCK_TT_LENS: Columns are 
!! \f$\ell,C_\ell^{dd},\partial{C_\ell^{dd}}/\partial{p_1},...,\partial{C_\ell^{dd}}/\partial{p_N}\f$ 
!! \arg \c DeflectionPSAndNoise.PLANCK_TT_LENS: Columns are lensing power spectra and minumum variance 
!! reconstruction noise by combining lensing estiamtors 
!! (same as in IndivLensNoiseSpectra in this case):
!! \f$ \ell, C_\ell^{Td}, C_\ell^{dd}, N_\ell^{dd}(TT)\f$  
!!
!!
!! <b> C. Temperature + polarization no lensing </b> <tt> (n_spectra = 3; include_lensing=F) </tt> 
!! 
!! Suppose my param file was called \c params.PLANCK_POL_NOLENS \n
!! In this case, 5 files are written out:
!! \arg \c FisherResults.PLANCK_POL_NOLENS 
!! \arg \c CMBPSAndNoise.PLANCK_POL_NOLENS: This time, the columns will be: 
!! \f$ \ell,C_\ell^{TT},C_\ell^{TE},C_\ell^{EE},N_\ell^{TT},N_\ell^{EE}\f$
!! \arg \c FisherMatrix.PLANCK_POL_NOLENS
!! \arg \c PSAndDerivs.PLANCK_POL_NOLENS
!! \arg \c PSEEAndDerivs.PLANCK_POL_NOLENS: EE fiducial Cls and derivatives.
!!
!! <b> D. Temperature + polarization with lensing </b> <tt> 
!! (n_spectra = 5; include_lensing=T) </tt> 
!! 
!! Suppose my param file was called \c params.PLANCK_TT_LENS \n
!! In this case, 7 files are written out:
!! \arg \c FisherResults.PLANCK_POL_LENS
!! \arg \c CMBPSAndNoise.PLANCK_POL_LENS          
!! \arg \c FisherMatrix.PLANCK_POL_LENS
!! \arg \c PSAndDerivs.PLANCK_POL_LENS
!! \arg \c PSEEAndDerivs.PLANCK_POL_LENS
!!
!! The above 4 are the same as the no lensing case. The next  are related to the lensing:
!! \arg \c IndivLensNoiseSpectra.PLANCK_POL_LENS: Noise bias in lensing reconstruction 
!! \f$ \ell, N_\ell^{dd}(TT),N_\ell^{dd}(TE),N_\ell^{dd}(TB),N_\ell^{dd}(EE),N_\ell^{dd}(EB)\f$  
!! \arg \c DeflectionPSAndDerivs.PLANCK_TT_LENS: Columns are 
!! \f$\ell,C_\ell^{dd},\partial{C_\ell^{dd}}/\partial{p_1},...,\partial{C_\ell^{dd}}/\partial{p_N}\f$ 
!! \arg \c DeflectionPSAndNoise.PLANCK_TT_LENS: Columns are lensing power spectra and minumum variance 
!! reconstruction noise by combining lensing estimators 
!! (i.e. minimum variance combination of lensing noise
!! spectra in  IndivLensNoiseSpectra file):
!! \f$ \ell,C_\ell^{Td}, C_\ell^{dd}, N_\ell^{dd}(mv)\f$  
!!
!! \section customize Things You May Want to Tweak
!! \arg Accuracy: In GetSpectra.f90 I recommend that you beef up the AccurayBoost and 
!! LAccuracyBoost parameters before getting your final results. The code will run \b much  slower,
!! but I have found it worthwhile, specially when varying curvature.
!! \arg Also, the pivot is now set at 0.05 /Mpc; You may want to change it to 0.002/Mpc in 
!! GetSpectra.f90

!<

!>Defines various global parameters and strcutures for the Fisher Matrix code

!>@author Sudeep Das, Princeton University

!>@version 1, 2008

module Fisher
  use precision
  use inifile

  character(len=80)  :: TAG
  character(len=80)  :: OUTDIR
  integer            :: verbosity
  logical            :: debugging
  integer, parameter :: NPARAMSMAX = 13
 
  
  integer, parameter :: I_OMEGA_B_H2 = 1
  integer, parameter :: I_OMEGA_C_H2 = 2
  integer, parameter :: I_NU_FRAC = 3
  integer, parameter :: I_THETA_A = 4
  integer, parameter :: I_OMEGA_K =5
  integer, parameter :: I_TAU = 6
  integer, parameter :: I_N_S = 7
  integer, parameter :: I_N_RUN = 8 
  integer, parameter :: I_A_S = 9
  integer, parameter :: I_W = 10
  integer, parameter :: I_HELIUM_FRACTION = 11
  integer, parameter :: I_NU_EFF = 12
  integer, parameter :: I_NU_MASSIVE = 13
 
  integer, parameter :: INDEX_TT = 1
  integer, parameter :: INDEX_TE = 2
  integer, parameter :: INDEX_EE = 3
  integer, parameter :: INDEX_BB = 4
  integer, parameter :: INDEX_Td = 5
  integer, parameter :: INDEX_dd = 6
  integer, parameter :: MAX_SPEC_INDEX = 6
 
  integer, parameter :: Maxx_l = 10000 !< Maximum possible l
  
  character(len=80),dimension(NPARAMSMAX),parameter:: paramNames = &
       &(/'Omegabh2', 'Omegach2','f_nu','ThetaA','Omegak','tau','ns',&
       &'nrun','As(1e-9)','w','YHe','N_nu_massless','N_nu_massive'/)
 
  !> Experimental Parameters
  type ExpParams
     character(len = ini_max_string_len) :: expFile !< not used currently
     real(dl)                 :: fsky !< fraction of sky covered
     integer                  :: nchannels !< number of channels
     real(dl),allocatable     :: frequency(:) !< frequecies of the channels (not really used)
     real(dl),allocatable     :: thetaFwhm(:) !< beam FWHM of each channel
     real(dl),allocatable     :: DeltaTemp(:) !< \f$\Delta T/T_0 \f$
     real(dl),allocatable     :: DeltaPol(:)  !< \f$\Delta P/T_0 \f$
  end type ExpParams
  
  !> A type holding the fiducial spectra, noise and derivatives
  type PowerSpectraAndDerivs
     real(dl), allocatable :: clFiducial(:,:)
     real(dl), allocatable :: clDeriv(:,:,:)
     real(dl), allocatable :: clNoise(:,:)
     logical               :: isAllocated
  end type PowerSpectraAndDerivs

 !> Spectra returned by CAMB and Noise
 type  CAMBandNoiseCLS
     real(dl) :: cl(2:Maxx_L,INDEX_TT:INDEX_DD)
     real(dl) :: clLensed(2:Maxx_L,INDEX_TT:INDEX_BB)
     real(dl) :: nl(2:Maxx_L,INDEX_TT:INDEX_DD)
  end type CAMBandNoiseCLS

  !> Information on parameters and step sizes going into
  !! the Fisher code.
  type FisherParams        

     real(dl)              :: fiducial(NPARAMSMAX)
     real(dl)              :: stepsize(NPARAMSMAX)
     real(dl)              :: error(NPARAMSMAX)
     logical               :: includeLensing
     integer               :: lMax
     integer               :: nSpectra
     integer, allocatable  :: SpectrumIndex(:)
     integer               :: nvary
     integer,allocatable   :: VaryingParamIndex(:)
     type(ExpParams)       :: pExp
     logical               :: isSet
  end type FisherParams
  
  !> The Fisher Matrix and Its inverse
  type FisherMatrix 
     real(dl), allocatable        :: theFisherMatrix(:,:)
     real(dl), allocatable        :: theInvFisherMatrix(:,:)
     logical                      :: isAllocated
     type(PowerSpectraAndDerivs)  :: PS
  end type FisherMatrix
  
end module Fisher

  
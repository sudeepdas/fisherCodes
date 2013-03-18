!> The main driver routine. 
!! @author Sudeep Das
!! @version 1, 2008
program FisherCombine
  use FileIO
  use Fisher
  use GetSpectra
  use FisherGenerate
  character(LEN=Ini_max_string_len):: InputFile
  
  type(FisherParams),allocatable   :: fp(:)
  type(FisherMatrix)   :: FM, FMtot

  
  integer              :: nExp , iExp
  character(len = 1024):: fullTag
  verbosity = 2
  InputFile = ''

  nExp = iargc()

  if (nExp == 0) STOP 'Parameter file required'
  
  write(*, '("You are combining", i2, " experments")') nExp 

  allocate(fp(1:nExp))

  tag = ''

  !Validate the consistency of the Fisher Matrices

   do iExp = 1, nExp
      if (iargc() /= 0)  call getarg(iExp,InputFile)
      if (InputFile == '') stop 'No parameter input file'
      call SetFisherParamsFromFile(InputFile,fp(iExp))
      
      if (iExp>1) then 
         call assert(fp(iExp-1)%nvary .eq.fp(iExp)%nvary ,'Number of parameters varied differ')
         call assert(all((fp(iExp-1)%VaryingParamIndex -fp(iExp)%VaryingParamIndex) .eq. 0),&
              &'The set of parameter being varied differ')
         call assert(all((fp(iExp-1)%Fiducial -fp(iExp)%Fiducial) .eq. 0),&
             &'Fiducial values differ')
         call assert(all((fp(iExp-1)%stepSize -fp(iExp)%stepSize) .eq. 0),&
             &'Some Step sizes differ')
      end if
   end do
  
  do iExp = 1, nExp
     if (iargc() /= 0)  call getarg(iExp,InputFile)
     if (InputFile == '') stop 'No parameter input file'

     
     tag = trim(Inputfile(index(InputFile,'.')+1:len(InputFile)))
     print *, "This experiment:", tag
     if (iExp .eq. 1) then 
        fullTag = trim(tag)
        call AllocFisherMatrix( Fp(iExp),FMtot)
        call FisherMatrixCreate(Fp(iExp),FM)
        FMtot = FM
        
     else
        fullTag = trim(fulltag)//'_'//trim(Inputfile(index(InputFile,'.')+1:len(InputFile)))
        
        call FisherMatrixCreate(Fp(iExp),FM)
        
        FMtot%theFisherMatrix = FMtot%theFisherMatrix + FM%theFisherMatrix
     
     end if
     
     call WritePSAndNoise(fp(iExp),FM%PS)
     call WritePS(fp(iExp),FM%PS)
     call FreeFisherMatrix(FM)
  end do

  call FisherMatrixInvert(Fp(nExp),FMtot)
  
  print*, "Results will be saved with tag: ", trim(fullTag)
  

  do i = 1, NPARAMSMAX
     write(*,'(a16,2x,4f10.6)') trim(paramNames(i)), fp(nExp)%fiducial(i), fp(nExp)%error(i) ,fp(nExp)%error(i)/ fp(nExp)%fiducial(i)*100. 
  end do
    

  call writeFisherResult(fp(nExp),fullTag)
  
  call writeFisherMatrix(fp(nExp),FMtot,fullTag)
 
end program FisherCombine

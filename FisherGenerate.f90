!> Evaluates the Fisher Matrix , Its inverse and ErrorBars
!!@author Sudeep Das, Princeton University
!!@version 1, 2008
module FisherGenerate
  use Fisher
  use GetSpectra
  use miscUtilsForCAMB
  implicit none

 contains
   !>Allocates memory for FisherMatrix
    subroutine AllocFisherMatrix(fp,FM) 
      type(FisherParams) :: fp
      type(FisherMatrix) :: FM

      call assert(fp%isSet,'Need to read in FisherParams first')
      allocate(FM%theFisherMatrix(fp%nvary,fp%nvary))
      allocate(FM%theInvFisherMatrix(fp%nvary,fp%nvary))
      FM%theFisherMatrix  = 0._dl
      FM%theInvFisherMatrix = 0._dl
      FM%IsAllocated = .true.
    end subroutine AllocFisherMatrix
    
    !> Frees memory
    subroutine FreeFisherMatrix(FM)
      type(FisherParams) :: fp
      type(FisherMatrix) :: FM

      
      deallocate(FM%theFisherMatrix)
      deallocate(FM%theInvFisherMatrix)
      
      FM%IsAllocated = .false.
    end subroutine FreeFisherMatrix
    
    !> Creates the FisherMatrix from values given in params file
    subroutine FisherMatrixCreate(fp,FM)
      type(FisherParams) :: fp
      type(FisherMatrix) :: FM
      type(PowerSpectraAndDerivs) :: PS 
      real(dl)                    :: covMatrix(2:fp%lmax,1:fp%nspectra,1:fp%nspectra)
      real(dl)                    :: invcovMatrix(2:fp%lmax,1:fp%nspectra,1:fp%nspectra)
      real(dl)                    :: lfactor(2:fp%lmax)
      integer                     :: p1,p2,spec1,spec2,il
      integer                     :: idd, iTd


      call AllocFisherMatrix(fp,FM)

      !Get Fiducial PowerSpec,  Derivatives and Noise
      call GetDerivsAndNoise(fp,PS)
      FM%PS = PS
      !Create the covariance matrix
      do il = 2,fp%lmax
         lfactor(il) = (2._dl/((2._dl*real(il,dl)+1)*fp%PExp%fsky))
      end do

      covMatrix = 0._dl


      if (not(fp%includeLensing)) then 
         !TT X  TT 

         covMatrix(2:fp%lmax, INDEX_TT,INDEX_TT) =&
              & (Ps%clFiducial(2:fp%lmax,INDEX_TT)+ PS%clNoise(2:fp%lmax,INDEX_TT))*&
              & (Ps%clFiducial(2:fp%lmax,INDEX_TT)+ PS%clNoise(2:fp%lmax,INDEX_TT))*&
              & lfactor(2:fp%lmax)


         if(fp%nspectra .ge. 3) then 
            !TT X TE
            covMatrix(2:fp%lmax, INDEX_TT,INDEX_TE) = &
                 & ((Ps%clFiducial(2:fp%lmax,INDEX_TT)+ PS%clNoise(2:fp%lmax,INDEX_TT)) * &
                 & (Ps%clFiducial(2:fp%lmax,INDEX_TE)+ PS%clNoise(2:fp%lmax,INDEX_TE))) * &
                 & lfactor(2:fp%lmax) 
            !TT X EE
            covMatrix(2:fp%lmax, INDEX_TT,INDEX_EE) = PS%clFiducial(2:fp%lmax,INDEX_TE)**2*&
                 & lfactor(2:fp%lmax) 

            !TE X TT
            covMatrix(:,INDEX_TE,INDEX_TT) = covMatrix(:,INDEX_TT, INDEX_TE)

            !TE X TE 
            covMatrix(2:fp%lmax, INDEX_TE,INDEX_TE) = &
                 & ((Ps%clFiducial(2:fp%lmax,INDEX_TE)+ PS%clNoise(2:fp%lmax,INDEX_TE)) * &
                 &  (Ps%clFiducial(2:fp%lmax,INDEX_TE)+ PS%clNoise(2:fp%lmax,INDEX_TE)) + &
                 &  (Ps%clFiducial(2:fp%lmax,INDEX_TT)+ PS%clNoise(2:fp%lmax,INDEX_TT)) * &
                 &  (Ps%clFiducial(2:fp%lmax,INDEX_EE)+ PS%clNoise(2:fp%lmax,INDEX_EE))) * &
                 & lfactor(2:fp%lmax)/2.

            !TE X EE
            covMatrix(2:fp%lmax, INDEX_TE,INDEX_EE) = &
                 & ((Ps%clFiducial(2:fp%lmax,INDEX_EE)+ PS%clNoise(2:fp%lmax,INDEX_EE)) * &
                 & (Ps%clFiducial(2:fp%lmax,INDEX_TE)+ PS%clNoise(2:fp%lmax,INDEX_TE))) * &
                 & lfactor(2:fp%lmax) 

            !EE X TT 
            covMatrix(:,INDEX_EE,INDEX_TT) = covMatrix(:,INDEX_TT,INDEX_EE)

            !EE X TE
            covMatrix(:,INDEX_EE,INDEX_TE) = covMatrix(:,INDEX_TE,INDEX_EE)

            !EE X EE
            covMatrix(2:fp%lmax, INDEX_EE,INDEX_EE) =&
                 & (Ps%clFiducial(2:fp%lmax,INDEX_EE)+ PS%clNoise(2:fp%lmax,INDEX_EE))*&
                 & (Ps%clFiducial(2:fp%lmax,INDEX_EE)+ PS%clNoise(2:fp%lmax,INDEX_EE))*&
                 & lfactor(2:fp%lmax)

         end if

         if (fp%nspectra .eq. 4) then !do BB,BB
            covMatrix(2:fp%lmax,INDEX_BB,INDEX_BB) = &
                 & (Ps%clFiducial(2:fp%lmax,INDEX_BB)+ PS%clNoise(2:fp%lmax,INDEX_BB))*&
                 & (Ps%clFiducial(2:fp%lmax,INDEX_BB)+ PS%clNoise(2:fp%lmax,INDEX_BB))*&
                 & lfactor(2:fp%lmax)
         end if
      end if


      if (fp%includeLensing) then 
         !Using Perotto et al (2006) but dropping the -ve terms from eqs 4.5 and 4.6
         !TT X  TT 

         covMatrix(2:fp%lmax, 1,1) =&
              & (Ps%clFiducial(2:fp%lmax,INDEX_TT)+ PS%clNoise(2:fp%lmax,INDEX_TT))*&
              & (Ps%clFiducial(2:fp%lmax,INDEX_TT)+ PS%clNoise(2:fp%lmax,INDEX_TT))*&
              & lfactor(2:fp%lmax)
         
         if (fp%nspectra .eq. 3) then 
            !TT X Td
            iTd = 2
            idd = 3
            covMatrix(2:fp%lmax,index_TT,iTd) = &
                 & (Ps%clFiducial(2:fp%lmax,iTd)* &
                 & (Ps%clFiducial(2:fp%lmax,INDEX_TT)+ PS%clNoise(2:fp%lmax,INDEX_TT)))*&
                 & lfactor(2:fp%lmax) 
            
            !TTxdd
            covMatrix(2:fp%lmax,index_TT,idd) = &
                 & (Ps%clFiducial(2:fp%lmax,iTd))**2*lfactor(2:fp%lmax)
            
            !TdXTT 
            covMatrix(2:fp%lmax,iTd,index_TT) = covMatrix(2:fp%lmax,index_TT,iTd)
            
            !TdXTd 
            covMatrix(2:fp%lmax,iTd,iTd)  = &
                 & ((Ps%clFiducial(2:fp%lmax,iTd))**2 + &
                 & (Ps%clFiducial(2:fp%lmax,INDEX_TT)+ PS%clNoise(2:fp%lmax,INDEX_TT))&
                 &*(Ps%clFiducial(2:fp%lmax,idd)+ PS%clNoise(2:fp%lmax,idd)))*&
                 & lfactor(2:fp%lmax)/2._dl
            
            !TdXdd 
            covMatrix(2:fp%lmax,iTd,idd) = &
                 &  Ps%clFiducial(2:fp%lmax,iTd)*&
                 & (Ps%clFiducial(2:fp%lmax,idd)+PS%clNoise(2:fp%lmax,idd))*lfactor(2:fp%lmax)
            
            !ddXTT 
            covMatrix(2:fp%lmax,idd,index_TT) = covMatrix(2:fp%lmax,index_TT,idd) 
            
            !ddXTd 
            covMatrix(2:fp%lmax,idd,iTd) = covMatrix(2:fp%lmax,iTd,idd) 
            
            !ddXdd 
            covMatrix(2:fp%lmax,idd,idd) = (Ps%clFiducial(2:fp%lmax,idd)+&
                 &PS%clNoise(2:fp%lmax,idd))**2*lfactor(2:fp%lmax)
        
         end if


         if (fp%nspectra .eq. 5) then
            iTd = 4
            idd = 5 
             !TT X TE
            covMatrix(2:fp%lmax, INDEX_TT,INDEX_TE) = &
                 & ((Ps%clFiducial(2:fp%lmax,INDEX_TT)+ PS%clNoise(2:fp%lmax,INDEX_TT)) * &
                 & (Ps%clFiducial(2:fp%lmax,INDEX_TE)+ PS%clNoise(2:fp%lmax,INDEX_TE))) * &
                 & lfactor(2:fp%lmax) 
            !TT X EE
            covMatrix(2:fp%lmax, INDEX_TT,INDEX_EE) = PS%clFiducial(2:fp%lmax,INDEX_TE)**2*&
                 & lfactor(2:fp%lmax) 

            !TE X TT
            covMatrix(:,INDEX_TE,INDEX_TT) = covMatrix(:,INDEX_TT, INDEX_TE)

            !TE X TE 
            covMatrix(2:fp%lmax, INDEX_TE,INDEX_TE) = &
                 & ((Ps%clFiducial(2:fp%lmax,INDEX_TE)+ PS%clNoise(2:fp%lmax,INDEX_TE)) * &
                 &  (Ps%clFiducial(2:fp%lmax,INDEX_TE)+ PS%clNoise(2:fp%lmax,INDEX_TE)) + &
                 &  (Ps%clFiducial(2:fp%lmax,INDEX_TT)+ PS%clNoise(2:fp%lmax,INDEX_TT)) * &
                 &  (Ps%clFiducial(2:fp%lmax,INDEX_EE)+ PS%clNoise(2:fp%lmax,INDEX_EE))) * &
                 & lfactor(2:fp%lmax)/2.

            !TE X EE
            covMatrix(2:fp%lmax, INDEX_TE,INDEX_EE) = &
                 & ((Ps%clFiducial(2:fp%lmax,INDEX_EE)+ PS%clNoise(2:fp%lmax,INDEX_EE)) * &
                 & (Ps%clFiducial(2:fp%lmax,INDEX_TE)+ PS%clNoise(2:fp%lmax,INDEX_TE))) * &
                 & lfactor(2:fp%lmax) 

            !EE X TT 
            covMatrix(:,INDEX_EE,INDEX_TT) = covMatrix(:,INDEX_TT,INDEX_EE)

            !EE X TE
            covMatrix(:,INDEX_EE,INDEX_TE) = covMatrix(:,INDEX_TE,INDEX_EE)

            !EE X EE
            covMatrix(2:fp%lmax, INDEX_EE,INDEX_EE) =&
                 & (Ps%clFiducial(2:fp%lmax,INDEX_EE)+ PS%clNoise(2:fp%lmax,INDEX_EE))*&
                 & (Ps%clFiducial(2:fp%lmax,INDEX_EE)+ PS%clNoise(2:fp%lmax,INDEX_EE))*&
                 & lfactor(2:fp%lmax)

            !TTxdd
            covMatrix(2:fp%lmax,index_TT,idd) = &
                 & (Ps%clFiducial(2:fp%lmax,iTd))**2*lfactor(2:fp%lmax)
            
            !TdXTT 
            covMatrix(2:fp%lmax,iTd,index_TT) = covMatrix(2:fp%lmax,index_TT,iTd)
            
            !TdXTd 
            covMatrix(2:fp%lmax,iTd,iTd)  = &
                 & ((Ps%clFiducial(2:fp%lmax,iTd))**2 + &
                 & (Ps%clFiducial(2:fp%lmax,INDEX_TT)+ PS%clNoise(2:fp%lmax,INDEX_TT))&
                 &*(Ps%clFiducial(2:fp%lmax,idd)+ PS%clNoise(2:fp%lmax,idd)))*&
                 & lfactor(2:fp%lmax)/2._dl
            
            !TdXdd 
            covMatrix(2:fp%lmax,iTd,idd) = &
                 &  Ps%clFiducial(2:fp%lmax,iTd)*&
                 & (Ps%clFiducial(2:fp%lmax,idd)+PS%clNoise(2:fp%lmax,idd))*lfactor(2:fp%lmax)
            
            !ddXTT 
            covMatrix(2:fp%lmax,idd,index_TT) = covMatrix(2:fp%lmax,index_TT,idd) 
            
            !ddXTd 
            covMatrix(2:fp%lmax,idd,iTd) = covMatrix(2:fp%lmax,iTd,idd) 
            
            !ddXdd 
            covMatrix(2:fp%lmax,idd,idd) = (Ps%clFiducial(2:fp%lmax,idd)+&
                 &PS%clNoise(2:fp%lmax,idd))**2*lfactor(2:fp%lmax)
         end if
         
      end if


      invCovMatrix = 0._dl
      !invert the covariance matrix
      if (fp%nspectra .eq. 1) then 
         invCovMatrix(2:fp%lmax,1,1) = 1._dl/CovMatrix(2:fp%lmax,1,1)
      else

         do il = 2,fp%lmax 

            call InvertMatrix(fp%nspectra,covMatrix(il,:,:),invCovMatrix(il,:,:))

            call  assert(abs(sum(matmul(invCovMatrix(il,:,:),covMatrix(il,:,:))) &
                 &- real(fp%nspectra,dl)) .le. 1.e-8,' FisherMatrixCreate: &
                 & Faliure to Invert covariance matrix ')        
         end do
      end if

      If (verbosity > 0) write(*,'("Computing Fisher Matrix ")')
      !Now compute the Fisher Matrix
      FM%theFisherMatrix = 0._dl
      do p1 = 1,fp%nvary
         do p2 = 1,fp%nvary
            do il = 2,fp%lmax
               do spec1 = 1,fp%nspectra
                  do spec2 = 1, fp%nspectra
                     FM%theFisherMatrix(p1,p2) = FM%theFisherMatrix(p1,p2) + &
                          & PS%clDeriv(il,spec1,p1)*invCovMatrix(il,spec1,spec2)*PS%clDeriv(il,spec2,p2)
                  end do
               end do
            end do
         end do
      end do
      !Compute the inverse Fisher Matrix
      !call InvertMatrix(fp%nvary,FM%theFisherMatrix,FM%theInvFisherMatrix)

      !FM%theInvFisherMatrix = MatInv(fp%nvary,FM%theFisherMatrix)

      !write(*,*) sum(matmul(FM%theInvFisherMatrix,FM%theFisherMatrix)), fp%nvary
      !print*,matmul(FM%theInvFisherMatrix,FM%theFisherMatrix)
       

      !do p1 = 1,fp%nvary
      !   fp%error(fp%VaryingParamIndex(p1)) = sqrt(FM%theInVFisherMatrix(p1,p1)) 
      !end do

    end subroutine FisherMatrixCreate

    !> Inverts the FisherMatrix
    subroutine FisherMatrixInvert(fp,FM)
      type(FisherMatrix) :: FM
      type(FisherParams) :: fp

      integer            :: p1
      
      call assert( fp%isSet,'Fisher Params are not set')
      call assert( FM%isAllocated, 'Create Fisher matrix first')
      !Compute the inverse Fisher Matrix
      If (verbosity > 0) write(*,'("Computing Inverse Fisher Matrix ")')

        !FM%theInvFisherMatrix = MatInv(fp%nvary,FM%theFisherMatrix)
      call InvertMatrix(fp%nvary,FM%theFisherMatrix,FM%theInvFisherMatrix)
      !call  assert(abs(sum(matmul(FM%theInvFisherMatrix,FM%theFisherMatrix))) &
      !     &- real(fp%nvary,dl) .le. 1.e-4,' FisherMatrixCreate: &
      !     & Faliure to Invert Fisher matrix ') 
      !print*, abs(sum(matmul(FM%theInvFisherMatrix,FM%theFisherMatrix))),real(fp%nvary,dl)
      do p1 = 1,fp%nvary
         fp%error(fp%VaryingParamIndex(p1)) = sqrt(FM%theInVFisherMatrix(p1,p1)) 
         if(fp%varyingParamIndex(p1) .eq. I_A_S) &
              &fp%error(fp%VaryingParamIndex(p1)) = fp%error(fp%VaryingParamIndex(p1))*1.d9 
      end do
    end subroutine FisherMatrixInvert

    !> Routine for matrix inversion using LU decomposition
    subroutine InvertMatrix(n, matrix,invmatrix)
      integer                  :: n 
      real(dl), dimension(n,n) :: matrix
      real(dl), dimension(n,n) :: invmatrix
      real(dl), dimension(n,n) :: temp
      
      real(dl), dimension(n)   :: colLU
      real(dl)                 :: dLU
      integer, dimension(n)    :: indx
      integer                  :: i 
      
      temp = matrix 
      call ludcmp(temp,indx,dLU)
      
      do i = 1,n
         colLU(:) = 0._dl
         colLU(i) = 1._dl
  
         call lubksb(temp,indx,colLU)
         !print*, colLU
         invMatrix(:,i) = colLU(:)
      end do
      
    end subroutine InvertMatrix

!---------------------------------------------------------------------------

!> fotran-90 routines for matrix inversion
    
    FUNCTION outerprod(a,b)
      REAL(dl), DIMENSION(:), INTENT(IN) :: a,b
      REAL(dl), DIMENSION(size(a),size(b)) :: outerprod
      outerprod = spread(a,dim=2,ncopies=size(b)) * &
           spread(b,dim=1,ncopies=size(a))
    END FUNCTION outerprod
    
    SUBROUTINE swap(a,b)
      REAL(DL), INTENT(INOUT) :: a(:),b(:)
      REAL(DL) :: dum(size(a))
      dum=a
      a=b
      b=dum
    END SUBROUTINE swap
    
    FUNCTION imaxloc(arr)
      REAL(DL), DIMENSION(:), INTENT(IN) :: arr
      INTEGER :: imaxloc
      INTEGER, DIMENSION(1) :: imax
      imax=maxloc(arr(:))
      imaxloc=imax(1)
    END FUNCTION imaxloc
    
    SUBROUTINE ludcmp(a,indx,d)
      !USE nrtype; USE nrutil, ONLY : assert_eq,imaxloc,nrerror,outerprod,swap
      IMPLICIT NONE
      REAL(DL), DIMENSION(:,:), INTENT(INOUT) :: a
      INTEGER, DIMENSION(:), INTENT(OUT) :: indx
      REAL(DL), INTENT(OUT) :: d
      REAL(DL), DIMENSION(size(a,1)) :: vv
      REAL(DL), PARAMETER :: TINY=1.0e-20_dl
      INTEGER :: j,n,imax
      !n=assert_eq(size(a,1),size(a,2),size(indx),'ludcmp')
      n = size(a,1)
      d=1.0
      vv=maxval(abs(a),dim=2)
      if (any(vv == 0.0)) STOP 'singular matrix in ludcmp'
      vv=1.0_dl/vv
      do j=1,n
         imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j)))
         if (j /= imax) then
            call swap(a(imax,:),a(j,:))
            d=-d
            vv(imax)=vv(j)
            
         end if
         indx(j)=imax
         if (a(j,j) == 0.0) a(j,j)=TINY
         a(j+1:n,j)=a(j+1:n,j)/a(j,j)
         a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
         
      end do
    END SUBROUTINE ludcmp
    
   SUBROUTINE lubksb(a,indx,b)
    
    IMPLICIT NONE
    REAL(DL), DIMENSION(:,:), INTENT(IN) :: a
    INTEGER, DIMENSION(:), INTENT(IN) :: indx
    REAL(DL), DIMENSION(:), INTENT(INOUT) :: b
    INTEGER :: i,n,ii,ll
    REAL(DL) :: summ
    !n=assert_eq(size(a,1),size(a,2),size(indx),'lubksb')
    n = size(a,1)
    ii=0
    do i=1,n
        ll=indx(i)
        summ=b(ll)
        b(ll)=b(i)
        if (ii /= 0) then
            summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
        else if (summ /= 0.0) then
            ii=i
        end if
        b(i)=summ
    end do
    do i=n,1,-1
        b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
    end do
    END SUBROUTINE lubksb
    
!!$
!!$!> code for inverting a matrix (uses lapack)
!!$    function MatInv(N,A)
!!$      IMPLICIT NONE
!!$      
!!$      integer, intent(in) :: N
!!$      INTEGER          :: INFO1, INFO2
!!$      INTEGER          :: IPIV(N)
!!$      double precision :: A(N,N), WORK(N),B(N,N)
!!$      double precision :: MatINv(N,N)
!!$      
!!$      EXTERNAL         DGETRF, DGETRI
!!$      
!!$      B = A !do not destroy A 
!!$      !*        Factorize 
!!$      CALL DGETRF(N,N,B,N,IPIV,INFO1)
!!$      !*
!!$      IF (INFO1.EQ.0) THEN
!!$         if (verbosity .ge. 3) print*,'LU decomposition done...&
!!$              & Now computing inverse'
!!$         
!!$         !*           Compute inverse           
!!$         CALL DGETRI(N,B,N,IPIV,WORK,N,INFO2)
!!$         if (INFO2 .ne. 0) then 
!!$            print*, 'MatInv> Program DGETRI returned exit code, INFO =',info2
!!$            STOP
!!$         end if
!!$      ELSE
!!$         
!!$         STOP' MatInv_Gen, The factor D is singular'
!!$      END IF
!!$      if (verbosity .ge. 3) print*,'Done...'
!!$      MatInv = B
!!$    end function MatInv

end module FisherGenerate

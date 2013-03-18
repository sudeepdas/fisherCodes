!Time-stamp: <2008-12-13 18:50:25 sudeep>
!>this module calculates the noise power spectra from the instrument and the 
!! noise for lensing reconstruction of the lensing deflection field, \f$ N_\ell^{dd}\f$
!! @author Sudeep Das
!! @version 1, 2008
!01/15/2006: adding polarization reconstruction. IN cldef.f90 added Integer 
!polcombi to control whether just TT reconstruction or Reconstruction using 
!Polarized CMB cross correlations is to be performed. 


!01/18/2006: Changed |l-L| from an integer to real. Interpolate and extrapolate
!to find C_|l-L| inside and outside the interval [lmin,lmax] respectively.
!see code lensing_noise_qromb.f90 for other details.

!01/18/2006: adding logical "doSharpGradFilter" to use the sharp cutoff at 
!high l's for the Gradient field as proposed in astro-ph/070127(Hu,DeDeo,Vale)

module LensingNoise
  use Fisher
  
  
  real, parameter      :: degree=Pi/180.
  real, parameter      :: arcmin=degree/60.
  
  integer, parameter   :: lspline = 100!lspline>24, lspinle~maxl/50 is good
  integer, parameter   :: lmin = 2
  integer, parameter   :: polcombimax=7 
  integer              :: polcombi !< no of polarization combinations  (1 or 7)
  logical, parameter   :: HuOkamoto2002 = .true. !< Does calcuations ala Hu and Okamoto http://arxiv.org/abs/astro-ph/0111606  .
  integer ,parameter   :: lens_index_TT=1,lens_index_TE=2,lens_index_TB=3,lens_index_ET=4,lens_index_EE=5,&
       & lens_index_EB=6, lens_index_BB=7
  
  public
  private::degree,arcmin, noiselspline,nlspline,lfunc,lspline,spline,splint,qromb,trapzd,polint

contains

  integer function lfunc(i,lmax)
    !-------- 
    ! The l values used for splining
    !-------- 
    integer :: i, l,lmax
    real :: x
    real :: facl
    facl = (real(lmax)/25.0)**(1./real(lspline-24)) !factor used in lfunc
    if(facl.le.1.) stop 'lspline too many'
    if (i .lt. 25) then
       l = i
    else

       l = int(25.0*facl**(i-24))
    end if

    lfunc = l

    return
  end function lfunc

 

  subroutine GetLensingNoise(lmax, polcombi, allcls)
    
    implicit none
    !calculates noise N_l^{dd}in quadratic estimation of deflection power 
    !spectrum,Cl^{dd}

    !assumes instrument noise and CAMBcls in allcls are set

    type(CAMBAndNoiseCls) :: allcls
    integer               :: lmax, polcombi
    real(dl)              :: cl(lmin:lmax,1:polcombiMAX),cl_mv(lmin:lmax)
    

    real(dl),dimension(lmin:lmax,1:polcombiMAX)::cltot,clunlensed,cllensed,cl_instrument
    real(dl),dimension(lmin:lmax,1:polcombiMAX)::GradFilter
    real(dl),dimension(lmin:lmax,1:3)::LensFilter


    integer::bigL,ell,il,polnow,ii,lcut,polcombi_noBB,polcombi_upper
    
    integer,parameter::lens_index_T=1,lens_index_E=2,lens_index_B=3
    real::clval
    real(dl),parameter::theta0=0._dl
    real(dl),parameter::theta1=2._dl*Pi

    real(dl)::value, mv 
    real :: larr(lspline-lmin+1), clarr(lspline-lmin+1), cl2arr(lspline-lmin+1)

    character(len=2),dimension(7)::PolX=(/"TT","TE","TB","ET","EE","EB","BB"/)

    

    cl_instrument(:,:)=0._dl
    !read in the instrument noise
    cl_instrument(lmin:lmax,lens_index_TT)=ALLCls%nl(lmin:lmax,INDEX_TT)
    cl_instrument(lmin:lmax,lens_index_EE)=ALLCls%nl(lmin:lmax,INDEX_EE)
    cl_instrument(lmin:lmax,lens_index_BB)=ALLCls%nl(lmin:lmax,INDEX_BB)
    
    

    clunlensed(:,:)=0._dl
    cllensed(:,:)=0._dl
    lensFilter=0._dl
    GradFilter=0._dl

   
    clunlensed(lmin:lmax,lens_index_TT) = allcls%cl(lmin:lmax,INDEX_TT) 
    cllensed(lmin:lmax,lens_index_TT) = allcls%cllensed(lmin:lmax,INDEX_TT)


    if(polcombi .gt. 1) then 

       clunlensed(lmin:lmax,lens_index_TE) = allcls%cl(lmin:lmax,INDEX_TE)
       clunlensed(lmin:lmax,lens_index_ET) = allcls%cl(lmin:lmax,INDEX_TE)
       clunlensed(lmin:lmax,lens_index_EE) = allcls%cl(lmin:lmax,INDEX_EE)

       cllensed(lmin:lmax,lens_index_TE)=&
            & allcls%cllensed(lmin:lmax,INDEX_TE)
       cllensed(lmin:lmax,lens_index_ET)=&
            & allcls%cllensed(lmin:lmax,INDEX_TE)
       cllensed(lmin:lmax,lens_index_EE)=&
            & allcls%cllensed(lmin:lmax,INDEX_EE)
       cllensed(lmin:lmax,lens_index_BB)=&
            & allcls%cllensed(lmin:lmax,INDEX_BB)

    end if
   
    cltot(:,:)=cllensed(:,:)+cl_instrument(:,:) !total cls
    
    !LensFilter and GradFilter are used in the Hu,Dedeo,Vale calculation 
    !Not used if HuAndOkmoto is true

    LensFilter(:,lens_index_T)=1._dl/Cltot(:,lens_index_TT) !Lensing Field Filter

    GradFilter(:,lens_index_TT)=clunlensed(:,lens_index_TT)/cltot(:,lens_index_TT) !Gradient 
    !Field Filter

    if(polcombi .gt. 1) then
       LensFilter(:,lens_index_E)=1._dl/Cltot(:,lens_index_EE)
       LensFilter(:,lens_index_B)=1._dl/Cltot(:,lens_index_BB)
       GradFilter(:,lens_index_TE)=clunlensed(:,lens_index_TE)/cltot(:,lens_index_TT)
       GradFilter(:,lens_index_TB)=clunlensed(:,lens_index_TE)/cltot(:,lens_index_TT)
       GradFilter(:,lens_index_ET)=clunlensed(:,lens_index_ET)/cltot(:,lens_index_EE)
       GradFilter(:,lens_index_EE)=clunlensed(:,lens_index_EE)/cltot(:,lens_index_EE)
       GradFilter(:,lens_index_EB)=clunlensed(:,lens_index_EE)/cltot(:,lens_index_EE)
    end if

    
    
   
    polcombi_noBB=polcombiMAX-1 !skipping BB estimator 
    
    polcombi_upper=polcombi_noBB
    iF(polcombi .eq. 1) polcombi_upper=1
    
    do polnow=1,polcombi_upper !loop over polarization combinations (TT,TE etc)
       
       do il=lmin,lspline

          larr(il-lmin+1)=lfunc(il,lmax)
          !print*,lfunc(il)
          bigL=larr(il-lmin+1)
          
          clval=0.
          do ell=lmin,lmax
            

             call qromb(noisekernel,theta0,theta1,value)
             clval=clval+value
          end do

          
          if (clval .lt. 1.E-22) then
             print*, 'underflow', clval,bigL
             clval=1.e-22
          end if
          
          
          clarr(il-lmin+1) = (real(bigL,dl)*real(bigL+1,dl))/(clval)

       end do





       !-------- Spline the logarithms (these will be more stable)


       call spline(log(larr), log(clarr), lspline-lmin+1, 1.e30, 1.e30, cl2arr)

       do il = lmin, lmax
          call splint(log(larr), log(clarr), cl2arr, lspline-lmin+1, log(real(il)), clval)
          cl(il,polnow) = exp(clval) 

       end do

       
       
       print*, 'done computing N_l_'//Polx(polnow)
    end do !end loop over polarizations

    cl_mv(lmin:lmax) = cl(lmin:lmax,lens_index_TT)

    If(polcombi .gt. 1) then 
       do il=lmin,lmax
          mv=sum(1._dl/real(cl(il,1:polcombi_noBB),dl))
          if(HuOkamoto2002) mv=sum(1._dl/real(cl(il,(/1,2,3,5,6/)),dl))
          cl_mv(il)=real(1._dl/mv,sp)
          
       end do
    end If
    !this ignores cross correlation between estimators
    
    allcls%nl(lmin:lmax,INDEX_DD) = cl_mv(lmin:lmax)
    
    !dump out the individual polarization noise spectra in a file
    open(unit=22,file= trim(OUTDIR)//'/IndivLensNoiseSpectra.'//trim(TAG))
    do il = lmin,lmax
       if (polcombi .gt. 1) then 
          write(22,'(i5,x,5(e20.8,x))') il, cl(il,(/1,2,3,5,6/))
       else
          write(22,'(i5,x,e20.8)') il, cl(il,1)
       end if
    end do
    close(22)
  contains  
    
    function noisekernel(theta)
      use precision  
      !this is the kernel for the theta integrsl
      real(dl), intent(in)::theta
      real(dl)::noisekernel

      real(dl):: Ldotl,LdotLminusl,flL(1:polcombiMAX),phi_l_Lminusl,elldotLminusl
      real(dl):: cos2phi,sin2phi,cosphi
      real(dl)::cl2unlensed(1:polcombiMAX),cl2tot(1:polcombiMAX) !interpolated values
      real(dl)::lensFilterl2(1:3),gradFilterl2(1:polcombiMAX)
      integer::llow,lhigh
      real(dl)::deltal
      real(dl)::kernel(1:polcombiMAX),mykernel
      real(dl)::modLminusl
      real::clvalue
      integer::i

      !returns the kernel for noise evaluation given the polnow (what is the 
      !polarization now?) and the value of bigL and small l 



      modLminusl=(sqrt(real(ell,dl)**2+real(bigL,dl)**2-2.*ell*bigL*cos(theta)))
      if(modLminusl .eq. 0) modlminusL=0.0001_dl
      !we shall use a non-integer value for |L-l| 
      !this is to avoid discontinuities in the integration kernel
      !that leads to instability and poor convergence of the integral
      !for the value of C_|L-l| we shall linearly interpolate between the 
      !C_l's for the two integral l's that bracket |L-l|
      !and extrapolate for |L-l| outside the interval.

      if(modLminusl.Lt. real(lmin,dl)) then !extrapolate

         deltal=real(lmin,dl)-modLminusl

         cl2unlensed(:)=clunlensed(lmin,:)&
              &+deltal*(clunlensed(lmin,:)-clunlensed(lmin+1,:))
         cl2tot(:)=cltot(lmin,:)&
              &+deltal*(cltot(lmin,:)-cltot(lmin+1,:))
         Lensfilterl2(:)=lensfilter(lmin,:)&
              &+deltal*(lensfilter(lmin,:)-lensfilter(lmin+1,:))
         gradfilterl2(:)=gradfilter(lmin,:)&
              &+deltal*(gradfilter(lmin,:)-gradfilter(lmin+1,:))

      else if(modLminusl .gt. real(lmax,dl)) then !extrapolate
         !smooth cutoff with a gaussian instead of abruptly
         !setting to zero.
         deltal=modLminusl-real(lmax,dl)
         cl2unlensed(:)=clunlensed(LMAX,:)&
              &+deltal*(clunlensed(lmax,:)-clunlensed(lmax-1,:))&
              &*exp(-deltal**2/1.d1)
         cl2tot(:)=cltot(lmax,:)&
              &+deltal*(cltot(lmax,:)-cltot(lmax-1,:))&
              &*exp(-deltal**2/1.d1)
         Lensfilterl2(:)=lensfilter(lmax,:)&
              &+deltal*(lensfilter(lmax,:)-lensfilter(lmax-1,:))&
              &*exp(-deltal**2/1.d1)

         gradfilterl2(:)=gradfilter(lmax,:)&
              &+deltal*(gradfilter(lmax,:)-gradfilter(lmax-1,:))&
              &*exp(-deltal**2/1.d1)
      else !interpolate

         llow=floor(modLminusl)
         lhigh=ceiling(modlminusl) !these are the l's bracketing l2=|L-l|

         !linearly interpolate the values at lhigh, llow to get value at l2


         deltal=modlminusL-real(llow,dl)

         cl2unlensed(:)=clunlensed(llow,:)&
              &-deltal*(clunlensed(llow,:)-clunlensed(lhigh,:))
         cl2tot(:)=cltot(llow,:)&
              &-deltal*(cltot(llow,:)-cltot(lhigh,:))
         Lensfilterl2(:)=lensfilter(llow,:)&
              &-deltal*(lensfilter(llow,:)-lensfilter(lhigh,:))
         gradfilterl2(:)=gradfilter(llow,:)&
              &-deltal*(gradfilter(llow,:)-gradfilter(lhigh,:))
      end if


      Ldotl=bigL*ell*cos(theta)
      LdotLminusl=real(bigL,dl)**2-Ldotl
      elldotLminusl=Ldotl-real(ell,dl)**2
      !phi_l_Lminusl=acos(ldotLminusl/(real(bigL,dl)*real(modLminusl,dl)))
      phi_l_Lminusl=acos(elldotLminusl/(real(ell,dl)*real(modLminusl,dl)))
      
      cos2phi=cos(2*phi_l_Lminusl)
      sin2phi=sin(2*phi_l_Lminusl)
      cosphi=cos(phi_l_Lminusl)
      
      select case (polnow) !set the kernel according to the polarization

      case (lens_index_TT)   

         flL(lens_index_TT)=clunlensed(ell,lens_index_TT)*Ldotl&
              &+cl2unlensed(lens_index_TT)*LdotLminusl
         mykernel=GradFilter(ell,lens_index_TT)&
              &*Lensfilterl2(lens_index_T)*flL(lens_index_TT)

         if(HuOkamoto2002) mykernel=flL(lens_index_TT)**2&
              &/(2.*cltot(ell,lens_index_TT)*cl2tot(lens_index_TT)) !eq 14.

      case(lens_index_TE)
         flL(lens_index_TE)=clunlensed(ell,lens_index_TE)*Ldotl*cos2phi&
              &+cl2unlensed(lens_index_TE)*LdotLminusl
         mykernel=GradFilter(ell,lens_index_TE)&
              &*Lensfilterl2(lens_index_E)*flL(lens_index_TE)*cos2phi

         If(HuOkamoto2002) then 
            flL(lens_index_TE) =clunlensed(ell,lens_index_TE)*Ldotl*cos2phi&
                 &+cl2unlensed(lens_index_TE)*LdotLminusl
            flL(lens_index_ET)=clunlensed(ell,lens_index_TE)*Ldotl&
                 &+cl2unlensed(lens_index_TE)*LdotLminusl*cos2phi

            mykernel=cltot(ell,lens_index_EE)*cl2tot(lens_index_TT)*flL(lens_index_TE)&
                 &-cltot(ell,lens_index_TE)*cl2tot(lens_index_TE)*flL(lens_index_ET)
            mykernel=mykernel&
                 &/(cltot(ell,lens_index_EE)*cl2tot(lens_index_TT)&
                 &*cltot(ell,lens_index_TT)*cl2tot(lens_index_EE)-&
                 &(cltot(ell,lens_index_TE)*cl2tot(lens_index_TE))**2)  !eq 13
            mykernel=mykernel*flL(lens_index_TE)
           
         end if

         
      case(lens_index_TB)
         flL(lens_index_TB)=clunlensed(ell,lens_index_TE)*Ldotl*sin2phi
         mykernel=GradFilter(ell,lens_index_TB)&
              &*Lensfilterl2(lens_index_B)*flL(lens_index_TB)*sin2phi

         If(HuOkamoto2002) then 
            mykernel=flL(lens_index_TB)**2&
                 &/(cltot(ell,lens_index_TT)*cl2tot(lens_index_BB)) !eq 15.
         end If
         
      case(lens_index_ET)
         flL(lens_index_ET)=clunlensed(ell,lens_index_TE)*Ldotl&
              &+cl2unlensed(lens_index_TE)*LdotLminusl*cos2phi
         mykernel=GradFilter(ell,lens_index_ET)&
              &*Lensfilterl2(lens_index_T)*flL(lens_index_ET) !correction
         !on Hu, Dedeo, Vale
         If(HuOkamoto2002) then 
            flL(lens_index_TE) =clunlensed(ell,lens_index_TE)*Ldotl*cos2phi&
                 &+cl2unlensed(lens_index_TE)*LdotLminusl
            flL(lens_index_ET)=clunlensed(ell,lens_index_TE)*Ldotl&
                 &+cl2unlensed(lens_index_TE)*LdotLminusl*cos2phi

            mykernel=cltot(ell,lens_index_EE)*cl2tot(lens_index_TT)*flL(lens_index_TE)&
                 &-cltot(ell,lens_index_TE)*cl2tot(lens_index_TE)*flL(lens_index_ET)
            mykernel=mykernel&
                 &/(cltot(ell,lens_index_EE)*cl2tot(lens_index_TT)&
                 &*cltot(ell,lens_index_TT)*cl2tot(lens_index_EE)-&
                 &(cltot(ell,lens_index_TE)*cl2tot(lens_index_TE))**2)  !eq 13
            mykernel=mykernel*flL(lens_index_TE)
         end if
         
         
      case(lens_index_EE)
         flL(lens_index_EE)=(clunlensed(ell,lens_index_EE)*Ldotl&
              &+cl2unlensed(lens_index_EE)*LdotLminusl)*cos2phi
         mykernel=GradFilter(ell,lens_index_EE)&
              &*Lensfilterl2(lens_index_E)*flL(lens_index_EE)*cos2phi

         
         If(HuOkamoto2002) then 
            mykernel=flL(lens_index_EE)**2&
                 &/(2.*cltot(ell,lens_index_EE)*cl2tot(lens_index_EE)) !eq 14.
         end If
         

      case(lens_index_EB)
         flL(lens_index_EB)=clunlensed(ell,lens_index_EE)*Ldotl*sin2phi
         mykernel=GradFilter(ell,lens_index_EB)&
              &*Lensfilterl2(lens_index_B)*flL(lens_index_EB)*sin2phi
         
         If(HuOkamoto2002) then 
            mykernel=flL(lens_index_EB)**2&
                 &/(cltot(ell,lens_index_EE)*cl2tot(lens_index_BB)) !eq 15.
         end If
         

      case(lens_index_BB)
         STOP 'cl_bb_unlensed is zero'
         mykernel=1e30 !some huge number \hat\kappa_BB is hopeless
      end select


      noisekernel=mykernel*ldotL*ell/(2*Pi)**2
      
      If(HuOkamoto2002) noisekernel=mykernel*ell/(2*Pi)**2
      
      !the following are some debug lines useful in comparing the
      !kernel with the IDL implementation noisekernel.pro

!!$      if(polnow.eq.lens_index_ET .and. ell .eq. 2 .and. bigL .eq. 3) then 
!!$         write(20,'(5e20.8)') theta, noisekernel,flL(polnow),gradfilter(ell,polnow),lensfilterl2(lens_index_E)
!!$         !print*,modLminusl
!!$      end if

      !the following would be the kernel from Hu 2001
      !noisekernel=flL(lens_index_TT)**2*ell/(2.*cltot(ell,lens_index_TT)*cltot(modLminusl,lens_index_TT))/(2*Pi)**2
      


      return
    end function noisekernel

    !-------- NR spline routines-------------
    subroutine spline(x,y,n,yp1,ypn,y2)
      integer n,NMAX
      real yp1,ypn
      real x(n), y(n), y2(n)
      parameter (NMAX=10000)	! Increased from 500 by Max
      integer i,k
      real p,qn,sig,un,u(NMAX)
      if (N.gt.NMAX) pause 'SPLINE NMAX DEATH ERROR' ! Added by Max
      !if (x(1).gt.x(n)) pause 'SPLINE WARNING: x NOT INCREASING' ! Added by Max
      if (yp1.gt..99e30) then
         y2(1)=0.
         u(1)=0.
      else
         y2(1)=-0.5
         u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do i=2,n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2.
         y2(i)=(sig-1.)/p
         u(i)=(6.*((y(i+1)-y(i))/(x(i+ &
              1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig &
              *u(i-1))/p
      end do
      if (ypn.gt..99e30) then
         qn=0.
         un=0.
      else
         qn=0.5
         un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do k=n-1,1,-1
         y2(k)=y2(k)*y2(k+1)+u(k)
      end do

    end subroutine spline


    subroutine splint(xa,ya,y2a,n,x,y)
      integer n
      real x,y
      real xa(n),y2a(n),ya(n)
      integer k,khi,klo
      real a,b,h

      !-------- If out of bounds, then linearly extrapolate
      ! else use the cubic spline
      if (x .lt. xa(1)) then
         y = ya(1) + (ya(1)-ya(2))*(x-xa(1))/(xa(1)-xa(2)) 
      else if (x .gt. xa(n)) then
         y = ya(n) + (ya(n)-ya(n-1))*(x-xa(n))/(xa(n)-xa(n-1))
      else
         klo=1
         khi=n
         do while (khi-klo.gt.1)
            k=(khi+klo)/2
            if(xa(k).gt.x)then
               khi=k
            else
               klo=k
            endif
         end do
         h=xa(khi)-xa(klo)
         if (h.eq.0.) pause 'bad xa input in splint'
         a=(xa(khi)-x)/h
         b=(x-xa(klo))/h
         y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h** &
              2)/6.
      end if

    end subroutine splint






  end subroutine GetLensingNoise

  subroutine qromb(func,a,b,ss)
    use precision
    implicit none

    integer,parameter:: JMAX=20,JMAXP=JMAX+1, K=5, KM=K-1
    real(dl):: a,b,ss,func
    real(dl),parameter::EPS=1.d-4
    !external func
    !   USES polint,trapzd
    integer::j

    real(dl):: dss,h(JMAXP),s(JMAXP)


    h(1)=1.d0
    do j=1,JMAX
       call trapzd(func,a,b,s(j),j)
       if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.d0,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
       endif
       s(j+1)=s(j)
       h(j+1)=0.25d0*h(j)
    enddo
    stop 'too many steps in qromb'


  end subroutine qromb


  subroutine trapzd(func,a,b,s,n)
    integer:: n
    real(dl):: a,b,s,func
    !external func
    integer:: it,j
    real(dl):: del,sum,tnm,x
    if (n.eq.1) then
       s=0.5d0*(b-a)*(func(a)+func(b))
    else
       it=2**(n-2)
       tnm=it
       del=(b-a)/tnm
       x=a+0.5d0*del
       sum=0.d0
       do j=1,it
          sum=sum+func(x)
          x=x+del
       enddo
       s=0.5d0*(s+(b-a)*sum/tnm)
    endif
    return
  end subroutine trapzd

  subroutine polint(xa,ya,n,x,y,dy)
    integer:: n
    real(dl)::dy,x,y,xa(n),ya(n)
    integer,parameter:: NMAX=10
    integer:: i,m,ns
    real(dl):: den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
    ns=1
    dif=abs(x-xa(1))
    do  i=1,n
       dift=abs(x-xa(i))
       if (dift.lt.dif) then
          ns=i
          dif=dift
       endif
       c(i)=ya(i)
       d(i)=ya(i)
    end do
    y=ya(ns)
    ns=ns-1
    do m=1,n-1
       do i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.d0)pause 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
       end do
       if (2*ns.lt.n-m)then
          dy=c(ns+1)
       else
          dy=d(ns)
          ns=ns-1
       endif
       y=y+dy
    enddo
    return
  end subroutine polint

  !a trial function to test rombint
  function trial(x)
    real(dl)::trial, x
    trial=x
    return
  end function trial

end module LensingNoise
    
  
     

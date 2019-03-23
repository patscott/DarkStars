! Extended precicison numerical integrator for DarkStars
!
! Pat Scott, Feb 2008; pat@physto.se
!--------------------------------------------------------------------------------     

      MODULE DkStrs_fint_ext

      use DkStrs_data
      use DkStrs_utils
      !This module has been disabled in the public release as it uses proprietry code.
      !use DkStrs_Romberg

      implicit none

      contains


      REAL(16) FUNCTION fint_ext(f,a,b,eps)
c_______________________________________________________________________
c  integrate function f between a and b
c  input
c    integration limits a and b
c  called by different routines
c  author: joakim edsjo (edsjo@physto.se) 96-05-16
c          2000-07-19 paolo gondolo added eps as argument
c          2007-07-27 pat scott added j > 6 condition to prevent
c           spurious early convergence, removing os initial definition.
c          2008-02-09 pat scott converted to work in and return extended
c           precision, and added switch out to Romberg integrator
c          2008-06-04 pat scott added switch out to Runge-Kutta integrator
c  based on paolo gondolos wxint.f routine.
c=======================================================================
      implicit none
      real*16 f,a,b,tot,eps,st,os,ost,del,summ,x,ss
      integer jmax,it,l,j,nfcn,jdid
      external f
      
      select case (integrator)
      case(1)
        !fint_ext = rkint(f,a,b,eps)
        write(*,*) 'Runge-Kutta integration has been disabled in the public release'
        write(*,*) 'because our implementation uses proprietary code.  Please plug'
        write(*,*) 'your favourite rk45 code into DkStrs_fint_ext if you would'
        write(*,*) 'like to use this option.'
        call die_quietly('DarkStars execution canceled.')
        return
      case(2)
        !call qromb16(f,a,b,eps,ss)
        !fint_ext = ss
        write(*,*) 'Romberg integration has been disabled in the public release'
        write(*,*) 'because our implementation uses proprietary code.  Please plug'
        write(*,*) 'your favourite qromb code into DkStrs_fint_ext if you would'
        write(*,*) 'like to use this option.'
        call die_quietly('DarkStars execution canceled.')
        return
      end select

      if (idtag.eq.'SPintegral') then
	    jmax=20 
	  else
	    jmax=25
      endif
	  fint_ext=0._16
      del=b-a
      ost=0.5_16*del*(f(a)+f(b))
      x=0.5_16*(b+a)
      st=0.5_16*(ost+del*f(x))
      ost=st
      it=1
      nfcn=3
      do j=3,jmax
        it=2*it
        del=0.5_16*del
        x=a+0.5_16*del
        summ=0._16
        do l=1,it
          summ=summ+f(x)
          nfcn=nfcn+1
          x=x+del
        enddo
        st=0.5_16*(st+del*summ)
        tot=(4._16*st-ost)/3._16
        jdid=j
        if (j.gt.6) then
           if (abs(tot-os).le.eps*abs(os)) then
              fint_ext=tot
              return
           endif
     	endif
        os=tot
        ost=st
      enddo
      if (idtag .ne. 'SPintegral' .and. idtag .ne. 'n_WIMPs_cent') then
        write(*,*) 'error in fint_ext: too many steps for model: ', idtag
        write(*,*) '  integral set to zero.'
      endif
      fint_ext=0._16

      END FUNCTION fint_ext

      
      REAL*16 FUNCTION rkint(f,a,b,eps)
      
      real*16 :: f,a,b,eps
      real*16 :: y(1), t, tout, relerr, abserr, emergency_abs_error, work(9)
      integer :: neqn, iflag, iwork(5)
      
      neqn = 1
      y(1) = 0._16
      t = a
      tout = b
      abserr = 0._16
      emergency_abs_error = 1.e-50_16
      relerr = eps
      iflag = 1      

      do while (iflag .lt. 6)
        !Uncomment this line if you want to plug in an rk45 integrator.
        !call rkf45(f,neqn,y,t,tout,relerr,abserr,iflag,work,iwork)
        select case (iflag)
        case(-2)
          iflag = 2
          abserr = 0._16
        case(2)
          rkint = y(1)
          return
        case(5)
          iflag = -2
          abserr = emergency_abs_error
        end select
      enddo

      write(*,*) 'Exit code: ', iflag
      call die_quietly('Error: Runge-Kutta integration failed in'//trim(idtag))

      END FUNCTION rkint


      END MODULE DkStrs_fint_ext

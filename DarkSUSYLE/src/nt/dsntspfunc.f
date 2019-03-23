
*************************
      real*8 function dsntspfunc(r)
      implicit none

      real*8 r,dsntsunmass,pi,gn
      parameter(pi=3.141592653589793234d0)
      parameter(gn=6.67259d-11) ! m^3 kg^-1 s^-1

      dsntspfunc=dsntsunmass(r)*gn/max(r,100.0d0)**2
c      write(*,*) 'dsntspfunc: ',r,dsntspfunc

      return
      end
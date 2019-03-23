
*************************
      real*8 function dscapspfunc(r)
      implicit none
      include 'dscapstar.h'

      real*8 r,dscapstarmass

      dscapspfunc=dscapstarmass(r)*gn/max(r,100.0d0)**2
c      write(*,*) 'dsntspfunc: ',r,dsntspfunc

      return
      end

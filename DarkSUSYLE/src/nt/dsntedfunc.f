
*************************
      real*8 function dsntedfunc(r)
      implicit none

      real*8 r,dsntearthdens,pi
      parameter(pi=3.141592653589793234d0)

      dsntedfunc=dsntearthdens(r)*1000.0d0*4.0d0*pi*r**2

      return
      end

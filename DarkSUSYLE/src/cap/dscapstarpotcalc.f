      subroutine dscapstarpotcalc

***********************************************************************
*** Calculate star's gravitational potential from a given
*** density distribution common block (as described in dscapstar.h).
*** Author: Joakim Edsjo
*** Date: 2003-11-25
*** Modified: 2007-06-26
*** Modified: 2007-07-04 Pat Scott (pat@physto.se) 
***********************************************************************

      implicit none
      include 'dsparam.h'
      include 'dscapstar.h'

      real*8 dscapstarpotint
      integer i
         
c...Calculate and tabulate the potential inside the star

      !write(*,*) 
     &!  'Tabulating potential inside the star...'
      do i=1,meshpoints
        starphi(i)=dscapstarpotint(starr(i)*r_star)
      enddo
      !write(*,*) '  ...done'

      return

      end


      


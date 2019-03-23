

      real*8 function dsai(x)
c   dsairy function
c   lb 990224
      implicit none
      real*8 x,pi,z,k13,dbskr3
      if (x.lt.0.d0) then
         write(*,*) 'dsairy fcn of negative argument - not implemented'
         stop
         else
           pi=3.141592653589793d0
           z=2.d0/3.d0*x**(1.5d0)
           k13=dbskr3(z,1)  ! cernlib function for k{1/3}
           dsai=1.d0/pi*sqrt(x/3.d0)*k13  ! num recipes (6.7.41)
         endif
      return
      end



      real*8 function dsbi(x)
c   dsairy function
c   lb 990224
      implicit none
      real*8 x,pi,z,i13,k13,dbsir3,dbskr3
      if (x.lt.0.d0) then
         write(*,*) 'dsairy fcn of negative argument - not implemented'
         stop
         else
           pi=3.141592653589793d0
           z=2.d0/3.d0*x**(1.5d0)
           i13=dbsir3(z,1)  ! cernlib function for i{1/3}
           k13=dbskr3(z,1)  ! cernlib function for k{1/3}
           dsbi=sqrt(x)*(2.d0*i13/sqrt(3.d0)+k13/pi)  ! num recipes (6.7.44)
         endif
      return
      end

      function erf(x)
      implicit none
      real*8 x,erfc,erf,t,z
      z = abs(x)
      t = 1.0d0/(1.0+0.5d0*z)
      erfc = t * exp(-z*z-1.26551223d0+t*(1.00002368d0+t*(0.37409196+
     &     t*(0.09678418d0+t*(-0.18628806d0+t*(0.27886807d0+t*
     &     (-1.13520398d0+t*(1.48851587d0+t*(-0.82215223d0+t*
     &     0.17087277d0)))))))))
      if (x.lt.0.d0) erfc = 2.0d0-erfc
      erf=1.0d0-erfc
      return
      end

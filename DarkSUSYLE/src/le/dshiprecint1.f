
      subroutine dshiprecint1(fun,foveru,lowlim,upplim,result)
      implicit none
      real*8 lowlim,upplim
      real*8 sla1,sla2,abserr,alist,blist,
     &  elist,epsabs,epsrel,result,rlist,foveru
      integer ier,iord,last,limit,neval
      dimension alist(1000),blist(1000),elist(1000),iord(1000),
     & rlist(1000)
      external fun,foveru
      epsabs=1.d-17
      epsrel=1.d-17
      limit=1000
      sla1=lowlim
      sla2=upplim
c compute the integral with the numerical integration
        call dsntdqagsec(fun,foveru,sla1,sla2,epsabs,epsrel,limit,
     &   result,
     &   abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      return
      end

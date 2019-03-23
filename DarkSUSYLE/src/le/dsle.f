      subroutine dsLEinit
      implicit none
      include 'dssusy.h'
      include 'dsversion.h'  ! set DarkSUSY version
      include 'dssubversion.h'
      include 'dsdir.h'      ! set DarkSUSY root directory

      integer i

c...Startup

      write(*,*)
      write(*,*) '    DarkSUSY found at '//dsinstall
      write(*,*)
     &  '    *********************************************************'
      write(*,*) 
     &  '    *** Welcome to DarkSUSY version                       ***'
      write(*,*) '    *** ',dsversion,'***'
      write(*,*) '    *** ',dssubversion,'***'
      write(*,*) 
     &  '    *********************************************************'
      write(*,*)
      write(*,*) '    Initializing DarkSUSYLE capture routines...'      
      
      call dscapsetup

      write(*,*) '    done.'
      write(*,*)

      return
      end


      real*8 function dsLEfint(f,a,b,eps)
c_______________________________________________________________________
c  integrate function f between a and b
c  input
c    integration limits a and b
c  called by different routines
c  author: joakim edsjo (edsjo@physto.se) 96-05-16
c          2000-07-19 paolo gondolo added eps as argument
c          2007-07-27 pat scott added j > 6 condition to prevent
c           spurious early convergence, removing os initial definition.
c          2008-03-28 pat scott added switch out to Romberg integrator
c  based on paolo gondolos wxint.f routine.
c=======================================================================
      implicit none
      include 'dsidtag.h'
      include 'dscapstar.h'
      real*8 f,a,b,tot,eps,st,os,ost,del,sum,x,ss
      integer jmax,it,l,j,nfcn,jdid
      external f
      parameter (jmax=25)
      
      if (integrator .eq. 2) then
        call qromb(f,a,b,eps,ss)
        dsLEfint = ss
        return
      endif      

      dsLEfint=0.d0
      del=b-a
      ost=0.5*del*(f(a)+f(b))
      x=0.5*(b+a)
      st=0.5*(ost+del*f(x))
      ost=st
      it=1
      nfcn=3
      do j=3,jmax
        it=2*it
        del=0.5*del
        x=a+0.5*del
        sum=0.0
        do l=1,it
          sum=sum+f(x)
          nfcn=nfcn+1
          x=x+del
        enddo
        st=0.5*(st+del*sum)
        tot=(4.0*st-ost)/3.0
        jdid=j
        if (j.gt.6) then
           if (abs(tot-os).le.eps*abs(os)) then
              dsLEfint=tot
              return
           endif
     	endif
        os=tot
        ost=st
      enddo
      write(*,*) 'error in dsLEfint: too many steps for model: ',
     &  idtag
      write(*,*) '  integral set to zero.'
      dsLEfint=0.0d0
      end

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
* ======================================================================
* nist guide to available math software.
* fullsource for module dqagse from package cmlib.
* retrieved from camsun on wed oct  8 08:26:30 1997.
* ======================================================================
      subroutine dsntdqagsec(f,foveru,
     &   a,b,epsabs,epsrel,limit,result,abserr,neval,
     1   ier,alist,blist,rlist,elist,iord,last)
c***begin prologue  dqagse
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a1
c***keywords  (end point) singularities,automatic integrator,
c             extrapolation,general-purpose,globally adaptive
c***author  piessens, robert, applied math. and progr. div. -
c             k. u. leuven
c           de doncker, elise, applied math. and progr. div. -
c             k. u. leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral i = integral of f over (a,b),
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.max(epsabs,epsrel*abs(i)).
c***description
c
c        computation of a definite integral
c        standard fortran subroutine
c        real*8 version
c
c        parameters
c         on entry
c            f      - real*8
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - real*8
c                     lower limit of integration
c
c            b      - real*8
c                     upper limit of integration
c
c            epsabs - real*8
c                     absolute accuracy requested
c            epsrel - real*8
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c            limit  - integer
c                     gives an upperbound on the number of subintervals
c                     in the partition of (a,b)
c
c         on return
c            result - real*8
c                     approximation to the integral
c
c            abserr - real*8
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine
c                             the estimates for integral and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                         = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more sub-
c                             divisions by increasing the value of limit
c                             (and taking the according dimension
c                             adjustments into account). however, if
c                             this yields no improvement it is advised
c                             to analyze the integrand in order to
c                             determine the integration difficulties. if
c                             the position of a local difficulty can be
c                             determined (e.g. singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling the
c                             integrator on the subranges. if possible,
c                             an appropriate special-purpose integrator
c                             should be used, which is designed for
c                             handling the type of difficulty involved.
c                         = 2 the occurrence of roundoff error is detec-
c                             ted, which prevents the requested
c                             tolerance from being achieved.
c                             the error may be under-estimated.
c                         = 3 extremely bad integrand behaviour
c                             occurs at some points of the integration
c                             interval.
c                         = 4 the algorithm does not converge.
c                             roundoff error is detected in the
c                             extrapolation table.
c                             it is presumed that the requested
c                             tolerance cannot be achieved, and that the
c                             returned result is the best which can be
c                             obtained.
c                         = 5 the integral is probably divergent, or
c                             slowly convergent. it must be noted that
c                             divergence can occur with any other value
c                             of ier.
c                         = 6 the input is invalid, because
c                             epsabs.le.0 and
c                             epsrel.lt.max(50*rel.mach.acc.,0.5d-28).
c                             result, abserr, neval, last, rlist(1),
c                             iord(1) and elist(1) are set to zero.
c                             alist(1) and blist(1) are set to a and b
c                             respectively.
c
c            alist  - real*8
c                     vector of dimension at least limit, the first
c                      last  elements of which are the left end points
c                     of the subintervals in the partition of the
c                     given integration range (a,b)
c
c            blist  - real*8
c                     vector of dimension at least limit, the first
c                      last  elements of which are the right end points
c                     of the subintervals in the partition of the given
c                     integration range (a,b)
c
c            rlist  - real*8
c                     vector of dimension at least limit, the first
c                      last  elements of which are the integral
c                     approximations on the subintervals
c
c            elist  - real*8
c                     vector of dimension at least limit, the first
c                      last  elements of which are the moduli of the
c                     absolute error estimates on the subintervals
c
c            iord   - integer
c                     vector of dimension at least limit, the first k
c                     elements of which are pointers to the
c                     error estimates over the subintervals,
c                     such that elist(iord(1)), ..., elist(iord(k))
c                     form a decreasing sequence, with k = last
c                     if last.le.(limit/2+2), and k = limit+1-last
c                     otherwise
c
c            last   - integer
c                     number of subintervals actually produced in the
c                     subdivision process
c***references  (none)
c***routines called  d1mach,dqelg,dqk21,dqpsrt
c***end prologue  dqagse
c
      real*8 a,abseps,abserr,alist,area,area1,area12,area2,a1,
     1  a2,b,blist,b1,b2,correc,dabs,defabs,defab1,defab2,d1mach,dmax1,
     2  dres,elist,epmach,epsabs,epsrel,erlarg,erlast,errbnd,errmax,
     3  error1,error2,erro12,errsum,ertest,f,oflow,resabs,reseps,result,
     4  res3la,rlist,rlist2,small,uflow,foveru
      integer id,ier,ierro,iord,iroff1,iroff2,iroff3,jupbnd,k,ksgn,
     1  ktmin,last,limit,maxerr,neval,nres,nrmax,numrl2
      logical extrap,noext
c
      dimension alist(limit),blist(limit),elist(limit),iord(limit),
     1 res3la(3),rlist(limit),rlist2(52)
c
      external f,foveru
c
c            the dimension of rlist2 is determined by the value of
c            limexp in subroutine dqelg (rlist2 should be of dimension
c            (limexp+2) at least).
c
c            list of major variables
c            -----------------------
c
c           alist     - list of left end points of all subintervals
c                       considered up to now
c           blist     - list of right end points of all subintervals
c                       considered up to now
c           rlist(i)  - approximation to the integral over
c                       (alist(i),blist(i))
c           rlist2    - array of dimension at least limexp+2 containing
c                       the part of the epsilon table which is still
c                       needed for further computations
c           elist(i)  - error estimate applying to rlist(i)
c           maxerr    - pointer to the interval with largest error
c                       estimate
c           errmax    - elist(maxerr)
c           erlast    - error on the interval currently subdivided
c                       (before that subdivision has taken place)
c           area      - sum of the integrals over the subintervals
c           errsum    - sum of the errors over the subintervals
c           errbnd    - requested accuracy max(epsabs,epsrel*
c                       abs(result))
c           *****1    - variable for the left interval
c           *****2    - variable for the right interval
c           last      - index for subdivision
c           nres      - number of calls to the extrapolation routine
c           numrl2    - number of elements currently in rlist2. if an
c                       appropriate approximation to the compounded
c                       integral has been obtained it is put in
c                       rlist2(numrl2) after numrl2 has been increased
c                       by one.
c           small     - length of the smallest interval considered up
c                       to now, multiplied by 1.5
c           erlarg    - sum of the errors over the intervals larger
c                       than the smallest interval considered up to now
c           extrap    - logical variable denoting that the routine is
c                       attempting to perform extrapolation i.e. before
c                       subdividing the smallest interval we try to
c                       decrease the value of erlarg.
c           noext     - logical variable denoting that extrapolation
c                       is no longer allowed (true value)
c
c            machine dependent constants
c            ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c           oflow is the largest positive magnitude.
c
c***first executable statement  dqagse
      epmach = d1mach(4)
c
c            test on validity of parameters
c            ------------------------------
      ier = 0
      neval = 0
      last = 0
      result = 0.0d+00
      abserr = 0.0d+00
      alist(1) = a
      blist(1) = b
      rlist(1) = 0.0d+00
      elist(1) = 0.0d+00
      if(epsabs.le.0.0d+00.and.epsrel.lt.dmax1(0.5d+02*epmach,0.5d-28))
     1   ier = 6
      if(ier.eq.6) go to 999
c
c           first approximation to the integral
c           -----------------------------------
c
      uflow = d1mach(1)
      oflow = d1mach(2)
      ierro = 0
      call dsntdqk21c(f,foveru,a,b,result,abserr,defabs,resabs)
c
c           test on accuracy.
c
      dres = dabs(result)
      errbnd = dmax1(epsabs,epsrel*dres)
      last = 1
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1
      if(abserr.le.1.0d+02*epmach*defabs.and.abserr.gt.errbnd) ier = 2
      if(limit.eq.1) ier = 1
      if(ier.ne.0.or.(abserr.le.errbnd.and.abserr.ne.resabs).or.
     1  abserr.eq.0.0d+00) go to 140
c
c           initialization
c           --------------
c
      rlist2(1) = result
      errmax = abserr
      maxerr = 1
      area = result
      errsum = abserr
      abserr = oflow
      nrmax = 1
      nres = 0
      numrl2 = 2
      ktmin = 0
      extrap = .false.
      noext = .false.
      iroff1 = 0
      iroff2 = 0
      iroff3 = 0
      ksgn = -1
      if(dres.ge.(0.1d+01-0.5d+02*epmach)*defabs) ksgn = 1
c
c           main do-loop
c           ------------
c
      do 90 last = 2,limit
c
c           bisect the subinterval with the nrmax-th largest error
c           estimate.
c
        a1 = alist(maxerr)
        b1 = 0.5d+00*(alist(maxerr)+blist(maxerr))
        a2 = b1
        b2 = blist(maxerr)
        erlast = errmax
        call dsntdqk21c(f,foveru,a1,b1,area1,error1,resabs,defab1)
        call dsntdqk21c(f,foveru,a2,b2,area2,error2,resabs,defab2)
c
c           improve previous approximations to integral
c           and error and test for accuracy.
c
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(defab1.eq.error1.or.defab2.eq.error2) go to 15
        if(dabs(rlist(maxerr)-area12).gt.0.1d-04*dabs(area12)
     1  .or.erro12.lt.0.99d+00*errmax) go to 10
        if(extrap) iroff2 = iroff2+1
        if(.not.extrap) iroff1 = iroff1+1
   10   if(last.gt.10.and.erro12.gt.errmax) iroff3 = iroff3+1
   15   rlist(maxerr) = area1
        rlist(last) = area2
        errbnd = dmax1(epsabs,epsrel*dabs(area))
c
c           test for roundoff error and eventually set error flag.
c
        if(iroff1+iroff2.ge.10.or.iroff3.ge.20) ier = 2
        if(iroff2.ge.5) ierro = 3
c
c           set error flag in the case that the number of subintervals
c           equals limit.
c
        if(last.eq.limit) ier = 1
c
c           set error flag in the case of bad integrand behaviour
c           at a point of the integration range.
c
        if(dmax1(dabs(a1),dabs(b2)).le.(0.1d+01+0.1d+03*epmach)*
     1  (dabs(a2)+0.1d+04*uflow)) ier = 4
c
c           append the newly-created intervals to the list.
c
        if(error2.gt.error1) go to 20
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 30
   20   alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
c
c           call subroutine dqpsrt to maintain the descending ordering
c           in the list of error estimates and select the subinterval
c           with nrmax-th largest error estimate (to be bisected next).
c
   30   call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
c ***jump out of do-loop
        if(errsum.le.errbnd) go to 115
c ***jump out of do-loop
        if(ier.ne.0) go to 100
        if(last.eq.2) go to 80
        if(noext) go to 90
        erlarg = erlarg-erlast
        if(dabs(b1-a1).gt.small) erlarg = erlarg+erro12
        if(extrap) go to 40
c
c           test whether the interval to be bisected next is the
c           smallest interval.
c
        if(dabs(blist(maxerr)-alist(maxerr)).gt.small) go to 90
        extrap = .true.
        nrmax = 2
   40   if(ierro.eq.3.or.erlarg.le.ertest) go to 60
c
c           the smallest interval has the largest error.
c           before bisecting decrease the sum of the errors over the
c           larger intervals (erlarg) and perform extrapolation.
c
        id = nrmax
        jupbnd = last
        if(last.gt.(2+limit/2)) jupbnd = limit+3-last
        do 50 k = id,jupbnd
          maxerr = iord(nrmax)
          errmax = elist(maxerr)
c ***jump out of do-loop
          if(dabs(blist(maxerr)-alist(maxerr)).gt.small) go to 90
          nrmax = nrmax+1
   50   continue
c
c           perform extrapolation.
c
   60   numrl2 = numrl2+1
        rlist2(numrl2) = area
        call dqelg(numrl2,rlist2,reseps,abseps,res3la,nres)
        ktmin = ktmin+1
        if(ktmin.gt.5.and.abserr.lt.0.1d-02*errsum) ier = 5
        if(abseps.ge.abserr) go to 70
        ktmin = 0
        abserr = abseps
        result = reseps
        correc = erlarg
        ertest = dmax1(epsabs,epsrel*dabs(reseps))
c ***jump out of do-loop
        if(abserr.le.ertest) go to 100
c
c           prepare bisection of the smallest interval.
c
   70   if(numrl2.eq.1) noext = .true.
        if(ier.eq.5) go to 100
        maxerr = iord(1)
        errmax = elist(maxerr)
        nrmax = 1
        extrap = .false.
        small = small*0.5d+00
        erlarg = errsum
        go to 90
   80   small = dabs(b-a)*0.375d+00
        erlarg = errsum
        ertest = errbnd
        rlist2(2) = area
   90 continue
c
c           set final result and error estimate.
c           ------------------------------------
c
  100 if(abserr.eq.oflow) go to 115
      if(ier+ierro.eq.0) go to 110
      if(ierro.eq.3) abserr = abserr+correc
      if(ier.eq.0) ier = 3
      if(result.ne.0.0d+00.and.area.ne.0.0d+00) go to 105
      if(abserr.gt.errsum) go to 115
      if(area.eq.0.0d+00) go to 130
      go to 110
  105 if(abserr/dabs(result).gt.errsum/dabs(area)) go to 115
c
c           test on divergence.
c
  110 if(ksgn.eq.(-1).and.dmax1(dabs(result),dabs(area)).le.
     1 defabs*0.1d-01) go to 130
      if(0.1d-01.gt.(result/area).or.(result/area).gt.0.1d+03
     1 .or.errsum.gt.dabs(area)) ier = 6
      go to 130
c
c           compute global integral sum.
c
  115 result = 0.0d+00
      do 120 k = 1,last
         result = result+rlist(k)
  120 continue
      abserr = errsum
  130 if(ier.gt.2) ier = ier-1
  140 neval = 42*last-21
  999 return
      end
      subroutine dsntdqk21c(f,foveru,a,b,result,abserr,resabs,resasc)
c***begin prologue  dqk21
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a2
c***keywords  21-point gauss-kronrod rules
c***author  piessens, robert, applied math. and progr. div. -
c             k. u. leuven
c           de doncker, elise, applied math. and progr. div. -
c             k. u. leuven
c***purpose  to compute i = integral of f over (a,b), with error
c                           estimate
c                       j = integral of abs(f) over (a,b)
c***description
c
c           integration rules
c           standard fortran subroutine
c           real*8 version
c
c           parameters
c            on entry
c              f      - real*8
c                       function subprogram defining the integrand
c                       function f(x). the actual name for f needs to be
c                       declared e x t e r n a l in the driver program.
c
c              a      - real*8
c                       lower limit of integration
c
c              b      - real*8
c                       upper limit of integration
c
c            on return
c              result - real*8
c                       approximation to the integral i
c                       result is computed by applying the 21-point
c                       kronrod rule (resk) obtained by optimal addition
c                       of abscissae to the 10-point gauss rule (resg).
c
c              abserr - real*8
c                       estimate of the modulus of the absolute error,
c                       which should not exceed abs(i-result)
c
c              resabs - real*8
c                       approximation to the integral j
c
c              resasc - real*8
c                       approximation to the integral of abs(f-i/(b-a))
c                       over (a,b)
c***references  (none)
c***routines called  d1mach
c***end prologue  dqk21
c
      real*8 a,absc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1,
     1  d1mach,epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc,
     2  resg,resk,reskh,result,uflow,wg,wgk,xgk,foveru
      integer j,jtw,jtwm1
      external f,foveru
c
      dimension fv1(10),fv2(10),wg(5),wgk(11),xgk(11)
c
c           the abscissae and weights are given for the interval (-1,1).
c           because of symmetry only the positive abscissae and their
c           corresponding weights are given.
c
c           xgk    - abscissae of the 21-point kronrod rule
c                    xgk(2), xgk(4), ...  abscissae of the 10-point
c                    gauss rule
c                    xgk(1), xgk(3), ...  abscissae which are optimally
c                    added to the 10-point gauss rule
c
c           wgk    - weights of the 21-point kronrod rule
c
c           wg     - weights of the 10-point gauss rule
c
c
c gauss quadrature weights and kronron quadrature abscissae and weights
c as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
c bell labs, nov. 1981.
c
      data wg  (  1) / 0.0666713443 0868813759 3568809893 332 d0 /
      data wg  (  2) / 0.1494513491 5058059314 5776339657 697 d0 /
      data wg  (  3) / 0.2190863625 1598204399 5534934228 163 d0 /
      data wg  (  4) / 0.2692667193 0999635509 1226921569 469 d0 /
      data wg  (  5) / 0.2955242247 1475287017 3892994651 338 d0 /
c
      data xgk (  1) / 0.9956571630 2580808073 5527280689 003 d0 /
      data xgk (  2) / 0.9739065285 1717172007 7964012084 452 d0 /
      data xgk (  3) / 0.9301574913 5570822600 1207180059 508 d0 /
      data xgk (  4) / 0.8650633666 8898451073 2096688423 493 d0 /
      data xgk (  5) / 0.7808177265 8641689706 3717578345 042 d0 /
      data xgk (  6) / 0.6794095682 9902440623 4327365114 874 d0 /
      data xgk (  7) / 0.5627571346 6860468333 9000099272 694 d0 /
      data xgk (  8) / 0.4333953941 2924719079 9265943165 784 d0 /
      data xgk (  9) / 0.2943928627 0146019813 1126603103 866 d0 /
      data xgk ( 10) / 0.1488743389 8163121088 4826001129 720 d0 /
      data xgk ( 11) / 0.0000000000 0000000000 0000000000 000 d0 /
c
      data wgk (  1) / 0.0116946388 6737187427 8064396062 192 d0 /
      data wgk (  2) / 0.0325581623 0796472747 8818972459 390 d0 /
      data wgk (  3) / 0.0547558965 7435199603 1381300244 580 d0 /
      data wgk (  4) / 0.0750396748 1091995276 7043140916 190 d0 /
      data wgk (  5) / 0.0931254545 8369760553 5065465083 366 d0 /
      data wgk (  6) / 0.1093871588 0229764189 9210590325 805 d0 /
      data wgk (  7) / 0.1234919762 6206585107 7958109831 074 d0 /
      data wgk (  8) / 0.1347092173 1147332592 8054001771 707 d0 /
      data wgk (  9) / 0.1427759385 7706008079 7094273138 717 d0 /
      data wgk ( 10) / 0.1477391049 0133849137 4841515972 068 d0 /
      data wgk ( 11) / 0.1494455540 0291690566 4936468389 821 d0 /
c
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           absc   - abscissa
c           fval*  - function value
c           resg   - result of the 10-point gauss formula
c           resk   - result of the 21-point kronrod formula
c           reskh  - approximation to the mean value of f over (a,b),
c                    i.e. to i/(b-a)
c
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  dqk21
      epmach = d1mach(4)
      uflow = d1mach(1)
c
      centr = 0.5d+00*(a+b)
      hlgth = 0.5d+00*(b-a)
      dhlgth = dabs(hlgth)
c
c           compute the 21-point kronrod approximation to
c           the integral, and estimate the absolute error.
c
      resg = 0.0d+00
      fc = f(centr,foveru)
      resk = wgk(11)*fc
      resabs = dabs(resk)
      do 10 j=1,5
        jtw = 2*j
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc,foveru)
        fval2 = f(centr+absc,foveru)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
   10 continue
      do 15 j = 1,5
        jtwm1 = 2*j-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc,foveru)
        fval2 = f(centr+absc,foveru)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
   15 continue
      reskh = resk*0.5d+00
      resasc = wgk(11)*dabs(fc-reskh)
      do 20 j=1,10
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = dabs((resk-resg)*hlgth)
      if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00)
     1  abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1
     1  ((epmach*0.5d+02)*resabs,abserr)
      return
      end
      !Pat's modified version of NR Ed 2 qromb - takes eps as input and is in double precision

      SUBROUTINE qromb(func,a,b,eps,ss)
      INTEGER JMAX,JMAXP,K,KM
      REAL*8 a,b,func,ss,eps
      EXTERNAL func
      PARAMETER (JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES polint,trapzd
      INTEGER j
      REAL*8 dss,h(JMAXP),s(JMAXP)
      h(1)=1.d0
      do 11 j=1,JMAX
        call trapzd(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.d0,ss,dss)
          if (abs(dss).le.eps*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25d0*h(j)
11    continue
      write(*,*) 'Too many steps in qromb for model: ',
     &  idtag
      write(*,*) '  integral set to zero.'
      ss=0.d0
      END
      !Pat's modified version of NR Ed 2 qromb - double precision

      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL*8 dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.d0)pause 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
      !Pat's modified version of NR Ed 2 qromb - double precision

      SUBROUTINE trapzd(func,a,b,s,n)
      INTEGER n
      REAL*8 a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL*8 del,sum,tnm,x
      if (n.eq.1) then
        s=0.5d0*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5d0*del
        sum=0.d0
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5d0*(s+(b-a)*sum/tnm)
      endif
      return
      END

      subroutine dshiprecint(fun,foveru,lowlim,upplim,result)
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
        call dsntdqagse(fun,foveru,sla1,sla2,epsabs,epsrel,limit,
     &   result,
     &   abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      return
      end

      subroutine dshiprecint2(fun,foveru,lowlim,upplim,result)
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
        call dsntdqagseb(fun,foveru,
     &   sla1,sla2,epsabs,epsrel,limit,result,
     &   abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      return
      end
* ======================================================================
* nist guide to available math software.
* fullsource for module dqagse from package cmlib.
* retrieved from camsun on wed oct  8 08:26:30 1997.
* ======================================================================
      subroutine dsntdqagse(f,foveru,
     &   a,b,epsabs,epsrel,limit,result,abserr,neval,
     1   ier,alist,blist,rlist,elist,iord,last)
c***begin prologue  dqagse
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a1
c***keywords  (end point) singularities,automatic integrator,
c             extrapolation,general-purpose,globally adaptive
c***author  piessens, robert, applied math. and progr. div. -
c             k. u. leuven
c           de doncker, elise, applied math. and progr. div. -
c             k. u. leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral i = integral of f over (a,b),
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.max(epsabs,epsrel*abs(i)).
c***description
c
c        computation of a definite integral
c        standard fortran subroutine
c        real*8 version
c
c        parameters
c         on entry
c            f      - real*8
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - real*8
c                     lower limit of integration
c
c            b      - real*8
c                     upper limit of integration
c
c            epsabs - real*8
c                     absolute accuracy requested
c            epsrel - real*8
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c            limit  - integer
c                     gives an upperbound on the number of subintervals
c                     in the partition of (a,b)
c
c         on return
c            result - real*8
c                     approximation to the integral
c
c            abserr - real*8
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine
c                             the estimates for integral and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                         = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more sub-
c                             divisions by increasing the value of limit
c                             (and taking the according dimension
c                             adjustments into account). however, if
c                             this yields no improvement it is advised
c                             to analyze the integrand in order to
c                             determine the integration difficulties. if
c                             the position of a local difficulty can be
c                             determined (e.g. singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling the
c                             integrator on the subranges. if possible,
c                             an appropriate special-purpose integrator
c                             should be used, which is designed for
c                             handling the type of difficulty involved.
c                         = 2 the occurrence of roundoff error is detec-
c                             ted, which prevents the requested
c                             tolerance from being achieved.
c                             the error may be under-estimated.
c                         = 3 extremely bad integrand behaviour
c                             occurs at some points of the integration
c                             interval.
c                         = 4 the algorithm does not converge.
c                             roundoff error is detected in the
c                             extrapolation table.
c                             it is presumed that the requested
c                             tolerance cannot be achieved, and that the
c                             returned result is the best which can be
c                             obtained.
c                         = 5 the integral is probably divergent, or
c                             slowly convergent. it must be noted that
c                             divergence can occur with any other value
c                             of ier.
c                         = 6 the input is invalid, because
c                             epsabs.le.0 and
c                             epsrel.lt.max(50*rel.mach.acc.,0.5d-28).
c                             result, abserr, neval, last, rlist(1),
c                             iord(1) and elist(1) are set to zero.
c                             alist(1) and blist(1) are set to a and b
c                             respectively.
c
c            alist  - real*8
c                     vector of dimension at least limit, the first
c                      last  elements of which are the left end points
c                     of the subintervals in the partition of the
c                     given integration range (a,b)
c
c            blist  - real*8
c                     vector of dimension at least limit, the first
c                      last  elements of which are the right end points
c                     of the subintervals in the partition of the given
c                     integration range (a,b)
c
c            rlist  - real*8
c                     vector of dimension at least limit, the first
c                      last  elements of which are the integral
c                     approximations on the subintervals
c
c            elist  - real*8
c                     vector of dimension at least limit, the first
c                      last  elements of which are the moduli of the
c                     absolute error estimates on the subintervals
c
c            iord   - integer
c                     vector of dimension at least limit, the first k
c                     elements of which are pointers to the
c                     error estimates over the subintervals,
c                     such that elist(iord(1)), ..., elist(iord(k))
c                     form a decreasing sequence, with k = last
c                     if last.le.(limit/2+2), and k = limit+1-last
c                     otherwise
c
c            last   - integer
c                     number of subintervals actually produced in the
c                     subdivision process
c***references  (none)
c***routines called  d1mach,dqelg,dqk21,dqpsrt
c***end prologue  dqagse
c
      real*8 a,abseps,abserr,alist,area,area1,area12,area2,a1,
     1  a2,b,blist,b1,b2,correc,dabs,defabs,defab1,defab2,d1mach,dmax1,
     2  dres,elist,epmach,epsabs,epsrel,erlarg,erlast,errbnd,errmax,
     3  error1,error2,erro12,errsum,ertest,f,oflow,resabs,reseps,result,
     4  res3la,rlist,rlist2,small,uflow,foveru
      integer id,ier,ierro,iord,iroff1,iroff2,iroff3,jupbnd,k,ksgn,
     1  ktmin,last,limit,maxerr,neval,nres,nrmax,numrl2
      logical extrap,noext
c
      dimension alist(limit),blist(limit),elist(limit),iord(limit),
     1 res3la(3),rlist(limit),rlist2(52)
c
      external f,foveru
c
c            the dimension of rlist2 is determined by the value of
c            limexp in subroutine dqelg (rlist2 should be of dimension
c            (limexp+2) at least).
c
c            list of major variables
c            -----------------------
c
c           alist     - list of left end points of all subintervals
c                       considered up to now
c           blist     - list of right end points of all subintervals
c                       considered up to now
c           rlist(i)  - approximation to the integral over
c                       (alist(i),blist(i))
c           rlist2    - array of dimension at least limexp+2 containing
c                       the part of the epsilon table which is still
c                       needed for further computations
c           elist(i)  - error estimate applying to rlist(i)
c           maxerr    - pointer to the interval with largest error
c                       estimate
c           errmax    - elist(maxerr)
c           erlast    - error on the interval currently subdivided
c                       (before that subdivision has taken place)
c           area      - sum of the integrals over the subintervals
c           errsum    - sum of the errors over the subintervals
c           errbnd    - requested accuracy max(epsabs,epsrel*
c                       abs(result))
c           *****1    - variable for the left interval
c           *****2    - variable for the right interval
c           last      - index for subdivision
c           nres      - number of calls to the extrapolation routine
c           numrl2    - number of elements currently in rlist2. if an
c                       appropriate approximation to the compounded
c                       integral has been obtained it is put in
c                       rlist2(numrl2) after numrl2 has been increased
c                       by one.
c           small     - length of the smallest interval considered up
c                       to now, multiplied by 1.5
c           erlarg    - sum of the errors over the intervals larger
c                       than the smallest interval considered up to now
c           extrap    - logical variable denoting that the routine is
c                       attempting to perform extrapolation i.e. before
c                       subdividing the smallest interval we try to
c                       decrease the value of erlarg.
c           noext     - logical variable denoting that extrapolation
c                       is no longer allowed (true value)
c
c            machine dependent constants
c            ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c           oflow is the largest positive magnitude.
c
c***first executable statement  dqagse
      epmach = d1mach(4)
c
c            test on validity of parameters
c            ------------------------------
      ier = 0
      neval = 0
      last = 0
      result = 0.0d+00
      abserr = 0.0d+00
      alist(1) = a
      blist(1) = b
      rlist(1) = 0.0d+00
      elist(1) = 0.0d+00
      if(epsabs.le.0.0d+00.and.epsrel.lt.dmax1(0.5d+02*epmach,0.5d-28))
     1   ier = 6
      if(ier.eq.6) go to 999
c
c           first approximation to the integral
c           -----------------------------------
c
      uflow = d1mach(1)
      oflow = d1mach(2)
      ierro = 0
      call dsntdqk21(f,foveru,a,b,result,abserr,defabs,resabs)
c
c           test on accuracy.
c
      dres = dabs(result)
      errbnd = dmax1(epsabs,epsrel*dres)
      last = 1
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1
      if(abserr.le.1.0d+02*epmach*defabs.and.abserr.gt.errbnd) ier = 2
      if(limit.eq.1) ier = 1
      if(ier.ne.0.or.(abserr.le.errbnd.and.abserr.ne.resabs).or.
     1  abserr.eq.0.0d+00) go to 140
c
c           initialization
c           --------------
c
      rlist2(1) = result
      errmax = abserr
      maxerr = 1
      area = result
      errsum = abserr
      abserr = oflow
      nrmax = 1
      nres = 0
      numrl2 = 2
      ktmin = 0
      extrap = .false.
      noext = .false.
      iroff1 = 0
      iroff2 = 0
      iroff3 = 0
      ksgn = -1
      if(dres.ge.(0.1d+01-0.5d+02*epmach)*defabs) ksgn = 1
c
c           main do-loop
c           ------------
c
      do 90 last = 2,limit
c
c           bisect the subinterval with the nrmax-th largest error
c           estimate.
c
        a1 = alist(maxerr)
        b1 = 0.5d+00*(alist(maxerr)+blist(maxerr))
        a2 = b1
        b2 = blist(maxerr)
        erlast = errmax
        call dsntdqk21(f,foveru,a1,b1,area1,error1,resabs,defab1)
        call dsntdqk21(f,foveru,a2,b2,area2,error2,resabs,defab2)
c
c           improve previous approximations to integral
c           and error and test for accuracy.
c
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(defab1.eq.error1.or.defab2.eq.error2) go to 15
        if(dabs(rlist(maxerr)-area12).gt.0.1d-04*dabs(area12)
     1  .or.erro12.lt.0.99d+00*errmax) go to 10
        if(extrap) iroff2 = iroff2+1
        if(.not.extrap) iroff1 = iroff1+1
   10   if(last.gt.10.and.erro12.gt.errmax) iroff3 = iroff3+1
   15   rlist(maxerr) = area1
        rlist(last) = area2
        errbnd = dmax1(epsabs,epsrel*dabs(area))
c
c           test for roundoff error and eventually set error flag.
c
        if(iroff1+iroff2.ge.10.or.iroff3.ge.20) ier = 2
        if(iroff2.ge.5) ierro = 3
c
c           set error flag in the case that the number of subintervals
c           equals limit.
c
        if(last.eq.limit) ier = 1
c
c           set error flag in the case of bad integrand behaviour
c           at a point of the integration range.
c
        if(dmax1(dabs(a1),dabs(b2)).le.(0.1d+01+0.1d+03*epmach)*
     1  (dabs(a2)+0.1d+04*uflow)) ier = 4
c
c           append the newly-created intervals to the list.
c
        if(error2.gt.error1) go to 20
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 30
   20   alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
c
c           call subroutine dqpsrt to maintain the descending ordering
c           in the list of error estimates and select the subinterval
c           with nrmax-th largest error estimate (to be bisected next).
c
   30   call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
c ***jump out of do-loop
        if(errsum.le.errbnd) go to 115
c ***jump out of do-loop
        if(ier.ne.0) go to 100
        if(last.eq.2) go to 80
        if(noext) go to 90
        erlarg = erlarg-erlast
        if(dabs(b1-a1).gt.small) erlarg = erlarg+erro12
        if(extrap) go to 40
c
c           test whether the interval to be bisected next is the
c           smallest interval.
c
        if(dabs(blist(maxerr)-alist(maxerr)).gt.small) go to 90
        extrap = .true.
        nrmax = 2
   40   if(ierro.eq.3.or.erlarg.le.ertest) go to 60
c
c           the smallest interval has the largest error.
c           before bisecting decrease the sum of the errors over the
c           larger intervals (erlarg) and perform extrapolation.
c
        id = nrmax
        jupbnd = last
        if(last.gt.(2+limit/2)) jupbnd = limit+3-last
        do 50 k = id,jupbnd
          maxerr = iord(nrmax)
          errmax = elist(maxerr)
c ***jump out of do-loop
          if(dabs(blist(maxerr)-alist(maxerr)).gt.small) go to 90
          nrmax = nrmax+1
   50   continue
c
c           perform extrapolation.
c
   60   numrl2 = numrl2+1
        rlist2(numrl2) = area
        call dqelg(numrl2,rlist2,reseps,abseps,res3la,nres)
        ktmin = ktmin+1
        if(ktmin.gt.5.and.abserr.lt.0.1d-02*errsum) ier = 5
        if(abseps.ge.abserr) go to 70
        ktmin = 0
        abserr = abseps
        result = reseps
        correc = erlarg
        ertest = dmax1(epsabs,epsrel*dabs(reseps))
c ***jump out of do-loop
        if(abserr.le.ertest) go to 100
c
c           prepare bisection of the smallest interval.
c
   70   if(numrl2.eq.1) noext = .true.
        if(ier.eq.5) go to 100
        maxerr = iord(1)
        errmax = elist(maxerr)
        nrmax = 1
        extrap = .false.
        small = small*0.5d+00
        erlarg = errsum
        go to 90
   80   small = dabs(b-a)*0.375d+00
        erlarg = errsum
        ertest = errbnd
        rlist2(2) = area
   90 continue
c
c           set final result and error estimate.
c           ------------------------------------
c
  100 if(abserr.eq.oflow) go to 115
      if(ier+ierro.eq.0) go to 110
      if(ierro.eq.3) abserr = abserr+correc
      if(ier.eq.0) ier = 3
      if(result.ne.0.0d+00.and.area.ne.0.0d+00) go to 105
      if(abserr.gt.errsum) go to 115
      if(area.eq.0.0d+00) go to 130
      go to 110
  105 if(abserr/dabs(result).gt.errsum/dabs(area)) go to 115
c
c           test on divergence.
c
  110 if(ksgn.eq.(-1).and.dmax1(dabs(result),dabs(area)).le.
     1 defabs*0.1d-01) go to 130
      if(0.1d-01.gt.(result/area).or.(result/area).gt.0.1d+03
     1 .or.errsum.gt.dabs(area)) ier = 6
      go to 130
c
c           compute global integral sum.
c
  115 result = 0.0d+00
      do 120 k = 1,last
         result = result+rlist(k)
  120 continue
      abserr = errsum
  130 if(ier.gt.2) ier = ier-1
  140 neval = 42*last-21
  999 return
      end
* ======================================================================
* nist guide to available math software.
* fullsource for module dqagse from package cmlib.
* retrieved from camsun on wed oct  8 08:26:30 1997.
* ======================================================================
      subroutine dsntdqagseb(f,foveru,
     &   a,b,epsabs,epsrel,limit,result,abserr,neval,
     1   ier,alist,blist,rlist,elist,iord,last)
c***begin prologue  dqagse
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a1
c***keywords  (end point) singularities,automatic integrator,
c             extrapolation,general-purpose,globally adaptive
c***author  piessens, robert, applied math. and progr. div. -
c             k. u. leuven
c           de doncker, elise, applied math. and progr. div. -
c             k. u. leuven
c***purpose  the routine calculates an approximation result to a given
c            definite integral i = integral of f over (a,b),
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.max(epsabs,epsrel*abs(i)).
c***description
c
c        computation of a definite integral
c        standard fortran subroutine
c        real*8 version
c
c        parameters
c         on entry
c            f      - real*8
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - real*8
c                     lower limit of integration
c
c            b      - real*8
c                     upper limit of integration
c
c            epsabs - real*8
c                     absolute accuracy requested
c            epsrel - real*8
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c            limit  - integer
c                     gives an upperbound on the number of subintervals
c                     in the partition of (a,b)
c
c         on return
c            result - real*8
c                     approximation to the integral
c
c            abserr - real*8
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine
c                             the estimates for integral and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                         = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more sub-
c                             divisions by increasing the value of limit
c                             (and taking the according dimension
c                             adjustments into account). however, if
c                             this yields no improvement it is advised
c                             to analyze the integrand in order to
c                             determine the integration difficulties. if
c                             the position of a local difficulty can be
c                             determined (e.g. singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling the
c                             integrator on the subranges. if possible,
c                             an appropriate special-purpose integrator
c                             should be used, which is designed for
c                             handling the type of difficulty involved.
c                         = 2 the occurrence of roundoff error is detec-
c                             ted, which prevents the requested
c                             tolerance from being achieved.
c                             the error may be under-estimated.
c                         = 3 extremely bad integrand behaviour
c                             occurs at some points of the integration
c                             interval.
c                         = 4 the algorithm does not converge.
c                             roundoff error is detected in the
c                             extrapolation table.
c                             it is presumed that the requested
c                             tolerance cannot be achieved, and that the
c                             returned result is the best which can be
c                             obtained.
c                         = 5 the integral is probably divergent, or
c                             slowly convergent. it must be noted that
c                             divergence can occur with any other value
c                             of ier.
c                         = 6 the input is invalid, because
c                             epsabs.le.0 and
c                             epsrel.lt.max(50*rel.mach.acc.,0.5d-28).
c                             result, abserr, neval, last, rlist(1),
c                             iord(1) and elist(1) are set to zero.
c                             alist(1) and blist(1) are set to a and b
c                             respectively.
c
c            alist  - real*8
c                     vector of dimension at least limit, the first
c                      last  elements of which are the left end points
c                     of the subintervals in the partition of the
c                     given integration range (a,b)
c
c            blist  - real*8
c                     vector of dimension at least limit, the first
c                      last  elements of which are the right end points
c                     of the subintervals in the partition of the given
c                     integration range (a,b)
c
c            rlist  - real*8
c                     vector of dimension at least limit, the first
c                      last  elements of which are the integral
c                     approximations on the subintervals
c
c            elist  - real*8
c                     vector of dimension at least limit, the first
c                      last  elements of which are the moduli of the
c                     absolute error estimates on the subintervals
c
c            iord   - integer
c                     vector of dimension at least limit, the first k
c                     elements of which are pointers to the
c                     error estimates over the subintervals,
c                     such that elist(iord(1)), ..., elist(iord(k))
c                     form a decreasing sequence, with k = last
c                     if last.le.(limit/2+2), and k = limit+1-last
c                     otherwise
c
c            last   - integer
c                     number of subintervals actually produced in the
c                     subdivision process
c***references  (none)
c***routines called  d1mach,dqelg,dqk21,dqpsrt
c***end prologue  dqagse
c
      real*8 a,abseps,abserr,alist,area,area1,area12,area2,a1,
     1  a2,b,blist,b1,b2,correc,dabs,defabs,defab1,defab2,d1mach,dmax1,
     2  dres,elist,epmach,epsabs,epsrel,erlarg,erlast,errbnd,errmax,
     3  error1,error2,erro12,errsum,ertest,f,oflow,resabs,reseps,result,
     4  res3la,rlist,rlist2,small,uflow,foveru
      integer id,ier,ierro,iord,iroff1,iroff2,iroff3,jupbnd,k,ksgn,
     1  ktmin,last,limit,maxerr,neval,nres,nrmax,numrl2
      logical extrap,noext
c
      dimension alist(limit),blist(limit),elist(limit),iord(limit),
     1 res3la(3),rlist(limit),rlist2(52)
c
      external f,foveru
c
c            the dimension of rlist2 is determined by the value of
c            limexp in subroutine dqelg (rlist2 should be of dimension
c            (limexp+2) at least).
c
c            list of major variables
c            -----------------------
c
c           alist     - list of left end points of all subintervals
c                       considered up to now
c           blist     - list of right end points of all subintervals
c                       considered up to now
c           rlist(i)  - approximation to the integral over
c                       (alist(i),blist(i))
c           rlist2    - array of dimension at least limexp+2 containing
c                       the part of the epsilon table which is still
c                       needed for further computations
c           elist(i)  - error estimate applying to rlist(i)
c           maxerr    - pointer to the interval with largest error
c                       estimate
c           errmax    - elist(maxerr)
c           erlast    - error on the interval currently subdivided
c                       (before that subdivision has taken place)
c           area      - sum of the integrals over the subintervals
c           errsum    - sum of the errors over the subintervals
c           errbnd    - requested accuracy max(epsabs,epsrel*
c                       abs(result))
c           *****1    - variable for the left interval
c           *****2    - variable for the right interval
c           last      - index for subdivision
c           nres      - number of calls to the extrapolation routine
c           numrl2    - number of elements currently in rlist2. if an
c                       appropriate approximation to the compounded
c                       integral has been obtained it is put in
c                       rlist2(numrl2) after numrl2 has been increased
c                       by one.
c           small     - length of the smallest interval considered up
c                       to now, multiplied by 1.5
c           erlarg    - sum of the errors over the intervals larger
c                       than the smallest interval considered up to now
c           extrap    - logical variable denoting that the routine is
c                       attempting to perform extrapolation i.e. before
c                       subdividing the smallest interval we try to
c                       decrease the value of erlarg.
c           noext     - logical variable denoting that extrapolation
c                       is no longer allowed (true value)
c
c            machine dependent constants
c            ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c           oflow is the largest positive magnitude.
c
c***first executable statement  dqagse
      epmach = d1mach(4)
c
c            test on validity of parameters
c            ------------------------------
      ier = 0
      neval = 0
      last = 0
      result = 0.0d+00
      abserr = 0.0d+00
      alist(1) = a
      blist(1) = b
      rlist(1) = 0.0d+00
      elist(1) = 0.0d+00
      if(epsabs.le.0.0d+00.and.epsrel.lt.dmax1(0.5d+02*epmach,0.5d-28))
     1   ier = 6
      if(ier.eq.6) go to 999
c
c           first approximation to the integral
c           -----------------------------------
c
      uflow = d1mach(1)
      oflow = d1mach(2)
      ierro = 0
      call dsntdqk21b(f,foveru,a,b,result,abserr,defabs,resabs)
c
c           test on accuracy.
c
      dres = dabs(result)
      errbnd = dmax1(epsabs,epsrel*dres)
      last = 1
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1
      if(abserr.le.1.0d+02*epmach*defabs.and.abserr.gt.errbnd) ier = 2
      if(limit.eq.1) ier = 1
      if(ier.ne.0.or.(abserr.le.errbnd.and.abserr.ne.resabs).or.
     1  abserr.eq.0.0d+00) go to 140
c
c           initialization
c           --------------
c
      rlist2(1) = result
      errmax = abserr
      maxerr = 1
      area = result
      errsum = abserr
      abserr = oflow
      nrmax = 1
      nres = 0
      numrl2 = 2
      ktmin = 0
      extrap = .false.
      noext = .false.
      iroff1 = 0
      iroff2 = 0
      iroff3 = 0
      ksgn = -1
      if(dres.ge.(0.1d+01-0.5d+02*epmach)*defabs) ksgn = 1
c
c           main do-loop
c           ------------
c
      do 90 last = 2,limit
c
c           bisect the subinterval with the nrmax-th largest error
c           estimate.
c
        a1 = alist(maxerr)
        b1 = 0.5d+00*(alist(maxerr)+blist(maxerr))
        a2 = b1
        b2 = blist(maxerr)
        erlast = errmax
        call dsntdqk21b(f,foveru,a1,b1,area1,error1,resabs,defab1)
        call dsntdqk21b(f,foveru,a2,b2,area2,error2,resabs,defab2)
c
c           improve previous approximations to integral
c           and error and test for accuracy.
c
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(defab1.eq.error1.or.defab2.eq.error2) go to 15
        if(dabs(rlist(maxerr)-area12).gt.0.1d-04*dabs(area12)
     1  .or.erro12.lt.0.99d+00*errmax) go to 10
        if(extrap) iroff2 = iroff2+1
        if(.not.extrap) iroff1 = iroff1+1
   10   if(last.gt.10.and.erro12.gt.errmax) iroff3 = iroff3+1
   15   rlist(maxerr) = area1
        rlist(last) = area2
        errbnd = dmax1(epsabs,epsrel*dabs(area))
c
c           test for roundoff error and eventually set error flag.
c
        if(iroff1+iroff2.ge.10.or.iroff3.ge.20) ier = 2
        if(iroff2.ge.5) ierro = 3
c
c           set error flag in the case that the number of subintervals
c           equals limit.
c
        if(last.eq.limit) ier = 1
c
c           set error flag in the case of bad integrand behaviour
c           at a point of the integration range.
c
        if(dmax1(dabs(a1),dabs(b2)).le.(0.1d+01+0.1d+03*epmach)*
     1  (dabs(a2)+0.1d+04*uflow)) ier = 4
c
c           append the newly-created intervals to the list.
c
        if(error2.gt.error1) go to 20
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 30
   20   alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
c
c           call subroutine dqpsrt to maintain the descending ordering
c           in the list of error estimates and select the subinterval
c           with nrmax-th largest error estimate (to be bisected next).
c
   30   call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
c ***jump out of do-loop
        if(errsum.le.errbnd) go to 115
c ***jump out of do-loop
        if(ier.ne.0) go to 100
        if(last.eq.2) go to 80
        if(noext) go to 90
        erlarg = erlarg-erlast
        if(dabs(b1-a1).gt.small) erlarg = erlarg+erro12
        if(extrap) go to 40
c
c           test whether the interval to be bisected next is the
c           smallest interval.
c
        if(dabs(blist(maxerr)-alist(maxerr)).gt.small) go to 90
        extrap = .true.
        nrmax = 2
   40   if(ierro.eq.3.or.erlarg.le.ertest) go to 60
c
c           the smallest interval has the largest error.
c           before bisecting decrease the sum of the errors over the
c           larger intervals (erlarg) and perform extrapolation.
c
        id = nrmax
        jupbnd = last
        if(last.gt.(2+limit/2)) jupbnd = limit+3-last
        do 50 k = id,jupbnd
          maxerr = iord(nrmax)
          errmax = elist(maxerr)
c ***jump out of do-loop
          if(dabs(blist(maxerr)-alist(maxerr)).gt.small) go to 90
          nrmax = nrmax+1
   50   continue
c
c           perform extrapolation.
c
   60   numrl2 = numrl2+1
        rlist2(numrl2) = area
        call dqelg(numrl2,rlist2,reseps,abseps,res3la,nres)
        ktmin = ktmin+1
        if(ktmin.gt.5.and.abserr.lt.0.1d-02*errsum) ier = 5
        if(abseps.ge.abserr) go to 70
        ktmin = 0
        abserr = abseps
        result = reseps
        correc = erlarg
        ertest = dmax1(epsabs,epsrel*dabs(reseps))
c ***jump out of do-loop
        if(abserr.le.ertest) go to 100
c
c           prepare bisection of the smallest interval.
c
   70   if(numrl2.eq.1) noext = .true.
        if(ier.eq.5) go to 100
        maxerr = iord(1)
        errmax = elist(maxerr)
        nrmax = 1
        extrap = .false.
        small = small*0.5d+00
        erlarg = errsum
        go to 90
   80   small = dabs(b-a)*0.375d+00
        erlarg = errsum
        ertest = errbnd
        rlist2(2) = area
   90 continue
c
c           set final result and error estimate.
c           ------------------------------------
c
  100 if(abserr.eq.oflow) go to 115
      if(ier+ierro.eq.0) go to 110
      if(ierro.eq.3) abserr = abserr+correc
      if(ier.eq.0) ier = 3
      if(result.ne.0.0d+00.and.area.ne.0.0d+00) go to 105
      if(abserr.gt.errsum) go to 115
      if(area.eq.0.0d+00) go to 130
      go to 110
  105 if(abserr/dabs(result).gt.errsum/dabs(area)) go to 115
c
c           test on divergence.
c
  110 if(ksgn.eq.(-1).and.dmax1(dabs(result),dabs(area)).le.
     1 defabs*0.1d-01) go to 130
      if(0.1d-01.gt.(result/area).or.(result/area).gt.0.1d+03
     1 .or.errsum.gt.dabs(area)) ier = 6
      go to 130
c
c           compute global integral sum.
c
  115 result = 0.0d+00
      do 120 k = 1,last
         result = result+rlist(k)
  120 continue
      abserr = errsum
  130 if(ier.gt.2) ier = ier-1
  140 neval = 42*last-21
  999 return
      end
      subroutine dsntdqk21(f,foveru,a,b,result,abserr,resabs,resasc)
c***begin prologue  dqk21
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a2
c***keywords  21-point gauss-kronrod rules
c***author  piessens, robert, applied math. and progr. div. -
c             k. u. leuven
c           de doncker, elise, applied math. and progr. div. -
c             k. u. leuven
c***purpose  to compute i = integral of f over (a,b), with error
c                           estimate
c                       j = integral of abs(f) over (a,b)
c***description
c
c           integration rules
c           standard fortran subroutine
c           real*8 version
c
c           parameters
c            on entry
c              f      - real*8
c                       function subprogram defining the integrand
c                       function f(x). the actual name for f needs to be
c                       declared e x t e r n a l in the driver program.
c
c              a      - real*8
c                       lower limit of integration
c
c              b      - real*8
c                       upper limit of integration
c
c            on return
c              result - real*8
c                       approximation to the integral i
c                       result is computed by applying the 21-point
c                       kronrod rule (resk) obtained by optimal addition
c                       of abscissae to the 10-point gauss rule (resg).
c
c              abserr - real*8
c                       estimate of the modulus of the absolute error,
c                       which should not exceed abs(i-result)
c
c              resabs - real*8
c                       approximation to the integral j
c
c              resasc - real*8
c                       approximation to the integral of abs(f-i/(b-a))
c                       over (a,b)
c***references  (none)
c***routines called  d1mach
c***end prologue  dqk21
c
      real*8 a,absc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1,
     1  d1mach,epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc,
     2  resg,resk,reskh,result,uflow,wg,wgk,xgk,foveru
      integer j,jtw,jtwm1
      external f,foveru
c
      dimension fv1(10),fv2(10),wg(5),wgk(11),xgk(11)
c
c           the abscissae and weights are given for the interval (-1,1).
c           because of symmetry only the positive abscissae and their
c           corresponding weights are given.
c
c           xgk    - abscissae of the 21-point kronrod rule
c                    xgk(2), xgk(4), ...  abscissae of the 10-point
c                    gauss rule
c                    xgk(1), xgk(3), ...  abscissae which are optimally
c                    added to the 10-point gauss rule
c
c           wgk    - weights of the 21-point kronrod rule
c
c           wg     - weights of the 10-point gauss rule
c
c
c gauss quadrature weights and kronron quadrature abscissae and weights
c as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
c bell labs, nov. 1981.
c
      data wg  (  1) / 0.0666713443 0868813759 3568809893 332 d0 /
      data wg  (  2) / 0.1494513491 5058059314 5776339657 697 d0 /
      data wg  (  3) / 0.2190863625 1598204399 5534934228 163 d0 /
      data wg  (  4) / 0.2692667193 0999635509 1226921569 469 d0 /
      data wg  (  5) / 0.2955242247 1475287017 3892994651 338 d0 /
c
      data xgk (  1) / 0.9956571630 2580808073 5527280689 003 d0 /
      data xgk (  2) / 0.9739065285 1717172007 7964012084 452 d0 /
      data xgk (  3) / 0.9301574913 5570822600 1207180059 508 d0 /
      data xgk (  4) / 0.8650633666 8898451073 2096688423 493 d0 /
      data xgk (  5) / 0.7808177265 8641689706 3717578345 042 d0 /
      data xgk (  6) / 0.6794095682 9902440623 4327365114 874 d0 /
      data xgk (  7) / 0.5627571346 6860468333 9000099272 694 d0 /
      data xgk (  8) / 0.4333953941 2924719079 9265943165 784 d0 /
      data xgk (  9) / 0.2943928627 0146019813 1126603103 866 d0 /
      data xgk ( 10) / 0.1488743389 8163121088 4826001129 720 d0 /
      data xgk ( 11) / 0.0000000000 0000000000 0000000000 000 d0 /
c
      data wgk (  1) / 0.0116946388 6737187427 8064396062 192 d0 /
      data wgk (  2) / 0.0325581623 0796472747 8818972459 390 d0 /
      data wgk (  3) / 0.0547558965 7435199603 1381300244 580 d0 /
      data wgk (  4) / 0.0750396748 1091995276 7043140916 190 d0 /
      data wgk (  5) / 0.0931254545 8369760553 5065465083 366 d0 /
      data wgk (  6) / 0.1093871588 0229764189 9210590325 805 d0 /
      data wgk (  7) / 0.1234919762 6206585107 7958109831 074 d0 /
      data wgk (  8) / 0.1347092173 1147332592 8054001771 707 d0 /
      data wgk (  9) / 0.1427759385 7706008079 7094273138 717 d0 /
      data wgk ( 10) / 0.1477391049 0133849137 4841515972 068 d0 /
      data wgk ( 11) / 0.1494455540 0291690566 4936468389 821 d0 /
c
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           absc   - abscissa
c           fval*  - function value
c           resg   - result of the 10-point gauss formula
c           resk   - result of the 21-point kronrod formula
c           reskh  - approximation to the mean value of f over (a,b),
c                    i.e. to i/(b-a)
c
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  dqk21
      epmach = d1mach(4)
      uflow = d1mach(1)
c
      centr = 0.5d+00*(a+b)
      hlgth = 0.5d+00*(b-a)
      dhlgth = dabs(hlgth)
c
c           compute the 21-point kronrod approximation to
c           the integral, and estimate the absolute error.
c
      resg = 0.0d+00
      fc = f(centr,foveru)
      resk = wgk(11)*fc
      resabs = dabs(resk)
      do 10 j=1,5
        jtw = 2*j
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc,foveru)
        fval2 = f(centr+absc,foveru)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
   10 continue
      do 15 j = 1,5
        jtwm1 = 2*j-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc,foveru)
        fval2 = f(centr+absc,foveru)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
   15 continue
      reskh = resk*0.5d+00
      resasc = wgk(11)*dabs(fc-reskh)
      do 20 j=1,10
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = dabs((resk-resg)*hlgth)
      if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00)
     1  abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1
     1  ((epmach*0.5d+02)*resabs,abserr)
      return
      end
      subroutine dsntdqk21b(f,foveru,a,b,result,abserr,resabs,resasc)
c***begin prologue  dqk21
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a1a2
c***keywords  21-point gauss-kronrod rules
c***author  piessens, robert, applied math. and progr. div. -
c             k. u. leuven
c           de doncker, elise, applied math. and progr. div. -
c             k. u. leuven
c***purpose  to compute i = integral of f over (a,b), with error
c                           estimate
c                       j = integral of abs(f) over (a,b)
c***description
c
c           integration rules
c           standard fortran subroutine
c           real*8 version
c
c           parameters
c            on entry
c              f      - real*8
c                       function subprogram defining the integrand
c                       function f(x). the actual name for f needs to be
c                       declared e x t e r n a l in the driver program.
c
c              a      - real*8
c                       lower limit of integration
c
c              b      - real*8
c                       upper limit of integration
c
c            on return
c              result - real*8
c                       approximation to the integral i
c                       result is computed by applying the 21-point
c                       kronrod rule (resk) obtained by optimal addition
c                       of abscissae to the 10-point gauss rule (resg).
c
c              abserr - real*8
c                       estimate of the modulus of the absolute error,
c                       which should not exceed abs(i-result)
c
c              resabs - real*8
c                       approximation to the integral j
c
c              resasc - real*8
c                       approximation to the integral of abs(f-i/(b-a))
c                       over (a,b)
c***references  (none)
c***routines called  d1mach
c***end prologue  dqk21
c
      real*8 a,absc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1,
     1  d1mach,epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc,
     2  resg,resk,reskh,result,uflow,wg,wgk,xgk,foveru
      integer j,jtw,jtwm1
      external f,foveru
c
      dimension fv1(10),fv2(10),wg(5),wgk(11),xgk(11)
c
c           the abscissae and weights are given for the interval (-1,1).
c           because of symmetry only the positive abscissae and their
c           corresponding weights are given.
c
c           xgk    - abscissae of the 21-point kronrod rule
c                    xgk(2), xgk(4), ...  abscissae of the 10-point
c                    gauss rule
c                    xgk(1), xgk(3), ...  abscissae which are optimally
c                    added to the 10-point gauss rule
c
c           wgk    - weights of the 21-point kronrod rule
c
c           wg     - weights of the 10-point gauss rule
c
c
c gauss quadrature weights and kronron quadrature abscissae and weights
c as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
c bell labs, nov. 1981.
c
      data wg  (  1) / 0.0666713443 0868813759 3568809893 332 d0 /
      data wg  (  2) / 0.1494513491 5058059314 5776339657 697 d0 /
      data wg  (  3) / 0.2190863625 1598204399 5534934228 163 d0 /
      data wg  (  4) / 0.2692667193 0999635509 1226921569 469 d0 /
      data wg  (  5) / 0.2955242247 1475287017 3892994651 338 d0 /
c
      data xgk (  1) / 0.9956571630 2580808073 5527280689 003 d0 /
      data xgk (  2) / 0.9739065285 1717172007 7964012084 452 d0 /
      data xgk (  3) / 0.9301574913 5570822600 1207180059 508 d0 /
      data xgk (  4) / 0.8650633666 8898451073 2096688423 493 d0 /
      data xgk (  5) / 0.7808177265 8641689706 3717578345 042 d0 /
      data xgk (  6) / 0.6794095682 9902440623 4327365114 874 d0 /
      data xgk (  7) / 0.5627571346 6860468333 9000099272 694 d0 /
      data xgk (  8) / 0.4333953941 2924719079 9265943165 784 d0 /
      data xgk (  9) / 0.2943928627 0146019813 1126603103 866 d0 /
      data xgk ( 10) / 0.1488743389 8163121088 4826001129 720 d0 /
      data xgk ( 11) / 0.0000000000 0000000000 0000000000 000 d0 /
c
      data wgk (  1) / 0.0116946388 6737187427 8064396062 192 d0 /
      data wgk (  2) / 0.0325581623 0796472747 8818972459 390 d0 /
      data wgk (  3) / 0.0547558965 7435199603 1381300244 580 d0 /
      data wgk (  4) / 0.0750396748 1091995276 7043140916 190 d0 /
      data wgk (  5) / 0.0931254545 8369760553 5065465083 366 d0 /
      data wgk (  6) / 0.1093871588 0229764189 9210590325 805 d0 /
      data wgk (  7) / 0.1234919762 6206585107 7958109831 074 d0 /
      data wgk (  8) / 0.1347092173 1147332592 8054001771 707 d0 /
      data wgk (  9) / 0.1427759385 7706008079 7094273138 717 d0 /
      data wgk ( 10) / 0.1477391049 0133849137 4841515972 068 d0 /
      data wgk ( 11) / 0.1494455540 0291690566 4936468389 821 d0 /
c
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           absc   - abscissa
c           fval*  - function value
c           resg   - result of the 10-point gauss formula
c           resk   - result of the 21-point kronrod formula
c           reskh  - approximation to the mean value of f over (a,b),
c                    i.e. to i/(b-a)
c
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  dqk21
      epmach = d1mach(4)
      uflow = d1mach(1)
c
      centr = 0.5d+00*(a+b)
      hlgth = 0.5d+00*(b-a)
      dhlgth = dabs(hlgth)
c
c           compute the 21-point kronrod approximation to
c           the integral, and estimate the absolute error.
c
      resg = 0.0d+00
      fc = f(centr,foveru)
      resk = wgk(11)*fc
      resabs = dabs(resk)
      do 10 j=1,5
        jtw = 2*j
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc,foveru)
        fval2 = f(centr+absc,foveru)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
   10 continue
      do 15 j = 1,5
        jtwm1 = 2*j-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc,foveru)
        fval2 = f(centr+absc,foveru)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
   15 continue
      reskh = resk*0.5d+00
      resasc = wgk(11)*dabs(fc-reskh)
      do 20 j=1,10
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = dabs((resk-resg)*hlgth)
      if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00)
     1  abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1
     1  ((epmach*0.5d+02)*resabs,abserr)
      return
      end
***********************************************************************
*** dsntsundens gives the density in the Sun as a function of radius
*** the radius should be given in m and the density is returned in
*** g/cm^3
*** Density and element mass fractions up to O16 are from the standard
*** solar model BP2000 of Bahcall, Pinsonneault and Basu,
*** ApJ 555 (2001) 990.
*** The mass fractions for heavier elements are from N. Grevesse and
*** A.J. Sauval, Space Science Reviews 85 (1998) 161 normalized such that
*** their total mass fractions matches that of the heavier elements in 
*** the BP2000 model.
***
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 2003-11-25
***********************************************************************

      real*8 function dsntsundens(r)
      implicit none

      include 'dssun.h'
      real*8 r,rpl
      integer i,j

c...Check if data file is loaded
      call dsntsunread


      if (r.ge.r_sun) then
        dsntsundens=0.0d0
        return
      endif

      if (r.eq.0.0d0) then
        dsntsundens=sdrho(1)
        return
      endif

      call dshunt(sdr,sdn,r/r_sun,j)
      if (j.lt.sdn) goto 20
      
      dsntsundens=0.0d0
      return

 20   rpl=(r-sdr(j)*r_sun)/(sdr(j+1)*r_sun-sdr(j)*r_sun)

      dsntsundens=sdrho(j)*(1.0d0-rpl)+sdrho(j+1)*rpl

      return

      end
***********************************************************************
*** dsntsunmfrac gives the mass fraction of element i (see dsntsunread.f
*** for definition of i) as a function of the solar radius r.
*** the radius should be given in m and returned is the mass fraction.
***
*** Element mass fractions up to O16 are from the standard
*** solar model BP2000 of Bahcall, Pinsonneault and Basu,
*** ApJ 555 (2001) 990.
*** The mass fractions for heavier elements are from N. Grevesse and
*** A.J. Sauval, Space Science Reviews 85 (1998) 161 normalized such that
*** their total mass fractions matches that of the heavier elements in 
*** the BP2000 model.
***
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 2003-11-26
***********************************************************************

      real*8 function dsntsunmfrac(r,itype)
      implicit none

      include 'dssun.h'
      real*8 r,rpl
      integer i,j,itype

c...Check if data file is loaded
      call dsntsunread

      if (r.gt.r_sun) then
        dsntsunmfrac=0.0d0
        return
      endif

      if (r.eq.0.0d0) then
        dsntsunmfrac=sdmfr(itype,1)
        return
      endif

      call dshunt(sdr,sdn,r/r_sun,j)
      if (j.lt.sdn) goto 20

      j=sdn-1

 20   rpl=(r-sdr(j)*r_sun)/(sdr(j+1)*r_sun-sdr(j)*r_sun)

      dsntsunmfrac=sdmfr(itype,j)*(1.0d0-rpl)
     &  +sdmfr(itype,j+1)*rpl

      return

      end
      subroutine dsntsunread

***********************************************************************
*** Reads in data about the solar model used and stores it in a
*** common block (as described in dssun.h).
*** Author: Joakim Edsjo
*** Date: 2003-11-25
*** Modified: 2004-01-28 (calculates potential instead of reading file)
***********************************************************************

      implicit none
      include 'dssun.h'
      include 'dssusy.h'
      include 'dsparam.h'

      real*8 dsntsunpotint,dsntsuncdensint
      integer i,fl,l,m

      logical sdread
      data sdread/.false./
      save sdread

      character*200 file,filene,scr

      real*8 tmp1,tmp2,tmp3,totfr,totfrheavy

c...If already initialized, return
      if (sdread) then
        return
      endif

      sdread=.true.

c...Determine which file we should read in
c...This is from the Standard solar model, BP2000
c      file=dsinstall//'share/DarkSUSY/'//'bp2000stdmodel.dat'
c...This is from the Standard solar model, BS05(OP)
      file=dsinstall//'share/DarkSUSY/'//'bs05op.dat'
c...This is form the alternative Standard solar model, BS05(AGS,OP)
c...with new heavy element measurements (fits worse with helioseismology
c...so we don't use it as a defualt).
c      file=dsinstall//'share/DarkSUSY/'//'bs05_agsop.dat'
c...The electron density is from the Standard solar model, BS05(OP)
      filene=dsinstall//'share/DarkSUSY/'//'nele_bs05op.dat'

c...delete possible spaces in file name
      fl=200
      do l=1,fl
 40     if (file(l:l).eq.' ') then
          fl=fl-1
          do m=l,fl
            file(m:m)=file(m+1:m+1)
          enddo
          if (fl.eq.l) goto 50
          goto 40
        endif
      enddo
 50   continue

c...delete possible spaces in file name
      fl=200
      do l=1,fl
 60     if (filene(l:l).eq.' ') then
          fl=fl-1
          do m=l,fl
            filene(m:m)=filene(m+1:m+1)
          enddo
          if (fl.eq.l) goto 70
          goto 60
        endif
      enddo
 70   continue


c...First define the abundances of heavy elements not listed in the
c...Bahcall et al standard solar model file
c...The number fractions are given by n_i = n_H log10(A-12)
c...The abundances listed here are from N. Grevesse and A.J. Sauval,
c...Space Science Reviews 85 (1998) 161.
      sdabund(7)  = 8.08d0  ! Ne
      sdabund(8)  = 6.33d0  ! Na
      sdabund(9)  = 7.58d0  ! Mg
      sdabund(10) = 6.47d0  ! Al
      sdabund(11) = 7.55d0  ! Si
      sdabund(12) = 7.33d0  ! S
      sdabund(13) = 6.40d0  ! Ar
      sdabund(14) = 6.36d0  ! Ca
      sdabund(15) = 7.50d0  ! Fe
      sdabund(16) = 6.25d0  ! Ni
      

c...Mass numbers
      sdaa(1)  = 1.0d0  ! H
      sdaa(2)  = 4.0d0  ! He4
      sdaa(3)  = 3.0d0  ! He3
      sdaa(4)  = 12.0d0 ! C12
      sdaa(5)  = 14.0d0 ! N14
      sdaa(6)  = 16.0d0 ! O8
      sdaa(7)  = 20.0d0 ! Ne
      sdaa(8)  = 23.0d0 ! Na
      sdaa(9)  = 24.0d0 ! Mg
      sdaa(10) = 27.0d0 ! Al
      sdaa(11) = 28.0d0 ! Si
      sdaa(12) = 32.0d0 ! S
      sdaa(13) = 40.0d0 ! Ar
      sdaa(14) = 40.0d0 ! Ca
      sdaa(15) = 56.0d0 ! Fe
      sdaa(16) = 59.0d0 ! Ni

c...Masses
      do i=1,16
        sdma(i)=sdaa(i)*(m_p+m_n)/2.0d0
      enddo

c...Now read in data from solar model data file
      write(*,*) 'dsntsunread: Opening file ',file
      open(unit=13,file=file,
     &  form='formatted',status='old')
      sdn=2  ! will add first line later

c...Read in scratch lines, NOTE: this might have to change with newer files
      do i=1,23
        read(13,'(A)',end=110) scr
      enddo

c...Calculate the total mass fraction (relative to H) of elements
c...heavier than O. This is just for normalization of the fractions
c...as a function of radius below
      totfrheavy=0.0d0
      do i=7,16
        totfrheavy=totfrheavy+sdma(i)*10**(sdabund(i)-12.0d0)
      enddo

c...Read in table and calculate remaining mass fractions from
c...average abundances sdabund
 100  read(13,*,end=110) sdm(sdn),sdr(sdn),tmp1,sdrho(sdn),tmp2,tmp3,
     &  sdmfr(1,sdn),sdmfr(2,sdn),sdmfr(3,sdn),sdmfr(4,sdn),
     &  sdmfr(5,sdn),sdmfr(6,sdn)
      totfr=sdmfr(1,sdn)+sdmfr(2,sdn)+sdmfr(3,sdn)+sdmfr(4,sdn)
     &  +sdmfr(5,sdn)+sdmfr(6,sdn)
c...Add the heavy elements (>O16)
      do i=7,16
        sdmfr(i,sdn)=(1.0d0-totfr)*sdma(i)*10**(sdabund(i)-12.0d0)
     &    /totfrheavy      
      enddo
      sdn=sdn+1
      if (sdn.gt.sdmax-1) then  ! need one more entry for last line
        goto 110
      endif
      goto 100
 110  continue
      sdn=sdn-1
      close(13)

c...Now add the first line with r/r_sun=0
      sdm(1)=0.0d0
      sdr(1)=0.0d0
      sdrho(1)=sdrho(2)
      do i=1,16
        sdmfr(i,1)=sdmfr(i,2)
      enddo

c...Now add the last line with r/r_sun=1
      sdn=sdn+1
      sdm(sdn)=1.0d0
      sdr(sdn)=1.0d0
      sdrho(sdn)=0.0d0
      do i=1,16
        sdmfr(i,sdn)=sdmfr(i,sdn-1)
      enddo


c...Electron density
c...Now read in data from solar model electron density file
      write(*,*) 'dsntsunread: Opening file ',filene
      open(unit=13,file=filene,
     &  form='formatted',status='old')
      sdnne=2  ! will add first line later

c...Read in scratch lines, NOTE: this might have to change with newer files
      do i=1,6
        read(13,'(A)',end=110) scr
      enddo

c...Read in table and calculate remaining mass fractions from
c...average abundances sdabund
 200  read(13,*,end=210) sdrne(sdnne),sdne(sdnne)
      sdnne=sdnne+1
      if (sdnne.gt.sdmax) then  ! need one more entry for last line
        write(*,*) 'ERROR in dsntsunread: array too small.'
        write(*,*) 'Increase sdmax.'
        goto 210
      endif
      goto 200
 210  continue
      sdnne=sdnne-1
      close(13)

c...Now add the first line with r/r_sun=0
      sdne(1)=sdne(2)
      sdrne(1)=0.0d0

c...Sun parameters
      m_sun = 1.9889d30   ! solar mass in kg
      r_sun = 6.9598d8    ! solar radius in m


c=======================================================================
c=======================================================================
c=======================================================================

c...Now calculate and tabulate the potential inside the Sun

      write(*,*) 'dsntsunread: tabulating potential inside the Sun...'
      do i=1,sdn
        sdphi(i)=dsntsunpotint(sdr(i)*r_sun)
      enddo
      write(*,*) '  ...done'

      write(*,*) 
     &  'dsntsunread: tabulating column density inside the Sun...'
      do i=1,sdn
        sdcdens(i,0)=dsntsuncdensint(sdr(i)*r_sun,'N')
        sdcdens(i,1)=dsntsuncdensint(sdr(i)*r_sun,'p')
        sdcdens(i,2)=dsntsuncdensint(sdr(i)*r_sun,'n')
      enddo
      write(*,*) '  ...done'
      cd_sun(0)=sdcdens(sdn,0)  ! total (p+n) column density g/cm^2
      cd_sun(1)=sdcdens(sdn,1)  ! total proton column density g/cm^2
      cd_sun(2)=sdcdens(sdn,2)  ! total neutron column density g/cm^2

c      write(*,*) 'Ip = ',cd_sun(1)/r_sun/100.0d0
c      write(*,*) 'In = ',cd_sun(2)/r_sun/100.0d0

      return

      end


      

**********************************************************************
*** dsntsuncdensint gives the column density in the Sun from the
*** centre out the tha given radius r (in meters)
*** if type = 'N', the total column density (up to that r) is calculated
***         = 'p', the column density on protons is calculated
***         = 'n', the column density on neutrons is calculated
*** in this routine, the actual integration is performed. for speed,
*** use dsntsuncdens instead which uses a tabulation of this result.
*** Author: joakim edsjo (edsjo@physto.se)
*** Date: November 24, 2005
**********************************************************************

      real*8 function dsntsuncdensint(r,type)
      implicit none
      include 'dssun.h'

      real*8 r,dsntsuncdfunc,dsf_int
      character*1 type
      external dsntsuncdfunc

      cdt=type ! transfer to common block
c...integrate sun density, the last factor of 100.0 comes from m->cm

      dsntsuncdensint=
     &  dsf_int(dsntsuncdfunc,0.0d0,min(r,r_sun),1.d-3)*100.0d0

      return
      end
***********************************************************************
*** dsntsuncdfunc returns the density of protons, neutrons or the total
*** density depending on the common block variable cdt. If
***   cdt='N': the total density is returned
***   cdt='p': the density in protons is returned
***   cdt='n': the density in neutrons is returned
*** the radius should be given in m and the density is returned in
*** g/cm^3.
*** This routine is used by dsntsuncdensint to calculate the column
*** density in the Sun.
***
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 2005-11-24
***********************************************************************

      real*8 function dsntsuncdfunc(r)
      implicit none

      include 'dssun.h'
      real*8 dsntsunmfrac,pfr,dsntsundens,r
      integer i,j

c...Check if data file is loaded
      call dsntsunread

      if (cdt.eq.'N') then
        dsntsuncdfunc=dsntsundens(r)
      elseif (cdt.eq.'p') then
        pfr=dsntsunmfrac(r,1)
        dsntsuncdfunc=dsntsundens(r)*(pfr + 0.5d0*(1.0d0-pfr))
      elseif (cdt.eq.'n') then
        pfr=dsntsunmfrac(r,1)
        dsntsuncdfunc=dsntsundens(r)*(0.5d0*(1.0d0-pfr))
      else
        write(*,*) 'ERROR in dsntsuncdfunc: wrong type: ',cdt
        stop
      endif

      return

      end
**********************************************************************
*** dsntsunpotint gives the gravitational potential inside and outside
*** of the sun as a function of the radius r (in meters).
*** in this routine, the actual integration is performed. for speed,
*** use dsntsunpot instead which uses a tabulation of this result.
*** author: joakim edsjo (edsjo@physto.se)
*** date: april 1, 1999
**********************************************************************

      real*8 function dsntsunpotint(r)
      implicit none
      include 'dssun.h'

      real*8 r,dsntspfunc,dsf_int,gn,dsntsunmass
      parameter(gn=6.67259d-11) ! m^3 kg^-1 s^-1
      external dsntspfunc

c...integrate sun density

      if (r.lt.r_sun) then
        dsntsunpotint=
     &    -gn * m_sun/r_sun  ! surface potential, -1.9069d11 m^2 s^-2
     &    -dsf_int(dsntspfunc,max(r,100.0d0),r_sun,1.d-3)
      else
        dsntsunpotint=-m_sun*gn/(max(r,100.0d0))
      endif

      return
      end

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
***********************************************************************
*** dsntsunmass gives the mass of the Sun as a function of radius
*** the radius should be given in m and the mass is given in kg
*** up to the specified radius.
*** Density and element mass fractions up to O16 are from the standard
*** solar model BP2000 of Bahcall, Pinsonneault and Basu,
*** ApJ 555 (2001) 990.
*** The mass fractions for heavier elements are from N. Grevesse and
*** A.J. Sauval, Space Science Reviews 85 (1998) 161 normalized such that
*** their total mass fractions matches that of the heavier elements in 
*** the BP2000 model.
***
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 2003-11-25
***********************************************************************

      real*8 function dsntsunmass(r)
      implicit none

      include 'dssun.h'
      real*8 r,rpl
      integer i,j

c...Check if data file is loaded
      call dsntsunread

      if (r.gt.r_sun) then
        dsntsunmass=m_sun
        return
      endif

      if (r.ge.0.0d0.and.r.lt.sdr(2)*r_sun) then ! do a better interpolation
        dsntsunmass=sdm(2)*m_sun*r**3/(sdr(2)*r_sun)**3
        return
      endif

      call dshunt(sdr,sdn,r/r_sun,j)
      if (j.lt.sdn) goto 20

      dsntsunmass=m_sun
      return

 20   rpl=(r-sdr(j)*r_sun)/(sdr(j+1)*r_sun-sdr(j)*r_sun)

      dsntsunmass=(sdm(j)*(1.0d0-rpl)+sdm(j+1)*rpl)*m_sun

      return

      end
************************************************************************
      subroutine dshunt(xx,n,x,indx)
*** returns the lowest index indx for which x>xx(indx).
*** if x<= xx(i) it returns 0, 
*** if x>xx(n) it returns indx=n
************************************************************************
      implicit none
      integer n,indx
      real*8 x,xx(n)
      integer inc,jhi,jm
      logical ascnd
c-----------------------------------------------------------------------
      ascnd=xx(n).gt.xx(1)
      if(indx.le.0.or.indx.gt.n)then
         indx=0
         jhi=n+1
         goto 3
      endif
      inc=1
      if(x.ge.xx(indx).eqv.ascnd)then
 1       jhi=indx+inc
         if(jhi.gt.n)then
            jhi=n+1
         else if(x.ge.xx(jhi).eqv.ascnd)then
            indx=jhi
            inc=inc+inc
            goto 1
         endif
      else
         jhi=indx
 2       indx=jhi-inc
         if(indx.lt.1)then
            indx=0
         else if(x.lt.xx(indx).eqv.ascnd)then
            jhi=indx
            inc=inc+inc
            goto 2
         endif
      endif
 3    if(jhi-indx.eq.1)return
      jm=(jhi+indx)/2
      if(x.gt.xx(jm).eqv.ascnd)then
         indx=jm
      else
         jhi=jm
      endif
      goto 3
      end
      real*8 function dsf_int(f,a,b,eps)
c_______________________________________________________________________
c  integrate function f between a and b
c  input
c    integration limits a and b
c  called by different routines
c  author: joakim edsjo (edsjo@physto.se) 96-05-16
c          2000-07-19 paolo gondolo added eps as argument 
c  based on paolo gondolos wxint.f routine.
c=======================================================================
      implicit none
      include 'dsidtag.h'
      real*8 f,a,b,tot,eps,st,os,ost,del,sum,x
      integer jmax,it,l,j,nfcn,jdid
      external f
      parameter (jmax=20)
      dsf_int=0.d0
      del=b-a
      ost=0.5*del*(f(a)+f(b))
      x=0.5*(b+a)
      st=0.5*(ost+del*f(x))
      os=(4.0*st-ost)/3.0
      ost=st
      it=1
      nfcn=3
      do j=3,jmax
        it=2*it
        del=0.5*del
        x=a+0.5*del
        sum=0.0
        do l=1,it
          sum=sum+f(x)
          nfcn=nfcn+1
          x=x+del
        enddo
        st=0.5*(st+del*sum)
        tot=(4.0*st-ost)/3.0
        jdid=j
        if (abs(tot-os).le.eps*abs(os)) then
           dsf_int=tot
           return
        endif
     	os=tot
        ost=st
      enddo
      write(*,*) 'error in dsf_int: too many steps for model: ',
     &  idtag
      write(*,*) '  integral set to zero.'
      dsf_int=0.0d0
      end
      SUBROUTINE dsspline(x,y,n,yp1,ypn,y2)
c spline routine, double precision
      implicit none
      INTEGER n,NMAX
      double precision yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=10000)
      INTEGER i,k
      double precision p,qn,sig,un,u(NMAX)
      if (yp1.gt..99d30) then
        y2(1)=0.d0
        u(1)=0.d0
      else
        y2(1)=-0.5d0
        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.d0
        y2(i)=(sig-1.d0)/p
        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))
     &         -(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))
     &         -sig*u(i-1))/p
11    continue
      if (ypn.gt..99d30) then
        qn=0.d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END
      SUBROUTINE dssplint(xa,ya,y2a,n,x,y)
c spline routine, double precision
      implicit none
      INTEGER n
      double precision x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      double precision a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.d0) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))
     &   *(h**2)/6.d0
      return
      END
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

      real*8 function d1mach(i)
      integer i
c
c  double-precision machine constants
c  d1mach( 1) = b**(emin-1), the smallest positive magnitude.
c  d1mach( 2) = b**emax*(1 - b**(-t)), the largest magnitude.
c  d1mach( 3) = b**(-t), the smallest relative spacing.
c  d1mach( 4) = b**(1-t), the largest relative spacing.
c  d1mach( 5) = log10(b)
c
      integer small(2)
      integer large(2)
      integer right(2)
      integer diver(2)
      integer log10(2)
      integer sc, cray1(38), j
      common /d9mach/ cray1
      save small, large, right, diver, log10, sc
      real*8 dmach(5)
      equivalence (dmach(1),small(1))
      equivalence (dmach(2),large(1))
      equivalence (dmach(3),right(1))
      equivalence (dmach(4),diver(1))
      equivalence (dmach(5),log10(1))
c  this version adapts automatically to most current machines.
c  r1mach can handle auto-double compiling, but this version of
c  d1mach does not, because we do not have quad constants for
c  many machines yet.
c  to compile on older machines, add a c in column 1
c  on the next line
      data sc/0/
c  and remove the c from column 1 in one of the sections below.
c  constants for even older machines can be obtained by
c          mail netlib@research.bell-labs.com
c          send old1mach from blas
c  please send corrections to dmg or ehg@bell-labs.com.
c
c     machine constants for the honeywell dps 8/70 series.
c      data small(1),small(2) / o402400000000, o000000000000 /
c      data large(1),large(2) / o376777777777, o777777777777 /
c      data right(1),right(2) / o604400000000, o000000000000 /
c      data diver(1),diver(2) / o606400000000, o000000000000 /
c      data log10(1),log10(2) / o776464202324, o117571775714 /, sc/987/
c
c     machine constants for pdp-11 fortrans supporting
c     32-bit integers.
c      data small(1),small(2) /    8388608,           0 /
c      data large(1),large(2) / 2147483647,          -1 /
c      data right(1),right(2) /  612368384,           0 /
c      data diver(1),diver(2) /  620756992,           0 /
c      data log10(1),log10(2) / 1067065498, -2063872008 /, sc/987/
c
c     machine constants for the univac 1100 series.
c      data small(1),small(2) / o000040000000, o000000000000 /
c      data large(1),large(2) / o377777777777, o777777777777 /
c      data right(1),right(2) / o170540000000, o000000000000 /
c      data diver(1),diver(2) / o170640000000, o000000000000 /
c      data log10(1),log10(2) / o177746420232, o411757177572 /, sc/987/
c
c     on first call, if no data uncommented, test machine types.
      if (sc .ne. 987) then
         dmach(1) = 1.d13
         if (      small(1) .eq. 1117925532
     *       .and. small(2) .eq. -448790528) then
*           *** ieee big endian ***
            small(1) = 1048576
            small(2) = 0
            large(1) = 2146435071
            large(2) = -1
            right(1) = 1017118720
            right(2) = 0
            diver(1) = 1018167296
            diver(2) = 0
            log10(1) = 1070810131
            log10(2) = 1352628735
         else if ( small(2) .eq. 1117925532
     *       .and. small(1) .eq. -448790528) then
*           *** ieee little endian ***
            small(2) = 1048576
            small(1) = 0
            large(2) = 2146435071
            large(1) = -1
            right(2) = 1017118720
            right(1) = 0
            diver(2) = 1018167296
            diver(1) = 0
            log10(2) = 1070810131
            log10(1) = 1352628735
         else if ( small(1) .eq. -2065213935
     *       .and. small(2) .eq. 10752) then
*               *** vax with d_floating ***
            small(1) = 128
            small(2) = 0
            large(1) = -32769
            large(2) = -1
            right(1) = 9344
            right(2) = 0
            diver(1) = 9472
            diver(2) = 0
            log10(1) = 546979738
            log10(2) = -805796613
         else if ( small(1) .eq. 1267827943
     *       .and. small(2) .eq. 704643072) then
*               *** ibm mainframe ***
            small(1) = 1048576
            small(2) = 0
            large(1) = 2147483647
            large(2) = -1
            right(1) = 856686592
            right(2) = 0
            diver(1) = 873463808
            diver(2) = 0
            log10(1) = 1091781651
            log10(2) = 1352628735
         else if ( small(1) .eq. 1120022684
     *       .and. small(2) .eq. -448790528) then
*           *** convex c-1 ***
            small(1) = 1048576
            small(2) = 0
            large(1) = 2147483647
            large(2) = -1
            right(1) = 1019215872
            right(2) = 0
            diver(1) = 1020264448
            diver(2) = 0
            log10(1) = 1072907283
            log10(2) = 1352628735
         else if ( small(1) .eq. 815547074
     *       .and. small(2) .eq. 58688) then
*           *** vax g-floating ***
            small(1) = 16
            small(2) = 0
            large(1) = -32769
            large(2) = -1
            right(1) = 15552
            right(2) = 0
            diver(1) = 15568
            diver(2) = 0
            log10(1) = 1142112243
            log10(2) = 2046775455
         else
            dmach(2) = 1.d27 + 1
            dmach(3) = 1.d27
            large(2) = large(2) - right(2)
            if (large(2) .eq. 64 .and. small(2) .eq. 0) then
               cray1(1) = 67291416
               do 10 j = 1, 20
                  cray1(j+1) = cray1(j) + cray1(j)
 10               continue
               cray1(22) = cray1(21) + 321322
               do 20 j = 22, 37
                  cray1(j+1) = cray1(j) + cray1(j)
 20               continue
               if (cray1(38) .eq. small(1)) then
*                  *** cray ***
                  call i1mcry(small(1), j, 8285, 8388608, 0)
                  small(2) = 0
                  call i1mcry(large(1), j, 24574, 16777215, 16777215)
                  call i1mcry(large(2), j, 0, 16777215, 16777214)
                  call i1mcry(right(1), j, 16291, 8388608, 0)
                  right(2) = 0
                  call i1mcry(diver(1), j, 16292, 8388608, 0)
                  diver(2) = 0
                  call i1mcry(log10(1), j, 16383, 10100890, 8715215)
                  call i1mcry(log10(2), j, 0, 16226447, 9001388)
               else
                  write(*,9000)
                  stop 779
                  end if
            else
               write(*,9000)
               stop 779
               end if
            end if
         sc = 987
         end if
*    sanity check
      if (dmach(4) .ge. 1.0d0) stop 778
      if (i .lt. 1 .or. i .gt. 5) then
         write(*,*) 'd1mach(i): i =',i,' is out of bounds.'
         stop
         end if
      d1mach = dmach(i)
      return
 9000 format(/' adjust d1mach by uncommenting data statements'/
     *' appropriate for your machine.')
* /* standard c source for d1mach -- remove the * in column 1 */
*#include <stdio.h>
*#include <float.h>
*#include <math.h>
*double d1mach_(long *i)
*{
*	switch(*i){
*	  case 1: return dbl_min;
*	  case 2: return dbl_max;
*	  case 3: return dbl_epsilon/flt_radix;
*	  case 4: return dbl_epsilon;
*	  case 5: return log10((double)flt_radix);
*	  }
*	fprintf(stderr, "invalid argument: d1mach(%ld)\n", *i);
*	exit(1); return 0; /* some compilers demand return values */
*}
      end
      subroutine i1mcry(a, a1, b, c, d)
**** special computation for old cray machines ****
      integer a, a1, b, c, d
      a1 = 16777216*b + c
      a = 16777216*a1 + d
      end
      subroutine dqpsrt(limit,last,maxerr,ermax,elist,iord,nrmax)
c***begin prologue  dqpsrt
c***refer to  dqage,dqagie,dqagpe,dqawse
c***routines called  (none)
c***revision date  810101   (yymmdd)
c***keywords  sequential sorting
c***author  piessens, robert, applied math. and progr. div. -
c             k. u. leuven
c           de doncker, elise, applied math. and progr. div. -
c             k. u. leuven
c***purpose  this routine maintains the descending ordering in the
c            list of the local error estimated resulting from the
c            interval subdivision process. at each call two error
c            estimates are inserted using the sequential search
c            method, top-down for the largest error estimate and
c            bottom-up for the smallest error estimate.
c***description
c
c           ordering routine
c           standard fortran subroutine
c           real*8 version
c
c           parameters (meaning at output)
c              limit  - integer
c                       maximum number of error estimates the list
c                       can contain
c
c              last   - integer
c                       number of error estimates currently in the list
c
c              maxerr - integer
c                       maxerr points to the nrmax-th largest error
c                       estimate currently in the list
c
c              ermax  - real*8
c                       nrmax-th largest error estimate
c                       ermax = elist(maxerr)
c
c              elist  - real*8
c                       vector of dimension last containing
c                       the error estimates
c
c              iord   - integer
c                       vector of dimension last, the first k elements
c                       of which contain pointers to the error
c                       estimates, such that
c                       elist(iord(1)),...,  elist(iord(k))
c                       form a decreasing sequence, with
c                       k = last if last.le.(limit/2+2), and
c                       k = limit+1-last otherwise
c
c              nrmax  - integer
c                       maxerr = iord(nrmax)
c***end prologue  dqpsrt
c
      real*8 elist,ermax,errmax,errmin
      integer i,ibeg,ido,iord,isucc,j,jbnd,jupbn,k,last,limit,maxerr,
     1  nrmax
      dimension elist(last),iord(last)
c
c           check whether the list contains more than
c           two error estimates.
c
c***first executable statement  dqpsrt
      if(last.gt.2) go to 10
      iord(1) = 1
      iord(2) = 2
      go to 90
c
c           this part of the routine is only executed if, due to a
c           difficult integrand, subdivision increased the error
c           estimate. in the normal case the insert procedure should
c           start after the nrmax-th largest error estimate.
c
   10 errmax = elist(maxerr)
      if(nrmax.eq.1) go to 30
      ido = nrmax-1
      do 20 i = 1,ido
        isucc = iord(nrmax-1)
c ***jump out of do-loop
        if(errmax.le.elist(isucc)) go to 30
        iord(nrmax) = isucc
        nrmax = nrmax-1
   20    continue
c
c           compute the number of elements in the list to be maintained
c           in descending order. this number depends on the number of
c           subdivisions still allowed.
c
   30 jupbn = last
      if(last.gt.(limit/2+2)) jupbn = limit+3-last
      errmin = elist(last)
c
c           insert errmax by traversing the list top-down,
c           starting comparison from the element elist(iord(nrmax+1)).
c
      jbnd = jupbn-1
      ibeg = nrmax+1
      if(ibeg.gt.jbnd) go to 50
      do 40 i=ibeg,jbnd
        isucc = iord(i)
c ***jump out of do-loop
        if(errmax.ge.elist(isucc)) go to 60
        iord(i-1) = isucc
   40 continue
   50 iord(jbnd) = maxerr
      iord(jupbn) = last
      go to 90
c
c           insert errmin by traversing the list bottom-up.
c
   60 iord(i-1) = maxerr
      k = jbnd
      do 70 j=i,jbnd
        isucc = iord(k)
c ***jump out of do-loop
        if(errmin.lt.elist(isucc)) go to 80
        iord(k+1) = isucc
        k = k-1
   70 continue
      iord(i) = last
      go to 90
   80 iord(k+1) = last
c
c           set maxerr and ermax.
c
   90 maxerr = iord(nrmax)
      ermax = elist(maxerr)
      return
      end
      subroutine dqelg(n,epstab,result,abserr,res3la,nres)
c***begin prologue  dqelg
c***refer to  dqagie,dqagoe,dqagpe,dqagse
c***routines called  d1mach
c***revision date  830518   (yymmdd)
c***keywords  convergence acceleration,epsilon algorithm,extrapolation
c***author  piessens, robert, applied math. and progr. div. -
c             k. u. leuven
c           de doncker, elise, applied math. and progr. div. -
c             k. u. leuven
c***purpose  the routine determines the limit of a given sequence of
c            approximations, by means of the epsilon algorithm of
c            p.wynn. an estimate of the absolute error is also given.
c            the condensed epsilon table is computed. only those
c            elements needed for the computation of the next diagonal
c            are preserved.
c***description
c
c           epsilon algorithm
c           standard fortran subroutine
c           real*8 version
c
c           parameters
c              n      - integer
c                       epstab(n) contains the new element in the
c                       first column of the epsilon table.
c
c              epstab - real*8
c                       vector of dimension 52 containing the elements
c                       of the two lower diagonals of the triangular
c                       epsilon table. the elements are numbered
c                       starting at the right-hand corner of the
c                       triangle.
c
c              result - real*8
c                       resulting approximation to the integral
c
c              abserr - real*8
c                       estimate of the absolute error computed from
c                       result and the 3 previous results
c
c              res3la - real*8
c                       vector of dimension 3 containing the last 3
c                       results
c
c              nres   - integer
c                       number of calls to the routine
c                       (should be zero at first call)
c***end prologue  dqelg
c
      real*8 abserr,dabs,delta1,delta2,delta3,dmax1,d1mach,
     1  epmach,epsinf,epstab,error,err1,err2,err3,e0,e1,e1abs,e2,e3,
     2  oflow,res,result,res3la,ss,tol1,tol2,tol3
      integer i,ib,ib2,ie,indx,k1,k2,k3,limexp,n,newelm,nres,num
      dimension epstab(52),res3la(3)
c
c           list of major variables
c           -----------------------
c
c           e0     - the 4 elements on which the computation of a new
c           e1       element in the epsilon table is based
c           e2
c           e3                 e0
c                        e3    e1    new
c                              e2
c           newelm - number of elements to be computed in the new
c                    diagonal
c           error  - error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
c           result - the element in the new diagonal with least value
c                    of error
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           oflow is the largest positive magnitude.
c           limexp is the maximum number of elements the epsilon
c           table can contain. if this number is reached, the upper
c           diagonal of the epsilon table is deleted.
c
c***first executable statement  dqelg
      epmach = d1mach(4)
      oflow = d1mach(2)
      nres = nres+1
      abserr = oflow
      result = epstab(n)
      if(n.lt.3) go to 100
      limexp = 50
      epstab(n+2) = epstab(n)
      newelm = (n-1)/2
      epstab(n) = oflow
      num = n
      k1 = n
      do 40 i = 1,newelm
        k2 = k1-1
        k3 = k1-2
        res = epstab(k1+2)
        e0 = epstab(k3)
        e1 = epstab(k2)
        e2 = res
        e1abs = dabs(e1)
        delta2 = e2-e1
        err2 = dabs(delta2)
        tol2 = dmax1(dabs(e2),e1abs)*epmach
        delta3 = e1-e0
        err3 = dabs(delta3)
        tol3 = dmax1(e1abs,dabs(e0))*epmach
        if(err2.gt.tol2.or.err3.gt.tol3) go to 10
c
c           if e0, e1 and e2 are equal to within machine
c           accuracy, convergence is assumed.
c           result = e2
c           abserr = abs(e1-e0)+abs(e2-e1)
c
        result = res
        abserr = err2+err3
c ***jump out of do-loop
        go to 100
   10   e3 = epstab(k1)
        epstab(k1) = e1
        delta1 = e1-e3
        err1 = dabs(delta1)
        tol1 = dmax1(e1abs,dabs(e3))*epmach
c
c           if two elements are very close to each other, omit
c           a part of the table by adjusting the value of n
c
        if(err1.le.tol1.or.err2.le.tol2.or.err3.le.tol3) go to 20
        ss = 0.1d+01/delta1+0.1d+01/delta2-0.1d+01/delta3
        epsinf = dabs(ss*e1)
c
c           test to detect irregular behaviour in the table, and
c           eventually omit a part of the table adjusting the value
c           of n.
c
        if(epsinf.gt.0.1d-03) go to 30
   20   n = i+i-1
c ***jump out of do-loop
        go to 50
c
c           compute a new element and eventually adjust
c           the value of result.
c
   30   res = e1+0.1d+01/ss
        epstab(k1) = res
        k1 = k1-2
        error = err2+dabs(res-e2)+err3
        if(error.gt.abserr) go to 40
        abserr = error
        result = res
   40 continue
c
c           shift the table.
c
   50 if(n.eq.limexp) n = 2*(limexp/2)-1
      ib = 1
      if((num/2)*2.eq.num) ib = 2
      ie = newelm+1
      do 60 i=1,ie
        ib2 = ib+2
        epstab(ib) = epstab(ib2)
        ib = ib2
   60 continue
      if(num.eq.n) go to 80
      indx = num-n+1
      do 70 i = 1,n
        epstab(i)= epstab(indx)
        indx = indx+1
   70 continue
   80 if(nres.ge.4) go to 90
      res3la(nres) = result
      abserr = oflow
      go to 100
c
c           compute error estimate
c
   90 abserr = dabs(result-res3la(3))+dabs(result-res3la(2))
     1  +dabs(result-res3la(1))
      res3la(1) = res3la(2)
      res3la(2) = res3la(3)
      res3la(3) = result
  100 abserr = dmax1(abserr,0.5d+01*epmach*dabs(result))
      return
      end

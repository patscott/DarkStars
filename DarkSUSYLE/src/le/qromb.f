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

      MODULE ez_vcool
      USE ez_vcool_data
      IMPLICIT NONE
      
      CONTAINS

      DOUBLE PRECISION FUNCTION Neutrino_Cooling(rho,T,lnrho,lnT,mue,xc,xo,xhe,nu_plasma,nu_brem,nu_pair,nu_photo)
      DOUBLE PRECISION, INTENT(IN) :: lnrho,lnT,rho,T,mue,xc,xo,xhe
      DOUBLE PRECISION, INTENT(OUT) :: nu_plasma,nu_brem,nu_pair,nu_photo
      DOUBLE PRECISION :: Anum(3), Z(3), nu_X(3)
      IF ( T .GE. 1D7 .AND. rho .GE. 1D0 ) THEN
         ! neutrino loss rates according to Itoh et.al., ApJS, 102, pg.411, 1996.
         nu_X(1)=xc
         nu_X(2)=xo
         nu_X(3)=xhe
         Anum = (/12D0,16D0,4D0/); Z = (/6D0,8D0,2D0/) ! Anum and Z are for the Itoh code
         call photo(lnrho,lnT,rho,T,mue,nu_photo)
         call plasma(lnrho,lnT,rho,T,mue,nu_plasma)
         call brem(lnrho,lnT,rho,T,nu_X,Anum,Z,nu_brem)
         call pair(lnrho,lnT,rho,T,mue,nu_pair)
         Neutrino_Cooling = -(nu_plasma + nu_brem + nu_pair + nu_photo)
      ELSE
         Neutrino_Cooling = 0D0
      END IF
      END FUNCTION Neutrino_Cooling

!  Neutrino emission rates from Itoh et.al., 
!  ApJS, 102, pg.411, 1996.

!  Q is in units of erg sec^{-1} cm^{-3}
!  epsilon=Q/\rho in units of erg sec^{-1} gram^{-1}
!  rho in units of gram cm^{-3}
!  mue is AMU per free electron
!  T in Kelvin

      subroutine init_neutrinorates
      LN_10 = LOG(10D0)       ! divide by LN_10 to convert ln to log10
      PI = 4D0*ATAN(1D0)
      ! call this during startup to create tables, etc.
      CV=0.5 + 2*sin2thetaw
      CAX=0.5
      CVprime=1.0-CV
      CAXprime=1.0-CAX
      call init_photo
      call init_plasma
      call init_brem
      call init_pair
      end subroutine init_neutrinorates

!  -----------------------------------------------------

      subroutine photo(lnrho,lnT,rho,T,mue,e_photo)
      DOUBLE PRECISION, INTENT(IN) :: lnrho,lnT,rho,T,mue
      DOUBLE PRECISION, INTENT(OUT) :: e_photo
      ! locals
      DOUBLE PRECISION :: pb(1:3),pc,tau,cc(1:3,0:3,0:6),dd(1:3,0:3,0:6)
      DOUBLE PRECISION :: aa(0:2),fphoto,prefactor,num,denom,denom2,lambda,xi,lqphoto,Qphoto,log10T
      DOUBLE PRECISION :: cos_tau, cos_1tau, cos_2tau, cos_3tau, cos_4tau, cos_5tau
      DOUBLE PRECISION :: sin_1tau, sin_2tau, sin_3tau, sin_4tau, sin_5tau
      integer kk
      common /nuphot/ cc,dd,pb
      if ( T .le. 1.e7 ) then
         e_photo=1.e-30
         return
      endif
      log10T = lnT / LN_10
      pc=0.5654 + min(1.0d0,log10T-7D0)
      if ( T.lt.1.e8 ) then
         tau = log10T-7D0; kk = 1
      else if ( T.lt.1.e9 ) then
         tau = log10T-8D0; kk = 2
      else
         tau = log10T-9D0; kk = 3
      endif
      cos_tau = cos(10.0*PI*tau)
      sin_1tau = sin(5.0*PI*tau/3.0)
      sin_2tau = sin(5.0*PI*2*tau/3.0)
      sin_3tau = sin(5.0*PI*tau)
      sin_4tau = sin(5.0*PI*4*tau/3.0)
      sin_5tau = sin(5.0*PI*5*tau/3.0)
      cos_1tau = cos(5.0*PI*tau/3.0)
      cos_2tau = cos(5.0*PI*2*tau/3.0)
      cos_3tau = cos(5.0*PI*tau)
      cos_4tau = cos(5.0*PI*4*tau/3.0)
      cos_5tau = cos(5.0*PI*4*tau/3.0)
      aa(0)=0.5*cc(kk,0,0) + 0.5*cc(kk,0,6)*cos_tau
      aa(0)=aa(0) + cc(kk,0,1)*cos_1tau + dd(kk,0,1)*sin_1tau
      aa(0)=aa(0) + cc(kk,0,2)*cos_2tau + dd(kk,0,2)*sin_2tau
      aa(0)=aa(0) + cc(kk,0,3)*cos_3tau + dd(kk,0,3)*sin_3tau
      aa(0)=aa(0) + cc(kk,0,4)*cos_4tau + dd(kk,0,4)*sin_4tau
      aa(0)=aa(0) + cc(kk,0,5)*cos_5tau + dd(kk,0,5)*sin_5tau
      aa(1)=0.5*cc(kk,1,0) + 0.5*cc(kk,1,6)*cos_tau
      aa(1)=aa(1) + cc(kk,1,1)*cos_1tau + dd(kk,1,1)*sin_1tau
      aa(1)=aa(1) + cc(kk,1,2)*cos_2tau + dd(kk,1,2)*sin_2tau
      aa(1)=aa(1) + cc(kk,1,3)*cos_3tau + dd(kk,1,3)*sin_3tau
      aa(1)=aa(1) + cc(kk,1,4)*cos_4tau + dd(kk,1,4)*sin_4tau
      aa(1)=aa(1) + cc(kk,1,5)*cos_5tau + dd(kk,1,5)*sin_5tau
      aa(2)=0.5*cc(kk,2,0) + 0.5*cc(kk,2,6)*cos_tau
      aa(2)=aa(2) + cc(kk,2,1)*cos_1tau + dd(kk,2,1)*sin_1tau
      aa(2)=aa(2) + cc(kk,2,2)*cos_2tau + dd(kk,2,2)*sin_2tau
      aa(2)=aa(2) + cc(kk,2,3)*cos_3tau + dd(kk,2,3)*sin_3tau
      aa(2)=aa(2) + cc(kk,2,4)*cos_4tau + dd(kk,2,4)*sin_4tau
      aa(2)=aa(2) + cc(kk,2,5)*cos_5tau + dd(kk,2,5)*sin_5tau
      lambda = T/5.9302e9
      xi = (rho/mue/1.e9)**(1.0/3.0) / lambda
      fphoto = (aa(0) + aa(1)*xi + aa(2)*xi**2)*exp(-pc*xi)/(xi**3 + pb(1)/lambda + pb(2)/lambda**2 + pb(3)/lambda**3)
      prefactor=0.5 * ( CV**2 + CAX**2 + nnu*(CVprime**2 + CAXprime**2) )
      num = CV**2-CAX**2 + nnu*(CVprime**2-CAXprime**2)
      denom = CV**2 + CAX**2 + nnu*(CVprime**2 + CAXprime**2)
      denom2 = (1.875e8*lambda + 1.653e8*lambda**2 + 8.499e8*lambda**3-1.604e8*lambda**4)
      lqphoto = 0.666*(1 + 2.045*lambda)**(-2.066)/(1 + (rho/mue)/denom2)
      Qphoto=prefactor * ( 1.0 - num*lqphoto/denom ) * rho/mue * lambda**5 * fphoto
      e_photo = max(1.d-30,Qphoto/rho)
      end subroutine photo

!       -----------------------------------------------------

      subroutine pair(lnrho,lnT,rho,T,mue,e_pair)
      DOUBLE PRECISION, INTENT(IN) :: lnrho,lnT,rho,T,mue
      DOUBLE PRECISION, INTENT(OUT) :: e_pair
      ! locals
      DOUBLE PRECISION aa(0:2),bb(1:2,1:3),cc(1:2),prefactor,num,denom,xi,pg,fpair,lqpair,Qpair
      DOUBLE PRECISION :: lambda, lambda2, lambda3, lambda4, lambda6, lambda8, lambda_sqrt
      common /nupair/ aa,bb,cc,prefactor,num,denom
      integer ii
      if (T .le. 1.e10) then
         ii = 1
      else 
         ii = 2
      endif
      lambda = T/5.9302e9
      lambda2 = lambda * lambda
      lambda3 = lambda * lambda2
      lambda4 = lambda2 * lambda2
      lambda6 = lambda2 * lambda4
      lambda8 = lambda4 * lambda4
      lambda_sqrt = sqrt(lambda)
      xi = (rho/mue/1.e9)**(1.0/3.0) / lambda
      pg = 1.0 - 13.04*lambda2 + 133.5*lambda4 + 1534*lambda6 + 918.6*lambda8
      fpair = (aa(0)+aa(1)*xi+aa(2)*xi**2)*exp(-cc(ii)*xi)/(xi**3+bb(ii,1)/lambda+bb(ii,2)/lambda2+bb(ii,3)/lambda3)
      lqpair = 1.0/(10.7480*lambda2+0.3967*lambda_sqrt+1.0050)/(1.0+(rho/mue)/(7.692e7*lambda3+9.715e6*lambda_sqrt))**0.3
      Qpair=prefactor * ( 1.0 + num*lqpair/denom ) * pg * exp(-2.0/lambda) * fpair
      e_pair = max(1.d-30,Qpair/rho)
      end subroutine pair
      
!       -----------------------------------------------------

      subroutine plasma(lnrho,lnT,rho,T,mue,e_plasma)
      DOUBLE PRECISION, INTENT(IN) :: lnrho,lnT,rho,T,mue
      DOUBLE PRECISION, INTENT(OUT) :: e_plasma
      ! locals
      DOUBLE PRECISION :: prefactor,gamma,fT,fL,x,y,fxy,QV,Q_plasma,lambda,log10T,logrho_mue
      DOUBLE PRECISION :: gamma_sqrt, gamma2
      common /nuplas/ prefactor
      gamma2 = 1.1095e11*rho/mue / T**2 / sqrt( 1.0 + (1.019e-6*rho/mue)**(2.0/3.0) ) 
      gamma = sqrt(gamma2)
      gamma_sqrt = sqrt(gamma)
      fT=2.4 + 0.6*gamma_sqrt + 0.51*gamma + 1.25*gamma*gamma_sqrt
      fL=(8.6*gamma2+1.35*gamma2*gamma_sqrt)/(225.0-17*gamma+gamma2)
      log10T = lnT / LN_10
      logrho_mue = log10(2.0*rho/mue)
      x=( 17.5 + logrho_mue - 3.0*log10T )/6.0
      y=( -24.5 + logrho_mue + 3.0*log10T )/6.0
      if (abs(x) .ge. 0.7 .or. y .le. 0.0) then
         fxy=1.0
      else
         fxy=1.05+exp(-(min(0.0d0,y-1.6+1.25*x)/(0.57-0.25*x))**2)*(0.39-1.25*x-0.35*sin(4.5*x)-0.3*exp(-(4.5*x+0.9)**2))
      endif
      lambda=T/5.9302e9
      QV = 3.e21 * lambda**9 * gamma2**3 * exp(-gamma) * ( fL + fT ) * fxy
      Q_plasma=prefactor*QV
      e_plasma=max(1.d-30,Q_plasma/rho)
      end subroutine plasma

!       -----------------------------------------------------

!  Only liquid metal phase in right now.

      subroutine brem(lnrho,lnT,rho,T,X,Anum,Z,e_brem)
      integer, parameter :: np = 3
      DOUBLE PRECISION, INTENT(IN) :: lnrho,lnT,rho,T,Anum(np),Z(np),X(np)
      DOUBLE PRECISION, INTENT(OUT) :: e_brem
      ! locals
      DOUBLE PRECISION :: T8,rho6,factor0,factor1,factor2,u,F1,F160,G1,G160,v,Fliquid(np),Gliquid(np)
      DOUBLE PRECISION :: Qliquid,w,Gamma,Gamma_m13,Gamma_m23
      DOUBLE PRECISION :: sin_u, sin_2u, sin_3u, sin_4u
      DOUBLE PRECISION :: cos_u, cos_2u, cos_3u, cos_4u, cos_5u
      integer :: ii
      DOUBLE PRECISION :: a(np,0:5),b(np,1:4),c(np),d(np),e(np,0:5)
      DOUBLE PRECISION :: f(np,1:4),g(np),h(np),i(np,0:5),j(np,1:4),k(np)
      DOUBLE PRECISION :: l(np),p(np,0:5),q(np,1:4),r(np),s(np),alpha(np,0:3),beta(np,0:3)
      common /nubrem/ a,b,c,d,e,f,g,h,i,j,k,l,p,q,r,s,alpha,beta
      T8=T/1.e8
      rho6=rho/1.e6
      factor0 = 0.5738*T8**6*rho
      factor1 = 0.5*(CV**2+CAX**2+nnu*(CVprime**2+CAXprime**2))
      factor2 = 0.5*(CV**2-CAX**2+nnu*(CVprime**2-CAXprime**2))
      u=2.0*PI*(lnrho/LN_10-3.0)/10.0
      sin_u = sin(u); sin_2u = sin(2*u); sin_3u = sin(3*u); sin_4u = sin(4*u)
      cos_u = cos(u); cos_2u = cos(2*u); cos_3u = cos(3*u); cos_4u = cos(4*u); cos_5u = cos(5*u)
      Qliquid=0.0
      do ii=1,np ! sum over number of partices
         F1 = 0.5*a(ii,0) + c(ii)*u+d(ii)
         F160 = 0.5*e(ii,0) + g(ii)*u+h(ii)
         G1 = 0.5*i(ii,0) + k(ii)*u+l(ii)
         G160 = 0.5*p(ii,0) + r(ii)*u+s(ii)
         F1 = F1 + a(ii,1)*cos_u + a(ii,2)*cos_2u + a(ii,3)*cos_3u
         F160 = F160 + e(ii,1)*cos_u + e(ii,2)*cos_2u + e(ii,3)*cos_3u
         G1 = G1 + i(ii,1)*cos_u + i(ii,2)*cos_2u + i(ii,3)*cos_3u
         G160 = G160 + p(ii,1)*cos_u + p(ii,2)*cos_2u + p(ii,3)*cos_3u
         F1 = F1 + a(ii,4)*cos_4u + a(ii,5)*cos_5u
         F160 = F160 + e(ii,4)*cos_4u + e(ii,5)*cos_5u
         G1 = G1 + i(ii,4)*cos_4u + i(ii,5)*cos_5u
         G160 = G160 + p(ii,4)*cos_4u + p(ii,5)*cos_5u
         F1 = F1 + b(ii,1)*sin_u + b(ii,2)*sin_2u + b(ii,3)*sin_3u + b(ii,4)*sin_4u
         F160 = F160 + f(ii,1)*sin_u + f(ii,2)*sin_2u + f(ii,3)*sin_3u + f(ii,4)*sin_4u
         G1 = G1 + j(ii,1)*sin_u + j(ii,2)*sin_2u + j(ii,3)*sin_3u + j(ii,4)*sin_4u
         G160 = G160 + q(ii,1)*sin_u + q(ii,2)*sin_2u + q(ii,3)*sin_3u + q(ii,4)*sin_4u
         Gamma = 0.2275*Z(ii)**2/T8*(rho6/Anum(ii))**(1.0/3.0)
         Gamma_m13 = Gamma**(-1.0/3.0)
         Gamma_m23 = Gamma_m13*Gamma_m13
         v = alpha(ii,0) + alpha(ii,1)*Gamma_m13 + alpha(ii,2)*Gamma_m23 + alpha(ii,3)/Gamma
         w = beta(ii,0) + beta(ii,1)*Gamma_m13 + beta(ii,2)*Gamma_m23 + beta(ii,3)/Gamma
         Fliquid(ii) = v*F1+(1.0-v)*F160
         Gliquid(ii) = w*G1+(1.0-w)*G160
         Qliquid = Qliquid + factor0*(X(ii)*Z(ii)**2/Anum(ii)*(factor1*Fliquid(ii)-factor2*Gliquid(ii)))
      end do
      e_brem=max(1.d-30,Qliquid/rho)
      end subroutine brem

!       -----------------------------------------------------

      subroutine init_photo
      DOUBLE PRECISION pb(1:3),cc(1:3,0:3,0:6),dd(1:3,0:3,0:6)
      common /nuphot/ cc,dd,pb
      pb(1)=6.290e-3
      pb(2)=7.483e-3
      pb(3)=3.061e-4
      cc(1,0,0) =  1.008e11
      cc(1,0,1) =  0.e0
      cc(1,0,2) =  0.e0
      cc(1,0,3) =  0.e0
      cc(1,0,4) =  0.e0
      cc(1,0,5) =  0.e0
      cc(1,0,6) =  0.e0
      cc(1,1,0) =  8.156e10
      cc(1,1,1) =  9.728e8
      cc(1,1,2) = -3.806e9
      cc(1,1,3) = -4.384e9
      cc(1,1,4) = -5.774e9
      cc(1,1,5) = -5.249e9
      cc(1,1,6) = -5.153e9
      cc(1,2,0) =  1.067e11
      cc(1,2,1) = -9.782e9
      cc(1,2,2) = -7.193e9
      cc(1,2,3) = -6.936e9
      cc(1,2,4) = -6.893e9
      cc(1,2,5) = -7.041e9
      cc(1,2,6) = -7.193e9
      dd(1,0,1) =  0.e0
      dd(1,0,2) =  0.e0
      dd(1,0,3) =  0.e0
      dd(1,0,4) =  0.e0
      dd(1,0,5) =  0.e0
      dd(1,1,1) = -1.879e10
      dd(1,1,2) = -9.667e9
      dd(1,1,3) = -5.602e9
      dd(1,1,4) = -3.370e9
      dd(1,1,5) = -1.825e9
      dd(1,2,1) = -2.919e10
      dd(1,2,2) = -1.185e10
      dd(1,2,3) = -7.270e9
      dd(1,2,4) = -4.222e9
      dd(1,2,5) = -1.560e9
      cc(2,0,0) =  9.889e10
      cc(2,0,1) = -4.524e8
      cc(2,0,2) = -6.088e6
      cc(2,0,3) =  4.269e7
      cc(2,0,4) =  5.172e7
      cc(2,0,5) =  4.910e7
      cc(2,0,6) =  4.388e7
      cc(2,1,0) =  1.813e11
      cc(2,1,1) = -7.556e9
      cc(2,1,2) = -3.304e9
      cc(2,1,3) = -1.031e9
      cc(2,1,4) = -1.764e9
      cc(2,1,5) = -1.851e9
      cc(2,1,6) = -1.928e9
      cc(2,2,0) =  9.750e10
      cc(2,2,1) =  3.484e10
      cc(2,2,2) =  5.199e9
      cc(2,2,3) = -1.695e9
      cc(2,2,4) = -2.865e9
      cc(2,2,5) = -3.395e9
      cc(2,2,6) = -3.418e9
      dd(2,0,1) = -1.135e8
      dd(2,0,2) =  1.256e8
      dd(2,0,3) =  5.149e7
      dd(2,0,4) =  3.436e7
      dd(2,0,5) =  1.005e7
      dd(2,1,1) =  1.652e9
      dd(2,1,2) = -3.119e9
      dd(2,1,3) = -1.839e9
      dd(2,1,4) = -1.458e9
      dd(2,1,5) = -8.956e8
      dd(2,2,1) = -1.548e10
      dd(2,2,2) = -9.338e9
      dd(2,2,3) = -5.899e9
      dd(2,2,4) = -3.035e9
      dd(2,2,5) = -1.598e9
      cc(3,0,0) =  9.581e10
      cc(3,0,1) =  4.107e8
      cc(3,0,2) =  2.305e8
      cc(3,0,3) =  2.236e8
      cc(3,0,4) =  1.580e8
      cc(3,0,5) =  2.165e8
      cc(3,0,6) =  1.721e8
      cc(3,1,0) =  1.459e12
      cc(3,1,1) =  1.314e11
      cc(3,1,2) = -1.169e11
      cc(3,1,3) = -1.765e11
      cc(3,1,4) = -1.867e11
      cc(3,1,5) = -1.983e11
      cc(3,1,6) = -1.896e11
      cc(3,2,0) =  2.424e11
      cc(3,2,1) = -3.669e9
      cc(3,2,2) = -8.691e9
      cc(3,2,3) = -7.967e9
      cc(3,2,4) = -7.932e9
      cc(3,2,5) = -7.987e9
      cc(3,2,6) = -8.333e9
      dd(3,0,1) =  4.724e8
      dd(3,0,2) =  2.976e8
      dd(3,0,3) =  2.242e8
      dd(3,0,4) =  7.937e7
      dd(3,0,5) =  4.859e7
      dd(3,1,1) = -7.094e11
      dd(3,1,2) = -3.697e11
      dd(3,1,3) = -2.189e11
      dd(3,1,4) = -1.273e11
      dd(3,1,5) = -5.705e10
      dd(3,2,1) = -2.254e10
      dd(3,2,2) = -1.551e10
      dd(3,2,3) = -7.793e9
      dd(3,2,4) = -4.489e9
      dd(3,2,5) = -2.185e9
      end subroutine init_photo

      subroutine init_pair
      DOUBLE PRECISION aa(0:2),bb(1:2,1:3),cc(1:2),prefactor,num,denom
      common /nupair/ aa,bb,cc,prefactor,num,denom
      prefactor=0.5 * ( CV**2 + CAX**2 + nnu*(CVprime**2 + CAXprime**2) )
      num = CV**2-CAX**2 + nnu*(CVprime**2-CAXprime**2)
      denom = CV**2 + CAX**2 + nnu*(CVprime**2 + CAXprime**2)
      aa(0)=6.002e19
      aa(1)=2.084e20
      aa(2)=1.872e21
      bb(1,1)=9.383e-1
      bb(1,2)=-4.141e-1
      bb(1,3)=4.829e-2
      bb(2,1)=1.2383
      bb(2,2)=-0.8141
      bb(2,3)=0.0
      cc(1)=5.5924
      cc(2)=4.9924
      end subroutine init_pair

      subroutine init_plasma
      DOUBLE PRECISION prefactor
      common /nuplas/ prefactor
      prefactor = CV**2 + nnu * CVprime**2
      end subroutine init_plasma

      subroutine init_brem
      integer, parameter :: np=3
      DOUBLE PRECISION :: a(np,0:5),b(np,1:4),c(np),d(np),e(np,0:5)
      DOUBLE PRECISION :: f(np,1:4),g(np),h(np),i(np,0:5),j(np,1:4),k(np)
      DOUBLE PRECISION :: l(np),p(np,0:5),q(np,1:4),r(np),s(np),alpha(np,0:3),beta(np,0:3)
      common /nubrem/ a,b,c,d,e,f,g,h,i,j,k,l,p,q,r,s,alpha,beta
!       bremsstrahlung fourier coefficients

!  first index labels species: 12C=1, 16O=2, 4He=3
!       Table 3
!       12C
        a(1,0)=0.17946
        a(1,1)=-0.05821
        a(1,2)=-0.01089
        a(1,3)=-0.01147
        a(1,4)=-0.00656
        a(1,5)=-0.00519
        b(1,1)=-0.04969
        b(1,2)=-0.01584
        b(1,3)=-0.00504
        b(1,4)=-0.00281
        c(1)=0.00945
        d(1)=0.34529
        e(1,0)=0.06781
        e(1,1)=-0.00944
        e(1,2)=-0.01289
        e(1,3)=-0.00589
        e(1,4)=-0.00404
        e(1,5)=-0.00330
        f(1,1)=-0.02213
        f(1,2)=-0.01136
        f(1,3)=-0.00467
        f(1,4)=-0.00131
        g(1)=-0.02342
        h(1)=0.24819
!       16O
        a(2,0)=0.20933
        a(2,1)=-0.06740
        a(2,2)=-0.01293
        a(2,3)=-0.01352
        a(2,4)=-0.00776
        a(2,5)=-0.00613
        b(2,1)=-0.05950
        b(2,2)=-0.01837
        b(2,3)=-0.00567
        b(2,4)=-0.00310
        c(2)=0.00952
        d(2)=0.35029
        e(2,0)=0.09304
        e(2,1)=-0.01656
        e(2,2)=-0.01489
        e(2,3)=-0.00778
        e(2,4)=-0.00520
        e(2,5)=-0.00418
        f(2,1)=-0.03076
        f(2,2)=-0.01390
        f(2,3)=-0.00522
        f(2,4)=-0.00161
        g(2)=-0.02513
        h(2)=.27480
!  4He
        a(3, 0 ) =  0.09037
        a(3, 1 ) = -0.03009
        a(3, 2 ) = -0.00564
        a(3, 3 ) = -0.00544
        a(3, 4 ) = -0.00290
        a(3, 5 ) = -0.00224
        b(3, 1 ) = -0.02148
        b(3, 2 ) = -0.00817
        b(3, 3 ) = -0.00300
        b(3, 4 ) = -0.00170
!        b(3, 5 ) =  0.
        c(3)      =  0.00671
        d(3)      =  0.28130
        e(3, 0 ) = -0.02006
        e(3, 1 ) =  0.01790
        e(3, 2 ) = -0.00783
        e(3, 3 ) = -0.00021
        e(3, 4 ) =  0.00024
        e(3, 5 ) = -0.00014
        f(3, 1 ) =  0.00538
        f(3, 2 ) = -0.00175
        f(3, 3 ) = -0.00346
        f(3, 4 ) = -0.00031
!        f(3, 5 ) =  0.
        g(3)      = -0.02199
        h(3)      =  0.17300

!       Table 4
!       12C
        i(1,0)=0.00766
        i(1,1)=-0.00710
        i(1,2)=-0.00028
        i(1,3)=0.00232
        i(1,4)=0.00044
        i(1,5)=0.00158
        j(1,1)=0.02300
        j(1,2)=-0.01078
        j(1,3)=0.00118
        j(1,4)=-0.00089
        k(1)=-0.01259
        l(1)=0.07917
        p(1,0)=-0.00769
        p(1,1)=0.00356
        p(1,2)=-0.00184
        p(1,3)=0.00146
        p(1,4)=0.00031
        p(1,5)=0.00069
        q(1,1)=0.01052
        q(1,2)=-0.00354
        q(1,3)=-0.00014
        q(1,4)=-0.00018
        r(1)=-0.00829
        s(1)=0.05211
!       16O
        i(2,0)=0.00951
        i(2,1)=-0.00838
        i(2,2)=-0.00011
        i(2,3)=0.00244
        i(2,4)=0.00046
        i(2,5)=0.00168
        j(2,1)=0.02455
        j(2,2)=-0.01167
        j(2,3)=0.00131
        j(2,4)=-0.00097
        k(2)=-0.01314
        l(2)=0.08263
        p(2,0)=-0.00700
        p(2,1)=0.00295
        p(2,2)=-0.00184
        p(2,3)=0.00166
        p(2,4)=0.00032
        p(2,5)=0.00082
        q(2,1)=0.01231
        q(2,2)=-0.00445
        q(2,3)=0.00002
        q(2,4)=-0.00026
        r(2)=-0.00921
        s(2)=0.05786
!  4He
        i(3, 0 ) =  0.00192
        i(3, 1 ) = -0.00301
        i(3, 2 ) = -0.00073
        i(3, 3 ) =  0.00182
        i(3, 4 ) =  0.00037
        i(3, 5 ) =  0.00116
        j(3, 1 ) =  0.01706
        j(3, 2 ) = -0.00753
        j(3, 3 ) =  0.00066
        j(3, 4 ) = -0.00060
!        j(3, 5 ) =  0.
        k(3)      = -0.01021
        l(3)      =  0.06417
        p(3, 0)  = -0.01112
        p(3, 1 ) =  0.00603
        p(3, 2 ) = -0.00149
        p(3, 3 ) =  0.00047
        p(3, 4 ) =  0.00040
        p(3, 5 ) =  0.00028
        q(3, 1 ) =  0.00422
        q(3, 2 ) = -0.00009
        q(3, 3 ) = -0.00066
        q(3, 4 ) = -0.00003
!        q(3, 5 ) =  0.
        r(3)      = -0.00561
        s(3)      =  0.03522

!       Table 6
!       12C
        alpha(1,0)=-0.05483
        alpha(1,1)=-0.01946
        alpha(1,2)=1.86310
        alpha(1,3)=-0.78873
        beta(1,0)=-0.06711
        beta(1,1)=0.06859
        beta(1,2)=1.74360
        beta(1,3)=-0.74498
!       16O
        alpha(2,0)=-0.06597
        alpha(2,1)=0.06048
        alpha(2,2)=1.74860
        alpha(2,3)=-0.74308
        beta(2,0)=-0.07356
        beta(2,1)=0.10865
        beta(2,2)=1.70150
        beta(2,3)= -0.73653
!  4He
        alpha(3,0) = -0.07980
        alpha(3,1) =  0.17057
        alpha(3,2) =  1.51980
        alpha(3,3) = -0.61058
        beta(3,0)  = -0.05881
        beta(3,1)  =  0.00165
        beta(3,2)  =  1.82700
        beta(3,3)  = -0.76993

      end subroutine init_brem

      END MODULE ez_vcool


      subroutine dsinit
      implicit none
      integer i
      include 'dssusy.h'
      include 'dsidtag.h'
      include 'dsandwcom.h'
      include 'dsrncom.h'
      include 'dsaccom.h'
      include 'dsprep.h'
      include 'dsascom.h'
      include 'dsversion.h'  ! set DarkSUSY version
      include 'dsdir.h'      ! set DarkSUSY root directory
      include 'dspbcom.h'
      real*8 dsmqpole1loop
      integer kk,ii
c
c... internal fixed-for-ever values go here
c
c  knu=(nue,numu,nutau)   kl=(e,mu,tau)    kqu=(u,c,t)    kqd=(d,s,b)
c  ksqu=(~u1,~c1,~t1,~u2,~c2,~t2)    ksqd=(~d1,~s1,~b1,~d2,~s2,~b2)
      data kse,       ksmu,         kstau /
     &     kse1,kse2, ksmu1,ksmu2,  kstau1,kstau2 /
      data ksu,       ksd,       ksc       /
     &     ksu1,ksu2, ksd1,ksd2, ksc1,ksc2 /
      data kss,       kst,       ksb   /
     &     kss1,kss2, kst1,kst2, ksb1,ksb2 /
      data kn,              kcha/
     &     kn1,kn2,kn3,kn4, kcha1,kcha2/
      data knu,               kl,          kqu,      kqd      /
     &     knue,knumu,knutau, ke,kmu,ktau, ku,kc,kt, kd,ks,kb /
      data ksnu                       /
     &     ksnue,ksnumu,ksnutau,0,0,0 /
      data ksl                                 /
     &     kse1,ksmu1,kstau1,kse2,ksmu2,kstau2 /
      data ksqu                          /
     &     ksu1,ksc1,kst1,ksu2,ksc2,kst2 /
      data ksqd              /
     &     ksd1,kss1,ksb1,ksd2,kss2,ksb2 /
      data pname/'error','nu_e','e','nu_mu','mu','nu_tau','tau','u',
     & 'd','c','s','t','b','gamma','w','z','gluon','h1','h2','a','h+',
     & 's-nu_1','s-l_1','s-l_2','s-nu_2','s-l_3','s-l_4','s-nu_3',
     & 's-l_5','s-l_6','s-qu_1','s-qu_2','s-qd_1','s-qd_2','s-qu_3',
     & 's-qu_4','s-qd_3','s-qd_4','s-qu_5','s-qu_6','s-qd_5','s-qd_6',
     & 'x0_1','x0_2','x0_3','x0_4','x+_1','x+_2','gluino','goldst0',
     & 'goldst+'/
c      data roption/'norun'/
c      data roption/'isasu'/
      data roption/'1loop'/

c...Degrees of freedom used for other coannihilation calculation
      data kdof/ 0,
     &     2,4,2,4,2, 4,12,12,12,12,
     &     12,12,2,6,3, 16,1,1,1,2,
     &     2,2,2,2,2, 2,2,2,2,6,
     &     6,6,6,6,6, 6,6,6,6,6,
     &     6,2,2,2,2, 4,4,16,1,2/

c...Startup
      write(*,*) 
     &  '*********************************************************'
      write(*,*) 
     &  '*** Welcome to DarkSUSY version                       ***'
      write(*,*) '*** ',dsversion,'***'
      write(*,*) 
     &  '*********************************************************'
      write(*,*) ' '
      write(*,*) 'Initializing DarkSUSY...'      



c... mathematical constants
      pi=3.141592653589793238d0
c... standard model charges
      wiso3(knue)   =+0.5d0
      wiso3(ke)     =-0.5d0
      wiso3(knumu)  =+0.5d0
      wiso3(kmu)    =-0.5d0
      wiso3(knutau) =+0.5d0
      wiso3(ktau)   =-0.5d0
      wiso3(ku)     =+0.5d0
      wiso3(kd)     =-0.5d0
      wiso3(kc)     =+0.5d0
      wiso3(ks)     =-0.5d0
      wiso3(kt)     =+0.5d0
      wiso3(kb)     =-0.5d0
      echarg(knue)  =0.d0
      echarg(ke)    =-1.d0
      echarg(knumu) =0.d0
      echarg(kmu)   =-1.d0
      echarg(knutau)=0.d0
      echarg(ktau)  =-1.d0
      echarg(ku)    =+2.d0/3.d0
      echarg(kd)    =-1.d0/3.d0
      echarg(kc)    =+2.d0/3.d0
      echarg(ks)    =-1.d0/3.d0
      echarg(kt)    =+2.d0/3.d0
      echarg(kb)    =-1.d0/3.d0
      ncolor(knue)  =1.d0
      ncolor(ke)    =1.d0
      ncolor(knumu) =1.d0
      ncolor(kmu)   =1.d0
      ncolor(knutau)=1.d0
      ncolor(ktau)  =1.d0
      ncolor(ku)    =3.d0
      ncolor(kd)    =3.d0
      ncolor(kc)    =3.d0
      ncolor(ks)    =3.d0
      ncolor(kt)    =3.d0
      ncolor(kb)    =3.d0
c
c... default values go here
c
c... gauge coupling constants at the z scale
      s2thw = 0.2319d0
      alphem = 0.0078186083d0 ! =1/127.9
c PDG 2002 value for alpha strong at the z scale
c now alph3 is set model by model at twice the lsp mass
      alph3mz = 0.1172d0

c... standard model masses 
      mass(kgamma) =  0.0d0
      mass(kgluon) =  0.0d0
      mass(kz)     = 91.187d0
      mass(kw)     = 80.33d0
      mass(knue)   =  0.0d0
      mass(ke)     =  0.000510999907d0
      mass(knumu)  =  0.0d0
      mass(kmu)    =  0.105658389d0
      mass(knutau) =  0.0d0
      mass(ktau)   =  1.777d0
      mass(ku)     =  0.0056d0
      mass(kd)     =  0.0099d0
c      mass(kc)     =  1.35d0
c for the c quark, we fix as input mc(mc) rather than the pole mass
c PDG 2002
      mcmc = 1.26d0
      mass(ks)     =  0.199d0
c      mass(kb)     =  5.0d0
c for the b quark, we fix as input mb(mb) rather than the pole mass
      mbmb = 4.2d0
c for the top quark, we fix the pole mass
c      call dssettopmass(175.d0)
c PDG 2002
      call dssettopmass(174.3d0)
c... standard model widths
      width(kz) = 2.490d0
      width(kw) = 2.07d0
c... cabibbo-kobayashi-maskawa mixing matrix
      ckms12=0.221d0     ! sin(theta_12)
      ckms23=0.040d0     ! sin(theta_23)
      ckms13=0.0035d0    ! sin(theta_13)
      ckmdelta=0.d0      ! delta no cp violation
c... program switches
c      higloop = 3        ! Carena-Espinosa-Quiros-Wagner
      higloop = 6        ! FeynHiggsFast
      prtlevel = 0
      neuloop = 1
      bsgqcd = 1
      msquarks = 0.d0
      msleptons = 0.d0
      incglue = .true.
      incgaga = .true.
      incgaz = .true.
      idtag = ' '
      luout = 6  ! unit where messages go
c... supersymmetric parameters
      mu=0.d0
      m1=0.d0
      m2=0.d0
      m3=0.d0
      ma=0.d0
      tanbe=0.d0
      do i=1,3
         mass2q(i)=0.d0
         mass2u(i)=0.d0
         mass2d(i)=0.d0
         mass2l(i)=0.d0
         mass2e(i)=0.d0
         asofte(i)=0.d0
         asoftu(i)=0.d0
         asoftd(i)=0.d0
      enddo
c... relic density
      omtype = 1
      omfast = 1
c... set-up defaults for modules
      call dsntset('default')
      call dshmset('default')
      call dspbset('default')
      call dsddset('si','default')
      call dsddset('sd','default')
      call dsepset('default')
      call dsacset('default')
c      call dsgcset('default',0.0d0)
      call dsanset('default')
      call dsgalpropset('default')

      call dsntsunread  ! initialize Sun potential

c...default for new model checks
      dsprepcalled=.false.

c... added by Piero Ullio 
c... set-up fermion family codes needed for coannihilations 
c... numbers are conventional, no physical meaning
      do i=0,50
        ivfam(i)=0
      enddo
      ivfam(knue)=11
      ivfam(ke)=12
      ivfam(ksnue)=11
      ivfam(kse1)=12
      ivfam(kse2)=12
      ivfam(knumu)=21
      ivfam(kmu)=22
      ivfam(ksnumu)=21
      ivfam(ksmu1)=22
      ivfam(ksmu2)=22
      ivfam(knutau)=31
      ivfam(ktau)=32
      ivfam(ksnutau)=31
      ivfam(kstau1)=32
      ivfam(kstau2)=32
      ivfam(ku)=41
      ivfam(kd)=42
      ivfam(ksu1)=41
      ivfam(ksu2)=41
      ivfam(ksd1)=42
      ivfam(ksd2)=42
      ivfam(kc)=51
      ivfam(ks)=52
      ivfam(ksc1)=51
      ivfam(ksc2)=51
      ivfam(kss1)=52
      ivfam(kss2)=52
      ivfam(kt)=61
      ivfam(kb)=62
      ivfam(kst1)=61
      ivfam(kst2)=61
      ivfam(ksb1)=62
      ivfam(ksb2)=62

c...Initialize Bessel routine zeros
      do kk=0,nbesselk
        do ii=1,nzerojk
          storage(kk,ii,1)=0.0d0
        enddo
      enddo

      write(*,*) 'Initialization of DarkSUSY complete.'

      return
      end



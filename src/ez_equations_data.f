      MODULE ez_equations_data
      IMPLICIT NONE
     
      ! composition equations for fractional abundances of H, He, C, O, and Ne
      INTEGER, PARAMETER :: Q_XH=1
      INTEGER, PARAMETER :: Q_XHE=3
      INTEGER, PARAMETER :: Q_XC=4
      INTEGER, PARAMETER :: Q_XO=2
      INTEGER, PARAMETER :: Q_XNE=5
      ! structure equations for P, R, T, L, and M (follow the composition equations in most blocks)
      INTEGER, PARAMETER :: Q_PQ=6
      INTEGER, PARAMETER :: Q_RQ=7
      INTEGER, PARAMETER :: Q_TQ=8
      INTEGER, PARAMETER :: Q_LQ=9
      INTEGER, PARAMETER :: Q_MQ=10
      ! surface boundary conditions for M, P, and T (follow the composition equations in the surface block)
      INTEGER, PARAMETER :: Q_BC_Msurf=6
      INTEGER, PARAMETER :: Q_BC_Psurf=7
      INTEGER, PARAMETER :: Q_BC_Tsurf=8
      ! center boundary conditions for M, L, and R (occur in a separate equation block)
      INTEGER, PARAMETER :: Q_BC_Mcntr=1
      INTEGER, PARAMETER :: Q_BC_Lcntr=2
      INTEGER, PARAMETER :: Q_BC_Rcntr=3
      ! mesh spacing equation
      INTEGER, PARAMETER :: Q_Q_dk=11
      INTEGER, PARAMETER :: Q_MAX=11, NUMQ = 11 ! number of equations
      
      DOUBLE PRECISION :: EQU(NUMQ)
      
      END MODULE ez_equations_data

      MODULE ez_equations
      USE ez_data
      USE ez_equations_data
      USE ez_shell_data
      IMPLICIT NONE 
      
      CONTAINS
      
      SUBROUTINE EQUNS1 ( Q, I, IEE )
      INTEGER, INTENT(IN) :: Q, I, IEE
      
      ! input Q is the equation block number
      !    (N_SURF_shell for surface shell, N_CNTR_shell for center shell, N_CNTR_shell+1 for center boundary conditions)
      ! input I is the variable number that has been modified (or 0 for none)
      ! input IEE identifies which shell holds the modified variable (or 0 if I = 0)
      !    IEE = 1, 2, or 3 when I is nonzero.  
      !
      ! Here are the equations for various values of Q, from the surface to the center.
      ! Q = N_SURF_shell: surface shell.  only one neighbor, so only 1st order equations possible.
         ! 5 1st order composition equations -- XH, XO, XHE, XC, XNE
         ! 3 surface boundary condition equations -- M, P, and T
      ! N_SURF_shell < Q < N_CNTR_shell: interior shells with neighbors on each side so can have 2nd order difference equations
         ! 5 2nd order composition equations -- XH, XO, XHE, XC, XNE
         ! 5 1st order structure equations -- P, R, T, L, M
         ! 1 mesh spacing equation -- equate this shell's spacing to the next outward shell's spacing
      ! Q = N_CNTR_shell: centermost shell.  only one neighbor, so only 1st order equations possible.
         ! 5 1st order composition equations -- XH, XO, XHE, XC, XNE
         ! 5 1st order structure equations -- P, R, T, L, M
         ! 1 mesh spacing equation -- equate this shell's spacing to the next outward shell's spacing
      ! Q = N_CNTR_shell+1: center boundary conditions.
         ! 3 equations -- M, R, and L
      !
      ! The shells for N_SURF_shell < Q <= N_CNTR_shell all have 11 equations to match their 11 variables.
      ! The surface shell only has 8 equations, but the 3 center BCs give us the required 11 equations
      ! to match the 11 variables of the surface shell.
      !
      ! In more detail, here's the dependencies of equations on variables spelled out explicitly.
      ! There are N_CNTR_shell blocks of variables, with 11 variables in each block.
      ! There are N_CNTR_shell+1 blocks of equations, with 8 in the first, 3 in the last, and 11 in each of the others.
      ! The first equation block uses variable blocks 1 and 2.
      ! The last equation block uses variable block N_CNTR_shell.
      ! The next-to-last equation block uses variable blocks N_CNTR_shell-1 and N_CNTR_shell.
      ! The interior equation blocks, 1 < Q < N_CNTR_shell, use variable blocks Q-1, Q, and Q+1.
      ! The Jacobian matrix for the problem gives the partial derivatives of the equations with respect to the variables.
      ! Because the equations are set up to depend only on neighboring shells, the Jacobian is a band diagonal matrix.
      !
      ! Notice that there are N_CNTR_shell different mesh spacing variables (one for each shell), but only
      ! N_CNTR_shell-1 mesh spacing equations (none for N_SURF_shell).  The equations don't do anything but link all the spacings
      ! together by equation one to the next.  The missing equation allows the values to float freely,
      ! while the N_CNTR_shell-1 equations force the spacing to be uniform.
      ! The mesh spacing, dQ/dk, is just the derivative of the mesh function, Q, with respect to shell number, k.
      ! The mesh function is a function of the local state that is designed to vary monotonically with significant
      ! variables.  To get a uniform mesh spacing, the system must "move" the mesh points to make the mesh function
      ! increase by the same amount from one shell to the next.  Since the mesh function depends on "important"
      ! variables, this means that the mesh will get denser in areas where there are rapid changes in these variables.
      ! This is the basic idea behind the automatic adjustment of the mesh to the ongoing evolution of the star.
      !
      ! Since we have a mesh spacing variable at each shell, we need to do a little extra to ensure that the spacing
      ! does in fact stay uniform.  There are equations that demand this, but the equations only get approximate
      ! solutions.  There could be additive errors such the the mesh spacing at the center of the star was different
      ! than the spacing at the surface.  To prevent this from happening, the hand of the programmer reaches in
      ! and forces the uniformity as part of the solution process.  The linear equation solver is encouraged to make
      ! the mesh spacing uniform by putting a large error penalty on any variations.  Then the proposed solution is
      ! "edited" by a post-processing step to set all of the spacing changes to the average change suggested by the
      ! equation solver.  This way, since the spacing starts uniform, it stays uniform as the star evolves. 
      !
      ! Note that EQUNS1 doesn't call any other subroutine to compute anything.  All it does is use function values
      ! previously computed by FUNCS1 and saved for use here.
      !
      DOUBLE PRECISION :: S12, S23, WTA, WTB, WTC, WTD
      INTEGER :: IC, IP, IPP, CQ, PQ, PPQ
      DOUBLE PRECISION :: PS, VX
      PS(VX) = 0.5D0*(VX + DABS(VX)) ! statement function PS(VX) = VX if VX >= 0; = 0 if VX < 0
      VX = 0 ! for the Absoft compiler.
      ! VAR(IC), VAR(IP), and VAR(IPP) hold values at current, previous and previous-previous meshpoints.
      ! VAR(IC) is nearest to surface; VAR(IPP) is nearest to center.
      EQU=0D0
      IF (Q .NE. N_SURF_SHELL) THEN
         CQ = Q-1; PQ = Q; PPQ = Q+1
      ELSE
         CQ = Q; PQ = Q+1
      END IF
      IC = FN_VALUE; IP = FN_VALUE; IPP = FN_VALUE ! FN_VALUE means use the unmodified value of the functions
      IF (IEE .EQ. 1) THEN ! variable in neighbor shell toward center modified
         IPP = I 
      ELSE IF (IEE .EQ. 2) THEN ! variable in middle shell modified
         IP = I
      ELSE IF (IEE .EQ. 3) THEN ! variable in neighbor shell toward surface modified
         IC = I
      END IF
      ! the composition equations
      IF ( Q .EQ. N_SURF_shell ) THEN ! surface shell
         ! the composition is unchanged between next-to-surface and surface
         ! IC values are for the surface shell.  IP are for the next shell below it.
         EQU(Q_XH) = FN(IC,CQ,F_XH) - FN(IP,PQ,F_XH)
         EQU(Q_XO) = FN(IC,CQ,F_XO) - FN(IP,PQ,F_XO)
         EQU(Q_XHE) = FN(IC,CQ,F_XHE) - FN(IP,PQ,F_XHE)
         EQU(Q_XC) = FN(IC,CQ,F_XC) - FN(IP,PQ,F_XC)
         EQU(Q_XNE) = FN(IC,CQ,F_XNE) - FN(IP,PQ,F_XNE)
      ELSEIF ( Q .EQ. N_CNTR_shell ) THEN ! centermost shell
         ! replace the second-order equation by a first-order one for the center
         ! IP values are for the centermost shell.  IC are for the next shell above it.
         S23 = -0.5D0*(FN(IP,PQ,F_SG)+FN(IC,CQ,F_SG)) 
         EQU(Q_XH) = S23*(FN(IC,CQ,F_XH) - FN(IP,PQ,F_XH)) + FN(IP,PQ,F_XH_dt)
         EQU(Q_XO) = S23*(FN(IC,CQ,F_XO) - FN(IP,PQ,F_XO)) + FN(IP,PQ,F_XO_dt)
         EQU(Q_XHE) = S23*(FN(IC,CQ,F_XHE) - FN(IP,PQ,F_XHE)) + FN(IP,PQ,F_XHE_dt)
         EQU(Q_XC) = S23*(FN(IC,CQ,F_XC) - FN(IP,PQ,F_XC)) + FN(IP,PQ,F_XC_dt)
         EQU(Q_XNE) = S23*(FN(IC,CQ,F_XNE) - FN(IP,PQ,F_XNE)) + FN(IP,PQ,F_XNE_dt)
      ELSE IF ( Q .LT. N_CNTR_shell ) THEN ! interior shells (have 2 neighbors, and 2nd order equations)
         ! for element i, the composition difference equation is the following:
         ! sigma(k+1/2)*(Xi(k+1)-Xi(k)) - sigma(k-1/2)*(Xi(k)-Xi(k-1)) =
         !        (d_Xi(k)/dt + Rnuc(k))*d_M/d_k +
         !        (Xi(k+1)-Xi(k))*[-d_M(k)/dt] - (Xi(k)-Xi(k-1))*[d_M(k-1)/dt]
         ! (sigma is the diffusion coefficient for convective mixing).
         ! When everything is moved to the right hand side, this becomes
         ! EQU = sigma(k+1/2)*(Xi(k+1)-Xi(k)) - sigma(k-1/2)*(Xi(k)-Xi(k-1)) -
         !       (d_Xi(k)/dt + Rnuc(k))*d_M/d_k -
         !       (Xi(k+1)-Xi(k))*[-d_M(k)/dt] + (Xi(k)-Xi(k-1))*[d_M(k-1)/dt]
         ! which is the same as
         ! EQU = (sigma(k+1/2)-[-d_M(k)/dt])*(Xi(k+1)-Xi(k)) -
         !       (sigma(k-1/2)-[d_M(k-1)/dt])*(Xi(k)-Xi(k-1)) -
         !       (d_Xi(k)/dt + Rnuc(k))*d_M/d_k
         ! The final term has already been computed by FUNCS1 and saved as XiT.
         ! sigma(k+1/2) = 0.5*(sigma(k)+sigma(k+1)) and similar for sigma(k-1/2).
         ! replace k by 2 to get the actual equation for this meshpoint
         ! recall that [x] = x for x >=0, = 0 otherwise, and is implemented by PS(x).
         ! S12 is (sigma(k-1/2)-[d_M(k-1)/dt])
         ! S23 is (sigma(k+1/2)-[-d_M(k)/dt])
         ! note that the S12 & S23 values are the same for all elements.
         ! the equation finally becomes
         ! EQU = S23*(Xi(IC)-Xi(IP)) - S12*(Xi(IP)-Xi(IPP)) - XiT(IP)
         ! IC values are for the neighbor nearer the surface.
         ! IPPs are for the neighbor nearer the center.
         ! IPs are for the shell in the middle.
         S12 = 0.5D0*(FN(IPP,PPQ,F_SG) + FN(IP,PQ,F_SG)) - PS(-FN(IP,PQ,F_M_dt))
         S23 = 0.5D0*(FN(IP,PQ,F_SG) + FN(IC,CQ,F_SG)) - PS(FN(IP,PQ,F_M_dt))
         EQU(Q_XH) = S23*(FN(IC,CQ,F_XH) - FN(IP,PQ,F_XH)) - S12*(FN(IP,PQ,F_XH) - FN(IPP,PPQ,F_XH)) - FN(IP,PQ,F_XH_dt)
         EQU(Q_XO) = S23*(FN(IC,CQ,F_XO) - FN(IP,PQ,F_XO)) - S12*(FN(IP,PQ,F_XO) - FN(IPP,PPQ,F_XO)) - FN(IP,PQ,F_XO_dt)
         EQU(Q_XHE) = S23*(FN(IC,CQ,F_XHE) - FN(IP,PQ,F_XHE)) - S12*(FN(IP,PQ,F_XHE) - FN(IPP,PPQ,F_XHE)) - FN(IP,PQ,F_XHE_dt)
         EQU(Q_XC) = S23*(FN(IC,CQ,F_XC) - FN(IP,PQ,F_XC)) - S12*(FN(IP,PQ,F_XC) - FN(IPP,PPQ,F_XC)) - FN(IP,PQ,F_XC_dt)
         EQU(Q_XNE) = S23*(FN(IC,CQ,F_XNE) - FN(IP,PQ,F_XNE)) - S12*(FN(IP,PQ,F_XNE) - FN(IPP,PPQ,F_XNE)) - FN(IP,PQ,F_XNE_dt)
      END IF
      ! the other equations
      IF ( Q .EQ. N_SURF_shell ) THEN ! surface boundary conditions
         ! just copy the values that were computed in FUNCS1 for mass loss, pressure, and temperature.
         EQU(Q_BC_Msurf) =  FN(IC,CQ,F_BC_Msurf)  ! this equation makes Mdot at surface match mass loss/gain if any.
         EQU(Q_BC_Psurf) =  FN(IC,CQ,F_BC_Psurf)  ! this is the photospheric pressure based on local gravity and opacity.
         EQU(Q_BC_Tsurf) =  FN(IC,CQ,F_BC_Tsurf)  ! this makes the temperature match the effective temperature from L and R.
      ELSEIF ( Q .GT. N_CNTR_shell ) THEN ! center boundary conditions
         ! The difference equations are designed to have the correct limiting behavior as r -> 0, 
         ! so the following are based directly on the interior versions described above.
         ! The variables for mass, luminosity, and radius (MQ, LQ, RQ) are all zero at the center.
         ! Since k increases inward, these all have negative d_Vx / d_k.
         ! So at the first meshpoint out from the center, Vx = -(d_Vx / d_k),
         ! and EQU = Vx + (d_Vx / d_k).
         EQU(Q_BC_Mcntr) = FN(IC,CQ,F_MQ) + FN(IC,CQ,F_MQ_dk) ! mass BC
         EQU(Q_BC_Lcntr) = FN(IC,CQ,F_LQ) + FN(IC,CQ,F_LQ1_dk) + FN(IC,CQ,F_LQ2_dk)*FN(IC,CQ,F_M_dt) ! luminosity BC
         EQU(Q_BC_Rcntr) = FN(IC,CQ,F_RQ) + FN(IC,CQ,F_RQ_dk) ! radius BC
      ELSE ! structure equations for all the shells except the centermost.  1st order equations using neighbor nearer to center.
         EQU(Q_Q_dk) = FN(IC,CQ,F_Q_dk) - FN(IP,PQ,F_Q_dk) ! the mesh spacing, Q_dk, is the same for all shells
         WTA = -0.5D0
         WTB = 0.5D0*(FN(IP,PQ,F_WT) + FN(IC,CQ,F_WT))
         WTC = WTA*WTB/(1D0 + WTB)
         WTD = -1 - WTC
         ! WTA, WTB, WTC, WTD determine the relative weighting of adjacent points for the shell average 
         EQU(Q_PQ) = FN(IC,CQ,F_PQ) - FN(IP,PQ,F_PQ) - WTC*FN(IP,PQ,F_PQ_dk) - WTD*FN(IC,CQ,F_PQ_dk)
            ! the pressure equation in terms of PQ rather than P
            ! PQ is C4 * log(P) + C5 * log(P + C9).   C4 = .05, C5 = 0.15
            ! The equations might have been ln(P(k+1))-ln(P(k)) = (d_ln(P) / d_k)(k+1/2), but instead of ln(P), 
            ! use PQ which is also the pressure term for the mesh spacing function MESH_Q.
            ! So the pressure equation is actually PQs(k+1)-PQs(k) = (d_PQ / d_k)(k+1/2) = PQ_dks(k+1/2).
            ! Instead of PQ_dks(k+1/2) = 0.5*(PQ_dks(k+1)+PQ_dks(k)), use WTC and WTD as weighting coefficients.
            ! The simple equation would be EQU = PQs(IC)-PQs(IP)+0.5*PQ_dks(IC)+0.5*PQ_dks(IP), but
            ! instead we have EQU = PQs(IC)-PQs(IP)-WTC*PQ_dks(IP)-WTD*PQ_dks(IC)
            ! with WTA=-0.5, WTB=0.5*(WTs(IP)+WTs(IC)); WTC=-0.5*WTB/(1 + WTB); WTD = -1 - WTC.
            ! if WTB is very large, we get WTC and WTD both -> -0.5
            ! if WTB is very small, we get WTC -> 0 and WTD -> -1 so PQ_dks(k+1/2) -> PQ_dks(k+1)
            ! so the average is biased towards the current mesh point.
         ! the M, T, and R equations are written in terms of computed functions of M, T, and R as follows:
         EQU(Q_MQ) = FN(IC,CQ,F_MQ) - FN(IP,PQ,F_MQ) - WTA*(FN(IC,CQ,F_MQ_dk) + FN(IP,PQ,F_MQ_dk))
            ! the mass equation in terms of VM rather than M.
            ! MQ is log(XMC / (XMC + M^(2/3))), XMC = C6 * MC^(2/3), C6 = .02,
            ! MC = mass in core to HPC, HPC = sqrt(P/(G*RHO^2)) at core
         EQU(Q_TQ) = FN(IC,CQ,F_TQ) - FN(IP,PQ,F_TQ) - WTC*FN(IP,PQ,F_TQ_dk) - WTD*FN(IC,CQ,F_TQ_dk)
            ! the temperature equation in terms of VT rather than T.
            ! TQ is C7 * log(T/(T+C10)).  C7 = .45
         EQU(Q_RQ) = FN(IC,CQ,F_RQ) - FN(IP,PQ,F_RQ) - WTD*FN(IP,PQ,F_RQ_dk) - WTC*FN(IC,CQ,F_RQ_dk)
            ! the radius equation in terms of VR rather than R.
            ! RQ is -C3 * log(R^2/C8 + 1).   C3 = .05, C8 = .0001
            ! RQ is R term for mesh spacing function MESH_Q and varies like R^2
            ! so the radius equation is RQs(k+1)-RQs(k) = (d_VR / d_k)(k+1/2) = RQ_dks(k+1/2)
         EQU(Q_LQ) = FN(IC,CQ,F_LQ) - FN(IP,PQ,F_LQ) - WTD*FN(IP,PQ,F_LQ1_dk) - WTC*FN(IC,CQ,F_LQ1_dk)
         EQU(Q_LQ) = EQU(Q_LQ) - FN(IP,PQ,F_LQ2_dk)*PS(-FN(IP,PQ,F_M_dt)) + FN(IC,CQ,F_LQ2_dk)*PS(FN(IC,CQ,F_M_dt))
            ! the luminosity equation in terms of VL rather than L.
            ! LQ is L * VML, VML is weight term that is about sqrt(2) at center and drops slowly to 1 at surface
      END IF
      END SUBROUTINE EQUNS1
      
      END MODULE ez_equations


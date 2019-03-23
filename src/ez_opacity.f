      MODULE ez_opacity
      USE ez_opacity_data
      IMPLICIT NONE
      
      CONTAINS

      DOUBLE PRECISION FUNCTION Opacity(TF, FR, XH1, XHE1)
      ! TF = log10(temperature)
      ! FR = log10(density)
      DOUBLE PRECISION, INTENT(IN) :: TF, FR, XH1, XHE1
      INTEGER :: JX
      DOUBLE PRECISION :: FKL, FKH, XF, XU, XT
      DOUBLE PRECISION :: DELT, DR, XTF, XFR, F1, F2, F3, F4
      INTEGER :: I1, I2
      ! Opacity tables from Alexander & Ferguson (1994; molecular), Itoh (1983;
      ! electron conduction) and Iglesias & Rogers (1992; the rest)
      XF = XH1 + XH1 + XHE1 + 1D0
      JX = MIN(KCSX,2)
      JX = JX - 2
      DO
         JX = JX + 1
         IF ( XF.LT.CSX(JX + 1) .AND. JX.LT.KCSX - 1 ) CYCLE
         IF ( XF.LT.CSX(JX) .OR. JX.LE.1 ) EXIT
         JX = JX - 2
      END DO
      XT = (XF-CSX(JX))/(CSX(JX+1)-CSX(JX)) ! XT is fractional distance from XF to CSX(JX)
      XU = 1D0 - XT ! XU is the fractional distance from XF to CSX(JX+1)
      ! Calculate a bicubic spline interpolation fit for the temperature
      ! and density opacity fit.  Do not stop if the input lies outside the
      ! array but rather use the value at the nearest edge-point.
      XTF = DMAX1(DMIN1(TF, TFM(MT)), TFM(1))
      XFR = DMAX1(DMIN1(FR, FRM(MR)), FRM(1))
      ! point (XTF,XFR) is nearest possible to target while still being in array range
      ! Find interval in which target point lies.
      I1 = 1 + (MT - 1)*(XTF - TFM(1))/(TFM(MT) - TFM(1))
      I2 = 1 + (MR - 1)*(XFR - FRM(1))/(FRM(MR) - FRM(1))
      DELT = TF - TFM(I1)
      DR = FR - FRM(I2)
      ! Evaluate the splines.
      F1 = DR*(FSPL(4,2,I1,I2,JX)+DR*(FSPL(4,3,I1,I2,JX)+DR*FSPL(4,4,I1,I2,JX)))
      F2 = DELT*(FSPL(4,1,I1,I2,JX)+F1)
      F3 = FSPL(3,1,I1,I2,JX)+DR*(FSPL(3,2,I1,I2,JX)+DR*(FSPL(3,3,I1,I2,JX)+DR*FSPL(3,4,I1,I2,JX)))+F2
      F4 = FSPL(2,1,I1,I2,JX)+DR*(FSPL(2,2,I1,I2,JX)+DR*(FSPL(2,3,I1,I2,JX)+DR*FSPL(2,4,I1,I2,JX)))+DELT*F3
      FKL = FSPL(1,1,I1,I2,JX)+DR*(FSPL(1,2,I1,I2,JX)+DR*(FSPL(1,3,I1,I2,JX)+DR*FSPL(1,4,I1,I2,JX)))+DELT*F4
      F1 = DR*(FSPL(4,2,I1,I2,JX+1)+DR*(FSPL(4,3,I1,I2,JX+1)+DR*FSPL(4,4,I1,I2,JX+1)))
      F2 = DELT*(FSPL(4,1,I1,I2,JX+1)+F1)
      F3 = FSPL(3,1,I1,I2,JX+1)+DR*(FSPL(3,2,I1,I2,JX+1)+DR*(FSPL(3,3,I1,I2,JX+1)+DR*FSPL(3,4,I1,I2,JX+1)))+F2
      F4 = FSPL(2,1,I1,I2,JX+1)+DR*(FSPL(2,2,I1,I2,JX+1)+DR*(FSPL(2,3,I1,I2,JX+1)+DR*FSPL(2,4,I1,I2,JX+1)))+DELT*F3
      FKH = FSPL(1,1,I1,I2,JX+1)+DR*(FSPL(1,2,I1,I2,JX+1)+DR*(FSPL(1,3,I1,I2,JX+1)+DR*FSPL(1,4,I1,I2,JX+1)))+DELT*F4
      Opacity = XT*10D0**FKH + XU*10D0**FKL ! interpolated opacity
      END FUNCTION Opacity

      SUBROUTINE OPSPLN
      DOUBLE PRECISION MAT(4,127)
      INTEGER JT, JR, JX, IT, IR, IQ, JQ, IC
      ! Calculate a bicubic spline interpolation for the temperature and density dependent opacity.
      ! This is based on routines from Princeton Splines, D. McCune found by L. Dray.
      DO JT = 1, MT
         TFM(JT) = 2.95D0 + 0.05D0*DFLOAT(JT) ! temperatures
      END DO
      DO JR = 1, MR
         FRM(JR) = 0.25D0*DFLOAT(JR) - 12.25D0 ! densities
      END DO
      DO JX = 1, 10 ! JX is the sub-composition
         DO JQ = 1, MT ! copy the opacity tables from CS
            DO IQ = 1, MR; FSPL(1, 1, JQ, IQ, JX) = CS(IQ, JQ, JX); END DO
         END DO
         ! Construct splines in the T direction.
         DO IR = 1, MR
            DO IT = 1, MT; MAT(1, IT) = FSPL(1, 1, IT, IR, JX); END DO
            CALL SPLINE ( MT, TFM, MAT ) ! SPLINE puts the data in MAT
            DO IT = 1, MT - 1
               FSPL(2, 1, IT, IR, JX) = MAT(2, IT)
               FSPL(3, 1, IT, IR, JX) = MAT(3, IT)
               FSPL(4, 1, IT, IR, JX) = MAT(4, IT)
            END DO
         END DO
         ! Construct splines in the rho direction.
         DO IT = 1, MT - 1
            ! Construct splines for each T coeff
            DO IC = 1, 4
               DO IR = 1, MR; MAT(1, IR) = FSPL(IC, 1, IT, IR, JX); END DO
               MAT(2, 1) = 0D0;  MAT(3, 1) = 0D0
               MAT(2, MR) = 0D0; MAT(3, MR) = 0D0
               CALL SPLINE ( MR, FRM, MAT )
               DO IR = 1, MR - 1
                  FSPL(IC, 2, IT, IR, JX) = MAT(2, IR)
                  FSPL(IC, 3, IT, IR, JX) = MAT(3, IR)
                  FSPL(IC, 4, IT, IR, JX) = MAT(4, IR)
               END DO
            END DO
         END DO
      END DO
      END SUBROUTINE OPSPLN

      SUBROUTINE SPLINE ( K, X, F )
      ! Calculate the coefficients of a 1-D cubic spline:
      ! Forsythe, Malcolm, Moler, Computer Methods for Mathematical
      ! Computations, Prentice-Hall, 1977, p.76
      ! input K         number of points
      ! input X         array of points
      ! in & out F      matrix of spline coefficients
      INTEGER, INTENT(IN) :: K
      DOUBLE PRECISION, INTENT(IN) :: X(*)
      DOUBLE PRECISION, INTENT(INOUT) :: F(4,*)
      INTEGER I, IB
      DOUBLE PRECISION T
      DO I = 2, 4; F(I,K) = 0D0; END DO
      ! Set up a tridiagonal system for A*y=B where y(i) are the second
      ! derivatives at the knots.
      ! f(2,i) are the diagonal elements of A
      ! f(4,i) are the off-diagonal elements of A
      ! f(3,i) are the B elements/3, and will become c/3 upon solution
      F(4,1) = X(2)-X(1)
      F(3,2) = (F(1,2) - F(1,1))/F(4,1)
      DO I = 2, K - 1
         F(4,I) = X(I+1) - X(I)
         F(2,I) = 2D0*(F(4,I-1) + F(4,I))
         F(3,I+1) = (F(1,I+1) - F(1,I))/F(4,I)
         F(3,I) = F(3,I+1) - F(3,I)
      END DO
      ! Boundaries.
      F(2,2) = F(4,1) + 2D0*F(4,2)
      F(3,2) = F(3,2)*F(4,2)/(F(4,1) + F(4,2))
      F(2,K-1) = 2D0*F(4,K-2) + F(4,K-1)
      F(3,K-1) = F(3,K-1)*F(4,K-2)/(F(4,K-1) + F(4,K-2))
      ! Forward elimination.
      T = F(4,2)/F(2,2)
      F(2,3) = F(2,3) - T*(F(4,2) - F(4,1))
      F(3,3) = F(3,3) - T*F(3,2)
      DO I = 4, K - 2
         T = F(4,I-1)/F(2,I-1)
         F(2,I) = F(2,I)-T*F(4,I-1)
         F(3,I) = F(3,I)-T*F(3,I-1)
      END DO
      T = (F(4,K-1) - F(4,K-1))/F(2,K-2)
      F(2,K-1) = F(2,K-1) - T*F(4,K-2)
      F(3,K-1) = F(3,K-1) - T*F(3,K-2)
      ! Back substitution.
      F(3,K-1) = F(3,K-1)/F(2,K-1)
      DO IB = 1, K - 4
         I = K - 1 - IB
         F(3,I) = (F(3,I) - F(4,I)*F(3,I+1))/F(2,I)
      END DO
      F(3,2) = (F(3,2) - (F(4,2) - F(4,1))*F(3,3))/F(2,2)
      ! Reset d array to step size.
      F(4,1) = X(2) - X(1)
      F(4,K-1) = X(K) - X(K-1)
      ! Set f(3,1) for not-a-knot.
      F(3,1) = (F(3,2)*(F(4,1) + F(4,2)) - F(3,3)*F(4,1))/F(4,2)
      F(3,K) = F(3,K-1) + (F(3,K-1) - F(3,K-2))*F(4,K-1)/F(4,K-2)
      ! Compute the polynomial coefficients.
      DO I = 1, K - 1
         F(2,I) = (F(1,I+1) - F(1,I))/F(4,I) - F(4,I)*(F(3,I+1) + 2D0*F(3,I))
         F(4,I) = (F(3,I+1) - F(3,I))/F(4,I)
         F(3,I) = 3D0*F(3,I)
         F(4,I) = F(4,I)
      END DO
      END SUBROUTINE SPLINE

      END MODULE ez_opacity


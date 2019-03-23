      MODULE ez_utils
      IMPLICIT NONE     
      CONTAINS

        DOUBLE PRECISION FUNCTION FIND0(XX1,YY1,XX2,YY2)
        DOUBLE PRECISION, INTENT(IN) :: XX1,YY1,XX2,YY2
        DOUBLE PRECISION :: a, b
        a = (XX1*YY2)-(XX2*YY1)
        b = YY2-YY1
        IF ((abs(a) .GT. abs(b)*1D30) .AND. ((YY1 .GE. 0D0 .AND. YY2 .LE. 0D0) .OR. (YY1 .LE. 0D0 .AND. YY2 .GE. 0D0))) THEN
            FIND0 = 0.5D0*(XX1+XX2)
        ELSE
            FIND0 = a/b
        END IF
        END FUNCTION FIND0
        
      SUBROUTINE CS_EVAL(ISX,ISY,II,xvalue,fvalue)
      ! cubic spline interpolation to evaluate a variable between meshpoints
      ! ISX and ISY are variable indices like SX_M or SX_R
      ! II is meshpoint near where want to do the evaluation
      ! xvalue is the value of the ISX parameter where want to evalutate the ISY parameter
      ! fvalue is the estimated value of ISY at location xvalue for ISX
      USE star_extras
      INTEGER, INTENT(IN) :: ISX, ISY, II
      DOUBLE PRECISION, INTENT(IN) :: xvalue
      DOUBLE PRECISION, INTENT(OUT) :: fvalue
      INTEGER I
      DOUBLE PRECISION CS_X(5), CS_A(5), CS_B(5), CS_C(5), CS_D(5) ! arrays for cubic spline interpolations
      I = MAX(II,SX_CNTR+2);
      CS_X(1) = SX(ISX,I-2); CS_X(2) = SX(ISX,I-1); CS_X(3) = SX(ISX,I);
      CS_X(4) = SX(ISX,I+1); CS_X(5) = SX(ISX,I+2) ! x values
      CS_A(1) = SX(ISY,I-2); CS_A(2) = SX(ISY,I-1); CS_A(3) = SX(ISY,I);
      CS_A(4) = SX(ISY,I+1); CS_A(5) = SX(ISY,I+2) ! function values
      CALL cspline_init(CS_X,CS_A,CS_B,CS_C,CS_D,5)
      CALL cspline_eval(fvalue,xvalue,CS_X,CS_A,CS_B,CS_C,CS_D,5) ! estimate value of function at xvalue
      END SUBROUTINE CS_EVAL
      
      SUBROUTINE CS_FINDMAX(ISX,ISY,II,xvalue,fvalue) ! cubic spline interpolation to find a max between meshpoints
      ! get max from near meshpoint II index in SX
      ! ISX and ISY are same as in CS_EVAL
      ! xvalue is the ISX location of the estimated max
      ! fvalue is the ISY value of the estimated max
      USE star_extras
      INTEGER, INTENT(IN) :: ISX, ISY, II
      DOUBLE PRECISION, INTENT(OUT) :: xvalue, fvalue 
      INTEGER I
      DOUBLE PRECISION xvalue1, fvalue1
      DOUBLE PRECISION CS_X(5), CS_A(5), CS_B(5), CS_C(5), CS_D(5) ! for cubic spline interpolations
      I = MAX(II,SX_CNTR+2);
      CS_X(1) = SX(ISX,I-2); CS_X(2) = SX(ISX,I-1); CS_X(3) = SX(ISX,I);
      CS_X(4) = SX(ISX,I+1); CS_X(5) = SX(ISX,I+2) ! x values
      CS_A(1) = SX(ISY,I-2); CS_A(2) = SX(ISY,I-1); CS_A(3) = SX(ISY,I);
      CS_A(4) = SX(ISY,I+1); CS_A(5) = SX(ISY,I+2) ! function values
      CALL cspline_init(CS_X,CS_A,CS_B,CS_C,CS_D,5)
      fvalue = SX(ISY,II); xvalue = SX(ISX,II) ! use meshpoint values if fail to find a max
      DO I=1,3 ! check the neighborhood
         CALL cspline_extreme(xvalue1,fvalue1,4-I,CS_X,CS_A,CS_B,CS_C,CS_D,5)
         IF (xvalue1 .GE. CS_X(1)) THEN
            xvalue = xvalue1; fvalue = fvalue1; EXIT
         END IF
      END DO
      END SUBROUTINE CS_FINDMAX
      
!     Cubic spline interpolation.

      SUBROUTINE cspline_init(x,a,b,c,d,n)
      ! initialize the second derivatives needed for cubic spline interpolation.
      INTEGER, INTENT(in) :: n
      DOUBLE PRECISION, INTENT(in) :: x(n), a(n)
      DOUBLE PRECISION, INTENT(out) :: b(n), c(n), d(n)
      ! locals
      DOUBLE PRECISION :: h(500), alpha(500), mu(500), z(500), l(500)
      INTEGER :: i, j
      IF (n .GT. 500) THEN
         stop 'arrays too big for cspline_init'
      END IF
      DO j=1,n-1
         h(j) = x(j+1)-x(j)
      END DO
      DO j=2,n-1
         alpha(j) = 3D0*(a(j+1)*h(j-1) - a(j)*(x(j+1)-x(j-1)) + a(j-1)*h(j))/(h(j-1)*h(j))
      END DO
      l(1) = 1D0
      mu(1) = 0D0
      z(1) = 0D0
      DO j=2,n-1
         l(j) = 2D0*(x(j+1)-x(j-1)) - h(j-1)*mu(j-1)
         mu(j) = h(j)/l(j)
         z(j) = (alpha(j)-h(j-1)*z(j-1))/l(j)
      END DO
      l(n) = 1D0
      z(n) = 0D0
      c(n) = 0D0
      DO i=1,n-1
         j = n-i
         c(j) = z(j) - mu(j)*c(j+1)
         b(j) = (a(j+1)-a(j))/h(j) - h(j)*(c(j+1)+2D0*c(j))/3D0
         d(j) = (c(j+1)-c(j))/(3D0*h(j))
      END DO
      END SUBROUTINE
      
      SUBROUTINE cspline_eval(fvalue,xvalue,x,a,b,c,d,n)
      INTEGER, INTENT(in) :: n
      DOUBLE PRECISION, INTENT(in) :: x(n), a(n), b(n), c(n), d(n), xvalue
      DOUBLE PRECISION, INTENT(out) :: fvalue
      INTEGER :: j
      DOUBLE PRECISION :: dx
      IF (xvalue .le. x(1)) THEN
         j = 1
      ELSE IF (xvalue .ge. x(n)) THEN
         j = n
      ELSE      ! for heavy use, improve this search and save j from last call as first guess.
         DO j=1,n-1 ! find j s.t. xvalue is between x(j) and x(j+1)
            IF (xvalue .ge. x(j) .and. xvalue .lt. x(j+1)) exit
         END DO
      END IF
      ! spline approximation to f(x) at x between x(j) and x(j+1) is
      !     S(x) == a(j) + b(j)*(x-x(j)) + c(j)*(x-x(j))^2 + d(j)*(x-x(j))^3
      dx = xvalue-x(j)
      fvalue = a(j) + dx*(b(j) + dx*(c(j) + dx*d(j)))
      END SUBROUTINE
      
      SUBROUTINE cspline_extreme(xe,fvalue,j,x,a,b,c,d,n)
      INTEGER, INTENT(in) :: j, n
      DOUBLE PRECISION, INTENT(in) :: x(n), a(n), b(n), c(n), d(n)
      DOUBLE PRECISION, INTENT(out) :: fvalue, xe
      ! locals
      DOUBLE PRECISION dx, dx1, dx2, c2, bd3
      IF (j .LT. 1 .OR. j .GE. n) THEN
         stop 'j out of range for cspline_extreme'
      END IF
      ! derivative is S'(x) == b(j) + 2*c(j)*(x-x(j)) + 3*d(j)*(x-x(j))^2
      c2 = c(j)*c(j); bd3 = 3D0*b(j)*d(j)
      IF (bd3 .GT. c2) THEN
         xe = x(1)-1D0; return
      END IF
      dx1 = (-c(j)+sqrt(c2-bd3))/(3D0*d(j))
      dx2 = (-c(j)-sqrt(c2-bd3))/(3D0*d(j))
      IF (dx1 .GE. 0 .AND. dx2 .GE. 0) THEN
         dx = min(dx1,dx2)
      ELSE IF (dx1 .LE. 0 .AND. dx2 .LE. 0) THEN
         dx = max(dx1,dx2)
      ELSE IF (dx1 .GT. dx2) THEN
         dx = dx1
      ELSE
         dx = dx2
      END IF
      xe = x(j)+dx
      IF (xe .LT. x(1) .OR. xe .GT. x(n)) THEN
         xe = x(1)-1D0; RETURN
      END IF
      fvalue = a(j) + dx*(b(j) + dx*(c(j) + dx*d(j)))
      END SUBROUTINE
      
      END MODULE ez_utils

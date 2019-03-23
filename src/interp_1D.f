      module interp_1D

      implicit none
      private
      public :: mkpmcub, mkpmcub_db
      contains
      
      
      
      subroutine minmod_db(z, nx, x, y)
         double precision, intent(OUT) :: z(nx)    
         integer, intent(IN) :: nx       ! length of vectors
         double precision, intent(IN) :: x(nx), y(nx)       
         z(1:nx) = 0.5 * (sign(1d0,x(1:nx)) + sign(1d0,y(1:nx))) * min(abs(x(1:nx)), abs(y(1:nx)))    
      end subroutine minmod_db
      
      
      subroutine median_db(z, nx, x1, x2, x3)
         double precision, intent(OUT) :: z(nx)    
         integer, intent(IN) :: nx       ! length of vectors
         double precision, intent(IN) :: x1(nx), x2(nx), x3(nx)
         call minmod_db(z(1:nx), nx, x2(1:nx) - x1(1:nx), x3(1:nx) - x1(1:nx))
         z(1:nx) = z(1:nx) + x1(1:nx)
      end subroutine median_db
            
      
      subroutine mkpmcub_db(x,nx,f)  ! make piecewise monotonic cubic interpolant
         integer, intent(IN) :: nx       ! length of x vector
         double precision, intent(IN)    :: x(nx)    ! junction points, strictly monotonic
         double precision, intent(INOUT) :: f(4,nx)  ! data & interpolation coefficients
      !
      !  on input:   f(1,i) = f(x(i))

      !  the interpolation formula is a simple cubic, 
      !     fi(x) = di + dx*(ci + dx*(bi + dx*ai)) with dx = x-xi
      
      !  and the slope is of course given by
      !     dfi_dx(x) = ci + 2*dx*(bi + 1.5*dx*ai)

      !  on output:  f(1,i) = di (unchanged from input)
      !              f(2,i) = ci
      !              f(3,i) = bi
      !              f(4,i) = ai
      !
         
         double precision :: e_mid_left, e_mid_right
         double precision, dimension(1) :: pp, qp, y
         double precision, dimension(nx) :: h, s_mid, d, s, d_mid, spL, spR, t, tmax
         
         ! intervals and divided differences
         h(1:nx-1) = x(2:nx) - x(1:nx-1)
         s_mid(1:nx-1) = (f(1,2:nx) - f(1,1:nx-1)) / h(1:nx-1) ! slope across interval       
         d(2:nx-1) = (s_mid(2:nx-1) - s_mid(1:nx-2)) / (x(3:nx) - x(1:nx-2)) ! curvature at point

         ! also need approx curvatures at endpoints -- eqn (5.4)
         e_mid_left = (d(3) - d(2)) / (x(4) - x(1))
         e_mid_right = (d(nx-1) - d(nx-2)) / (x(nx) - x(nx-3))
         d(1) = d(2) - e_mid_left * (x(3) - x(1))
         d(nx) = d(nx-1) + e_mid_right * (x(nx) - x(nx-2))
                        
         ! s(i) = candidate slope at x(i) using slopes of adjacent intervals  eqn (2.8)
         call minmod_db(s(2:nx-1), nx-2, s_mid(2:nx-1), s_mid(1:nx-2))

         ! d_mid(i) = curvature at x(i+1/2) based on curvatures at x(i) and x(i+1)  eqn (3.4)
         call minmod_db(d_mid(1:nx-1), nx-1, d(1:nx-1), d(2:nx))
         
         ! spL(i) = slope at x(i) of parabola from left   eqn for p'_(i-1/2)(xi) = in Theorem 1
         spL(2:nx) = s_mid(1:nx-1) + d_mid(1:nx-1)*h(1:nx-1)

         ! spR(i) = slope at x(i) of parabola from right   eqn for p'_(i+1/2)(xi) = in Theorem 1
         spR(1:nx-1) = s_mid(1:nx-1) - d_mid(1:nx-1)*h(1:nx-1)
         
         ! M3-A method (A stands for average)  eqn (3.18)
         f(2,2:nx-1) = 0.5*(spL(2:nx-1)+spR(2:nx-1)) ! average of left and right parabola slopes
         
         ! for M3-quartic method, replace the average by the following:
         ! initial f(2,3:nx-2) = slope of quartic at middle of 5 points
         ! updated f(2,3:nx-2) = median(f(2,3:nx-2), spL(2,3:nx-2), spR(2,3:nx-2))   eqn (3.21)
         
         ! t(i) = candidate slope at x(i) using slopes of adjacent parabolas  eqn (3.13)
         call minmod_db(t(2:nx-1), nx-2, spL(2:nx-1), spR(2:nx-1))
         
         ! tmax(i) = max monotonicity preserving slope at x(i)  eqn (3.19a)
         tmax(2:nx-1) = sign(1d0,t(2:nx-1))*max(3*abs(s(2:nx-1)),1.5*abs(t(2:nx-1)))
         
         ! finally, apply the monotonicity preserving bounds  eqn (3.19b)
         call minmod_db(f(2,2:nx-1), nx-4, f(2,2:nx-1), tmax(2:nx-1))

         ! slope at i=1
         pp(1) = s_mid(1) - d(2)*h(1) ! slope at x(1) of parabola (F1,F2,F3) eqn (5.1)
         qp(1) = pp(1) + e_mid_left*h(1)*(x(3)-x(1)) ! slope at x(1) of cubic (F1,F2,F3,F4)  eqn (4.5a)
         call median_db(f(2,1), 1, s_mid(1), pp(1), qp(1)) ! eqn (5.7)
         y(1) = 3*s_mid(1)
         call minmod_db(f(2,1), 1, f(2,1), y) ! eqn (5.2)
         
         ! slope at i=nx
         pp(1) = s_mid(nx-1) + d(nx-1)*h(nx-1) ! slope at x(nx) of parabola(F(nx-2),F(nx-1),F(nx)) eqn (5.1)
         qp(1) = pp(1) + e_mid_right*(x(nx-2)-x(nx))*h(nx-1) ! slope at x(nx) of cubic (F(nx-3),F(nx-2),F(nx-1),F(nx))  eqn (4.5c)
         call median_db(f(2,nx), 1, s_mid(nx-1), pp(1), qp(1)) ! eqn (5.7)
         y(1) = 3*s_mid(nx-1)
         call minmod_db(f(2,nx), 1, f(2,nx), y) ! eqn (5.2)
   
         ! 2nd and 3rd derivatives
         f(3,1:nx-1) = (3*s_mid(1:nx-1) - 2*f(2,1:nx-1) - f(2,2:nx)) / h(1:nx-1)       
         f(4,1:nx-1) = (f(2,1:nx-1) + f(2,2:nx) - 2*s_mid(1:nx-1)) / h(1:nx-1)**2
      
      end subroutine mkpmcub_db
      
      
      subroutine minmod(z, nx, x, y)
         real, intent(OUT) :: z(nx)    
         integer, intent(IN) :: nx       ! length of vectors
         real, intent(IN) :: x(nx), y(nx)       
         z(1:nx) = 0.5 * (sign(1.0,x(1:nx)) + sign(1.0,y(1:nx))) * min(abs(x(1:nx)), abs(y(1:nx)))    
      end subroutine minmod
      
      
      subroutine median(z, nx, x1, x2, x3)
         real, intent(OUT) :: z(nx)    
         integer, intent(IN) :: nx       ! length of vectors
         real, intent(IN) :: x1(nx), x2(nx), x3(nx)
         call minmod(z(1:nx), nx, x2(1:nx) - x1(1:nx), x3(1:nx) - x1(1:nx))
         z(1:nx) = z(1:nx) + x1(1:nx)
      end subroutine median
            
      
      subroutine mkpmcub(x,nx,f)  ! make piecewise monotonic cubic interpolant
         integer, intent(IN) :: nx       ! length of x vector
         real, intent(IN)    :: x(nx)    ! junction points, strictly monotonic
         real, intent(INOUT) :: f(4,nx)  ! data & interpolation coefficients
      !
      !  on input:   f(1,i) = f(x(i))

      !  the interpolation formula is a simple cubic, 
      !     fi(x) = di + dx*(ci + dx*(bi + dx*ai)) with dx = x-xi
      
      !  and the slope is of course given by
      !     dfi_dx(x) = ci + 2*dx*(bi + 1.5*dx*ai)

      !  on output:  f(1,i) = di (unchanged from input)
      !              f(2,i) = ci
      !              f(3,i) = bi
      !              f(4,i) = ai
      !
         
         real :: e_mid_left, e_mid_right
         real, dimension(1) :: pp, qp, y
         real, dimension(nx) :: h, s_mid, d, s, d_mid, spL, spR, t, tmax
         
         ! intervals and divided differences
         h(1:nx-1) = x(2:nx) - x(1:nx-1)
         s_mid(1:nx-1) = (f(1,2:nx) - f(1,1:nx-1)) / h(1:nx-1) ! slope across interval       
         d(2:nx-1) = (s_mid(2:nx-1) - s_mid(1:nx-2)) / (x(3:nx) - x(1:nx-2)) ! curvature at point

         ! also need approx curvatures at endpoints -- eqn (5.4)
         e_mid_left = (d(3) - d(2)) / (x(4) - x(1))
         e_mid_right = (d(nx-1) - d(nx-2)) / (x(nx) - x(nx-3))
         d(1) = d(2) - e_mid_left * (x(3) - x(1))
         d(nx) = d(nx-1) + e_mid_right * (x(nx) - x(nx-2))
                        
         ! s(i) = candidate slope at x(i) using slopes of adjacent intervals  eqn (2.8)
         call minmod(s(2:nx-1), nx-2, s_mid(2:nx-1), s_mid(1:nx-2))

         ! d_mid(i) = curvature at x(i+1/2) based on curvatures at x(i) and x(i+1)  eqn (3.4)
         call minmod(d_mid(1:nx-1), nx-1, d(1:nx-1), d(2:nx))
         
         ! spL(i) = slope at x(i) of parabola from left   eqn for p'_(i-1/2)(xi) = in Theorem 1
         spL(2:nx) = s_mid(1:nx-1) + d_mid(1:nx-1)*h(1:nx-1)

         ! spR(i) = slope at x(i) of parabola from right   eqn for p'_(i+1/2)(xi) = in Theorem 1
         spR(1:nx-1) = s_mid(1:nx-1) - d_mid(1:nx-1)*h(1:nx-1)
         
         ! M3-A method (A stands for average)  eqn (3.18)
         f(2,2:nx-1) = 0.5*(spL(2:nx-1)+spR(2:nx-1)) ! average of left and right parabola slopes
         
         ! for M3-quartic method, replace the average by the following:
         ! initial f(2,3:nx-2) = slope of quartic at middle of 5 points
         ! updated f(2,3:nx-2) = median(f(2,3:nx-2), spL(2,3:nx-2), spR(2,3:nx-2))   eqn (3.21)
         
         ! t(i) = candidate slope at x(i) using slopes of adjacent parabolas  eqn (3.13)
         call minmod(t(2:nx-1), nx-2, spL(2:nx-1), spR(2:nx-1))
         
         ! tmax(i) = max monotonicity preserving slope at x(i)  eqn (3.19a)
         tmax(2:nx-1) = sign(1.0,t(2:nx-1))*max(3*abs(s(2:nx-1)),1.5*abs(t(2:nx-1)))
         
         ! finally, apply the monotonicity preserving bounds  eqn (3.19b)
         call minmod(f(2,2:nx-1), nx-4, f(2,2:nx-1), tmax(2:nx-1))

         ! slope at i=1
         pp(1) = s_mid(1) - d(2)*h(1) ! slope at x(1) of parabola (F1,F2,F3) eqn (5.1)
         qp(1) = pp(1) + e_mid_left*h(1)*(x(3)-x(1)) ! slope at x(1) of cubic (F1,F2,F3,F4)  eqn (4.5a)
         call median(f(2,1), 1, s_mid(1), pp(1), qp(1)) ! eqn (5.7)
         y(1) = 3*s_mid(1)
         call minmod(f(2,1), 1, f(2,1), y) ! eqn (5.2)
         
         ! slope at i=nx
         pp(1) = s_mid(nx-1) + d(nx-1)*h(nx-1) ! slope at x(nx) of parabola(F(nx-2),F(nx-1),F(nx)) eqn (5.1)
         qp(1) = pp(1) + e_mid_right*(x(nx-2)-x(nx))*h(nx-1) ! slope at x(nx) of cubic (F(nx-3),F(nx-2),F(nx-1),F(nx))  eqn (4.5c)
         call median(f(2,nx), 1, s_mid(nx-1), pp(1), qp(1)) ! eqn (5.7)
         y = 3*s_mid(nx-1)
         call minmod(f(2,nx), 1, f(2,nx), y) ! eqn (5.2)
   
         ! 2nd and 3rd derivatives
         f(3,1:nx-1) = (3*s_mid(1:nx-1) - 2*f(2,1:nx-1) - f(2,2:nx)) / h(1:nx-1)       
         f(4,1:nx-1) = (f(2,1:nx-1) + f(2,2:nx) - 2*s_mid(1:nx-1)) / h(1:nx-1)**2

      end subroutine mkpmcub


      end module interp_1D

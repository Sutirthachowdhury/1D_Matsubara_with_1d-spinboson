subroutine mapverlet(hel,x0,p0)

   use global
   implicit none

   real*8  hel(nstate,nstate)
   real*8  x0(nstate),p0(nstate)
   
   INTEGER(kind=8) :: i, j
   REAL(kind=8) :: x_dot(2), p_dot(2), dtt

   !======================================================!                                                                                                                                              

   dtt = 0.5*dt / 10.D0

      DO j=1, 10
         
         ! half velocity step                                                                                                                                                                            
         !p_dot(1) = - ((2.0*(nstate+2))/(nstate+4.0))*SUM(hel(1,:) * x0(:))
         !p_dot(2) = - ((2.0*(nstate+2))/(nstate+4.0))*SUM(hel(2,:) * x0(:))

         p_dot(1) = - SUM(hel(1,:) * x0(:))
         p_dot(2) = - SUM(hel(2,:) * x0(:))

         p0(:) = p0(:) + 0.5D0 * p_dot(:) * dtt

         ! full position step                                                                                                                                                                            
         !x_dot(1) = ((2.0*(nstate+2))/(nstate+4.0))*SUM(hel(1,:) * p0(:))
         !x_dot(2) = ((2.0*(nstate+2))/(nstate+4.0))*SUM(hel(2,:) * p0(:))

         x_dot(1) = SUM(hel(1,:) * p0(:))
         x_dot(2) = SUM(hel(2,:) * p0(:))

         x0(:) = x0(:) + x_dot(:) * dtt

         ! second half velocity step                                                                                                                                                                     
         !p_dot(1) = - ((2.0*(nstate+2))/(nstate+4.0))*SUM(hel(1,:) * x0(:))
         !p_dot(2) = - ((2.0*(nstate+2))/(nstate+4.0))*SUM(hel(2,:) * x0(:))

         p_dot(1) = - SUM(hel(1,:) * x0(:))
         p_dot(2) = - SUM(hel(2,:) * x0(:))

         p0(:) = p0(:) + 0.5D0 * p_dot(:) * dtt

      END DO

      !((2.0*(nstate+2))/(nstate+4.0))



 end subroutine mapverlet


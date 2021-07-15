subroutine forces_t0(Q_mats,f_mats)

    use global

    implicit none

    real*8, intent(in) :: Q_mats(ndof,lb_m:ub_m)    
    real*8, intent(out) :: f_mats(ndof,lb_m:ub_m)
    integer :: i,j,k,l,n
    real*8 :: matsfreq(lb_m:ub_m)

    !------- matsubara frequency

    do n = lb_m,ub_m
      matsfreq(n) = (2.*pi*n)/beta
    enddo  

    !!! Compute the quartic force term and spring forces
  
    f_mats(:,:) = 0.0d0

    do i = lb_m,ub_m
      f_mats(:,i) = - matsfreq(i)**2 * Q_mats(:,i)
    enddo
  
    do i = lb_m,ub_m
      do j = lb_m,ub_m
        do k = lb_m,ub_m
          do l = lb_m,ub_m
            f_mats(:,i) = f_mats(:,i) - Q_mats(:,j)*Q_mats(:,k)*Q_mats(:,l)*a_4(i,j,k,l)*4.
          enddo 
        enddo 
      enddo 
    enddo 

    
  
  return

end subroutine forces_t0  
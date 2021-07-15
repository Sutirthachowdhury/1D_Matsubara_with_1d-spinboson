subroutine bigphase_sub(Q_mats,P_mats,phase)

    use global

    implicit none

    real*8, intent(in) :: Q_mats(1,lb_m:ub_m),P_mats(1,lb_m:ub_m)
    real*8, intent(out) :: phase

    real*8 :: matsfreq(lb_m:ub_m),theta
    integer :: i

          !------------- calculation of matsubara phase term----------
    do i = lb_m,ub_m
        matsfreq(i) = (2.0*pi*i)/beta
    enddo  

      theta = 0.
      do i=lb_m,ub_m
        theta = theta + matsfreq(i)*Q_mats(1,-i)*P_mats(1,i)
      enddo

    phase = cos(beta*theta)


end subroutine bigphase_sub
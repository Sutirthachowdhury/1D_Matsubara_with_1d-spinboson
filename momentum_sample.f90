subroutine momentum_sample(P_mats)

    use global

    implicit none

    real*8, intent(out) :: P_mats(ndof,lb_m:ub_m)
    integer :: i,j,k,l
    real*8 :: gran

    do j = 1,ndof
        do i = lb_m,ub_m
            P_mats(j,i) = sqrt(mnuc/beta)*gran()
        enddo
    enddo   

end subroutine momentum_sample
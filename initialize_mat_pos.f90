subroutine initialize_mat_pos(Q_mats)

    use global

    real*8, intent(out) :: Q_mats(ndof,lb_m:ub_m)
    real*8 :: gran,drp,rand,rnum

    integer :: i,j,k

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Initialize positions 
    do j = 1,ndof
        do i=lb_m,ub_m
            Q_mats(j,i) = 0.
        enddo
    enddo   

    !! Initialize positions (from gran)
    !drp = 1.0/sqrt(betap)

    !do j = 1,ndof
    !    do i =lb_m,ub_m
    !        Q_mats(j,i) = gran()*drp    
    !    enddo
    !enddo   

     !! Initialize positions (from random distributions)

    !do j = 1,ndof
    !    do i = lb_m,ub_m
    !        call random_number(rand)
    !        rnum=rand
    !        Q_mats(j,i)=(rnum-0.5)*2
    !    enddo
    !enddo



end subroutine initialize_mat_pos
subroutine evolve_t0(Q_mats,P_mats)

    use global

    implicit none
    real*8, intent(inout) :: Q_mats(ndof,lb_m:ub_m),P_mats(ndof,lb_m:ub_m)
    real*8 :: f_mats(ndof,lb_m:ub_m)

    integer :: i,j,k,l,m

    call forces_t0(Q_mats,f_mats)

    !!! Step 1: Half velocity update
    do m = 1,ndof
        do i = lb_m,ub_m
            P_mats(m,i) = P_mats(m,i) + (dt/2.)*f_mats(m,i)
        enddo
    enddo
    
    !!! Step 2: Full position update    
    do m = 1,ndof
        do i = lb_m,ub_m
            Q_mats(m,i) = Q_mats(m,i) + dt*P_mats(m,i)/(mnuc)
        enddo    
    enddo

    call forces_t0(Q_mats,f_mats)

      !!! Step 1: Half velocity update
    do m = 1,ndof
        do i = lb_m,ub_m
            P_mats(m,i) = P_mats(m,i) + (dt/2.)*f_mats(m,i)
        enddo
    enddo


end subroutine evolve_t0
subroutine initialize_mapp(x0_mc,p0_mc)

    use global

    real*8, intent(out) :: x0_mc(nstate,nb),p0_mc(nstate,nb)

    real*8 :: gran,sigma
    real*8 :: Rfocus(nstate),rand,theta
    integer :: i,j,k,init_state 

    x0_mc = 0.0
    p0_mc = 0.0

    init_state = 1

    ! unoccupied (jeremy estimator)
    Rfocus(:)=1.0
    ! unoccupied (nandini estimator)
    !Rfocus(:)=0.707

   !occupied(jeremy estimator)
    Rfocus(init_state)=1.732
    ! occupied (nandini)
    !Rfocus(init_state)=1.398

    do i= 1,nb
 
   do j = 1,nstate 
      call random_number(rand)
      theta = 2*pi*rand

      x0_mc(j,i) =  Rfocus(j)*COS(theta) !gran()*d_pos
      p0_mc(j,i) =  Rfocus(j)*SIN(theta) !gran()*d_mom

   enddo
   
enddo
    
    !sigma = sqrt(0.5)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Initialize positions 
    !do j = 1,nstate
    !    do i=1,nb
    !        x0_mc(j,i) = gran()*sigma
    !        p0_mc(j,i) = gran()*sigma
    !    enddo
    !enddo   

end subroutine initialize_mapp
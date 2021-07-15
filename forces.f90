subroutine forces(Q_mats,x0,p0,f_mats)

    use global

    implicit none


    real*8, intent(in) :: Q_mats(ndof,lb_m:ub_m)
    real*8, intent(inout) :: x0(nstate,nb),p0(nstate,nb)    
    real*8, intent(out) :: f_mats(ndof,lb_m:ub_m)
    integer :: i,j,k,l,m
    real*8 :: pot

    ! mapping hamiltonian variables
    real*8  hel(nstate,nstate,nb)
    real*8  dhel(nstate,nstate,ndof,nb)
  

    real*8 :: q(ndof,nb),f_local(ndof,nb)

    f_local = 0.0d0
    
    !--------- going back to bead representation-----------
    call back_transform_matrix(Q_mats,q)

   
    do k = 1,nb

      call potential(q(:,k),hel(:,:,k))

      !----- half step mapverlet-------------------
      call mapverlet(hel(:,:,k),x0(:,k),p0(:,k))

      !---- state-indep force-------------
      do m=1,ndof
        f_local(m,k) = -k1*q(m,k)
      enddo
      !----------------------------

    !------------ state dep fprce------------
      call grad_of_potential(q(:,k),dhel(:,:,:,k))

      do m=1,ndof
        do i = 1,nstate
          do j = 1,nstate
          
              f_local(m,k) = f_local(m,k) - 0.5*dhel(i,j,m,k)*(x0(i,k)*x0(j,k)+p0(i,k)*p0(j,k)-del(i,j)) 
                !+ ((nstate)/(2*(nstate+4.0)))*dhel(i,j,m,k)*(x0(i,k)*x0(j,k)+p0(i,k)*p0(j,k)-0.5*del(i,j))
       
          end do
        end do
      end do
   !----------------------------------------------------

  enddo 


  call transform_matrix(f_local,f_mats)

  
  return

end subroutine forces   
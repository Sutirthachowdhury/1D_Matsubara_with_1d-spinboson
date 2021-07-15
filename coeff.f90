subroutine coeff(pot,coeff1)

    use global, only : nb,nstate,ndof,betap,eigen,EVECTOR 
  
    implicit none
  
    real*8, intent(in) :: pot(nstate,nstate)
    real*8, intent(out) :: coeff1(nstate,nstate)
  
    integer :: i,j
  
    call DIAG (eigen,EVECTOR,pot)
  
  !  write(*,*) "eigen=",eigen
  !  write(*,*) "evec=",EVECTOR (1,1), EVECTOR (1,2)
  !  write(*,*) "evec=",EVECTOR (2,1), EVECTOR (2,2)
  
  
    coeff1 = 0.
  
    do i=1,nstate
       do j=1,nstate
          coeff1(i,j) = sum(EVECTOR(i,:)*dexp(-betap*eigen(:))*EVECTOR(j,:))
       enddo
    enddo
  
  !  do i=1,nstate
  !     write(*,*) (coeff1 (i,j),j=1,nstate)
  !  end do
  
  !========================================================================================!
  ! Ananth-Miller approximated version (JCP2010) based on Chandler's short time approximation
  !  do i=1,nstate
  !     coeff2(i,i)=0.5*pot(i,i)/dble(nb)*dexp(-0.5*betap*pot(i,i))
  !     coeff1(i,i)=dexp(-0.5*betap*pot(i,i))
  !     do j=1,nstate
  !        if(j.ne.i) coeff1(i,j)=-0.5*betap*pot(i,j)*&
  !             dexp(-0.5*betap*pot(j,j))
  !        if(j.ne.i) coeff2(i,j)=0.5*pot(i,j)/dble(nb)*&
  !             dexp(-0.5*betap*pot(j,j))*(-0.5*betap*pot(j,j)+1.0d0)
  !     end do
  !  end do
  !=========================================================================================!
  
    return
  
  end subroutine coeff
  
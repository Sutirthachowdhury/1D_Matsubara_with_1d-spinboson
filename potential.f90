subroutine potential(r,hel)

    use global
   
    implicit none
  
    real*8, intent(in) :: r(ndof)
    real*8, intent(out) :: hel(nstate,nstate)
  
    integer :: i,j
  
  ! bare electronic hamiltonian 

    hel(1,1) = alpha
    hel(2,2) = alpha
    hel(1,2) = 0.5*delta
    hel(2,1) = hel(1,2)
  
      ! cq = k (in Thoss JCP 2013), qs = x (in Thoss 2013 JCP).
  
     !diagonal term only   
     do j=1,nstate
        do i=1,ndof
           hel(j,j) = hel(j,j) + cq(j,i)*r(i) 
        end do
     end do
  
    return
  
  end subroutine potential
  
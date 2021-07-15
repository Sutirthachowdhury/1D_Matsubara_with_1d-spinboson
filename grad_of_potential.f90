subroutine grad_of_potential(qs,dhel)
    use global

    implicit none
  
    integer i,j,l
    
    real*8, intent(in) :: qs(ndof)
    
    real*8, intent(out) :: dhel(nstate,nstate,ndof)
  
    dhel=0.
  
    do j=1,nstate
       do i=1,ndof
          dhel(j,j,i) =  cq(j,i)
       enddo
    enddo
    
    return
  end subroutine grad_of_potential
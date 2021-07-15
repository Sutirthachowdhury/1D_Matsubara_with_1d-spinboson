!==========================================================
subroutine init_pot()
    !==========================================================
    ! Compute the transformation tensor for exact evaluation
    ! of the potential using trapezoidal rule
    !==========================================================
    
    use global
    implicit none

    !real*8, intent(out) :: a_4(lb_m:ub_m,lb_m:ub_m,lb_m:ub_m:lb_m:ub_m)

    integer, parameter :: ngrid = 2048 
    real, parameter :: tol = 1e-15
    integer :: i,j,k,l,m
    real :: xgrid(0:ngrid)
    real :: xmin,xmax,dx 
    real :: temp,fct,cte 
    real :: nm_t 
    
    !!! Build the integration grid
    xmax = 1.0
    xmin = 0.
    dx = (xmax-xmin)/ngrid
    
    do i= 0,ngrid
      xgrid(i) = xmin + i*dx
    enddo
    
    !write(*,*) xgrid(2048)
    !write(6,*) '# Potential integrals:'
       
    !!! Compute the quartic potential term
  
     cte=1.0/4.0
      
      do i = lb_m, ub_m
        do j = lb_m, ub_m
          do k = lb_m, ub_m
            do l = lb_m, ub_m
              
                temp = 0.
                do m = 0,ngrid
                    fct = 1.
                    if((m.eq.0).or.(m.eq.ngrid)) fct = 0.5  !trapezoidal rule 
                    temp = temp + fct * nm_t(xgrid(m),i) * nm_t(xgrid(m),j) * nm_t(xgrid(m),k) * nm_t(xgrid(m),l) * dx
                enddo
              if(abs(temp).lt.tol) temp=0.
              a_4(i,j,k,l) = cte * temp 
              !write(6,*)'# quartic: ',i,j,k,l,a_4(i,j,k,l)
            enddo
          enddo
        enddo
      enddo
    
    
    !write(6,*)
    
    return
    end subroutine init_pot
    !==========================================================
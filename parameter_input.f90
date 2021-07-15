subroutine parameter_input
    use global 
    
    implicit none

    integer i
    
    !=========================================
    open(10,file='params.in',status='old')
    read(10,*) nstep
    read(10,*) ntraj
    read(10,*) nmc
    read(10,*) beta
    read(10,*) mnuc
    read(10,*) k1
    read(10,*) kappa
    read(10,*) alpha
    read(10,*) delta
    read(10,*) dt
    read(10,*) ourseed
    close(10)
    !==========================================

    allocate(step(nstate+1))
    allocate(a_4(lb_m:ub_m,lb_m:ub_m,lb_m:ub_m,lb_m:ub_m))
    allocate(eigen(nstate),EVECTOR(nstate,nstate))  

    betap = beta/dble(nb)
    print *, 'betap', betap


    open(10,file='steps.in',status='old')
    read(10,*) step(1), step(2), step(3)
    close(10)

    !--- freq of HO-----
    omega(1) = 1.0    
    !--------------------
          
    !grad part
    cq(1,1) = sqrt(2.0)*kappa
    cq(2,1) = -sqrt(2.0)*kappa
   
    ! kronicker delta for mapping hamiltonian (we will use later)---- 

     del=0.

     do i=1,nstate
        del(i,i)=1.
     end do

    !----------------------------------- 

end subroutine parameter_input

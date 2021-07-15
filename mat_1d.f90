program mat_1d

    use global
    implicit none

    integer :: istep,iup,idn,isample,ieq     
    integer,allocatable :: seed(:)
    integer seed_dimension,time


    integer i,j,k,l,m,n,ncount

    !======= premitive bead variables=============

    !real*8, allocatable :: xe_mc(:,:),xe(:,:) ! nuclear bead positions 
    !real*8, allocatable :: ve(:,:) ! nuclear beads momentum
    !=======================================


    real*8, allocatable:: Q_mats(:,:),P_mats(:,:),f_mats(:,:) !"Matsubara" positions, momentum and force
    real*8, allocatable:: x0_mc(:,:),p0_mc(:,:) ! mapping positons and momentum (in bead rep.)
    real*8, allocatable:: Q_mats_dyn(:,:),P_mats_dyn(:,:),f_mats_dyn(:,:),x0(:,:),p0(:,:) ! this variables will use in dynamics
    real*8, allocatable :: matsfreq(:)
    real*8, allocatable :: pop0(:)
    real*8, allocatable :: pop(:,:) ! this is population

    real*8 :: q(ndof,nb)

    real*8 :: dp,gran,rand,theta 
   

    complex*16 :: pref ! this is the weighting term from non-adiabatic part
    real*8 :: phase ! this is the real part of e^{i*beta*theta}

    real*8 :: partition ! weighting the cosine term

    ! getting the inputs
    call parameter_input
 
    !-------- random seed-------------------
    CALL RANDOM_SEED(size=seed_dimension)
    ALLOCATE (seed(seed_dimension))
    do i=1,seed_dimension
      seed(i) =time()+3*i-1 + ourseed
    end do
    CALL RANDOM_SEED(PUT=seed)
    !-----------------------------------------

    !allocate(xe_mc(ndof,nb),xe(ndof,nb),ve(ndof,nb))
    allocate(Q_mats(ndof,lb_m:ub_m),P_mats(ndof,lb_m:ub_m),f_mats(ndof,lb_m:ub_m))
    allocate(Q_mats_dyn(ndof,lb_m:ub_m),P_mats_dyn(ndof,lb_m:ub_m),f_mats_dyn(ndof,lb_m:ub_m))
    allocate(x0(nstate,nb),p0(nstate,nb))
    allocate(matsfreq(lb_m:ub_m))
    allocate(x0_mc(nstate,nb),p0_mc(nstate,nb))
    
    allocate(pop0(nstate))
    allocate(pop(nstate,nstep))

    !--- initialization of correlation function 
    pop = 0.0
    pop0 = 0.0
    !----------------------------

    !------- initializing the partition function----

    partition = 0.0

    !---- initialize the matsubara pos----
    call initialize_mat_pos(Q_mats) ! for matsubara positions
    !===========================================================

    do j = 1,ntraj
    
      call monte_carlo(Q_mats)
      call initialize_mapp(x0_mc,p0_mc)
       
       !store sample configuration
       Q_mats_dyn = Q_mats
       x0=x0_mc
       p0=p0_mc
  
      call momentum_sample(P_mats)

      !------------- calculation of mat phase----------

      call bigphase_sub(Q_mats_dyn,P_mats,phase)

      !---------- partition function calculation ---------

      partition = partition + phase

      !---------------------------------------

      !----------- population at t =0----------------------
      do k=1,nstate
        do i=1,nb
           pop0(k) =  pop0(k) &
                + (phase)*((0.5*(x0(k,i)**2+p0(k,i)**2-1.0))/real(nb))
        enddo
     enddo
      !================dynamics step ==========================================
     
      do istep = 1,nstep

        call evolve(Q_mats_dyn,P_mats,x0,p0)
        
        do k=1,nstate
          do i=1,nb
             pop(k,istep) =  pop(k,istep) &
                  + (phase)*((0.5*(x0(k,i)**2+p0(k,i)**2-1.0))/real(nb)) 

          enddo
       enddo


      enddo 


    enddo 

    pop = pop/real(partition)
    pop0 = pop0/real(partition)

    write(3001,*) pop0

    do ncount=1,nstep
      write(200,222) ncount*dt,(pop(k,ncount),k=1,nstate)
    end do

222      format(60(e13.6,2x))
    stop
end program mat_1d
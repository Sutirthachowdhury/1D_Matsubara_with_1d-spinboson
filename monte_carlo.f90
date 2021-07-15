subroutine monte_carlo(Q_mats)

    use global
    implicit none

    real*8, allocatable :: Q_new(:,:)
    real*8, intent(inout) :: Q_mats(ndof,lb_m:ub_m)

    integer :: tot_mode,imode
    integer :: init_nmc
    integer :: i,j,iup,idn,k,init_mc,l

    real*8 :: iacc_x,iatt_x,rnum,rand
    real*8 :: wt,wt1,wtt
    real*8 :: energy,energy_old,del_e
 
   
    allocate(Q_new(ndof,lb_m:ub_m))

     !=========MONTE CARLO INTEGRATION RUN==========   
 
    !initialize all zeroes

      iatt_x=0.d0
      iacc_x=0.d0
    
  do init_mc=1,nmc

    !-----randomly pick a mode to move-----                                 
    call random_number(rand)
    tot_mode = int(rand*mats_mode)
    imode = int(tot_mode-2)   ! here "-1" refers the lowest mats modes, in this case as M=3 so lowest mat mode is "-1"
    !--------------------------------
 
    !move individual modes, generating trial moves
        
      Q_new = Q_mats
    
      do j=1,ndof
        call random_number(rand)
        rnum = rand
        Q_new(j,imode) = Q_mats(j,imode) + (rnum-0.5d0)*step(1)  ! here need to check with gran as well
      
      end do

    call sampling_u0(Q_new,wtt)
    
    energy = wtt 

    iatt_x = iatt_x + 1.0
    
    !checking against original weight
     
      call sampling_u0(Q_mats,wt1)
      energy_old = wt1 

    ! calculating the energy difference
    del_e = (energy-energy_old)

    !accepting/rejecting step
    wt = dexp(-beta*del_e)

    call random_number(rand)
    rnum=rand

   if(rnum.lt.wt) then
      
         Q_mats(:,imode) = Q_new(:,imode)
         iacc_x=iacc_x+1.d0
         !write(*,*) int_mc,imode,"nuc"
     
   end if
          
  enddo 

  open(20,file='final_config.dat')
  write(20,*) Q_mats
 
  close(20)

end subroutine monte_carlo
!============================================
subroutine sampling_u0(Q_new,U0)

  use global
  implicit none
 
  real*8, intent(in) :: Q_new(ndof,lb_m:ub_m)
  real*8, intent(out) :: U0
 
  real*8 :: q(ndof,nb),pot
  integer :: l
 
   !------------- transform the modes onto bead representation ----------------
  call back_transform_matrix(Q_new,q)
 
  !--- compute the potential in bead representation------
 
 pot = 0.
 
 do l = 1,nb
  pot = pot + 0.5*k1*q(1,l)**2   !state-indepedent part is hard-coded here, one can do it more generally
 enddo 
 
 pot = pot/real(nb)
 
 U0 = pot
  
 end subroutine sampling_u0
!---------------------------------------------------------------
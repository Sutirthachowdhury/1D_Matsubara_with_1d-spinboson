subroutine estimator(Q_mats,x0_mc,p0_mc,theta)

    use global
    implicit none
  
    real*8, intent(in) :: Q_mats(ndof,lb_m:ub_m),x0_mc(nstate,nb),p0_mc(nstate,nb)
    complex*16, intent(out) :: theta
  
    
    integer :: i,iup,idn,j,k,l
    real*8 :: identity(nstate,nstate)
    real*8 :: coeff1(nstate,nstate)
    real*8 :: inv1(nstate,nstate)
    real*8 :: hel(nstate,nstate)
    complex*16 :: weight(nstate,nstate)
    complex*16 :: cmat(nstate,nstate)
    complex*16 :: full_mat(nstate,nstate)
    real*8 :: q(ndof,nb)
  
    !definition of initial gamma matrix--------
    cmat(1,1) = cmplx(1.,0.)
    cmat(1,2) = cmplx(0.,0.)
    cmat(2,1) = cmat(1,2)
    cmat(2,2) = cmat(1,1)
    !---------------------------------------
  
    full_mat(1,1) = cmplx(1.,0.)
    full_mat(1,2) = cmplx(0.,0.)
    full_mat(2,1) = full_mat(1,2)
    full_mat(2,2) = full_mat(1,1)
  
    !---- transform back to bead coordinate
    call back_transform_matrix(Q_mats,q)
  
    do i=1,nb
  
      iup=i+1
      if(iup.gt.nb) iup=1
      idn=i-1
      if(idn.lt.1) idn=nb
   
      call potential(q(:,i),hel)
      call coeff(hel,coeff1)
      inv1=coeff1
  
       !initialization of weight
     weight=cmplx(0.,0.)
   
     do j=1,nstate
       do k=1,nstate
          !for prefactor
          weight(j,k) = ((x0_mc(j,i)+eye*p0_mc(j,i))*(x0_mc(k,i)-eye*p0_mc(k,i)))&
               -(0.5d0*del(j,k))
       enddo
    enddo
  
    full_mat = MATMUL(weight,inv1)
  
  !full prefactor
    
   cmat = MATMUL(cmat,full_mat)
  
  enddo
  
  !----- trace----------
  theta=cmplx(0.,0.)
  
  do l = 1,nstate
     theta = theta + cmat(l,l)
  enddo
  !------------------------
  
  end subroutine estimator
  !==========================================================================
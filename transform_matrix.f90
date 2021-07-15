subroutine transform_matrix (f_local,f_mats)
    use global
    implicit none
    integer ::  i,j,k,l,m,n
 
    !     ------------------------------------------------------------------                                                                                     
 !     Orthogonal Transformation matrix for free ring-polymer cartesian to normal modes                                                                          
 !     ------------------------------------------------------------------                                                                                        
 
    real*8, intent(in) :: f_local(ndof,nb)
    real*8, intent(out) :: f_mats(ndof,lb_m:ub_m)
 
    real*8 :: cmat(nb,lb_m:ub_m)
    real*8 :: fnew(ndof,lb_m:ub_m)
    real*8 :: pibyn
 
       pibyn = dacos(-1.d0)/nb
 
   do j = 1,nb
 
          do l = lb_m,ub_m
 
             if(l.eq.0)then
                cmat(j,l) = 1.0  !1.0/sqrt(real(nb))                                                                                                             
 
 
             else if(l.ge.lb_m.AND.l.lt.0) then
                cmat(j,l) = dsqrt(2.0d0)*dsin(2.0*pibyn*j*l)
 
             else if(l.gt.0.AND.l.le.ub_m) then
                cmat(j,l) = dsqrt(2.0d0)*dcos(2.0*pibyn*j*l)
 
             endif
 
          enddo
 
       enddo
      !========================================================================                                                                                  
 
       fnew = 0.0d0
 
       do j = 1,nb
          do k = lb_m,ub_m
             do m = 1,ndof
                fnew(m,k) = fnew(m,k) + f_local(m,j)*cmat(j,k)
        
 
             enddo
          enddo
       enddo
 
 
    f_mats = fnew/real(nb)
       !P_norm = pnew/real(nb)
       
 
       return
     end subroutine transform_matrix
     !==============================================================               
 

      subroutine realft (data,m,n,mode)
      implicit real*8  (a-h,o-z)
      integer*8 plana,planb
!      integer*8 k,j
!     ------------------------------------------------------------------
!     FFT of m real arrays (if mode = 1) or complex Hermitian
!     arrays in real storage (if mode = -1), using -lfftw3.
!     Works equally well with f77 and ifc.
!     ------------------------------------------------------------------
!
      dimension data(m,n)
      parameter (nmax = 1024)
      dimension copy(nmax)
      data np /0/
      save copy,scale,plana,planb,np
!
      if (n .ne. np) then
         if (n .gt. nmax) stop 'realft 1'
         scale = dsqrt(1.d0/n)
         call dfftw_plan_r2r_1d (plana,n,copy,copy,0,64)
         call dfftw_plan_r2r_1d (planb,n,copy,copy,1,64)
         np = n
      endif
      do k = 1,m
         do j = 1,n
            copy(j) = data(k,j)
         enddo
         if (mode .eq. 1) then
            call dfftw_execute (plana)
         else if (mode .eq. -1) then
            call dfftw_execute (planb)
         else
            stop 'realft 2'
         endif
         do j = 1,n
            data(k,j) = scale*copy(j)
         enddo
      enddo
      return
    end 


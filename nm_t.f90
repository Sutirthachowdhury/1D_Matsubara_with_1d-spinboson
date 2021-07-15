FUNCTION nm_t(x,j)
!==========================================================
! Normal mode transformation matrix 
! Note: The NM is computed in the N\to\infty limit and 
!       includes a sqrt(N) factor. In this notation,
!       the 'j' index refers to a Matsubara mode 
!       whereas 'x' refers to i/N.
!==========================================================

use global, only : pi,ub_m,lb_m
implicit none
integer :: j
real :: x,nm_t

if (j.eq.0) then
  nm_t = 1.
elseif (j.gt.0) then
  nm_t = sqrt(2.)*sin(2.*pi*x*j)
else
  nm_t = sqrt(2.)*cos(2.*pi*x*j)
endif

END FUNCTION nm_t
!==========================================================
!***************************************************************************
!*
!*  TRIDAGR
!*
!***************************************************************************
subroutine tridagr(ar,br,cr,knownr,potr)
implicit none
include 'runscf.h'
!***************************************************************************
!*
!  tridagr solves for potr from the linear tridiagonal system of equations
!  with br being the diagonal elements of the matrix, ar and cr are the
!  off-diagonal elements and knownr is the right hand side.  The code
!  comes from section 2.4, page 43 of Numerical Recipes in Fortran, 2nd ed.
!*
!***************************************************************************
!*
!*  Subroutine Arguments

real, dimension(numr) :: ar, cr

real, dimension(numr,numphi) :: br

real, dimension(numr,numz,numphi) :: knownr, potr

!*
!***************************************************************************
!*
!*  Local Variables

real, dimension(numr,numphi) :: bet, gam  

integer :: J, K, L

!*
!***************************************************************************
!  initialize the local variables
do L = 1, numphi
   do J = 1, numr
      gam(J,L) = 0.0
      bet(J,L) = 0.0
   enddo
enddo

! setup
do L = 1, numphi
   bet(2,L) = br(2,L)
enddo
do L = 1, numphi
   do K = zlwb, zupb
      potr(2,K,L) = knownr(2,K,L) / bet(2,L)
   enddo
enddo

!  decomposition and forward substitution
do L = 1, numphi
   do J = 3, numr-1
      gam(J,L) = cr(J-1) / bet(J-1,L)
      bet(J,L) = br(J,L) - ar(J)*gam(J,L)
   enddo
enddo
do L = 1, numphi
   do K = zlwb, zupb
      do J = 3, numr-1
         potr(J,K,L) = (knownr(J,K,L)-ar(J)*potr(J-1,K,L)) / bet(J,L)
      enddo
   enddo
enddo

! back subsitution
do L = 1, numphi
   do K = zlwb, zupb
      do J = numr-2, 2, -1
         potr(J,K,L) = potr(J,K,L) - gam(J+1,L)*potr(J+1,K,L)
      enddo
   enddo
enddo

return
end subroutine tridagr

subroutine binary_scf(model_number, initial_model_type, ra, rb, rc, rhom1, rhom2, frac, pin, eps, &
                     qfinal)
implicit none
include 'runscf.h'
!**************************************************************************************************
!
!  subroutine arguments
!

integer, intent(in) :: model_number
integer, intent(in) :: initial_model_type
integer, intent(in) :: ra
integer, intent(in) :: rb
integer, intent(in) :: rc
real, intent(in) :: rhom1
real, intent(in) :: rhom2
real, intent(in) :: frac
real, intent(in) :: pin
real, intent(in) :: eps
integer, intent(out) :: qfinal

!
!**************************************************************************************************
!
!   global varaibles
!

real, dimension(numr,numz,numphi) :: pot, rho
common /poisson/ pot, rho

real, dimension(numr) :: r, rhf, rinv, rhfinv
real, dimension(numz) :: zhf
real, dimension(numphi) :: phi
common /grid/ rhf, r, rhfinv, rinv, zhf, phi

real, dimension(numr) :: rhf_g, r_g, rhfinv_g, rinv_g
real, dimension(numz) :: zhf_g
common /global_grid/ rhf_g,r_g,rhfinv_g,rinv_g,zhf_g

real, dimension(numphi) :: cosine
real, dimension(numphi) :: sine
common /trig/ cosine, sine

real :: dr, dz, dphi, drinv, dzinv, dphiinv
common /coord_differentials/ dr, dz, dphi, drinv, dzinv, dphiinv

!
!**************************************************************************************************
!
!    local variables
!

real, dimension(numr, numz, numphi) :: h, pot_it, pot_old, temp

real, dimension(numr, numphi) :: psi

real, dimension(maxit) :: c1, c2, mass1, mass2, omsq, hm1, hm2

real :: cnvgom, cnvgc1, cnvgc2, cnvgh1, cnvgh2

real :: dpot, dpsi, psitmp1, psitmp2, pottmp1, pottmp2

real :: ret1, ret2

real :: xavg1, xavg2, com, separation

real :: virial_error, virial_error1, virial_error2

real :: virial_error_prev

real :: volume_factor

real :: gamma

real :: time1, time2

integer :: za, phia
integer :: zb, phib
integer :: zc, phic

real :: temp_hm1, temp_hm2

real, dimension(3) :: rhm1, rhm2, temp_rhm1, temp_rhm2

integer :: phi1, phi2, phi3, phi4

integer :: I, J, K, L, Q

!
!**************************************************************************************************

call cpu_time(time1)

do L = 1, numphi
   do K = 1, numz
      do J = 1, numr
         h(J,K,L) = 0.0
      enddo
   enddo
enddo

do L = 1, numphi
   do J = 1, numr
      psi(J,L) = 0.0
   enddo
enddo

qfinal = 1

phi1 = int(numphi / 4.0) - 1
phi2 = int(numphi / 4.0) + 1
phi3 = int(3.0 * numphi / 4.0) - 1
phi4 = int(3.0 * numphi / 4.0) + 1

phia = 1
phib = 1
phic = numphi / 2 + 1
za = 2
zb = 2
zc = 2

volume_factor = 2.0 * dr * dz * dphi

gamma = 1.0 + 1.0 / pin

mass1 = 0.0
mass2 = 0.0
omsq = 0.0
hm1 = 0.0
hm2 = 0.0
c1 = 0.0
c2 = 0.0

virial_error = 1.0

! calculate the initial total mass
do K = philwb, phiupb
   do J = zlwb, zupb
      do I = rlwb, rupb
         temp(I,J,K) = rhf(I) * rho(I,J,K)
      enddo
   enddo
enddo
call binary_sum(temp, ret1, ret2)
mass1(1) = volume_factor * ret1
mass2(1) = volume_factor * ret2

! calculate the initial center of mass
do K = philwb, phiupb
   do J = zlwb, zupb
      do I = rlwb, rupb
         temp(I,J,K) = rhf(I) * rhf(I) * cosine(K) * rho(I,J,K)
      enddo
   enddo
enddo
call binary_sum(temp, ret1, ret2)
xavg1 = volume_factor * ret1 / mass1(1)
xavg2 = volume_factor * ret2 / mass2(1)
separation = xavg1 - xavg2
com = separation * mass2(1) / ( mass1(1) + mass2(1) )
com = xavg1 - com

open(unit=13,file='iteration_log',form='formatted',status='unknown',position='append')
write(13,*) mass1(1), mass2(1), xavg1, xavg2, com, separation

do Q = 2, maxit-1                                   ! START OF THE ITERATION CYCLE

   ! solve the Poisson equation for the current density field
   if ( Q == 2 ) then
      call potential_solver(0)
      do L = 1, numphi
         do K = 1, numz
            do J = 1, numr
               pot_old(J,K,L) = pot(J,K,L)
            enddo
         enddo
      enddo
   else
      call potential_solver(1)
   endif
   do L = 1, numphi
      do K = 1, numz
         do J = 1, numr
            pot_it(J,K,L) = (1.0 - frac) * pot(J,K,L) + frac * pot_old(J,K,L)
            pot_old(J,K,L) = pot(J,K,L)
         enddo
      enddo
   enddo

   ! compute the form factor for the centrifugal potential
   do K = philwb, phiupb
      do I = rlwb-1, rupb+1
         psi(I,K) = - 0.5 * ( (rhf(I)*cosine(K) - com)**2 + rhf(I)*rhf(I)*sine(K)*sine(K) )
      enddo
   enddo

   ! calculate the angular frequency
   pottmp1 = 0.5 * ( pot_it(ra,za,phia) + pot_it(ra-1,za,phia) )
   psitmp1 = 0.5 * ( psi(ra,phia) + psi(ra-1,phia) )

   pottmp2 = 0.5 * ( pot_it(rb,zb,phib) + pot_it(rb+1,zb,phib) )
   psitmp2 = 0.5 * ( psi(rb,phib) + psi(rb+1,phib) )

   dpot = pottmp1 - pottmp2
   dpsi = psitmp2 - psitmp1
   omsq(Q) = dpot / dpsi

   ! calculate the two integration constants
   pottmp1 = 0.5 * ( pot_it(rb,zb,phib) + pot_it(rb+1,zb,phib) )
   psitmp1 = 0.5 * ( psi(rb,phib) + psi(rb+1,phib) )

   pottmp2 = 0.5 * ( pot_it(rc,zc,phic) + pot_it(rc+1,zc,phic) )
   psitmp2 = 0.5 * ( psi(rc,phic) + psi(rc+1,phic) )

   c1(Q) = pottmp1 + omsq(Q) * psitmp1
   c2(Q) = pottmp2 + omsq(Q) * psitmp2

   ! now compute the new  enthalpy field from the potential and SCF constants
   temp_hm1 = 0.0
   temp_rhm1 = 0.0
   temp_hm2 = 0.0
   temp_rhm2 = 0.0
   do K = philwb, phi1
      do J = zlwb, zupb
         do I = rlwb, rupb
            h(I,J,K) = c1(Q) - pot_it(I,J,K) - omsq(Q) * psi(I,K)
            if ( h(I,J,K) > temp_hm1 ) then
               temp_hm1 = h(I,J,K)
               rhm1(1) = rhf(I)
               rhm1(2) = zhf(J)
               rhm1(3) = phi(K)
            endif
         enddo
      enddo
   enddo
   do K = phi2, phi3
      do J = zlwb, zupb
         do I = rlwb, rupb
            h(I,J,K) = c2(Q) - pot_it(I,J,K) - omsq(Q) * psi(I,K)
            if ( h(I,J,K) > temp_hm2 ) then
               temp_hm2 = h(I,J,K)
               rhm2(1) = rhf(I)
               rhm2(2) = zhf(J)
               rhm2(3) = phi(K)
            endif
         enddo
      enddo
   enddo
   do K = phi4, phiupb
      do J = zlwb,  zupb
         do I = rlwb, rupb
            h(I,J,K) = c1(Q) - pot_it(I,J,K) - omsq(Q) * psi(I,K)
            if ( h(I,J,K) > temp_hm1 ) then
               temp_hm1 = h(I,J,K)
               rhm1(1) = rhf(I)
               rhm1(2) = zhf(J)
               rhm1(3) = phi(K)
            endif
         enddo
      enddo
   enddo
   hm1(Q) = temp_hm1
   hm2(Q) = temp_hm2

   ! calculate the new density field from the enthalpy
   do K = 1, phi1
      do J = zlwb, zupb
         do I = rlwb, rupb
            if ( h(I,J,K) > 0.0 ) then
               rho(I,J,K) = rhom1 * ( h(I,J,K) / hm1(Q) )**pin
            else
               rho(I,J,K) = 0.0
            endif
         enddo
      enddo
   enddo
   do K = phi2, phi3
      do J = zlwb, zupb
         do I = rlwb, rupb
            if ( h(I,J,K) > 0.0 ) then
               rho(I,J,K) = rhom2 * ( h(I,J,K) / hm2(Q) )**pin
            else
               rho(I,J,K) = 0.0
            endif
         enddo
      enddo
   enddo
   do K = phi4, phiupb
      do J = zlwb, zupb
         do I = rlwb, rupb
            if ( h(I,J,K) > 0.0 ) then
               rho(I,J,K) = rhom1 * ( h(I,J,K) / hm1(Q) )**pin
            else
               rho(I,J,K) = 0.0
            endif
         enddo
      enddo
   enddo

   ! zero out the density field between the axis and the inner boundary points
   do K = philwb, phi1
      do J = zlwb, zupb
         do I = rlwb, rupb
            if ( rhf(I) <= rhf(rb) ) then
               rho(I,J,K) = 0.0
            endif
         enddo
      enddo
   enddo
   do K = phi2, phi3
      do J = zlwb, zupb
         do I = rlwb, rupb
            if ( rhf(I) <= rhf(rc) ) then
               rho(I,J,K) = 0.0
            endif
         enddo
      enddo
   enddo
   do K = phi4, phiupb
      do J = zlwb, zupb
         do I = rlwb, rupb
            if ( rhf(I) <= rhf(rb) ) then
               rho(I,J,K) = 0.0
            endif
         enddo
      enddo
   enddo

   ! impose the equatorial boundary condition
   do K = philwb, phiupb
      do J = rlwb, rupb
         rho(J,zlwb-1,K) = rho(J,zlwb,K)
      enddo
   enddo

   ! impose the axial boundary condition
   do L = 1, numphi_by_two
      do K = zlwb, zupb
         rho(rlwb-1,K,L)               = rho(rlwb,K,L+numphi_by_two)
         rho(rlwb-1,K,L+numphi_by_two) = rho(rlwb,K,L)
      enddo
   enddo

   ! calculate the total mass for each star
   do K = philwb, phiupb
      do J = zlwb, zupb
         do I = rlwb, rupb
            temp(I,J,K) = rhf(I) * rho(I,J,K)
         enddo
      enddo
   enddo
   call binary_sum(temp, ret1, ret2)
   mass1(Q) = volume_factor * ret1 
   mass2(Q) = volume_factor * ret2

   ! calculate the center of mass for each star
   do K = philwb, phiupb
      do J  = zlwb, zupb
         do I = rlwb, rupb
            temp(I,J,K) = rhf(I) * rhf(I) * cosine(K) * rho(I,J,K)
         enddo
      enddo
   enddo
   call binary_sum(temp, ret1, ret2)
   xavg1 = volume_factor * ret1 / mass1(Q)
   xavg2 = volume_factor * ret2 / mass2(Q)
   separation = xavg1 - xavg2
   com = separation * mass2(Q) / (mass1(Q) + mass2(Q) )
   com = xavg1 - com
 
   ! has the solution converged?
   cnvgom = abs( (omsq(Q) - omsq(Q-1)) / omsq(Q) )
   cnvgc1 = abs( (c1(Q) - c1(Q-1)) / c1(Q) )
   cnvgc2 = abs( (c2(Q) - c2(Q-1)) / c2(Q) )
   cnvgh1 = abs( (hm1(Q) - hm1(Q-1)) / hm1(Q) )
   cnvgh2 = abs( (hm2(Q) - hm2(Q-1)) / hm2(Q) )

   virial_error_prev = virial_error
   call compute_virial_error(psi, h, sqrt(omsq(Q)), pin, volume_factor, virial_error1, &
                             virial_error2, virial_error)


   write(13,*) Q, mass1(Q), mass2(Q), xavg1, xavg2, com, omsq(Q), c1(Q), c2(Q), &
               hm1(Q), hm2(Q), cnvgom, cnvgc1, cnvgc2, cnvgh1, cnvgh2, virial_error1, &
               virial_error2, virial_error

   if ( cnvgom < eps .and. cnvgc1 < eps .and. cnvgc2 < eps .and. &
        cnvgh1 < eps .and. cnvgh2 < eps ) then
      exit
   endif

   if ( virial_error > virial_error_prev .and. Q > 10  ) then
      exit
   endif

enddo                                               ! END OF THE ITERATION CYCLE

if ( q >= maxit ) then
   qfinal = maxit
else
   qfinal = q + 1
endif

! syncchronize the potential with the converged density distribution
call potential_solver(1)

! and synchronize the centrifugal potential
do K = philwb, phiupb
   do I = rlwb-1, rupb+1
      psi(I,K) = - 0.5 * ( (rhf(I)*cosine(K) - com)**2 + rhf(I)*rhf(I)*sine(K)*sine(K) )
   enddo
enddo

! now calculate the final angular frequency
pottmp1 = 0.5 * ( pot_it(ra,za,phia) + pot_it(ra-1,za,phia) )
psitmp1 = 0.5 * ( psi(ra,phia) + psi(ra-1,phia) )

pottmp2 = 0.5 * ( pot_it(rb,zb,phib) + pot_it(rb+1,zb,phib) )
psitmp2 = 0.5 * ( psi(rb,phib) + psi(rb+1,phib) )

dpot = pottmp1 - pottmp2
dpsi = psitmp2 - psitmp1
omsq(qfinal) = dpot / dpsi

! and calculate the two integration constants
pottmp1 = 0.5 * ( pot_it(rb,zb,phib) + pot_it(rb+1,zb,phib) )
psitmp1 = 0.5 * ( psi(rb,phib) + psi(rb+1,phib) )

pottmp2 = 0.5 * ( pot_it(rc,zc,phic) + pot_it(rc+1,zc,phic) )
psitmp2 = 0.5 * ( psi(rc,phic) + psi(rc+1,phic) )

c1(qfinal) = pottmp1 + omsq(qfinal) * psitmp1
c2(qfinal) = pottmp2 + omsq(qfinal) * psitmp2

hm1(qfinal) = hm1(qfinal-1)
hm2(qfinal) = hm2(qfinal-1)
mass1(qfinal) = mass1(qfinal-1)
mass2(qfinal) = mass2(qfinal-1)

call binary_output(c1, c2, omsq, hm1, hm2, mass1, mass2, psi, h, qfinal,     &
                   initial_model_type, model_number, ra, za, phia, rb, zb,   &
                   phib, rc, zc, phic, rhm1, rhm2, pin, rhom1, rhom2, xavg1, &
                   xavg2, separation, com, volume_factor, eps)

call output(1000, 'rho', rho)

call cpu_time(time2)

write(13,*) 'Model: ', model_number, ' done in time: ', time2 - time1
close(13)

end subroutine binary_scf

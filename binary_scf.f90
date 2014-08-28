subroutine binary_scf(model_number, initial_model_type, ra, rb, rc, rhom1, rhom2, frac, pin, eps, &
                     qfinal)
implicit none
include 'runscf.h'
!include 'mpif.h'
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

real, dimension(numr_dd,numz_dd,numphi) :: pot, rho
common /poisson/ pot, rho

real, dimension(numr_dd) :: r, rhf, rinv, rhfinv
real, dimension(numz_dd) :: zhf
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

logical :: iam_on_top, iam_on_bottom, iam_on_axis,           &
           iam_on_edge, iam_root
integer :: column_num, row_num
integer :: iam, down_neighbor, up_neighbor,                  &
           in_neighbor, out_neighbor, root,                  &
           REAL_SIZE, INT_SIZE, numprocs
integer, dimension(numr_procs,numz_procs) :: pe_grid
common /processor_grid/ iam, numprocs, iam_on_top,           &
                        iam_on_bottom, iam_on_axis,          &
                        iam_on_edge, down_neighbor,          &
                        up_neighbor, in_neighbor,            &
                        out_neighbor, root, column_num,      &
                        row_num, pe_grid, iam_root,          &
                        REAL_SIZE, INT_SIZE

!
!**************************************************************************************************
!
!    local variables
!

real, dimension(numr_dd, numz_dd, numphi) :: h, pot_it, pot_old, temp

real, dimension(numr_dd, numphi) :: psi

real, dimension(maxit) :: c1, c2, mass1, mass2, omsq, hm1, hm2

real :: cnvgom, cnvgc1, cnvgc2, cnvgh1, cnvgh2

real :: dpot, dpsi, psitmp1, psitmp2, pottmp1, pottmp2

real :: ret1, ret2

real :: global_ret1, global_ret2

real :: xavg1, xavg2, com, separation

real :: virial_error, virial_error1, virial_error2

real :: virial_error_prev

real :: volume_factor

real :: gamma

real :: time1, time2

integer :: za, phia
integer :: zb, phib
integer :: zc, phic

logical :: temp_processor

integer :: ra_pe, rb_pe, rc_pe

logical :: i_have_ra, i_have_rb, i_have_rc

integer :: ra_local_index, rb_local_index, rc_local_index

real :: temp_hm1, temp_hm2

real, dimension(3) :: rhm1, rhm2, temp_rhm1, temp_rhm2

integer :: phi1, phi2, phi3, phi4

integer :: I, J, K, L, Q

integer :: ierror

!integer, dimension(MPI_STATUS_SIZE) :: istatus 

!
!**************************************************************************************************

call cpu_time(time1)

h = 0.0

psi = 0.0

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

i_have_ra = .false.
ra_local_index = 0
ra_pe = 0

if ( iam_on_bottom ) then
   do I = rlwb, rupb
      if ( abs(rhf_g(ra) - rhf(I)) < 0.1 * dr ) then
         i_have_ra = .true.
         ra_local_index = I
      endif
   enddo
endif
!if ( iam_root ) then 
!   temp_processor = .false.
!   if ( i_have_ra ) ra_pe = iam
!   do I = 1, numprocs-1
!      call mpi_recv(temp_processor, 1, MPI_LOGICAL, I, 100+I, MPI_COMM_WORLD, istatus, ierror)
!      if ( temp_processor ) ra_pe = I
!   enddo
!else
!   call mpi_send(i_have_ra, 1, MPI_LOGICAL, root, 100+iam, MPI_COMM_WORLD, ierror)
!endif
!call mpi_bcast(ra_pe, 1, INT_SIZE, root, MPI_COMM_WORLD, ierror)
!call mpi_bcast(ra_local_index, 1, INT_SIZE, ra_pe, MPI_COMM_WORLD, ierror)

i_have_rb = .false.
if ( iam_on_bottom ) then
   do I = rlwb, rupb
      if ( abs(rhf_g(rb) - rhf(I)) < 0.1 * dr ) then
         i_have_rb = .true.
         rb_local_index = I
      endif
   enddo
endif
!if ( iam_root ) then
!   temp_processor = .false.
!   if ( i_have_rb ) rb_pe = iam
!   do I = 1, numprocs-1
!      call mpi_recv(temp_processor, 1, MPI_LOGICAL, I, 100+I, MPI_COMM_WORLD, istatus, ierror)
!      if ( temp_processor ) rb_pe = 1
!   enddo
!else
!   call mpi_send(i_have_rb, 1, MPI_LOGICAL, root, 100+iam, MPI_COMM_WORLD, ierror)
!endif
!call mpi_bcast(rb_pe, 1, INT_SIZE, root, MPI_COMM_WORLD, ierror)
!call mpi_bcast(rb_local_index, 1, INT_SIZE, rb_pe, MPI_COMM_WORLD, ierror)


i_have_rc = .false.
if ( iam_on_bottom ) then
   do I = rlwb, rupb
         if ( abs(rhf_g(rc) - rhf(I)) < 0.1 * dr ) then
          i_have_rc = .true.
          rc_local_index = I
        endif
   enddo
endif
!if ( iam_root ) then
!   temp_processor = .false.
!   if ( i_have_rc ) rc_pe = iam
!   do I = 1, numprocs-1
!     call mpi_recv(temp_processor, 1, MPI_LOGICAL, I, 100+I, MPI_COMM_WORLD, istatus, ierror)
!     if ( temp_processor ) rc_pe = I
!   enddo
!else
!   call mpi_send(i_have_rc, 1, MPI_LOGICAL, root, 100+iam, MPI_COMM_WORLD, ierror)
!endif
!call mpi_bcast(rc_pe, 1, INT_SIZE, root, MPI_COMM_WORLD, ierror)
!call mpi_bcast(rc_local_index, 1, INT_SIZE, rc_pe, MPI_COMM_WORLD, ierror)

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
!call mpi_reduce(ret1, global_ret1, 1, REAL_SIZE, MPI_SUM, root, MPI_COMM_WORLD, ierror)
!call mpi_reduce(ret2, global_ret2, 1, REAL_SIZE, MPI_SUM, root, MPI_COMM_WORLD, ierror)
!call mpi_bcast(global_ret1, 1, REAL_SIZE, root, MPI_COMM_WORLD, ierror)
!call mpi_bcast(global_ret2, 1, REAL_SIZE, root, MPI_COMM_WORLD, ierror)
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
!call mpi_reduce(ret1, global_ret1, 1, REAL_SIZE, MPI_SUM, root, MPI_COMM_WORLD, ierror)
!call mpi_reduce(ret2, global_ret2, 1, REAL_SIZE, MPI_SUM, root, MPI_COMM_WORLD, ierror)
!call mpi_bcast(global_ret1, 1, REAL_SIZE, root, MPI_COMM_WORLD, ierror)
!call mpi_bcast(global_ret2, 1, REAL_SIZE, root, MPI_COMM_WORLD, ierror)
xavg1 = volume_factor * ret1 / mass1(1)
xavg2 = volume_factor * ret2 / mass2(1)
separation = xavg1 - xavg2
com = separation * mass2(1) / ( mass1(1) + mass2(1) )
com = xavg1 - com

if ( iam_root ) then
   open(unit=13,file='iteration_log',form='formatted',status='unknown',position='append')
   write(13,*) iam, 1, mass1(1), mass2(1), xavg1, xavg2, com, separation
endif

do Q = 2, maxit-1                                   ! START OF THE ITERATION CYCLE

   ! solve the Poisson equation for the current density field
   if ( Q == 2 ) then
      call potential_solver(0)
      pot_old = pot
   else
      call potential_solver(1)
   endif
   pot_it = (1.0 - frac) * pot + frac * pot_old
   pot_old = pot

   ! compute the form factor for the centrifugal potential
   do K = philwb, phiupb
      do I = rlwb-1, rupb+1
         psi(I,K) = - 0.5 * ( (rhf(I)*cosine(K) - com)**2 + rhf(I)*rhf(I)*sine(K)*sine(K) )
      enddo
   enddo

   ! calculate the angular frequency
   if ( i_have_ra ) then
      pottmp1 = 0.5 * ( pot_it(ra_local_index,za,phia) + pot_it(ra_local_index-1,za,phia) )
      psitmp1 = 0.5 * ( psi(ra_local_index,phia) + psi(ra_local_index-1,phia) )
   endif
   !call mpi_bcast(pottmp1, 1, REAL_SIZE, ra_pe, MPI_COMM_WORLD, ierror)
   !call mpi_bcast(psitmp1, 1, REAL_SIZE, ra_pe, MPI_COMM_WORLD, ierror)
   if ( i_have_rb ) then
      pottmp2 = 0.5 * ( pot_it(rb_local_index,zb,phib) + pot_it(rb_local_index+1,zb,phib) )
      psitmp2 = 0.5 * ( psi(rb_local_index,phib) + psi(rb_local_index+1,phib) )
   endif
   !call mpi_bcast(pottmp2, 1,  REAL_SIZE, rb_pe, MPI_COMM_WORLD, ierror)
   !call mpi_bcast(psitmp2, 1,  REAL_SIZE, rb_pe, MPI_COMM_WORLD, ierror)
   dpot = pottmp1 - pottmp2
   dpsi = psitmp2 - psitmp1
   omsq(Q) = dpot / dpsi

   ! calculate the two integration constants
   if ( i_have_rb ) then
      pottmp1 = 0.5 * ( pot_it(rb_local_index,zb,phib) + pot_it(rb_local_index+1,zb,phib) )
      psitmp1 = 0.5 * ( psi(rb_local_index,phib) + psi(rb_local_index+1,phib) )
   endif
   !call mpi_bcast(pottmp1, 1, REAL_SIZE, rb_pe, MPI_COMM_WORLD, ierror)
   !call mpi_bcast(psitmp1, 1, REAL_SIZE, rb_pe, MPI_COMM_WORLD, ierror)
   if ( i_have_rc ) then
      pottmp2 = 0.5 * ( pot_it(rc_local_index,zc,phic) + pot_it(rc_local_index+1,zc,phic) )
      psitmp2 = 0.5 * ( psi(rc_local_index,phic) + psi(rc_local_index+1,phic) )
   endif
   !call mpi_bcast(pottmp2, 1, REAL_SIZE, rc_pe, MPI_COMM_WORLD, ierror)
   !call mpi_bcast(psitmp2, 1, REAL_SIZE, rc_pe, MPI_COMM_WORLD, ierror)
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
   !if ( iam_root ) then
      hm1(Q) = temp_hm1
      hm2(Q) = temp_hm2
      !temp_rhm1 = rhm1
      !temp_rhm2 = rhm2
      !do I = 1, numprocs-1
         !call mpi_recv(temp_hm1, 1, REAL_SIZE, I, 100+I, MPI_COMM_WORLD, istatus, ierror)
         !call mpi_recv(temp_hm2, 1, REAL_SIZE, I, 200+I, MPI_COMM_WORLD, istatus, ierror)
         !call mpi_recv(temp_rhm1, 3, REAL_SIZE, I, 300+I, MPI_COMM_WORLD, istatus, ierror)
         !call mpi_recv(temp_rhm2, 3, REAL_SIZE, I, 400+I, MPI_COMM_WORLD, istatus, ierror)
         !if ( temp_hm1 > hm1(Q) ) then
            !hm1(Q) = temp_hm1
            !rhm1 = temp_rhm1
         !endif
         !if ( temp_hm2 > hm2(Q) ) then
            !hm2(Q) = temp_hm2
            !rhm2 = temp_rhm2
         !endif
      !enddo
   !else
      !call mpi_send(temp_hm1, 1, REAL_SIZE, root, 100+iam, MPI_COMM_WORLD, ierror)
      !call mpi_send(temp_hm2, 1, REAL_SIZE, root, 200+iam, MPI_COMM_WORLD, ierror)
      !call mpi_send(rhm1, 3, REAL_SIZE, root, 300+iam, MPI_COMM_WORLD, ierror)
      !call mpi_send(rhm2, 3, REAL_SIZE, root, 400+iam, MPI_COMM_WORLD, ierror)
   !endif
   !call mpi_bcast(hm1(Q), 1, REAL_SIZE, root, MPI_COMM_WORLD, ierror)
   !call mpi_bcast(hm2(Q), 1, REAL_SIZE, root, MPI_COMM_WORLD, ierror)
   !call mpi_bcast(rhm1, 3, REAL_SIZE, root, MPI_COMM_WORLD, ierror)
   !call mpi_bcast(rhm2, 3, REAL_SIZE, root, MPI_COMM_WORLD, ierror)

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
            if ( rhf(I) <= rhf_g(rb) ) then
               rho(I,J,K) = 0.0
            endif
         enddo
      enddo
   enddo
   do K = phi2, phi3
      do J = zlwb, zupb
         do I = rlwb, rupb
            if ( rhf(I) <= rhf_g(rc) ) then
               rho(I,J,K) = 0.0
            endif
         enddo
      enddo
   enddo
   do K = phi4, phiupb
      do J = zlwb, zupb
         do I = rlwb, rupb
            if ( rhf(I) <= rhf_g(rb) ) then
               rho(I,J,K) = 0.0
            endif
         enddo
      enddo
   enddo

   ! impose the equatorial boundary condition
   if ( iam_on_bottom ) then
      do K = philwb, phiupb
         do I = rlwb, rupb
            rho(I,zlwb-1,K) = rho(I,zlwb,K)
         enddo
      enddo
   endif

   ! impose the axial boundary condition
   if ( iam_on_axis ) then
      rho(rlwb-1,:,:) = cshift(rho(rlwb,:,:),dim=2,shift=numphi/2)
   endif

   ! calculate the total mass for each star
   do K = philwb, phiupb
      do J = zlwb, zupb
         do I = rlwb, rupb
            temp(I,J,K) = rhf(I) * rho(I,J,K)
         enddo
      enddo
   enddo
   call binary_sum(temp, ret1, ret2)
   !call mpi_reduce(ret1, global_ret1, 1, REAL_SIZE, MPI_SUM, root, MPI_COMM_WORLD, ierror)
   !call mpi_reduce(ret2, global_ret2, 1, REAL_SIZE, MPI_SUM, root, MPI_COMM_WORLD, ierror)
   !call mpi_bcast(global_ret1, 1, REAL_SIZE, root, MPI_COMM_WORLD, ierror)
   !call mpi_bcast(global_ret2, 1, REAL_SIZE, root, MPI_COMM_WORLD, ierror)
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
   !call mpi_reduce(ret1, global_ret1, 1, REAL_SIZE, MPI_SUM, root, MPI_COMM_WORLD, ierror)
   !call mpi_reduce(ret2, global_ret2, 1, REAL_SIZE, MPI_SUM, root, MPI_COMM_WORLD, ierror)
   !call mpi_bcast(global_ret1, 1, REAL_SIZE, root, MPI_COMM_WORLD, ierror)
   !call mpi_bcast(global_ret2, 1, REAL_SIZE, root, MPI_COMM_WORLD, ierror)
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


   if ( iam_root ) then
      write(13,*) Q, mass1(Q), mass2(Q), xavg1, xavg2, com, omsq(Q), c1(Q), c2(Q), &
                  hm1(Q), hm2(Q), cnvgom, cnvgc1, cnvgc2, cnvgh1, cnvgh2, virial_error1, &
                  virial_error2, virial_error
   endif

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
if ( i_have_ra ) then
   pottmp1 = 0.5 * ( pot_it(ra_local_index,za,phia) + pot_it(ra_local_index-1,za,phia) )
   psitmp1 = 0.5 * ( psi(ra_local_index,phia) + psi(ra_local_index-1,phia) )
endif
!call mpi_bcast(pottmp1, 1, REAL_SIZE, ra_pe, MPI_COMM_WORLD, ierror)
!call mpi_bcast(psitmp1, 1, REAL_SIZE, ra_pe, MPI_COMM_WORLD, ierror)
if ( i_have_rb ) then
   pottmp2 = 0.5 * ( pot_it(rb_local_index,zb,phib) + pot_it(rb_local_index+1,zb,phib) )
   psitmp2 = 0.5 * ( psi(rb_local_index,phib) + psi(rb_local_index+1,phib) )
endif
!call mpi_bcast(pottmp2, 1,  REAL_SIZE, rb_pe, MPI_COMM_WORLD, ierror)
!call mpi_bcast(psitmp2, 1,  REAL_SIZE, rb_pe, MPI_COMM_WORLD, ierror)
dpot = pottmp1 - pottmp2
dpsi = psitmp2 - psitmp1
omsq(qfinal) = dpot / dpsi

! and calculate the two integration constants
if ( i_have_rb ) then
   pottmp1 = 0.5 * ( pot_it(rb_local_index,zb,phib) + pot_it(rb_local_index+1,zb,phib) )
   psitmp1 = 0.5 * ( psi(rb_local_index,phib) + psi(rb_local_index+1,phib) )
endif
!call mpi_bcast(pottmp1, 1, REAL_SIZE, rb_pe, MPI_COMM_WORLD, ierror)
!call mpi_bcast(psitmp1, 1, REAL_SIZE, rb_pe, MPI_COMM_WORLD, ierror)
if ( i_have_rc ) then
   pottmp2 = 0.5 * ( pot_it(rc_local_index,zc,phic) + pot_it(rc_local_index+1,zc,phic) )
   psitmp2 = 0.5 * ( psi(rc_local_index,phic) + psi(rc_local_index+1,phic) )
endif
!call mpi_bcast(pottmp2, 1, REAL_SIZE, rc_pe, MPI_COMM_WORLD, ierror)
!call mpi_bcast(psitmp2, 1, REAL_SIZE, rc_pe, MPI_COMM_WORLD, ierror)
c1(qfinal) = pottmp1 + omsq(qfinal) * psitmp1
c2(qfinal) = pottmp2 + omsq(qfinal) * psitmp2

hm1(qfinal) = hm1(qfinal-1)
hm2(qfinal) = hm2(qfinal-1)
mass1(qfinal) = mass1(qfinal-1)
mass2(qfinal) = mass2(qfinal-1)

if ( iam_root ) then
   write(13,*) iam, 'Model: ', model_number, ' done in time: ', time2 - time1
endif

call binary_output(c1, c2, omsq, hm1, hm2, mass1, mass2, psi, h, qfinal,     &
                   initial_model_type, model_number, ra, za, phia, rb, zb,   &
                   phib, rc, zc, phic, rhm1, rhm2, pin, rhom1, rhom2, xavg1, &
                   xavg2, separation, com, volume_factor, eps)

call cpu_time(time1)

call output(1000, 'density.bin', rho)

call cpu_time(time2)

if ( iam_root ) then
   write(13,*) iam, 'Model: ', model_number, ' disk I/O done in time: ', time2 - time1
   close(13)
endif

end subroutine binary_scf

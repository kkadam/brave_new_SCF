    integer, parameter :: numr = 258
    integer, parameter :: numz = 258
    integer, parameter :: numphi = 512

      ! logical, parameter :: have_green_funcs = .false.

       integer, parameter :: numr_procs = 1
       integer, parameter :: numz_procs = 1

       integer, parameter :: numr_dd = ( (numr - 2)/numr_procs ) + 2     !numr

       integer, parameter :: numz_dd = ( (numz - 2)/numz_procs ) + 2     !numz 

       integer, parameter :: rlwb = 2, rupb = numr_dd - 1       !2, numr-1

       integer, parameter :: zlwb = 2, zupb = numz_dd - 1       !2, numz-1

       integer, parameter :: philwb = 1, phiupb = numphi

       integer, parameter :: numr_dd_z = (numr-2)/numz_procs + 2

       integer, parameter :: numphi_dd = numphi/numr_procs

       integer, parameter :: numphi_by_two = numphi / 2

       integer, parameter :: numphi_by_four = numphi / 4

       real, parameter :: numphiinv = 1.0 / numphi

       integer, parameter :: maxit = 100

! restrictions on the above parameters:
!
!  numphi must be a power of two
!
!  (numr - 2) must be evenly divisible by numr_procs
!
!  (numz - 2) must be evenly divisible by numz_procs
!
!  numphi has to be evenly divisible by numr_procs for the
!         data redistribution in the Poisson solver
!
!  (numr-2) must be evenly divisible by numz_procs for the
!           data redistribution in the Poisson solver
!
!  numr_procs must be greater than or equal to two
!
!  numz_procs must be greater than or equal to two
!
!  numr_procs must be an even number
!
!  numz_procs must be an even number

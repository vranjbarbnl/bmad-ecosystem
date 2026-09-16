program test_dm_interp
  !! Unit tests for dm_interp_mod.f90, exercised entirely through its
  !! public API (no bmad dependency -- compiles standalone, same as
  !! dm_cli). Run via test_dm_bmad_module.py, which also covers the two
  !! things that can't live in this process: the theta<6 hard stop (which
  !! calls `stop 1` and would kill this whole test run) and a large-N
  !! statistical convergence smoke check.
  !!
  !! Every expected value below is either:
  !!   (a) an EXACT grid point of the real tables (dm_table_fccee_clean.npz,
  !!       spin_table_fccee.npz, dm_cdf_table_fccee.npz) -- trilinear
  !!       interpolation reproduces stored grid values exactly (weight 0
  !!       or 1 on every axis), so these are tight-tolerance checks, not
  !!       "close enough" ones; or
  !!   (b) a value derived purely from other already-checked quantities
  !!       plus hand-verified Knuth-algorithm arithmetic (the Poisson
  !!       tests), so no external "golden number" needs to be trusted.
  use dm_interp
  implicit none

  integer :: n_pass = 0, n_fail = 0

  ! Grid-point constants, pulled directly from the source .npz files
  ! (see the conversation this test was written from for the exact
  ! python one-liners used to extract them).
  real(8), parameter :: E1 = 10.0d0, A1 = 0.05d0, T1 = 0.0d0
  real(8), parameter :: DE1 = 1.0d-6
  ! Only the (3,3) entry is nonzero, so row-major vs column-major fill
  ! order doesn't matter here.
  real(8), parameter :: M1(3,3) = reshape( &
    [0.0d0, 0.0d0, 0.0d0,  0.0d0, 0.0d0, 0.0d0,  0.0d0, 0.0d0, -1.0d0], [3,3])

  real(8), parameter :: E2 = 33.59818286283781d0, A2 = 0.2419483326789812d0, T2 = 10.526315789473683d0
  real(8), parameter :: DE2 = 0.011684608507279112d0
  real(8), parameter :: MEAN_T2 = 0.8735352963211952d0
  real(8), parameter :: TQ2_AT_HALF = 0.9124057337411259d0  ! t_quantile(E2,A2,T2, q=0.5)

  real(8), parameter :: TOL_EXACT = 1.0d-9

  ! State for the scripted_rng test double (see below).
  real(8), allocatable :: scripted_values(:)
  integer :: scripted_idx = 0

  call test_energy_grid_point_edge()
  call test_energy_grid_point_interior()
  call test_energy_out_of_range()
  call test_spin_grid_point_edge()
  call test_spin_out_of_range_identity()
  call test_dm_map_physical_bound()
  call test_cdf_in_range()
  call test_poisson_zero_mean_never_calls_rng()
  call test_poisson_known_sequence()
  call test_stochastic_energy_kick_scripted()

  write(*,*)
  write(*,'(A,I0,A,I0,A)') 'RESULTS: ', n_pass, ' passed, ', n_fail, ' failed'
  if (n_fail > 0) stop 1

contains

  subroutine check(cond, name)
    logical, intent(in) :: cond
    character(len=*), intent(in) :: name
    if (cond) then
      n_pass = n_pass + 1
      write(*,'(A,A)') 'PASS: ', name
    else
      n_fail = n_fail + 1
      write(*,'(A,A)') 'FAIL: ', name
    end if
  end subroutine check

  subroutine check_close(actual, expected, tol, name)
    real(8), intent(in) :: actual, expected, tol
    character(len=*), intent(in) :: name
    call check(abs(actual - expected) <= tol, name)
    if (abs(actual - expected) > tol) then
      write(*,'(A,ES16.8,A,ES16.8)') '   actual=', actual, ' expected=', expected
    end if
  end subroutine check_close

  !---------------------------------------------------------------------

  subroutine test_energy_grid_point_edge()
    real(8) :: dE
    call dm_energy_kick(E1, A1, T1, dE)
    call check_close(dE, DE1, TOL_EXACT, 'energy: exact at grid corner (E=10,a0=0.05,theta=0)')
  end subroutine test_energy_grid_point_edge

  subroutine test_energy_grid_point_interior()
    real(8) :: dE
    call dm_energy_kick(E2, A2, T2, dE)
    call check_close(dE, DE2, TOL_EXACT, 'energy: exact at interior grid point')
  end subroutine test_energy_grid_point_interior

  subroutine test_energy_out_of_range()
    real(8) :: dE
    call dm_energy_kick(5.0d0, 0.02d0, 25.0d0, dE)
    call check_close(dE, 0.0d0, TOL_EXACT, 'energy: returns 0 outside grid (below E_min, below a0_min, above theta_max)')
  end subroutine test_energy_out_of_range

  subroutine test_spin_grid_point_edge()
    real(8) :: M(3,3)
    call dm_spin_transfer_matrix(E1, A1, T1, M)
    call check(maxval(abs(M - M1)) <= TOL_EXACT, 'spin: exact matrix at grid corner (E=10,a0=0.05,theta=0)')
  end subroutine test_spin_grid_point_edge

  subroutine test_spin_out_of_range_identity()
    real(8) :: M(3,3), identity(3,3)
    identity = 0.0d0
    identity(1,1) = 1.0d0; identity(2,2) = 1.0d0; identity(3,3) = 1.0d0
    call dm_spin_transfer_matrix(5.0d0, 0.02d0, 25.0d0, M)
    call check(maxval(abs(M - identity)) <= TOL_EXACT, 'spin: identity matrix outside grid')
  end subroutine test_spin_out_of_range_identity

  subroutine test_dm_map_physical_bound()
    !! |s_out| <= 1 must hold even for an unphysical (un-normalized)
    !! s_in, and at points off the exact grid (real interpolation, not
    !! just the exact-grid-point shortcut).
    real(8) :: dE, s_out(3)
    real(8), parameter :: test_points(3,4) = reshape( &
      [45.6d0, 0.2d0, 10.0d0,   20.0d0, 0.5d0, 15.0d0, &
       80.0d0, 0.8d0, 18.0d0,   12.0d0, 0.07d0, 7.0d0], [3,4])
    real(8), parameter :: s_ins(3,3) = reshape( &
      [1.0d0, 1.0d0, 1.0d0,   0.0d0, 0.0d0, 1.0d0,   1.0d0, 0.0d0, 0.0d0], [3,3])
    integer :: i, j
    logical :: all_ok

    all_ok = .true.
    do i = 1, 4
      do j = 1, 3
        call dm_map(test_points(1,i), test_points(2,i), test_points(3,i), s_ins(:,j), dE, s_out)
        if (sqrt(sum(s_out**2)) > 1.0d0 + TOL_EXACT) all_ok = .false.
        if (dE < 0.0d0) all_ok = .false.
      end do
    end do
    call check(all_ok, 'dm_map: |s_out| <= 1 and dE_over_E >= 0 across sample points/spins')
  end subroutine test_dm_map_physical_bound

  subroutine test_cdf_in_range()
    call check(dm_cdf_in_range(45.6d0, 0.2d0, 10.0d0), &
      'cdf_in_range: true at FCC-ee nominal point')
    call check(.not. dm_cdf_in_range(45.6d0, 0.2d0, 3.0d0), &
      'cdf_in_range: false below min theta (3 deg)')
    call check(.not. dm_cdf_in_range(45.6d0, 0.95d0, 10.0d0), &
      'cdf_in_range: false above max usable a0 (0.95 > 0.854)')
    call check(.not. dm_cdf_in_range(200.0d0, 0.2d0, 10.0d0), &
      'cdf_in_range: false above E grid max (200 GeV)')
  end subroutine test_cdf_in_range

  subroutine test_poisson_zero_mean_never_calls_rng()
    integer :: n
    call dm_sample_poisson(0.0d0, poison_rng, n)
    call check(n == 0, 'poisson: mean=0 gives n=0 without calling the RNG')
  end subroutine test_poisson_zero_mean_never_calls_rng

  subroutine test_poisson_known_sequence()
    !! mean_n=5, every uniform draw = 0.5 exactly. Hand-verified via
    !! Knuth's algorithm: L=exp(-5)=6.737947e-3. Running product of 0.5's
    !! is 0.5,0.25,0.125,0.0625,0.03125,0.015625,0.0078125,0.00390625 --
    !! first term <= L is the 8th (0.00390625 < 6.737947e-3), so the loop
    !! incremented n on draws 1-7 and exits on draw 8 without a further
    !! increment. Expected n = 7.
    integer :: n
    call scripted_rng_reset([0.5d0, 0.5d0, 0.5d0, 0.5d0, 0.5d0, 0.5d0, 0.5d0, 0.5d0])
    call dm_sample_poisson(5.0d0, scripted_rng, n)
    call check(n == 7, 'poisson: mean=5, all draws=0.5 gives n=7 (hand-verified)')
  end subroutine test_poisson_known_sequence

  subroutine test_stochastic_energy_kick_scripted()
    !! Force exactly N=1 photon at the E2/A2/T2 grid point via scripted
    !! draws, then check the resulting dE_over_E equals TQ2_AT_HALF --
    !! this exercises mean_t interpolation, the Poisson gate, and
    !! t_quantile interpolation together, all through the public API.
    !!
    !! mean_n at this exact grid point = DE2 / MEAN_T2 (both already
    !! independently checked above / known from the source table), so
    !! L = exp(-mean_n) is computable here rather than re-asserted as a
    !! separate magic number.
    real(8) :: mean_n, L, dE
    integer :: n_photons

    mean_n = DE2 / MEAN_T2
    L = exp(-mean_n)

    ! draw1 continues the Poisson loop (must exceed L), draw2 exits it
    ! with n=1 (must be <= L/draw1), draw3 is consumed as the photon's
    ! energy-fraction quantile q=0.5.
    call check(0.999d0 > L, 'sanity: scripted draw1=0.999 exceeds L (else this test is miscalibrated)')
    call check(0.5d0 <= L / 0.999d0, 'sanity: scripted draw2=0.5 stays <= L/draw1')

    call scripted_rng_reset([0.999d0, 0.5d0, 0.5d0])
    call dm_sample_energy_kick(E2, A2, T2, scripted_rng, dE, n_photons)

    call check(n_photons == 1, 'stochastic: scripted draws force exactly n_photons=1')
    call check_close(dE, TQ2_AT_HALF, TOL_EXACT, &
      'stochastic: dE_over_E with n=1 equals t_quantile(q=0.5) at that grid point')
  end subroutine test_stochastic_energy_kick_scripted

  !---------------------------------------------------------------------
  ! Fake RNGs (test doubles)
  !---------------------------------------------------------------------

  subroutine poison_rng(r)
    !! Fails the test loudly if ever called -- used to assert a code
    !! path that must short-circuit before touching the RNG at all.
    real(8), intent(out) :: r
    write(*,*) 'FAIL: poison_rng was called (should have short-circuited)'
    n_fail = n_fail + 1
    r = 0.5d0
  end subroutine poison_rng

  subroutine scripted_rng_reset(values)
    real(8), intent(in) :: values(:)
    if (allocated(scripted_values)) deallocate(scripted_values)
    allocate(scripted_values(size(values)))
    scripted_values = values
    scripted_idx = 0
  end subroutine scripted_rng_reset

  subroutine scripted_rng(r)
    real(8), intent(out) :: r
    scripted_idx = scripted_idx + 1
    if (scripted_idx > size(scripted_values)) then
      write(*,*) 'FAIL: scripted_rng exhausted (test script needs more values)'
      n_fail = n_fail + 1
      r = 0.5d0
      return
    end if
    r = scripted_values(scripted_idx)
  end subroutine scripted_rng

end program test_dm_interp

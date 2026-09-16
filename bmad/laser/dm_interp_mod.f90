module dm_interp
  !! Density-matrix (beyond-LCFA) lookup-table interpolation for Bmad.
  !!
  !! Reads the two tables exported by DensityMatrix/export_dm_tables_fortran.py
  !! (from dm_table_fccee_clean.npz and spin_table_fccee.npz) and provides:
  !!   - dm_energy_kick(E_GeV, a0, theta_deg, dE_over_E)
  !!   - dm_spin_transfer_matrix(E_GeV, a0, theta_deg, M)
  !!   - dm_map(E_GeV, a0, theta_deg, s_in, dE_over_E, s_out)   [combined, mirrors lcfa_map]
  !!
  !! IMPORTANT -- unit convention: all lookups here take the electron energy
  !! in **GeV**, matching the table's native axis (E_GeV), NOT eV. The
  !! Python density-matrix code has a well-documented history of silent
  !! unit-conversion bugs (see UNITS_FIXED.md / CURRENT_STATUS_SUMMARY.md in
  !! DensityMatrix/) -- callers must convert explicitly, e.g.
  !!   E_GeV = E_eV * 1.0d-9
  !! Do not add an implicit eV->GeV conversion inside this module; keep the
  !! conversion visible at the call site.
  !!
  !! Interpolation method: trilinear (linear in each of E_GeV, a0, theta_deg).
  !! This is a DELIBERATE SIMPLIFICATION relative to the Python reference:
  !!   - dm_lookup.DensityMatrixLookup uses cubic RegularGridInterpolator
  !!     for the energy table.
  !!   - spin_lookup.SpinTransferLookup uses linear RegularGridInterpolator
  !!     per matrix element (so the spin table matches Python exactly in
  !!     interpolation *method*, just reimplemented).
  !! A tensor-product cubic spline was not ported to avoid a much larger,
  !! harder-to-verify Fortran implementation; the resulting extra
  !! interpolation error should be checked against Python with
  !! compare_dm_fortran_python.py before this is trusted for production
  !! (see that script for the actual measured discrepancy).
  !!
  !! Out-of-range behavior matches the Python lookups:
  !!   - energy: returns dE_over_E = 0 if (E_GeV,a0,theta_deg) is outside
  !!     the grid on any axis (matches dm_lookup.py's fill_value=0.0).
  !!   - spin: returns the identity matrix if outside the grid on any axis
  !!     (matches spin_lookup.py's fill_value=1.0 diagonal / 0.0 off-diag).
  !!
  !! --- Stochastic photon-number/energy sampling (dm_sample_energy_kick) ---
  !!
  !! Reads a THIRD table (from export_dm_cdf_table_fortran.py, sourced from
  !! dm_cdf_table_fccee.npz) holding mean_t(E,a0,theta) and the photon
  !! energy-fraction quantile function t_quantile(E,a0,theta,q). These are
  !! the only two quantities exported from that table: emission_weight/
  !! p_emit are DELIBERATELY EXCLUDED because they carry a confirmed ~4x
  !! absolute-normalization bug at theta>=6 deg (and blow up to ~1e14 at
  !! theta<6 deg) -- see CDF_TABLE_README.md's own "absolute p_emit
  !! normalization" caveat, now quantified. mean_t and t_quantile are
  !! normalized moments of the same underlying (buggy-in-scale) spectrum,
  !! so the bad overall prefactor cancels out of them; they were verified
  !! finite and well-behaved for theta>=6 deg, a0<=0.854 before export.
  !!
  !! The photon count is instead derived from the ALREADY-VALIDATED
  !! mean-kick energy table:
  !!     mean_N(E,a0,theta) = dE_over_E(E,a0,theta) / mean_t(E,a0,theta)
  !! then Poisson-sampled, and each sampled photon's energy fraction drawn
  !! from t_quantile via inverse-CDF (uniform draw -> quantile lookup),
  !! following the same pattern Bmad's own bend_photon_init uses for
  !! synchrotron photon energies (bmad/photon/photon_init_mod.f90):
  !! draw a uniform, invert a CDF. Unlike bend_photon_init's closed-form
  !! inverse, ours interpolates a tabulated quantile grid.
  !!
  !! This construction guarantees E[sampled dE/E] = dE_over_E exactly (by
  !! construction: mean_N * mean_t = dE_over_E), so the stochastic path
  !! reproduces the validated deterministic mean-kick table in expectation
  !! while adding genuine photon-count and photon-energy fluctuations --
  !! important for realistic single-photon-driven loss tails, which a
  !! mean-only kick suppresses by construction.
  !!
  !! Restricted to theta_deg >= DM_STOCHASTIC_MIN_THETA_DEG and
  !! a0 <= DM_STOCHASTIC_MAX_A0 (both read from the table file's own
  !! metadata line, not hardcoded): calling dm_sample_energy_kick outside
  !! that range is a HARD STOP, not a silent fallback -- these are known
  !! numerically-broken regions of the underlying calculation, not just
  !! unlikely-to-be-used parameter combinations, and this project's history
  !! (see UNITS_FIXED.md, CURRENT_STATUS_SUMMARY.md) is full of silent
  !! wrong-physics bugs that a loud failure would have caught immediately.
  !!
  !! No RNG is baked in: callers pass a `uniform_rng` procedure conforming
  !! to the `uniform_rng_if` abstract interface below, so this module stays
  !! standalone-testable (dm_cli binds it to the Fortran intrinsic
  !! random_number) while production Bmad code binds it to
  !! random_mod's ran_uniform for consistency with Bmad's own RNG
  !! seeding/engine selection.
  !!
  implicit none
  private

  abstract interface
    subroutine uniform_rng_if(r)
      real(8), intent(out) :: r
    end subroutine uniform_rng_if
  end interface

  type :: grid3
    integer :: nE = 0, nA0 = 0, nTheta = 0
    real(8), allocatable :: E_GeV(:), a0(:), theta_deg(:)
  end type grid3

  type :: energy_table_type
    type(grid3) :: grid
    real(8), allocatable :: dE_over_E(:,:,:)  ! (nE, nA0, nTheta)
  end type energy_table_type

  type :: spin_table_type
    type(grid3) :: grid
    real(8), allocatable :: M(:,:,:,:,:)      ! (nE, nA0, nTheta, 3, 3)
  end type spin_table_type

  type :: cdf_table_type
    type(grid3) :: grid
    integer :: nQ = 0
    real(8), allocatable :: q_grid(:)
    real(8), allocatable :: mean_t(:,:,:)         ! (nE, nA0, nTheta)
    real(8), allocatable :: t_quantile(:,:,:,:)   ! (nE, nA0, nTheta, nQ)
    real(8) :: min_theta_deg = 0.0d0
    real(8) :: max_a0 = 0.0d0
  end type cdf_table_type

  type(energy_table_type), save :: energy_table
  type(spin_table_type), save :: spin_table
  type(cdf_table_type), save :: cdf_table
  logical, save :: energy_loaded = .false.
  logical, save :: spin_loaded = .false.
  logical, save :: cdf_loaded = .false.

  character(len=*), parameter :: default_energy_file = "dm_energy_table.dat"
  character(len=*), parameter :: default_spin_file   = "dm_spin_table.dat"
  character(len=*), parameter :: default_cdf_file    = "dm_cdf_table.dat"

  public :: load_dm_energy_table, load_dm_spin_table, load_dm_cdf_table
  public :: dm_energy_kick, dm_spin_transfer_matrix, dm_map
  public :: dm_in_range
  public :: uniform_rng_if
  public :: dm_sample_poisson, dm_sample_energy_kick, dm_cdf_in_range

contains

  !---------------------------------------------------------------------
  ! Loading
  !---------------------------------------------------------------------

  subroutine read_grid_header(unit, grid)
    integer, intent(in) :: unit
    type(grid3), intent(out) :: grid
    character(len=1024) :: line
    integer :: i

    read(unit, '(A)') line  ! discard leading "#" comment line
    read(unit, *) grid%nE, grid%nA0, grid%nTheta

    allocate(grid%E_GeV(grid%nE))
    allocate(grid%a0(grid%nA0))
    allocate(grid%theta_deg(grid%nTheta))

    do i = 1, grid%nE
      read(unit, *) grid%E_GeV(i)
    end do
    do i = 1, grid%nA0
      read(unit, *) grid%a0(i)
    end do
    do i = 1, grid%nTheta
      read(unit, *) grid%theta_deg(i)
    end do
  end subroutine read_grid_header

  subroutine load_dm_energy_table(filename)
    !! Load the energy-loss table. No-op if already loaded.
    character(len=*), intent(in), optional :: filename
    character(len=256) :: fname
    integer :: unit, ios, i, j, k

    if (energy_loaded) return
    fname = default_energy_file
    if (present(filename)) fname = filename

    open(newunit=unit, file=fname, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,*) "ERROR: cannot open DM energy table: ", trim(fname)
      stop 1
    end if

    call read_grid_header(unit, energy_table%grid)
    allocate(energy_table%dE_over_E(energy_table%grid%nE, energy_table%grid%nA0, &
                                     energy_table%grid%nTheta))

    do i = 1, energy_table%grid%nE
      do j = 1, energy_table%grid%nA0
        do k = 1, energy_table%grid%nTheta
          read(unit, *) energy_table%dE_over_E(i, j, k)
        end do
      end do
    end do
    close(unit)
    energy_loaded = .true.
  end subroutine load_dm_energy_table

  subroutine load_dm_spin_table(filename)
    !! Load the spin transfer-matrix table. No-op if already loaded.
    character(len=*), intent(in), optional :: filename
    character(len=256) :: fname
    integer :: unit, ios, i, j, k
    real(8) :: row(9)

    if (spin_loaded) return
    fname = default_spin_file
    if (present(filename)) fname = filename

    open(newunit=unit, file=fname, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,*) "ERROR: cannot open DM spin table: ", trim(fname)
      stop 1
    end if

    call read_grid_header(unit, spin_table%grid)
    allocate(spin_table%M(spin_table%grid%nE, spin_table%grid%nA0, &
                           spin_table%grid%nTheta, 3, 3))

    do i = 1, spin_table%grid%nE
      do j = 1, spin_table%grid%nA0
        do k = 1, spin_table%grid%nTheta
          read(unit, *) row
          ! row is M11 M12 M13 M21 M22 M23 M31 M32 M33 (row-major)
          spin_table%M(i, j, k, 1, 1:3) = row(1:3)
          spin_table%M(i, j, k, 2, 1:3) = row(4:6)
          spin_table%M(i, j, k, 3, 1:3) = row(7:9)
        end do
      end do
    end do
    close(unit)
    spin_loaded = .true.
  end subroutine load_dm_spin_table

  subroutine load_dm_cdf_table(filename)
    !! Load the stochastic CDF table (mean_t, t_quantile). No-op if
    !! already loaded.
    character(len=*), intent(in), optional :: filename
    character(len=256) :: fname
    character(len=1024) :: line
    integer :: unit, ios, i, j, k
    integer :: nE, nA0, nTheta, nQ

    if (cdf_loaded) return
    fname = default_cdf_file
    if (present(filename)) fname = filename

    open(newunit=unit, file=fname, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,*) "ERROR: cannot open DM CDF table: ", trim(fname)
      stop 1
    end if

    read(unit, '(A)') line  ! discard leading "#" comment line
    read(unit, *) nE, nA0, nTheta, nQ
    read(unit, *) cdf_table%min_theta_deg, cdf_table%max_a0

    cdf_table%grid%nE = nE
    cdf_table%grid%nA0 = nA0
    cdf_table%grid%nTheta = nTheta
    cdf_table%nQ = nQ

    allocate(cdf_table%grid%E_GeV(nE))
    allocate(cdf_table%grid%a0(nA0))
    allocate(cdf_table%grid%theta_deg(nTheta))
    allocate(cdf_table%q_grid(nQ))
    allocate(cdf_table%mean_t(nE, nA0, nTheta))
    allocate(cdf_table%t_quantile(nE, nA0, nTheta, nQ))

    do i = 1, nE
      read(unit, *) cdf_table%grid%E_GeV(i)
    end do
    do i = 1, nA0
      read(unit, *) cdf_table%grid%a0(i)
    end do
    do i = 1, nTheta
      read(unit, *) cdf_table%grid%theta_deg(i)
    end do
    do i = 1, nQ
      read(unit, *) cdf_table%q_grid(i)
    end do

    do i = 1, nE
      do j = 1, nA0
        do k = 1, nTheta
          read(unit, *) cdf_table%mean_t(i, j, k)
        end do
      end do
    end do

    do i = 1, nE
      do j = 1, nA0
        do k = 1, nTheta
          read(unit, *) cdf_table%t_quantile(i, j, k, :)
        end do
      end do
    end do

    close(unit)
    cdf_loaded = .true.
  end subroutine load_dm_cdf_table

  !---------------------------------------------------------------------
  ! Bracketing / trilinear interpolation helpers
  !---------------------------------------------------------------------

  subroutine bracket(x, n, xv, idx, w)
    !! Find idx such that x(idx) <= xv <= x(idx+1) (x assumed sorted
    !! ascending), and the linear weight w in [0,1] for that bracket.
    !! Does NOT clamp xv to the table range -- callers must check
    !! dm_in_range() first if they want the fill_value behavior instead
    !! of edge extrapolation.
    real(8), intent(in) :: x(:)
    integer, intent(in) :: n
    real(8), intent(in) :: xv
    integer, intent(out) :: idx
    real(8), intent(out) :: w

    idx = 1
    do while (idx < n - 1 .and. x(idx + 1) < xv)
      idx = idx + 1
    end do

    if (x(idx + 1) == x(idx)) then
      w = 0.0d0
    else
      w = (xv - x(idx)) / (x(idx + 1) - x(idx))
    end if
    w = max(0.0d0, min(1.0d0, w))
  end subroutine bracket

  logical function in_range_1d(x, n, xv)
    real(8), intent(in) :: x(:)
    integer, intent(in) :: n
    real(8), intent(in) :: xv
    in_range_1d = (xv >= x(1) .and. xv <= x(n))
  end function in_range_1d

  logical function dm_in_range(grid_kind, E_GeV, a0, theta_deg)
    !! Check whether (E_GeV, a0, theta_deg) lies within the grid used by
    !! the requested table. grid_kind = 'energy' or 'spin'.
    character(len=*), intent(in) :: grid_kind
    real(8), intent(in) :: E_GeV, a0, theta_deg
    type(grid3) :: g

    if (grid_kind == 'energy') then
      if (.not. energy_loaded) call load_dm_energy_table()
      g = energy_table%grid
    else
      if (.not. spin_loaded) call load_dm_spin_table()
      g = spin_table%grid
    end if

    dm_in_range = in_range_1d(g%E_GeV, g%nE, E_GeV) .and. &
                  in_range_1d(g%a0, g%nA0, a0) .and. &
                  in_range_1d(g%theta_deg, g%nTheta, theta_deg)
  end function dm_in_range

  real(8) function trilinear(data, iE, iA0, iTheta, wE, wA0, wTheta)
    real(8), intent(in) :: data(:,:,:)
    integer, intent(in) :: iE, iA0, iTheta
    real(8), intent(in) :: wE, wA0, wTheta
    real(8) :: c00, c01, c10, c11, c0, c1

    c00 = data(iE, iA0, iTheta)     * (1.0d0 - wTheta) + data(iE, iA0, iTheta+1)     * wTheta
    c01 = data(iE, iA0+1, iTheta)   * (1.0d0 - wTheta) + data(iE, iA0+1, iTheta+1)   * wTheta
    c10 = data(iE+1, iA0, iTheta)   * (1.0d0 - wTheta) + data(iE+1, iA0, iTheta+1)   * wTheta
    c11 = data(iE+1, iA0+1, iTheta) * (1.0d0 - wTheta) + data(iE+1, iA0+1, iTheta+1) * wTheta

    c0 = c00 * (1.0d0 - wA0) + c01 * wA0
    c1 = c10 * (1.0d0 - wA0) + c11 * wA0

    trilinear = c0 * (1.0d0 - wE) + c1 * wE
  end function trilinear

  !---------------------------------------------------------------------
  ! Public lookups
  !---------------------------------------------------------------------

  subroutine dm_energy_kick(E_GeV, a0, theta_deg, dE_over_E, filename)
    !! Fractional energy loss dE/E. Returns 0 if outside the table's grid
    !! range (matches dm_lookup.py's fill_value=0.0).
    real(8), intent(in) :: E_GeV, a0, theta_deg
    real(8), intent(out) :: dE_over_E
    character(len=*), intent(in), optional :: filename
    integer :: iE, iA0, iTheta
    real(8) :: wE, wA0, wTheta

    if (.not. energy_loaded) call load_dm_energy_table(filename)

    if (.not. dm_in_range('energy', E_GeV, a0, theta_deg)) then
      dE_over_E = 0.0d0
      return
    end if

    call bracket(energy_table%grid%E_GeV, energy_table%grid%nE, E_GeV, iE, wE)
    call bracket(energy_table%grid%a0, energy_table%grid%nA0, a0, iA0, wA0)
    call bracket(energy_table%grid%theta_deg, energy_table%grid%nTheta, theta_deg, iTheta, wTheta)

    dE_over_E = trilinear(energy_table%dE_over_E, iE, iA0, iTheta, wE, wA0, wTheta)
  end subroutine dm_energy_kick

  subroutine dm_spin_transfer_matrix(E_GeV, a0, theta_deg, M, filename)
    !! 3x3 spin transfer matrix such that s_final = M @ s_initial.
    !! Returns the identity matrix if outside the table's grid range
    !! (matches spin_lookup.py's fill_value=1.0 diag / 0.0 off-diag).
    real(8), intent(in) :: E_GeV, a0, theta_deg
    real(8), intent(out) :: M(3,3)
    character(len=*), intent(in), optional :: filename
    integer :: iE, iA0, iTheta, r, c
    real(8) :: wE, wA0, wTheta

    if (.not. spin_loaded) call load_dm_spin_table(filename)

    if (.not. dm_in_range('spin', E_GeV, a0, theta_deg)) then
      M = 0.0d0
      M(1,1) = 1.0d0
      M(2,2) = 1.0d0
      M(3,3) = 1.0d0
      return
    end if

    call bracket(spin_table%grid%E_GeV, spin_table%grid%nE, E_GeV, iE, wE)
    call bracket(spin_table%grid%a0, spin_table%grid%nA0, a0, iA0, wA0)
    call bracket(spin_table%grid%theta_deg, spin_table%grid%nTheta, theta_deg, iTheta, wTheta)

    do r = 1, 3
      do c = 1, 3
        M(r, c) = trilinear(spin_table%M(:,:,:,r,c), iE, iA0, iTheta, wE, wA0, wTheta)
      end do
    end do
  end subroutine dm_spin_transfer_matrix

  subroutine dm_map(E_GeV, a0, theta_deg, s_in, dE_over_E, s_out, energy_file, spin_file)
    !! Combined convenience call mirroring lcfa_map's interface: given the
    !! electron energy [GeV], laser strength a0, crossing angle [deg], and
    !! incoming spin vector, return the fractional energy loss and the
    !! outgoing spin vector. Enforces |s_out| <= 1 (matches
    !! spin_lookup.evolve_spin's safety renormalization).
    real(8), intent(in) :: E_GeV, a0, theta_deg
    real(8), intent(in) :: s_in(3)
    real(8), intent(out) :: dE_over_E
    real(8), intent(out) :: s_out(3)
    character(len=*), intent(in), optional :: energy_file, spin_file
    real(8) :: M(3,3)
    real(8) :: s_mag

    call dm_energy_kick(E_GeV, a0, theta_deg, dE_over_E, energy_file)
    call dm_spin_transfer_matrix(E_GeV, a0, theta_deg, M, spin_file)

    s_out = matmul(M, s_in)

    s_mag = sqrt(sum(s_out**2))
    if (s_mag > 1.0d0) s_out = s_out / s_mag
  end subroutine dm_map

  !---------------------------------------------------------------------
  ! Stochastic photon-number / photon-energy sampling
  !---------------------------------------------------------------------

  logical function dm_cdf_in_range(E_GeV, a0, theta_deg, filename)
    !! Whether (E_GeV, a0, theta_deg) is in the CDF table's *usable* range
    !! -- i.e. within the full grid AND at/above min_theta_deg AND at/below
    !! max_a0 (both trimmed at export time to exclude confirmed-bad
    !! numerical regions; see load_dm_cdf_table / the module header).
    real(8), intent(in) :: E_GeV, a0, theta_deg
    character(len=*), intent(in), optional :: filename

    if (.not. cdf_loaded) call load_dm_cdf_table(filename)

    dm_cdf_in_range = in_range_1d(cdf_table%grid%E_GeV, cdf_table%grid%nE, E_GeV) .and. &
                       theta_deg >= cdf_table%min_theta_deg .and. &
                       theta_deg <= cdf_table%grid%theta_deg(cdf_table%grid%nTheta) .and. &
                       a0 >= cdf_table%grid%a0(1) .and. &
                       a0 <= cdf_table%max_a0
  end function dm_cdf_in_range

  real(8) function cdf_mean_t(E_GeV, a0, theta_deg)
    !! Trilinear-interpolated mean photon energy fraction <t> at
    !! (E_GeV, a0, theta_deg). Caller must have already range-checked.
    real(8), intent(in) :: E_GeV, a0, theta_deg
    integer :: iE, iA0, iTheta
    real(8) :: wE, wA0, wTheta

    call bracket(cdf_table%grid%E_GeV, cdf_table%grid%nE, E_GeV, iE, wE)
    call bracket(cdf_table%grid%a0, cdf_table%grid%nA0, a0, iA0, wA0)
    call bracket(cdf_table%grid%theta_deg, cdf_table%grid%nTheta, theta_deg, iTheta, wTheta)

    cdf_mean_t = trilinear(cdf_table%mean_t, iE, iA0, iTheta, wE, wA0, wTheta)
  end function cdf_mean_t

  real(8) function cdf_t_quantile(E_GeV, a0, theta_deg, q)
    !! Inverse-CDF: photon energy fraction t such that q of the emission
    !! probability lies below it, interpolated quadrilinearly over
    !! (E_GeV, a0, theta_deg, q). Caller must have already range-checked
    !! (E_GeV, a0, theta_deg); q is clamped to [0,1].
    real(8), intent(in) :: E_GeV, a0, theta_deg, q
    integer :: iE, iA0, iTheta, iQ
    real(8) :: wE, wA0, wTheta, wQ, q_clamped
    real(8) :: t_lo, t_hi

    call bracket(cdf_table%grid%E_GeV, cdf_table%grid%nE, E_GeV, iE, wE)
    call bracket(cdf_table%grid%a0, cdf_table%grid%nA0, a0, iA0, wA0)
    call bracket(cdf_table%grid%theta_deg, cdf_table%grid%nTheta, theta_deg, iTheta, wTheta)

    q_clamped = max(0.0d0, min(1.0d0, q))
    call bracket(cdf_table%q_grid, cdf_table%nQ, q_clamped, iQ, wQ)

    t_lo = trilinear(cdf_table%t_quantile(:,:,:,iQ),   iE, iA0, iTheta, wE, wA0, wTheta)
    t_hi = trilinear(cdf_table%t_quantile(:,:,:,iQ+1), iE, iA0, iTheta, wE, wA0, wTheta)

    cdf_t_quantile = t_lo * (1.0d0 - wQ) + t_hi * wQ
  end function cdf_t_quantile

  subroutine dm_sample_poisson(mean_n, uniform_rng, n)
    !! Knuth's algorithm: sample n ~ Poisson(mean_n) using repeated draws
    !! from uniform_rng. O(mean_n) draws -- fine here since mean_n is
    !! always small in the table's usable range (confirmed < 0.3 for the
    !! full FCC-ee-relevant theta>=6 deg, a0<=0.854 grid). mean_n <= 0
    !! always returns n = 0.
    real(8), intent(in) :: mean_n
    procedure(uniform_rng_if) :: uniform_rng
    integer, intent(out) :: n
    real(8) :: l, p, u

    if (mean_n <= 0.0d0) then
      n = 0
      return
    end if

    l = exp(-mean_n)
    n = 0
    p = 1.0d0
    do
      call uniform_rng(u)
      p = p * u
      if (p <= l) exit
      n = n + 1
    end do
  end subroutine dm_sample_poisson

  subroutine dm_sample_energy_kick(E_GeV, a0, theta_deg, uniform_rng, dE_over_E, n_photons, &
                                     overlap_scale, energy_file, cdf_file)
    !! Stochastically sample the fractional energy loss for one crossing:
    !!   mean_N = dE_over_E(E,a0,theta) / mean_t(E,a0,theta)   [* overlap_scale]
    !!   N ~ Poisson(mean_N)
    !!   dE_over_E = sum of N independent draws of t_quantile(uniform())
    !! E[dE_over_E] equals the deterministic dm_energy_kick() value exactly
    !! (times overlap_scale), by construction -- see module header.
    !!
    !! overlap_scale (default 1.0) multiplies mean_N *before* the Poisson
    !! draw, not the sampled total afterward: a partial-overlap crossing
    !! should have a genuinely lower *chance* of any photon at all, not a
    !! uniformly shrunk version of whatever the discrete outcome was --
    !! rescaling after sampling would reintroduce the same
    !! outlier-suppression this whole feature exists to avoid.
    !!
    !! HARD STOPS if (E_GeV, a0, theta_deg) is outside dm_cdf_in_range --
    !! see that function and the module header for why this is not a
    !! silent fallback.
    real(8), intent(in) :: E_GeV, a0, theta_deg
    procedure(uniform_rng_if) :: uniform_rng
    real(8), intent(out) :: dE_over_E
    integer, intent(out) :: n_photons
    real(8), intent(in), optional :: overlap_scale
    character(len=*), intent(in), optional :: energy_file, cdf_file
    real(8) :: dE_mean, mean_t, mean_n, scale, q, t
    integer :: i

    if (.not. energy_loaded) call load_dm_energy_table(energy_file)
    if (.not. cdf_loaded) call load_dm_cdf_table(cdf_file)

    if (.not. dm_cdf_in_range(E_GeV, a0, theta_deg)) then
      write(*,'(A,3ES14.6)') "FATAL (dm_interp): dm_sample_energy_kick called outside the " // &
        "CDF table's usable range (E_GeV, a0, theta_deg) = ", E_GeV, a0, theta_deg
      write(*,'(A,ES14.6,A,ES14.6,A)') "  Usable range requires theta_deg >= ", &
        cdf_table%min_theta_deg, " deg and a0 <= ", cdf_table%max_a0, &
        " (both numerically-broken regions trimmed at export time -- see " // &
        "export_dm_cdf_table_fortran.py)."
      stop 1
    end if

    scale = 1.0d0
    if (present(overlap_scale)) scale = overlap_scale

    call dm_energy_kick(E_GeV, a0, theta_deg, dE_mean, energy_file)
    mean_t = cdf_mean_t(E_GeV, a0, theta_deg)

    if (mean_t <= 0.0d0 .or. dE_mean <= 0.0d0) then
      dE_over_E = 0.0d0
      n_photons = 0
      return
    end if

    mean_n = (dE_mean / mean_t) * scale

    call dm_sample_poisson(mean_n, uniform_rng, n_photons)

    dE_over_E = 0.0d0
    do i = 1, n_photons
      call uniform_rng(q)
      t = cdf_t_quantile(E_GeV, a0, theta_deg, q)
      dE_over_E = dE_over_E + t
    end do
  end subroutine dm_sample_energy_kick

end module dm_interp

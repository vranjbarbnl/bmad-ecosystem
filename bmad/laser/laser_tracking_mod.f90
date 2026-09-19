module laser_tracking_mod
  !! Laser tracking module for Bmad using the full density-matrix (DM)
  !! lookup tables, beyond the LCFA approximation.
  !!
  !! *** THIS VERSION IS DM-ONLY -- THERE IS NO LCFA FALLBACK. ***
  !! It replaces the production LCFA-based laser_tracking_mod.f90 that
  !! lives on the `main` branch of ~/bmad-dev/bmad-ecosystem (the one
  !! actually used for the published LCFA results) -- that checkout is
  !! untouched. This is the `density-matrix` branch of the separate
  !! bmad-ecosystem-dm checkout.
  !!
  !! Hook entry point: laser_track1_preprocess (same name/signature as the
  !! LCFA version, so bmad/custom/track1_preprocess.f90 needs no changes).
  !! Applies a full 3D spin-transfer kick (all of start_orb%spin(1:3), not
  !! just Xi_y) and an energy kick to start_orb%vec(6), using
  !! dm_interp%dm_map. RF/SR remain handled by Bmad, same as the LCFA
  !! version.
  !!
  !! Element identification: any element with a defined custom attribute
  !! "LASER_XI" is treated as a laser kick element (same convention as the
  !! LCFA version; LASER_XI is the normalized vector potential a0).
  !!
  !! Attribute differences from the LCFA version:
  !!   - LASER_THETA_DEG is now REQUIRED (not just optional) on every
  !!     laser element. The DM lookup table is indexed on
  !!     (E_GeV, a0, theta_deg) explicitly -- unlike the LCFA path, there
  !!     is no head-on-collision assumption baked into the chi formula to
  !!     fall back on. A missing LASER_THETA_DEG stops the run rather than
  !!     silently assuming a geometry (this codebase's history is full of
  !!     silent unit/parameter bugs that took days to track down -- see
  !!     DensityMatrix/CURRENT_STATUS_SUMMARY.md and UNITS_FIXED.md; a
  !!     wrong-but-silent table axis is exactly that failure mode).
  !!
  !! Units: dm_interp's dm_map takes E_GeV (not eV) -- the eV->GeV
  !! conversion is done explicitly below and only there.
  !!
  !! Interpolation is trilinear (see dm_interp_mod.f90 for why, and
  !! DensityMatrix/compare_dm_fortran_python.py for the measured
  !! trilinear-vs-Python-cubic discrepancy on the energy table, accepted
  !! as small compared to the DM-vs-LCFA physics uncertainty).
  !!
  !! *** ENERGY RECOIL IS STOCHASTIC (photon-number + photon-energy MC);
  !! SPIN IS NOT. *** The energy kick samples a Poisson-distributed number
  !! of photons and an independent energy fraction per photon (see
  !! dm_interp_mod.f90's dm_sample_energy_kick for the full derivation --
  !! in short: the Poisson mean is taken from the validated mean-kick
  !! table divided by the CDF table's mean_t, sidestepping a confirmed
  !! ~4x normalization bug in that table's own absolute rate). This is
  !! what actually lets rare hard-photon events drive realistic losses,
  !! which the old deterministic mean-kick-every-crossing approach
  !! suppresses by construction.
  !!
  !! Spin is deliberately left as the SAME unconditional mean
  !! dm_spin_transfer_matrix() kick applied every crossing (unchanged from
  !! the non-stochastic version of this file) -- NOT conditioned on how
  !! many photons the energy sampler drew. Gating it on "photon fired"
  !! was considered and rejected: the existing matrix already represents
  !! the quantum-mechanically-averaged spin change over the WHOLE crossing
  !! (including the no-emission branch), so applying it only when N>=1
  !! would silently shrink the mean depolarization by roughly the mean
  !! photon number (often <<1 here), a large unintended deviation from the
  !! validated SPIN_TRANSPORT_VALIDATED.md numbers -- and rescaling by
  !! 1/P(N>=1) to compensate would extrapolate a linear small-signal
  !! transfer matrix far outside where it was derived. Correlating spin
  !! outliers with the specific sampled photon energy requires a genuinely
  !! new (E,a0,theta,t)-conditioned spin table, which does not exist yet
  !! and is tracked as separate follow-on work.
  !!
  !! RNG: uses random_mod's ran_uniform (via the bmad_uniform_rng wrapper
  !! below), for consistency with Bmad's own RNG seeding/engine selection,
  !! not a private RNG stream.
  !!
  !! Runtime requirement: dm_energy_table.dat, dm_spin_table.dat, AND
  !! dm_cdf_table.dat (from DensityMatrix/export_dm_tables_fortran.py and
  !! export_dm_cdf_table_fortran.py) must be present in the working
  !! directory the executable is run from, same convention as
  !! lcfa_kernels.csv was for the LCFA version.
  !!
  !! The CDF table is only usable for theta_deg >= 6 and a0 <= 0.854 (both
  !! read from the table file itself, not hardcoded here) -- calling a
  !! laser element outside that range is a HARD STOP inside
  !! dm_sample_energy_kick, by design (see dm_interp_mod.f90).
  use bmad
  use random_mod, only: ran_uniform
  use dm_interp
  implicit none
  private

  ! Debug control flag - set to .false. to disable all debug output
  logical, parameter :: LASER_DEBUG = .false.

  type :: laser_state
    logical :: init = .false.
    logical :: on = .true.
    integer :: turns_on = 0
    integer :: cooldown_left = 0
    integer :: step_idx = 0
  end type laser_state

  type(laser_state), allocatable, save :: state(:)
  logical, save :: state_ready = .false.

  public :: laser_track1_preprocess

contains

  subroutine ensure_state(ix_ele_needed)
    integer, intent(in) :: ix_ele_needed
    integer :: n_old, n_new
    if (.not. allocated(state)) then
      n_new = max(1, ix_ele_needed+10)
      allocate(state(0:n_new))
      state = laser_state()
      state_ready = .true.
    else if (size(state) <= ix_ele_needed) then
      n_old = size(state)
      n_new = max(ix_ele_needed+10, int(1.2d0*n_old))
      call extend_state(n_new)
    end if
  end subroutine ensure_state

  subroutine extend_state(n_new)
    integer, intent(in) :: n_new
    type(laser_state), allocatable :: tmp(:)
    integer :: n_old
    n_old = size(state)
    allocate(tmp(0:n_new))
    tmp = laser_state()
    tmp(0:n_old-1) = state
    call move_alloc(tmp, state)
  end subroutine extend_state

  subroutine bmad_uniform_rng(r)
    !! Binds dm_interp's uniform_rng_if to Bmad's own RNG (random_mod's
    !! ran_uniform), so photon sampling here uses the same RNG
    !! seeding/engine selection as the rest of the tracking run rather
    !! than a private stream.
    real(8), intent(out) :: r
    call ran_uniform(r)
  end subroutine bmad_uniform_rng

  subroutine laser_track1_preprocess(start_orb, ele, param, err_flag, finished, radiation_included, track)
    use attribute_mod, only: attribute_name
    type(coord_struct), intent(inout) :: start_orb
    type(ele_struct), intent(inout) :: ele
    type(lat_param_struct), intent(in) :: param
    logical, intent(out) :: err_flag
    logical, intent(out) :: finished
    logical, intent(inout) :: radiation_included
    logical, intent(in) :: track

    logical :: has_laser, has_theta
    real(8) :: laser_a0, sigma_t, phi_span, laser_z_center, theta_deg
    integer :: steps, every, off_after, cool_turns
    logical :: scale_overlap
    integer :: turn_idx, idx, step_idx
    real(8) :: delta, E_tot, E_GeV, dE_over_E
    real(8) :: overlap_scale, t_particle, z_relative
    real(8) :: s_in(3), s_out(3), s_out_scaled(3), M(3,3)
    integer :: n_photons

    err_flag = .false.
    finished = .false.

    idx = ele%ix_ele
    call ensure_state(idx)

    ! LASER_XI is the normalized vector potential a0, same convention as
    ! the LCFA version (chi_e = gamma * a0 * ... there; here a0 is a
    ! direct DM table axis).
    has_laser = fetch_attr(ele, "LASER_XI", laser_a0)
    if (.not. has_laser) then
      return  ! not a laser element
    end if

    ! LASER_XI <= 0 means the laser is genuinely off this turn (e.g. during
    ! the OFF phase of a gating ramper) -- a0=0 is below the DM tables' grid
    ! (min a0=0.05), so a genuine "no interaction" case would otherwise hit
    ! dm_sample_energy_kick's hard-stop for being out of range. The LCFA
    ! version had the equivalent short-circuit via its `if (chi <= 0) return`
    ! check; this is the DM equivalent, checked before requiring
    ! LASER_THETA_DEG since a genuinely-off laser doesn't need one this turn.
    if (laser_a0 <= 0.0d0) then
      finished = .true.
      return
    end if

    ! REQUIRED for the DM lookup (see module header) -- no silent default.
    has_theta = fetch_attr(ele, "LASER_THETA_DEG", theta_deg)
    if (.not. has_theta) then
      write(*,'(A,I0,A)') 'FATAL (laser_tracking_mod, DM-only): element ', idx, &
        ' has LASER_XI set but no LASER_THETA_DEG. The density-matrix ' // &
        'lookup table is indexed on (E_GeV, a0, theta_deg) -- there is ' // &
        'no head-on-collision fallback to silently assume. Set ' // &
        'LASER_THETA_DEG [deg] on this element and rerun.'
      err_flag = .true.
      stop 1
    end if

    ! Fetch optional attributes with defaults (unchanged from LCFA version)
    call fetch_attr_default(ele, "LASER_SIGMA_T_S", 1.0d-12, sigma_t)
    call fetch_attr_default(ele, "LASER_PHI_SPAN", 0.0d0, phi_span)
    call fetch_attr_default(ele, "LASER_Z_CENTER", 0.0d0, laser_z_center)
    call fetch_attr_default_int(ele, "LASER_STEPS", 1, steps)
    call fetch_attr_default_int(ele, "LASER_EVERY", 1, every)
    call fetch_attr_default_int(ele, "LASER_OFF_AFTER", 0, off_after)
    call fetch_attr_default_int(ele, "LASER_COOL_TURNS", 0, cool_turns)
    call fetch_attr_default_logical(ele, "LASER_SCALE_OVERLAP", .false., scale_overlap)

    ! Basic gating and sweep state
    ! Note: turn tracking not available in track1_preprocess hook, using static turn_idx = 0
    turn_idx = 0
    call update_state(state(idx), turn_idx, every, off_after, cool_turns, steps)
    if (.not. state(idx)%on) then
      finished = .true.
      return
    end if
    step_idx = state(idx)%step_idx

    ! Compute overlap scaling based on particle's longitudinal position
    ! (unchanged from LCFA version)
    if (scale_overlap) then
      z_relative = start_orb%vec(5) - laser_z_center
      t_particle = z_relative / c_light
      overlap_scale = exp(-0.5d0 * (t_particle*t_particle) / max(1.d-30, sigma_t*sigma_t))
    else
      overlap_scale = 1.0d0
    end if

    ! Electron energy in GeV -- dm_interp's dm_map takes E_GeV explicitly,
    ! never eV. Convert exactly once, here.
    delta = start_orb%vec(6)
    E_tot = (1.0d0 + delta) * ele%value(p0c$)  ! assumes p0c ~ E0 for electrons, in eV
    E_GeV = E_tot * 1.0d-9

    if (E_GeV <= 0.0d0) then
      finished = .true.
      return
    end if

    ! Full 3D spin vector in (all components, unlike the LCFA path which
    ! only tracks the vertical component Xi_y in spin(2)).
    s_in = start_orb%spin(1:3)

    if (LASER_DEBUG .and. idx == ele%ix_ele .and. abs(start_orb%vec(1)) < 1.d-3) then
      write(*,'(A,I0,A,ES12.4,A,ES12.4,A,ES12.4)') 'DM LASER DEBUG: Element ', idx, &
        ', E_GeV=', E_GeV, ', a0=', laser_a0, ', theta_deg=', theta_deg
      write(*,'(A,3ES12.4)') '  spin_in = ', s_in
    end if

    ! Energy: stochastic. overlap_scale multiplies the Poisson MEAN before
    ! sampling (inside dm_sample_energy_kick), not the sampled total
    ! afterward -- see module header for why post-hoc rescaling would
    ! reintroduce outlier suppression.
    call dm_sample_energy_kick(E_GeV, laser_a0, theta_deg, bmad_uniform_rng, &
                                 dE_over_E, n_photons, overlap_scale=overlap_scale)

    ! Spin: deterministic mean kick, unconditional on n_photons -- see
    ! module header for why this is NOT gated on n_photons >= 1.
    call dm_spin_transfer_matrix(E_GeV, laser_a0, theta_deg, M)
    s_out = matmul(M, s_in)
    if (sqrt(sum(s_out**2)) > 1.0d0) s_out = s_out / sqrt(sum(s_out**2))
    s_out_scaled = s_in + (s_out - s_in) * overlap_scale

    if (LASER_DEBUG .and. idx == ele%ix_ele .and. abs(start_orb%vec(1)) < 1.d-3) then
      write(*,'(A,I0,A,ES12.4,A,3ES12.4,A,F8.4)') '  n_photons=', n_photons, &
        ', dE_over_E=', dE_over_E, ', spin_out=', s_out_scaled, ', overlap_scale=', overlap_scale
    end if

    ! Apply energy kick: reduce delta by dE/E.
    start_orb%vec(6) = start_orb%vec(6) - dE_over_E
    start_orb%spin(1:3) = s_out_scaled

    finished = .true.
  end subroutine laser_track1_preprocess

  subroutine update_state(st, turn_idx, every, off_after, cool_turns, steps)
    type(laser_state), intent(inout) :: st
    integer, intent(in) :: turn_idx, every, off_after, cool_turns, steps
    integer :: eff_every
    if (.not. st%init) then
      st%init = .true.
      st%on = .true.
      st%turns_on = 0
      st%cooldown_left = 0
      st%step_idx = 0
    end if

    eff_every = max(1, every)

    if (.not. st%on) then
      if (cool_turns > 0) then
        st%cooldown_left = max(0, st%cooldown_left - 1)
        if (st%cooldown_left == 0) st%on = .true.
      else
        st%on = .true.
      end if
      st%turns_on = 0
    else
      st%turns_on = st%turns_on + 1
      if (off_after > 0 .and. st%turns_on >= off_after) then
        st%on = .false.
        st%cooldown_left = max(0, cool_turns)
      end if
    end if

    if (steps > 0) then
      st%step_idx = modulo(turn_idx / eff_every, max(1, steps))
    else
      st%step_idx = 0
    end if
  end subroutine update_state

  logical function fetch_attr(ele, name, val)
    type(ele_struct), intent(in) :: ele
    character(len=*), intent(in) :: name
    real(8), intent(out) :: val
    type(all_pointer_struct) :: a_ptr
    logical :: err_flag
    val = 0.0d0
    call pointer_to_attribute(ele, name, .false., a_ptr, err_flag)
    if (err_flag) then
      fetch_attr = .false.
    else
      val = a_ptr%r
      fetch_attr = .true.
    end if
  end function fetch_attr

  subroutine fetch_attr_default(ele, name, default, val)
    type(ele_struct), intent(in) :: ele
    character(len=*), intent(in) :: name
    real(8), intent(in) :: default
    real(8), intent(out) :: val
    if (.not. fetch_attr(ele, name, val)) val = default
  end subroutine fetch_attr_default

  subroutine fetch_attr_default_int(ele, name, default, val)
    type(ele_struct), intent(in) :: ele
    character(len=*), intent(in) :: name
    integer, intent(in) :: default
    integer, intent(out) :: val
    real(8) :: tmp
    if (fetch_attr(ele, name, tmp)) then
      val = nint(tmp)
    else
      val = default
    end if
  end subroutine fetch_attr_default_int

  subroutine fetch_attr_default_logical(ele, name, default, val)
    type(ele_struct), intent(in) :: ele
    character(len=*), intent(in) :: name
    logical, intent(in) :: default
    logical, intent(out) :: val
    real(8) :: tmp
    if (fetch_attr(ele, name, tmp)) then
      val = (tmp /= 0.0d0)
    else
      val = default
    end if
  end subroutine fetch_attr_default_logical

end module laser_tracking_mod

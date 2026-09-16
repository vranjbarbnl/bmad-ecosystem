program dm_cli
  !! Simple CLI to evaluate the density-matrix lookup tables, for
  !! cross-checking against the Python reference (dm_lookup.py /
  !! spin_lookup.py) via compare_dm_fortran_python.py.
  !!
  !! --stochastic --n-samples N runs dm_sample_energy_kick N times using
  !! the Fortran intrinsic RNG (via intrinsic_uniform_rng below) and
  !! reports summary statistics, for validate_dm_stochastic.py to check
  !! against the theoretical Poisson-compound-sum moments and against
  !! the deterministic dm_energy_kick() mean.
  use dm_interp
  implicit none

  real(8) :: E_GeV, a0, theta_deg
  real(8) :: s_in(3), s_out(3), dE_over_E
  real(8) :: M(3,3)
  character(len=256) :: energy_file, spin_file, cdf_file
  logical :: have_energy_file, have_spin_file, have_cdf_file
  logical :: stochastic
  integer :: n_samples, seed_arg
  integer :: r, c

  call parse_args(E_GeV, a0, theta_deg, s_in, energy_file, have_energy_file, &
                   spin_file, have_spin_file, cdf_file, have_cdf_file, &
                   stochastic, n_samples, seed_arg)

  if (stochastic) then
    call run_stochastic(E_GeV, a0, theta_deg, n_samples, seed_arg, &
                          energy_file, have_energy_file, cdf_file, have_cdf_file)
    stop
  end if

  if (have_energy_file) then
    call dm_energy_kick(E_GeV, a0, theta_deg, dE_over_E, energy_file)
  else
    call dm_energy_kick(E_GeV, a0, theta_deg, dE_over_E)
  end if

  if (have_spin_file) then
    call dm_spin_transfer_matrix(E_GeV, a0, theta_deg, M, spin_file)
  else
    call dm_spin_transfer_matrix(E_GeV, a0, theta_deg, M)
  end if

  s_out = matmul(M, s_in)
  if (sqrt(sum(s_out**2)) > 1.0d0) s_out = s_out / sqrt(sum(s_out**2))

  write(*,'(A,1X,ES14.6)') "E_GeV", E_GeV
  write(*,'(A,1X,ES14.6)') "a0", a0
  write(*,'(A,1X,ES14.6)') "theta_deg", theta_deg
  write(*,'(A,1X,ES14.6)') "dE_over_E", dE_over_E
  do r = 1, 3
    do c = 1, 3
      write(*,'(A,I0,I0,1X,ES14.6)') "M", r, c, M(r,c)
    end do
  end do
  write(*,'(A,1X,ES14.6)') "s_out_x", s_out(1)
  write(*,'(A,1X,ES14.6)') "s_out_y", s_out(2)
  write(*,'(A,1X,ES14.6)') "s_out_z", s_out(3)

contains

  subroutine intrinsic_uniform_rng(r)
    !! Binds dm_interp's uniform_rng_if to the Fortran intrinsic RNG, for
    !! standalone testing (production Bmad code binds it to
    !! random_mod's ran_uniform instead -- see laser_tracking_mod.f90).
    real(8), intent(out) :: r
    call random_number(r)
  end subroutine intrinsic_uniform_rng

  subroutine run_stochastic(E_GeV, a0, theta_deg, n_samples, seed_arg, &
                              energy_file, have_energy_file, cdf_file, have_cdf_file)
    real(8), intent(in) :: E_GeV, a0, theta_deg
    integer, intent(in) :: n_samples, seed_arg
    character(len=256), intent(in) :: energy_file, cdf_file
    logical, intent(in) :: have_energy_file, have_cdf_file
    real(8) :: dE_over_E_i, dE_mean_deterministic
    integer :: n_photons_i, i
    integer, allocatable :: seed(:)
    integer :: seed_size
    real(8) :: sum_dE, sum_dE2, sum_n, sum_n2
    integer :: max_n_photons

    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = seed_arg
    call random_seed(put=seed)

    if (have_energy_file) then
      call dm_energy_kick(E_GeV, a0, theta_deg, dE_mean_deterministic, energy_file)
    else
      call dm_energy_kick(E_GeV, a0, theta_deg, dE_mean_deterministic)
    end if

    sum_dE = 0.0d0; sum_dE2 = 0.0d0; sum_n = 0.0d0; sum_n2 = 0.0d0
    max_n_photons = 0

    do i = 1, n_samples
      if (have_cdf_file) then
        call dm_sample_energy_kick(E_GeV, a0, theta_deg, intrinsic_uniform_rng, &
                                     dE_over_E_i, n_photons_i, cdf_file=cdf_file)
      else
        call dm_sample_energy_kick(E_GeV, a0, theta_deg, intrinsic_uniform_rng, &
                                     dE_over_E_i, n_photons_i)
      end if
      sum_dE = sum_dE + dE_over_E_i
      sum_dE2 = sum_dE2 + dE_over_E_i**2
      sum_n = sum_n + n_photons_i
      sum_n2 = sum_n2 + real(n_photons_i, 8)**2
      max_n_photons = max(max_n_photons, n_photons_i)
    end do

    write(*,'(A,1X,ES14.6)') "E_GeV", E_GeV
    write(*,'(A,1X,ES14.6)') "a0", a0
    write(*,'(A,1X,ES14.6)') "theta_deg", theta_deg
    write(*,'(A,1X,I0)') "n_samples", n_samples
    write(*,'(A,1X,ES14.6)') "dE_over_E_deterministic", dE_mean_deterministic
    write(*,'(A,1X,ES14.6)') "dE_over_E_sampled_mean", sum_dE / n_samples
    write(*,'(A,1X,ES14.6)') "dE_over_E_sampled_var", sum_dE2/n_samples - (sum_dE/n_samples)**2
    write(*,'(A,1X,ES14.6)') "n_photons_mean", sum_n / n_samples
    write(*,'(A,1X,ES14.6)') "n_photons_var", sum_n2/n_samples - (sum_n/n_samples)**2
    write(*,'(A,1X,I0)') "n_photons_max", max_n_photons
  end subroutine run_stochastic

  subroutine parse_args(E_GeV, a0, theta_deg, s_in, energy_file, have_energy_file, &
                          spin_file, have_spin_file, cdf_file, have_cdf_file, &
                          stochastic, n_samples, seed_arg)
    real(8), intent(out) :: E_GeV, a0, theta_deg, s_in(3)
    character(len=256), intent(out) :: energy_file, spin_file, cdf_file
    logical, intent(out) :: have_energy_file, have_spin_file, have_cdf_file
    logical, intent(out) :: stochastic
    integer, intent(out) :: n_samples, seed_arg
    character(len=256) :: arg
    integer :: i

    E_GeV = -1.0d0
    a0 = -1.0d0
    theta_deg = 0.0d0
    s_in = (/ 0.0d0, 0.0d0, 0.0d0 /)
    energy_file = ""
    spin_file = ""
    cdf_file = ""
    have_energy_file = .false.
    have_spin_file = .false.
    have_cdf_file = .false.
    stochastic = .false.
    n_samples = 10000
    seed_arg = 12345

    i = 1
    do while (i <= command_argument_count())
      call get_command_argument(i, arg)
      select case (arg)
      case ("--E_GeV")
        call get_command_argument(i+1, arg); read(arg, *) E_GeV; i = i + 1
      case ("--a0")
        call get_command_argument(i+1, arg); read(arg, *) a0; i = i + 1
      case ("--theta_deg")
        call get_command_argument(i+1, arg); read(arg, *) theta_deg; i = i + 1
      case ("--sx")
        call get_command_argument(i+1, arg); read(arg, *) s_in(1); i = i + 1
      case ("--sy")
        call get_command_argument(i+1, arg); read(arg, *) s_in(2); i = i + 1
      case ("--sz")
        call get_command_argument(i+1, arg); read(arg, *) s_in(3); i = i + 1
      case ("--energy-table")
        call get_command_argument(i+1, energy_file); have_energy_file = .true.; i = i + 1
      case ("--spin-table")
        call get_command_argument(i+1, spin_file); have_spin_file = .true.; i = i + 1
      case ("--cdf-table")
        call get_command_argument(i+1, cdf_file); have_cdf_file = .true.; i = i + 1
      case ("--stochastic")
        stochastic = .true.
      case ("--n-samples")
        call get_command_argument(i+1, arg); read(arg, *) n_samples; i = i + 1
      case ("--seed")
        call get_command_argument(i+1, arg); read(arg, *) seed_arg; i = i + 1
      case ("-h", "--help")
        call usage()
      end select
      i = i + 1
    end do

    if (E_GeV <= 0.0d0 .or. a0 <= 0.0d0) call usage()
  end subroutine parse_args

  subroutine usage()
    write(*,*) "Usage: dm_cli --E_GeV <E> --a0 <a0> --theta_deg <theta> " // &
               "[--sx <sx> --sy <sy> --sz <sz>] " // &
               "[--energy-table dm_energy_table.dat] [--spin-table dm_spin_table.dat]"
    write(*,*) "       dm_cli --E_GeV <E> --a0 <a0> --theta_deg <theta> --stochastic " // &
               "[--n-samples 10000] [--seed 12345] " // &
               "[--energy-table dm_energy_table.dat] [--cdf-table dm_cdf_table.dat]"
    stop 1
  end subroutine usage

end program dm_cli

module mod_orbital
  use class_orbital, only: orbital
  implicit none
  private

  public :: init_orbital, finish_orbital

  public :: medium

  public :: norbital_compat
  public :: energy_ionize_compat, energy_kinetic_compat
  public :: energy_singlet_compat, energy_triplet_compat
  public :: number_electrons_compat
  public :: number_ionize_compat, gvalue_ionize_compat
  public :: number_singlet_compat, gvalue_singlet_compat
  public :: number_triplet_compat, gvalue_triplet_compat

  type(orbital) :: medium

  integer :: norbital_compat
  
  real, allocatable :: energy_ionize_compat(:), energy_kinetic_compat(:)
  real, allocatable :: energy_singlet_compat(:), energy_triplet_compat(:)
  real, allocatable :: number_electrons_compat(:)

  real, allocatable :: number_ionize_compat(:), gvalue_ionize_compat(:)
  real, allocatable :: number_singlet_compat(:), gvalue_singlet_compat(:)
  real, allocatable :: number_triplet_compat(:), gvalue_triplet_compat(:)

  logical :: initialized = .false.
  logical :: debug = .false.
  
contains

  subroutine init_orbital
    use mod_file_utils, only: get_unused_unit, unit_stdin, unit_stdout
    integer :: unit
    character (len=100) :: file_orbital
    integer :: norbital
    real :: number_density
    namelist /param_orbital/ file_orbital
    namelist /orbital_dim/ norbital
    namelist /param_medium/ number_density

    if (initialized) return

    read(unit_stdin, param_orbital)
    write(unit_stdout, param_orbital)
    read(unit_stdin, param_medium)
    write(unit_stdout, param_medium)
    
    call get_unused_unit(unit)
    open(unit,file=trim(file_orbital))
    read(unit,orbital_dim)
    close(unit)
    write(6,orbital_dim)
    if (debug) write(6,*) norbital

    medium = orbital(norbital, file_orbital)

    medium%number_density = number_density

    ! for backward compatibility
    norbital_compat = norbital
    if (.not.allocated(energy_ionize_compat)) allocate(energy_ionize_compat(norbital_compat))
    if (.not.allocated(energy_kinetic_compat)) allocate(energy_kinetic_compat(norbital_compat))
    if (.not.allocated(energy_singlet_compat)) allocate(energy_singlet_compat(norbital_compat))
    if (.not.allocated(energy_triplet_compat)) allocate(energy_triplet_compat(norbital_compat))
    if (.not.allocated(number_electrons_compat)) allocate(number_electrons_compat(norbital_compat))
    energy_ionize_compat = medium%energy_ionize
    energy_kinetic_compat = medium%energy_kinetic
    energy_singlet_compat = medium%energy_singlet
    energy_triplet_compat = medium%energy_triplet
    number_electrons_compat = medium%number_electrons
    
    if (.not.allocated(number_ionize_compat)) allocate(number_ionize_compat(norbital_compat))
    if (.not.allocated(gvalue_ionize_compat)) allocate(gvalue_ionize_compat(norbital_compat))
    if (.not.allocated(number_singlet_compat)) allocate(number_singlet_compat(norbital_compat))
    if (.not.allocated(gvalue_singlet_compat)) allocate(gvalue_singlet_compat(norbital_compat))
    if (.not.allocated(number_triplet_compat)) allocate(number_triplet_compat(norbital_compat))
    if (.not.allocated(gvalue_triplet_compat)) allocate(gvalue_triplet_compat(norbital_compat))
    
    initialized = .true.
  end subroutine init_orbital
  
  subroutine finish_orbital
    if (.not.initialized) return
    if (allocated(energy_ionize_compat)) deallocate(energy_ionize_compat)
    if (allocated(energy_kinetic_compat)) deallocate(energy_kinetic_compat)
    if (allocated(energy_singlet_compat)) deallocate(energy_singlet_compat)
    if (allocated(energy_triplet_compat)) deallocate(energy_triplet_compat)
    if (allocated(number_electrons_compat)) deallocate(number_electrons_compat)

    if (allocated(number_ionize_compat)) deallocate(number_ionize_compat)
    if (allocated(gvalue_ionize_compat)) deallocate(gvalue_ionize_compat)
    if (allocated(number_singlet_compat)) deallocate(number_singlet_compat)
    if (allocated(gvalue_singlet_compat)) deallocate(gvalue_singlet_compat)
    if (allocated(number_triplet_compat)) deallocate(number_triplet_compat)
    if (allocated(gvalue_triplet_compat)) deallocate(gvalue_triplet_compat)
    initialized = .false.
  end subroutine finish_orbital
  
end module mod_orbital

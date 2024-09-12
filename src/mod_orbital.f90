module mod_orbital
  use class_orbital, only: orbital
  implicit none
  private

  public :: init_orbital, finish_orbital

  public :: medium

  type(orbital) :: medium

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
    
    call get_unused_unit(unit)
    open(unit, file=trim(file_orbital))
    read(unit, orbital_dim)
    read(unit, param_medium)
    close(unit)
    write(unit_stdout, orbital_dim)
    write(unit_stdout, param_medium)
    if (debug) write(unit_stdout, *) norbital

    medium = orbital(norbital, file_orbital)

    medium%number_density = number_density

    initialized = .true.
  end subroutine init_orbital
  
  subroutine finish_orbital
    if (.not.initialized) return
    initialized = .false.
  end subroutine finish_orbital
  
end module mod_orbital

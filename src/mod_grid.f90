module mod_grid
  use class_grid, only: grid
  implicit none
  private
  public :: init_grid, finish_grid

  public :: egrid
  
  type(grid) :: egrid
  
  logical :: initialized = .false.
  logical :: debug = .false.
  
contains

  subroutine init_grid
    use mod_file_utils, only: unit_stdin, unit_stdout

    ! number of division in energy grid
    ! nediv means T of i+nediv is a half of T of i
    integer :: nediv
    ! egrid_max, egrid_min are the maximum, minimum energy
    real :: egrid_max, egrid_min
    
    namelist /grid_param/ nediv, egrid_max, egrid_min

    if (initialized) return

    read(unit_stdin, grid_param)
    write(unit_stdout, grid_param)

    egrid = grid(nediv, egrid_max, egrid_min)
    
    initialized = .true.
  end subroutine init_grid

  subroutine finish_grid
    if (.not.initialized) return
    initialized = .false.
  end subroutine finish_grid
  
end module mod_grid

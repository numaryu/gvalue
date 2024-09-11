module mod_grid
  use class_grid, only: grid
  implicit none
  private
  public :: init_grid, finish_grid

  public :: egrid_compat
  public :: negrid_compat
  public :: egrid_max_compat
  public :: mesh_compat

  public :: egrid
  
  type(grid) :: egrid
  
  logical :: initialized = .false.
  logical :: debug = .false.
  
  ! values for backward compatibility
  ! degrid_compat = (1/2)**(1/nediv) = T_(i+1)/T_(i) with nediv=40
  real :: degrid_compat = 0.5**(1.0/40)
  real :: egrid_max_compat = 1.e5, egrid_min_compat = 1.e-2
  integer, parameter :: negrid_compat = 1500
  
  real :: egrid_compat(negrid_compat)
  
contains

  subroutine init_grid
    use mod_file_utils, only: unit_stdin, unit_stdout
    integer :: i

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
    
    do i=1,negrid_compat
       egrid_compat(i) = egrid_max_compat*degrid_compat**(i-1)
    end do

    write(6,*) 'grid diff = ',sum(abs(egrid_compat(1:min(negrid_compat,egrid%number))-egrid%val(1:min(negrid_compat,egrid%number))))
    write(6,*) 'grid diff = ',negrid_compat,egrid%number
    initialized = .true.
  end subroutine init_grid

  subroutine finish_grid
    if (.not.initialized) return
    initialized = .false.
  end subroutine finish_grid

  integer function mesh_compat(tx,tp)
    !
    !     Functtion of numbering electron energy
    !
    real, intent(in) :: tx, tp
    mesh_compat = ifix(40.0*(alog(tp/tx)/alog(2.0))+1)
    return
  end function mesh_compat
  
end module mod_grid

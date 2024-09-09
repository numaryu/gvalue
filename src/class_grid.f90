module class_grid
  implicit none
  private

  public :: grid
  
  ! log-spaced grid type in descending order from val_max to val_min
  type grid
     private
     ! number of grids
     integer, public :: number = 0
     ! nhalf means T of i+nediv is a half of T of i
     integer :: nhalf = 0
     ! ratio of neiboring grids Ti/Ti+1
     real :: div = 0.
     ! grids
     real, pointer, public :: val(:) => null()
     ! maximum and minimumx of grids
     real, public :: val_max = 0., val_min = 0.

   contains
     final :: destroy
     procedure :: grid_number
  end type grid

  ! declare constructor
  interface grid
     module procedure init_grid
  end interface grid

contains

  type(grid) function init_grid(nhalf, val_max, val_min)
    integer, intent(in) :: nhalf
    real, intent(in) :: val_max, val_min
    integer :: i
    init_grid%nhalf = nhalf
    init_grid%val_max = val_max
    init_grid%val_min = val_min

    init_grid%div = 0.5**(1.0/nhalf)
    init_grid%number = grid_number(init_grid, val_min)

    if (.not.associated(init_grid%val)) allocate(init_grid%val(init_grid%number))
    do i=1,init_grid%number
       init_grid%val(i) = val_max*init_grid%div**(i-1)
    end do
  end function init_grid

  subroutine destroy(self)
    type(grid), intent(in out) :: self
    if (associated(self%val)) nullify(self%val)
  end subroutine destroy
  
  integer function grid_number(self, val)
    class(grid) :: self
    real, intent(in) :: val
    grid_number = int(self%nhalf*(log(self%val_max/val)/log(2.))+1)
    return
  end function grid_number

end module class_grid

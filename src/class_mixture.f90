module class_mixture
  implicit none
  private

  public :: mixture

  type mixture
     private

     ! mixture of stopping power S(energy)
     real, pointer, public :: stop_power_mixture(:) => null()

     ! mixture of degradation y(energy)
     real, pointer, public :: degradation_mixture(:) => null()

     ! mixture of mean free path(energy)
     real, pointer, public :: mean_free_path_mixture(:) => null()

   contains
     final :: destroy
  end type mixture

  interface mixture
     module procedure init_mixture
  end interface mixture

contains

  type(mixture) function init_mixture(ngrid, ngeneration)
    integer, intent(in) :: ngrid, ngeneration

    if (.not.associated(init_mixture%stop_power_mixture)) allocate(init_mixture%stop_power_mixture(ngrid))
    init_mixture%stop_power_mixture=0.
    if (.not.associated(init_mixture%degradation_mixture)) allocate(init_mixture%degradation_mixture(ngrid))
    if (.not.associated(init_mixture%mean_free_path_mixture)) &
         allocate(init_mixture%mean_free_path_mixture(ngrid))
  end function init_mixture

  subroutine destroy(self)
    type(mixture), intent(in out) :: self
    if (associated(self%stop_power_mixture)) nullify(self%stop_power_mixture)
    if (associated(self%degradation_mixture)) nullify(self%degradation_mixture)
    if (associated(self%mean_free_path_mixture)) nullify(self%mean_free_path_mixture)
  end subroutine destroy

end module class_mixture

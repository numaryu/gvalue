module class_medium
  implicit none
  private

  public :: medium

  type medium
     private
     ! name of medium
     character (len=100), public :: name = 'unknown'
     ! number density
     real, public :: number_density

     ! number of orbital
     integer, public :: number = 0
     ! ionization Ii, kinetic Ei, singlet Esi, triplet Eti energies (per orbital)
     real, pointer, public :: energy_ionize(:) => null(), energy_kinetic(:) => null()
     real, pointer, public :: energy_singlet(:) => null(), energy_triplet(:) => null()
     ! number of electrons (per orbital)
     integer, pointer, public :: number_electrons(:) => null()

     ! yield Ni and G value of processes (per orbital)
     real, pointer, public :: yield_ionize(:) => null(), gvalue_ionize(:) => null()
     real, pointer, public :: yield_singlet(:) => null(), gvalue_singlet(:) => null()
     real, pointer, public :: yield_triplet(:) => null(), gvalue_triplet(:) => null()

     ! stopping power S(energy), Si(energy) (per orbital)
     real, pointer :: stop_power_orbital(:,:) => null()
     real, pointer, public :: stop_power(:) => null()

     ! degradation y(generation, energy)
     real, pointer, public :: degradation_gen(:,:) => null()

     ! total cross section Qion(energy), Qsngl(energy), Qtrpl(energy)
     !                     Qion,i(energy), Qsngl,i(energy), Qtripl,i(energy) (per orbital)
     real, pointer, public :: total_cross_section_ionize_orbital(:,:) => null()
     real, pointer, public :: total_cross_section_singlet_orbital(:,:) => null()
     real, pointer, public :: total_cross_section_triplet_orbital(:,:) => null()
     real, pointer, public :: total_cross_section_ionize(:) => null()
     real, pointer, public :: total_cross_section_singlet(:) => null()
     real, pointer, public :: total_cross_section_triplet(:) => null()
     real, pointer, public :: total_cross_section_total(:) => null()

     ! platzman energy*y(energy)*Q(energy) (per orbital)
     real, pointer, public :: platzman_ionize_orbital(:,:) => null()
     real, pointer, public :: platzman_singlet_orbital(:,:) => null()
     real, pointer, public :: platzman_triplet_orbital(:,:) => null()

     ! range
     real, pointer, public :: range_ionize(:) => null()
     real, pointer, public :: range_singlet(:) => null()
     real, pointer, public :: range_triplet(:) => null()

   contains
     final :: destroy
     procedure, public :: init_orbital_vars
     procedure :: calculate_stopping_power
     procedure :: calculate_degradation
     procedure :: calculate_yield
  end type medium

  ! declare constructor
  interface medium
     module procedure init_medium
  end interface medium

contains

  type(medium) function init_medium(name, file_medium, number_density)
    character (len=*), intent(in) :: name
    character (len=*), intent(in) :: file_medium
    real, intent(in) :: number_density

    init_medium%name = name
    init_medium%number_density = number_density

    call read_params(init_medium, file_medium)

    if (.not.associated(init_medium%yield_ionize)) allocate(init_medium%yield_ionize(init_medium%number))
    if (.not.associated(init_medium%gvalue_ionize)) allocate(init_medium%gvalue_ionize(init_medium%number))
    if (.not.associated(init_medium%yield_singlet)) allocate(init_medium%yield_singlet(init_medium%number))
    if (.not.associated(init_medium%gvalue_singlet)) allocate(init_medium%gvalue_singlet(init_medium%number))
    if (.not.associated(init_medium%yield_triplet)) allocate(init_medium%yield_triplet(init_medium%number))
    if (.not.associated(init_medium%gvalue_triplet)) allocate(init_medium%gvalue_triplet(init_medium%number))

  end function init_medium

  subroutine destroy(self)
    type(medium), intent(in out) :: self
    if (associated(self%energy_ionize)) nullify(self%energy_ionize)
    if (associated(self%energy_kinetic)) nullify(self%energy_kinetic)
    if (associated(self%energy_singlet)) nullify(self%energy_singlet)
    if (associated(self%energy_triplet)) nullify(self%energy_triplet)
    if (associated(self%number_electrons)) nullify(self%number_electrons)

    if (associated(self%yield_ionize)) nullify(self%yield_ionize)
    if (associated(self%gvalue_ionize)) nullify(self%gvalue_ionize)
    if (associated(self%yield_singlet)) nullify(self%yield_singlet)
    if (associated(self%gvalue_singlet)) nullify(self%gvalue_singlet)
    if (associated(self%yield_triplet)) nullify(self%yield_triplet)
    if (associated(self%gvalue_triplet)) nullify(self%gvalue_triplet)

    call finish_orbital_vars(self)
  end subroutine destroy

  subroutine read_params(self, file)
    use mod_file_utils, only: get_unused_unit
    class(medium) :: self
    integer, parameter :: norbital_max = 1000
    character (len=*), intent(in) :: file

    integer :: unit
    
    real :: energy_ionize(norbital_max) = 0.
    real :: energy_kinetic(norbital_max) = 0.
    real :: energy_singlet(norbital_max) = 0.
    real :: energy_triplet(norbital_max) = 0.
    integer :: number_electrons(norbital_max) = 0
    
    namelist /params_per_orbitals/ energy_ionize, energy_kinetic, &
         energy_singlet, energy_triplet, number_electrons
    
    call get_unused_unit(unit)
    open(unit,file=trim(file))
    read(unit,params_per_orbitals)
    close(unit)

    self%number = count(number_electrons /= 0)

    if (.not.associated(self%energy_ionize)) allocate(self%energy_ionize(self%number))
    if (.not.associated(self%energy_kinetic)) allocate(self%energy_kinetic(self%number))
    if (.not.associated(self%energy_singlet)) allocate(self%energy_singlet(self%number))
    if (.not.associated(self%energy_triplet)) allocate(self%energy_triplet(self%number))
    if (.not.associated(self%number_electrons)) allocate(self%number_electrons(self%number))

    self%energy_ionize(1:self%number) = energy_ionize(1:self%number)
    self%energy_kinetic(1:self%number) = energy_kinetic(1:self%number)
    self%energy_singlet(1:self%number) = energy_singlet(1:self%number)
    self%energy_triplet(1:self%number) = energy_triplet(1:self%number)
    self%number_electrons(1:self%number) = number_electrons(1:self%number)
  end subroutine read_params

  subroutine init_orbital_vars(self, ngrid, ngeneration)
    class(medium) :: self
    integer, intent(in) :: ngrid, ngeneration
    if (.not.associated(self%stop_power_orbital)) allocate(self%stop_power_orbital(self%number, ngrid))
    self%stop_power_orbital = 0.0
    if (.not.associated(self%stop_power)) allocate(self%stop_power(ngrid))
!!!    self%stop_power = 0.0 ! not necessary

    if (.not.associated(self%degradation_gen)) allocate(self%degradation_gen(ngeneration, ngrid))
    self%degradation_gen = 0.0

    if (.not.associated(self%total_cross_section_ionize_orbital)) &
         allocate(self%total_cross_section_ionize_orbital(self%number, ngrid))
    if (.not.associated(self%total_cross_section_singlet_orbital)) &
         allocate(self%total_cross_section_singlet_orbital(self%number, ngrid))
    if (.not.associated(self%total_cross_section_triplet_orbital)) &
         allocate(self%total_cross_section_triplet_orbital(self%number, ngrid))
    self%total_cross_section_ionize_orbital = 0.0
    self%total_cross_section_singlet_orbital = 0.0
    self%total_cross_section_triplet_orbital = 0.0
    
    if (.not.associated(self%total_cross_section_ionize)) &
         allocate(self%total_cross_section_ionize(ngrid))
    if (.not.associated(self%total_cross_section_singlet)) &
         allocate(self%total_cross_section_singlet(ngrid))
    if (.not.associated(self%total_cross_section_triplet)) &
         allocate(self%total_cross_section_triplet(ngrid))
    if (.not.associated(self%total_cross_section_total)) &
         allocate(self%total_cross_section_total(ngrid))
!!!    self%total_cross_section_ionize = 0.0 ! not necessary
!!!    self%total_cross_section_singlet = 0.0 ! not necessary
!!!    self%total_cross_section_triplet = 0.0 ! not necessary

    if (.not.associated(self%platzman_ionize_orbital)) &
         allocate(self%platzman_ionize_orbital(self%number, ngrid))
    if (.not.associated(self%platzman_singlet_orbital)) &
         allocate(self%platzman_singlet_orbital(self%number, ngrid))
    if (.not.associated(self%platzman_triplet_orbital)) &
         allocate(self%platzman_triplet_orbital(self%number, ngrid))
    self%platzman_ionize_orbital = 0.0
    self%platzman_singlet_orbital = 0.0
    self%platzman_triplet_orbital = 0.0

    if (.not.associated(self%range_ionize)) allocate(self%range_ionize(ngrid))
    if (.not.associated(self%range_singlet)) allocate(self%range_singlet(ngrid))
    if (.not.associated(self%range_triplet)) allocate(self%range_triplet(ngrid))
    self%range_ionize = 0.0
    self%range_singlet = 0.0
    self%range_triplet = 0.0

  end subroutine init_orbital_vars

  subroutine finish_orbital_vars(self)
    class(medium) :: self
    if (associated(self%stop_power_orbital)) nullify(self%stop_power_orbital)
    if (associated(self%stop_power)) nullify(self%stop_power)

    if (associated(self%degradation_gen)) nullify(self%degradation_gen)

    if (associated(self%total_cross_section_ionize_orbital)) &
         nullify(self%total_cross_section_ionize_orbital)
    if (associated(self%total_cross_section_singlet_orbital)) &
         nullify(self%total_cross_section_singlet_orbital)
    if (associated(self%total_cross_section_triplet_orbital)) &
         nullify(self%total_cross_section_triplet_orbital)

    if (associated(self%total_cross_section_ionize)) nullify(self%total_cross_section_ionize)
    if (associated(self%total_cross_section_singlet)) nullify(self%total_cross_section_singlet)
    if (associated(self%total_cross_section_triplet)) nullify(self%total_cross_section_triplet)
    if (associated(self%total_cross_section_total)) nullify(self%total_cross_section_total)

    if (associated(self%platzman_ionize_orbital)) nullify(self%platzman_ionize_orbital)
    if (associated(self%platzman_singlet_orbital)) nullify(self%platzman_singlet_orbital)
    if (associated(self%platzman_triplet_orbital)) nullify(self%platzman_triplet_orbital)

    if (associated(self%range_ionize)) nullify(self%range_ionize)
    if (associated(self%range_singlet)) nullify(self%range_singlet)
    if (associated(self%range_triplet)) nullify(self%range_triplet)
  end subroutine finish_orbital_vars

  subroutine calculate_stopping_power(self, egrid)
    use class_grid, only: grid
    class(medium) :: self
    class(grid), intent(in) :: egrid
    integer :: io, ie
    integer :: ne1, ne2, ne3

    do io = 1, self%number ! orbital index

       if (self%energy_ionize(io) < egrid%val_max) then
          ne1 = egrid%grid_number(self%energy_ionize(io))

          ! energy >=  energy_ionize
          do ie = 1, ne1
             self%stop_power_orbital(io, ie) = &
                  integrate_E_sigma(self, "total", io, &
                  egrid%val(ie), self%energy_singlet(io), 0.5*(egrid%val(ie) + self%energy_ionize(io))) &
                  + 0.5*integrate_E_sigma(self, "exchange", io, &
                  egrid%val(ie), self%energy_triplet(io), self%energy_singlet(io))
          end do
       else
          ne1 = 0
       end if

       if (self%energy_singlet(io) < egrid%val_max) then
          ne2 = egrid%grid_number(self%energy_singlet(io))

          ! energy_ionize > energy >=  energy_singlet
          do ie = ne1+1, ne2
             self%stop_power_orbital(io, ie) = &
                  integrate_E_sigma(self, "total", io, &
                  egrid%val(ie), self%energy_singlet(io), egrid%val(ie)) + &
                  0.5*integrate_E_sigma(self, "exchange", io, &
                  egrid%val(ie), self%energy_triplet(io), self%energy_singlet(io))
          end do
       else
          ne2 = 0
       end if

       if (self%energy_triplet(io) < egrid%val_max) then
          ne3 = egrid%grid_number(self%energy_triplet(io))

          ! energy_singlet > energy >=  energy_triplet
          do ie = ne2+1, ne3
             self%stop_power_orbital(io, ie) = 0.5*integrate_E_sigma(self, "exchange", io, &
                  egrid%val(ie), self%energy_triplet(io), egrid%val(ie))
          end do
       end if

    end do

    self%stop_power(:) = sum(self%stop_power_orbital(:,:), dim=1)
    
  end subroutine calculate_stopping_power

  real function integrate_E_sigma(self, type, io, energy, range_min, range_max)
    class(medium) :: self
    character (len=*) :: type
    integer, intent(in) :: io
    real, intent(in) :: energy
    real, intent(in) :: range_min, range_max
    real :: val
    real :: energy_ionize, energy_kinetic

    abstract interface
       real function indefinite_integral(x, energy, energy_ionize, energy_kinetic)
         real, intent(in) :: x
         real, intent(in) :: energy, energy_ionize, energy_kinetic
       end function indefinite_integral
    end interface
      
    procedure (indefinite_integral), pointer :: indefinite_E_sigma => null()
      
    select case(type)
    case ("total")
       indefinite_E_sigma => indefinite_E_sigma_Etot
    case ("direct")
       indefinite_E_sigma => indefinite_E_sigma_Edir
    case ("exchange")
       indefinite_E_sigma => indefinite_E_sigma_Eexc
    end select
    
    energy_ionize = self%energy_ionize(io)
    energy_kinetic = self%energy_kinetic(io)

    val = indefinite_E_sigma(range_max, energy, energy_ionize, energy_kinetic) &
         - indefinite_E_sigma(range_min, energy, energy_ionize, energy_kinetic)
    integrate_E_sigma = val * self%number_electrons(io) * self%number_density

    return

  contains
    
    real function indefinite_E_sigma_Etot(x, energy, energy_ionize, energy_kinetic)
      real, intent(in) :: x
      real, intent(in) :: energy, energy_ionize, energy_kinetic
      real :: val
      val = indefinite_E_sigma_Edir(x, energy, energy_ionize, energy_kinetic) + &
           indefinite_E_sigma_Eexc(x, energy, energy_ionize, energy_kinetic)
      indefinite_E_sigma_Etot = val
      return 
    end function indefinite_E_sigma_Etot
    
    real function indefinite_E_sigma_Edir(x, energy, energy_ionize, energy_kinetic)
      use mod_constants, only: bb
      real, intent(in) :: x
      real, intent(in) :: energy, energy_ionize, energy_kinetic
      real :: val

      val = bb/(energy + energy_ionize + energy_kinetic) * &
           ( log(x) - 4./3.*energy_kinetic/x )

      indefinite_E_sigma_Edir = val
      return
    end function indefinite_E_sigma_Edir

    real function indefinite_E_sigma_Eexc(x, energy, energy_ionize, energy_kinetic)
      use mod_constants, only: bb
      real, intent(in) :: x
      real, intent(in) :: energy, energy_ionize, energy_kinetic
      real :: val

      val = bb/(energy + energy_ionize + energy_kinetic) * &
           ( log(energy + energy_ionize - x) &
           + (energy + energy_ionize)/(energy + energy_ionize - x) &
           + 4./3.*energy_kinetic* &
           ( - 1./(energy + energy_ionize - x) + &
           + 0.5*(energy + energy_ionize)/(energy + energy_ionize - x)**2 )&
           )

      indefinite_E_sigma_Eexc = val
      return
    end function indefinite_E_sigma_Eexc

  end function integrate_E_sigma

  subroutine calculate_degradation(self, egrid, mediamix, igen)
    use mod_constants, only: bb
    use class_grid, only: grid
    use class_mixture, only: mixture
    class(medium) :: self
    class(grid), intent(in) :: egrid
    class(mixture), intent(in) :: mediamix
    integer, intent(in) :: igen
    real :: energy1, denergy1, energy2, denergy2
    real :: energy_ionize, energy_kinetic
    integer :: io
    integer :: i, j
    real :: energy1_max, energy2_max
    integer :: nenergy_max
    integer :: nenergy1_max, nenergy2_max

    nenergy_max = egrid%grid_number(minval(self%energy_triplet))

    do io = 1, self%number
       energy_ionize = self%energy_ionize(io)
       energy_kinetic = self%energy_kinetic(io)

       energy1_max = (egrid%val_max - (2**(igen-2)-1)*energy_ionize)/(2**(igen-2))
       if (igen == 2) energy1_max = egrid%val_max
       energy2_max = (egrid%val_max - (2**(igen-1)-1)*energy_ionize)/(2**(igen-1))
       if (energy2_max < minval(self%energy_triplet)) cycle

       nenergy1_max = egrid%grid_number(energy1_max)
       nenergy2_max = egrid%grid_number(energy2_max)

       ! i: index of T2
       do i = nenergy2_max + 1, nenergy_max
          energy2 = egrid%val(i)
          denergy2 = egrid%val(i-1) - egrid%val(i)
          ! j: index of T1
          do j = nenergy1_max + 1, nenergy_max
             energy1 = egrid%val(j)
             denergy1 = egrid%val(j-1) - egrid%val(j)
             if (energy1 > 2.0*energy2 + energy_ionize) then
                self%degradation_gen(igen, i) = self%degradation_gen(igen, i) &
                     + self%number_electrons(io) * self%number_density &
                     * self%degradation_gen(igen-1, j) &
                     * denergy1 * denergy2 &
                     * bb/(energy1 + energy_ionize + energy_kinetic) * &
                     ( 1.0/(energy2 + energy_ionize)**2 &
                     + 4.0/3.0 * energy_kinetic/(energy2 + energy_ionize)**3 &
                     + 4.0/3.0 * energy_kinetic/(energy1 - energy2)**3 &
                     + 1.0/(energy1 - energy2)**2 &
                     )
             else
                exit
             end if
          end do
       end do
    end do

    do i = 2, nenergy_max
       self%degradation_gen(igen, i) = self%degradation_gen(igen, i-1) &
            + self%degradation_gen(igen,i)
    end do
    where (mediamix%stop_power_mixture /= 0.)
       self%degradation_gen(igen, :) = self%degradation_gen(igen, :)/mediamix%stop_power_mixture(:)
    end where

    return

  end subroutine calculate_degradation

  subroutine calculate_yield(self, egrid, mediamix)
    use class_grid, only: grid
    use class_mixture, only: mixture
    class(medium) :: self
    class(grid), intent(in) :: egrid
    class(mixture), intent(in) :: mediamix
    integer :: io, ie
    integer :: ne1, ne2, ne3
    integer :: ne1_max = 1, ne2_max = 1, ne3_max = 1
    real :: t_over_s(egrid%number)

    do io = 1, self%number

       if (self%energy_ionize(io) < egrid%val_max) then
          ne1 = egrid%grid_number(self%energy_ionize(io))
          if (ne1 > ne1_max) ne1_max = ne1

          ! energy >=  energy_ionize
          do ie = 1, ne1
             self%total_cross_section_ionize_orbital(io, ie) = &
                  integrate_sigma(self, "direct", io, &
                  egrid%val(ie), self%energy_ionize(io), egrid%val(ie))

             self%total_cross_section_singlet_orbital(io, ie) = &
                  integrate_sigma(self, "direct", io, &
                  egrid%val(ie), self%energy_singlet(io), self%energy_ionize(io)) &
                  + 0.5*integrate_sigma(self, "exchange", io, &
                  egrid%val(ie), self%energy_singlet(io), self%energy_ionize(io))

             self%total_cross_section_triplet_orbital(io, ie) = &
                  0.5*integrate_sigma(self, "exchange", io, &
                  egrid%val(ie), self%energy_triplet(io), self%energy_ionize(io) )
          end do
       else
          ne1 = 0
       end if

       if (self%energy_singlet(io) < egrid%val_max) then
          ne2 = egrid%grid_number(self%energy_singlet(io))
          if (ne2 > ne2_max) ne2_max = ne2

          ! energy_ionize > energy >=  energy_singlet
          do ie = ne1+1, ne2
             self%total_cross_section_singlet_orbital(io, ie) = &
                  integrate_sigma(self, "direct", io, &
                  egrid%val(ie), self%energy_singlet(io), egrid%val(ie)) &
                  + 0.5*integrate_sigma(self, "exchange", io, &
                  egrid%val(ie), self%energy_singlet(io), egrid%val(ie))

             self%total_cross_section_triplet_orbital(io, ie) = &
                  0.5*integrate_sigma(self, "exchange", io, &
                  egrid%val(ie), self%energy_triplet(io), egrid%val(ie))
          end do
       else
          ne2 = 0
       end if

       if (self%energy_triplet(io) < egrid%val_max) then
          ne3 = egrid%grid_number(self%energy_triplet(io))
          if (ne3 > ne3_max) ne3_max = ne3

          ! Singlet_energy > energy >=  energy_triplet
          do ie = ne2+1, ne3
             self%total_cross_section_triplet_orbital(io, ie) = &
                  0.5*integrate_sigma(self, "exchange", io, &
                  egrid%val(ie), self%energy_triplet(io), egrid%val(ie))
          end do
       end if

    end do

    self%total_cross_section_ionize(:) = sum(self%total_cross_section_ionize_orbital(:,:), dim=1)
    self%total_cross_section_singlet(:) = sum(self%total_cross_section_singlet_orbital(:,:), dim=1)
    self%total_cross_section_triplet(:) = sum(self%total_cross_section_triplet_orbital(:,:), dim=1)

    do io = 1, self%number
       self%platzman_ionize_orbital(io, :) = egrid%val * &
            mediamix%degradation_mixture * &
            self%total_cross_section_ionize_orbital(io, :)
       self%platzman_singlet_orbital(io, :) = egrid%val * &
            mediamix%degradation_mixture * &
            self%total_cross_section_singlet_orbital(io, :)
       self%platzman_triplet_orbital(io, :) = egrid%val * &
            mediamix%degradation_mixture * &
            self%total_cross_section_triplet_orbital(io, :)

       ! integrate over energy
       ! minus sign is because the energy grid is in the descending order
       self%yield_ionize(io) = (-log(egrid%div)) * &
            ( sum(self%platzman_ionize_orbital(io, :)) &
            - 0.5*self%platzman_ionize_orbital(io, 1) )

       if (self%energy_singlet(io) < minval(self%energy_ionize)) then
          self%yield_singlet(io) = (-log(egrid%div)) * &
               ( sum(self%platzman_singlet_orbital(io, :)) &
               - 0.5*self%platzman_singlet_orbital(io, 1) )
       else
          self%yield_ionize(io) = self%yield_ionize(io) + (-log(egrid%div)) * &
               ( sum(self%platzman_singlet_orbital(io, :)) &
               - 0.5*self%platzman_singlet_orbital(io, 1) )
          self%yield_singlet(io) = 0.
       end if

       if (self%energy_triplet(io) < minval(self%energy_ionize)) then
          self%yield_triplet(io) = (-log(egrid%div)) * &
               ( sum(self%platzman_triplet_orbital(io, :)) &
               - 0.5*self%platzman_triplet_orbital(io, 1) )
       else
          self%yield_ionize(io) = self%yield_ionize(io) + (-log(egrid%div)) * &
               ( sum(self%platzman_triplet_orbital(io, :)) &
               - 0.5*self%platzman_triplet_orbital(io, 1) )
          self%yield_triplet(io) = 0.
       end if
    end do

    self%total_cross_section_total = self%total_cross_section_ionize + &
         self%total_cross_section_singlet + self%total_cross_section_triplet

    where (self%stop_power /= 0.)
       t_over_s = egrid%val(:) / self%stop_power(:) * (-log(egrid%div))
    end where

    do ie = ne3_max, 1, -1
       self%range_triplet(ie) = self%range_triplet(ie+1) + t_over_s(ie)
    end do
    do ie = ne2_max, 1, -1
       self%range_singlet(ie) = self%range_singlet(ie+1) + t_over_s(ie)
    end do
    do ie = ne1_max, 1, -1
       self%range_ionize(ie) = self%range_ionize(ie+1) + t_over_s(ie)
    end do

    self%yield_ionize = self%number_electrons * self%number_density * self%yield_ionize
    self%yield_singlet = self%number_electrons * self%number_density * self%yield_singlet
    self%yield_triplet = self%number_electrons * self%number_density * self%yield_triplet

    self%gvalue_ionize = 100. * self%yield_ionize/egrid%val_max
    self%gvalue_singlet = 100. * self%yield_singlet/egrid%val_max
    self%gvalue_triplet = 100. * self%yield_triplet/egrid%val_max

  end subroutine calculate_yield
  
  real function integrate_sigma(self, type, io, energy, range_min, range_max)
    class(medium) :: self
    character (len=*) :: type
    integer, intent(in) :: io
    real, intent(in) :: energy
    real, intent(in) :: range_min, range_max
    real :: val
    real :: energy_ionize, energy_kinetic

    abstract interface
       real function indefinite_integral(x, energy, energy_ionize, energy_kinetic)
         real, intent(in) :: x
         real, intent(in) :: energy, energy_ionize, energy_kinetic
       end function indefinite_integral
    end interface

    procedure (indefinite_integral), pointer :: indefinite_sigma => null()

    select case(type)
    case ("total")
       indefinite_sigma => indefinite_sigma_Etot
    case ("direct")
       indefinite_sigma => indefinite_sigma_Edir
    case ("exchange")
       indefinite_sigma => indefinite_sigma_Eexc
    end select

    energy_ionize = self%energy_ionize(io)
    energy_kinetic = self%energy_kinetic(io)

    val = indefinite_sigma(range_max, energy, energy_ionize, energy_kinetic) &
         - indefinite_sigma(range_min, energy, energy_ionize, energy_kinetic)
    integrate_sigma = val

    return

  contains

    real function indefinite_sigma_Etot(x, energy, energy_ionize, energy_kinetic)
      real, intent(in) :: x
      real, intent(in) :: energy, energy_ionize, energy_kinetic
      real :: val
      val = indefinite_sigma_Edir(x, energy, energy_ionize, energy_kinetic) + &
           indefinite_sigma_Eexc(x, energy, energy_ionize, energy_kinetic)
      indefinite_sigma_Etot = val
      return
    end function indefinite_sigma_Etot

    real function indefinite_sigma_Edir(x, energy, energy_ionize, energy_kinetic)
      use mod_constants, only: bb
      real, intent(in) :: x
      real, intent(in) :: energy, energy_ionize, energy_kinetic
      real :: val

      val = bb/(energy + energy_ionize + energy_kinetic) * &
           ( -1./x - 2./3.*energy_kinetic/x**2 )

      indefinite_sigma_Edir = val
      return
    end function indefinite_sigma_Edir

    real function indefinite_sigma_Eexc(x, energy, energy_ionize, energy_kinetic)
      use mod_constants, only: bb
      real, intent(in) :: x
      real, intent(in) :: energy, energy_ionize, energy_kinetic
      real :: val

      val = bb/(energy + energy_ionize + energy_kinetic) * &
           ( 1./(energy + energy_ionize - x) &
           + 2./3.*energy_kinetic/(energy + energy_ionize -x)**2 )

      indefinite_sigma_Eexc = val
      return
    end function indefinite_sigma_Eexc

  end function integrate_sigma

end module class_medium

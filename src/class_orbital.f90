module class_orbital
  implicit none
  private

  public :: orbital

  integer, parameter :: ngen_degradation = 6
  
  type orbital
     private
     ! number of orbital
     integer, public :: number = 0
     ! ionization Ii, kinetic Ei, singlet Esi, triplet Eti energies (per orbital)
     real, pointer, public :: energy_ionize(:) => null(), energy_kinetic(:) => null()
     real, pointer, public :: energy_singlet(:) => null(), energy_triplet(:) => null()
     ! number of electrons (per orbital)
     integer, pointer, public :: number_electrons(:) => null()
     ! number density
     real, public :: number_density

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

     ! platzman energy*y(energy)*Q(energy) (per orbital)
     real, pointer, public :: platzman_ionize_orbital(:,:) => null()
     real, pointer, public :: platzman_singlet_orbital(:,:) => null()
     real, pointer, public :: platzman_triplet_orbital(:,:) => null()

     ! mean free path
     real, pointer, public :: mean_free_path(:) => null()

   contains
     final :: destroy
     procedure, public :: init_orbital_vars
     procedure :: calculate_stopping_power
     procedure :: calculate_degradation
     procedure :: calculate_yield
     procedure :: print_results
  end type orbital

  ! declare constructor
  interface orbital
     module procedure init_orbital
  end interface orbital

contains

  type(orbital) function init_orbital(norbital,file_orbital)
    integer, intent(in) :: norbital
    character (len=*), intent(in) :: file_orbital
    
    init_orbital%number = norbital

    if (.not.associated(init_orbital%energy_ionize)) allocate(init_orbital%energy_ionize(norbital))
    if (.not.associated(init_orbital%energy_kinetic)) allocate(init_orbital%energy_kinetic(norbital))
    if (.not.associated(init_orbital%energy_singlet)) allocate(init_orbital%energy_singlet(norbital))
    if (.not.associated(init_orbital%energy_triplet)) allocate(init_orbital%energy_triplet(norbital))
    if (.not.associated(init_orbital%number_electrons)) allocate(init_orbital%number_electrons(norbital))

    if (.not.associated(init_orbital%yield_ionize)) allocate(init_orbital%yield_ionize(norbital))
    if (.not.associated(init_orbital%gvalue_ionize)) allocate(init_orbital%gvalue_ionize(norbital))
    if (.not.associated(init_orbital%yield_singlet)) allocate(init_orbital%yield_singlet(norbital))
    if (.not.associated(init_orbital%gvalue_singlet)) allocate(init_orbital%gvalue_singlet(norbital))
    if (.not.associated(init_orbital%yield_triplet)) allocate(init_orbital%yield_triplet(norbital))
    if (.not.associated(init_orbital%gvalue_triplet)) allocate(init_orbital%gvalue_triplet(norbital))

    call read_params(init_orbital, norbital, file_orbital)
  end function init_orbital

  subroutine destroy(self)
    type(orbital), intent(in out) :: self
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

  subroutine read_params(self, norbital, file)
    use mod_file_utils, only: get_unused_unit
    class(orbital) :: self
    integer, intent(in) :: norbital
    character (len=*), intent(in) :: file

    integer :: unit
    
    real :: energy_ionize(norbital)
    real :: energy_kinetic(norbital)
    real :: energy_singlet(norbital)
    real :: energy_triplet(norbital)
    integer :: number_electrons(norbital)
    
    namelist /orbital_params/ energy_ionize, energy_kinetic, &
         energy_singlet, energy_triplet, number_electrons
    
    call get_unused_unit(unit)
    open(unit,file=trim(file))
    read(unit,orbital_params)
    close(unit)

    self%energy_ionize(1:norbital) = energy_ionize(1:norbital)
    self%energy_kinetic(1:norbital) = energy_kinetic(1:norbital)
    self%energy_singlet(1:norbital) = energy_singlet(1:norbital)
    self%energy_triplet(1:norbital) = energy_triplet(1:norbital)
    self%number_electrons(1:norbital) = number_electrons(1:norbital)
  end subroutine read_params

  subroutine init_orbital_vars(self, ngrid)
    class(orbital) :: self
    integer, intent(in) :: ngrid
    if (.not.associated(self%stop_power_orbital)) allocate(self%stop_power_orbital(self%number, ngrid))
    self%stop_power_orbital = 0.0
    if (.not.associated(self%stop_power)) allocate(self%stop_power(ngrid))
!!!    self%stop_power = 0.0 ! not necessary

    if (.not.associated(self%degradation_gen)) allocate(self%degradation_gen(ngen_degradation, ngrid))
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

    if (.not.associated(self%mean_free_path)) allocate(self%mean_free_path(ngrid))
    self%mean_free_path = 0.0

  end subroutine init_orbital_vars

  subroutine finish_orbital_vars(self)
    class(orbital) :: self
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

    if (associated(self%platzman_ionize_orbital)) nullify(self%platzman_ionize_orbital)
    if (associated(self%platzman_singlet_orbital)) nullify(self%platzman_singlet_orbital)
    if (associated(self%platzman_triplet_orbital)) nullify(self%platzman_triplet_orbital)

    if (associated(self%mean_free_path)) nullify(self%mean_free_path)
  end subroutine finish_orbital_vars

  subroutine calculate_stopping_power(self)
    use mod_grid, only: egrid
    class(orbital) :: self
    integer :: io, ie
    integer :: ne1, ne2, ne3

    do io = 1, self%number ! orbital index

       ne1 = egrid%grid_number(self%energy_ionize(io))
       ne2 = egrid%grid_number(self%energy_singlet(io))
       ne3 = egrid%grid_number(self%energy_triplet(io))

       ! energy >=  energy_ionize
       do ie = 1, ne1
          self%stop_power_orbital(io, ie) = &
               integrate_E_sigma(self, "total", io, &
               egrid%val(ie), self%energy_singlet(io), 0.5*(egrid%val(ie) + self%energy_ionize(io))) &
               + 0.5*integrate_E_sigma(self, "exchange", io, &
               egrid%val(ie), self%energy_triplet(io), self%energy_singlet(io))
       end do

       ! energy_ionize > energy >=  energy_singlet
       do ie = ne1+1, ne2
          self%stop_power_orbital(io, ie) = &
               integrate_E_sigma(self, "total", io, &
               egrid%val(ie), self%energy_singlet(io), egrid%val(ie)) + &
               0.5*integrate_E_sigma(self, "exchange", io, &
               egrid%val(ie), self%energy_triplet(io), self%energy_singlet(io))
       end do

       ! energy_singlet > energy >=  energy_triplet
       do ie = ne2+1, ne3
          self%stop_power_orbital(io, ie) = 0.5*integrate_E_sigma(self, "exchange", io, &
               egrid%val(ie), self%energy_triplet(io), egrid%val(ie))
       end do
    end do 

    self%stop_power(:) = sum(self%stop_power_orbital(:,:), dim=1)
    
  end subroutine calculate_stopping_power

  real function integrate_E_sigma(self, type, io, energy, range_min, range_max)
    class(orbital) :: self
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
      use mod_constants, only: bb => bb_compat
      real, intent(in) :: x
      real, intent(in) :: energy, energy_ionize, energy_kinetic
      real :: val

      val = bb/(energy + energy_ionize + energy_kinetic) * &
           ( log(x) - 4./3.*energy_kinetic/x )

      indefinite_E_sigma_Edir = val
      return
    end function indefinite_E_sigma_Edir

    real function indefinite_E_sigma_Eexc(x, energy, energy_ionize, energy_kinetic)
      use mod_constants, only: bb => bb_compat
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

  recursive subroutine calculate_degradation(self, ngen_in)
    use mod_constants, only: bb => bb_compat
    use mod_grid, only: egrid
    class(orbital) :: self
    integer, intent(in), optional :: ngen_in
    integer :: ngen
    real :: energy1, denergy1, energy2, denergy2
    real :: energy_ionize, energy_kinetic
    integer :: io
    integer :: i, j
    real :: energy1_max, energy2_max
    integer :: nenergy_max
    integer :: nenergy1_max, nenergy2_max

    if (.not.present(ngen_in)) then
       ngen = ngen_degradation
    else
       ngen = ngen_in
    end if
    
    nenergy_max = egrid%grid_number(minval(self%energy_triplet))

    if (ngen == 1) then

       where(self%stop_power /= 0.)
          self%degradation_gen(1, :) = 1./self%stop_power(:)
       end where

    else

       call calculate_degradation(self, ngen-1)
       
       do io = 1, self%number
          energy_ionize = self%energy_ionize(io)
          energy_kinetic = self%energy_kinetic(io)

          energy1_max = (egrid%val_max - (2**(ngen-2)-1)*energy_ionize)/(2**(ngen-2))
          if (ngen == 2) energy1_max = egrid%val_max
          energy2_max = (egrid%val_max - (2**(ngen-1)-1)*energy_ionize)/(2**(ngen-1))
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
                   self%degradation_gen(ngen, i) = self%degradation_gen(ngen, i) &
                        + self%number_electrons(io) * self%number_density &
                        * self%degradation_gen(ngen-1, j) &
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

       do i = nenergy1_max + 1, nenergy_max
          self%degradation_gen(ngen, i) = self%degradation_gen(ngen, i-1) &
               + self%degradation_gen(ngen,i)
       end do
       where (self%stop_power /= 0.)
          self%degradation_gen(ngen, :) = self%degradation_gen(ngen, :)/self%stop_power(:)
       end where
    end if

    return

  end subroutine calculate_degradation

  subroutine calculate_yield(self)
    use mod_grid, only: egrid
    class(orbital) :: self
    integer :: io, ie
    integer :: ne1, ne2, ne3

    do io = 1, self%number

       ne1 = egrid%grid_number(self%energy_ionize(io))
       ne2 = egrid%grid_number(self%energy_singlet(io))
       ne3 = egrid%grid_number(self%energy_triplet(io))

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

       ! singlet_energy > energy >=  energy_triplet
       do ie = ne2+1, ne3
          self%total_cross_section_triplet_orbital(io, ie) = &
               0.5*integrate_sigma(self, "exchange", io, &
               egrid%val(ie), self%energy_triplet(io), egrid%val(ie))
       end do

    end do

    self%total_cross_section_ionize(:) = sum(self%total_cross_section_ionize_orbital(:,:), dim=1)
    self%total_cross_section_singlet(:) = sum(self%total_cross_section_singlet_orbital(:,:), dim=1)
    self%total_cross_section_triplet(:) = sum(self%total_cross_section_triplet_orbital(:,:), dim=1)

    do io = 1, self%number
       self%platzman_ionize_orbital(io, :) = egrid%val * &
            sum(self%degradation_gen, dim=1) * &
            self%total_cross_section_ionize_orbital(io, :)
       self%platzman_singlet_orbital(io, :) = egrid%val * &
            sum(self%degradation_gen, dim=1) * &
            self%total_cross_section_singlet_orbital(io, :)
       self%platzman_triplet_orbital(io, :) = egrid%val * &
            sum(self%degradation_gen, dim=1) * &
            self%total_cross_section_triplet_orbital(io, :)

       self%yield_ionize(io) = sum(self%platzman_ionize_orbital(io, :)) &
            - 0.5*self%platzman_ionize_orbital(io, 1)
       self%yield_singlet(io) = sum(self%platzman_singlet_orbital(io, :)) &
            - 0.5*self%platzman_singlet_orbital(io, 1)
       self%yield_triplet(io) = sum(self%platzman_triplet_orbital(io, :)) &
            - 0.5*self%platzman_triplet_orbital(io, 1)
    end do

    self%mean_free_path = self%number_density*( &
         self%total_cross_section_ionize + &
         self%total_cross_section_singlet + &
         self%total_cross_section_triplet )
    where (self%mean_free_path /= 0.)
       self%mean_free_path = 1./self%mean_free_path
    end where

    ! integrate over energy
    ! minus sign is because the energy grid is in the descending order
    self%yield_ionize = self%number_electrons * self%number_density * self%yield_ionize * (-log(egrid%div))
    self%yield_singlet = self%number_electrons * self%number_density * self%yield_singlet * (-log(egrid%div))
    self%yield_triplet = self%number_electrons * self%number_density * self%yield_triplet * (-log(egrid%div))

    self%gvalue_ionize = 100. * self%yield_ionize/egrid%val_max
    self%gvalue_singlet = 100. * self%yield_singlet/egrid%val_max
    self%gvalue_triplet = 100. * self%yield_triplet/egrid%val_max

  end subroutine calculate_yield
  
  real function integrate_sigma(self, type, io, energy, range_min, range_max)
    class(orbital) :: self
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
      use mod_constants, only: bb => bb_compat
      real, intent(in) :: x
      real, intent(in) :: energy, energy_ionize, energy_kinetic
      real :: val

      val = bb/(energy + energy_ionize + energy_kinetic) * &
           ( -1./x - 2./3.*energy_kinetic/x**2 )

      indefinite_sigma_Edir = val
      return
    end function indefinite_sigma_Edir

    real function indefinite_sigma_Eexc(x, energy, energy_ionize, energy_kinetic)
      use mod_constants, only: bb => bb_compat
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

  subroutine print_results(self)
    use mod_file_utils, only: get_unused_unit, unit_stdout
    use mod_grid, only: egrid
    class(orbital) :: self
    integer :: io, ie
    integer :: unit
    character (len=100) :: file

    file = 'results.dat'

    call get_unused_unit(unit)
    open(unit, file=trim(file))

    write(unit,'("#")')
    write(unit,'("#",a)') ' Yield:'
    write(unit,'("#",4(a20))') 'Orbital', 'ionize', 'singlet', 'triplet'
    do io = 1, self%number
       write(unit,'("#",15x,i5,3(8x,f12.4))') &
            io, self%yield_ionize(io), self%yield_singlet(io), self%yield_triplet(io)
    end do
    if (self%number > 1) then
       write(unit,'("#",a20,3(8x,f12.4))') &
            'sum', sum(self%yield_ionize), sum(self%yield_singlet), sum(self%yield_triplet)
    end if

    write(unit,'("#")')
    write(unit,'("#",a)') ' G-value:'
    write(unit,'("#",4(a20))') 'Orbital', 'ionize', 'singlet', 'triplet'
    do io = 1, self%number
       write(unit,'("#",15x,i5,3(8x,f12.4))') &
            io, self%gvalue_ionize(io), self%gvalue_singlet(io), self%gvalue_triplet(io)
    end do
    if (self%number > 1) then
       write(unit,'("#",a20,3(8x,f12.4))') &
            'sum', sum(self%gvalue_ionize), sum(self%gvalue_singlet), sum(self%gvalue_triplet)
    end if
    write(unit,'("#")')

    write(unit,'("#",10(1x,a20))') 'Energy', 'Stopping Power', 'Degradation (sum)', &
         'Total Cross Sec. i', 'Total Cross Sec. s', 'Total Cross Sec. t', &
         'Platzman i', 'Platzman s', 'Platzman t', 'Mean Free Path'

    do ie = 1, egrid%number
       write(unit,'(1x,10(1x,e20.12))') egrid%val(ie), self%stop_power(ie), &
            sum(self%degradation_gen(:, ie)), &
            self%total_cross_section_ionize(ie), &
            self%total_cross_section_singlet(ie), &
            self%total_cross_section_triplet(ie), &
            sum(self%platzman_ionize_orbital(:, ie)), &
            sum(self%platzman_singlet_orbital(:, ie)), &
            sum(self%platzman_triplet_orbital(:, ie)), &
            self%mean_free_path(ie)
    end do
    close(unit)

    unit = unit_stdout

    write(unit,*)
    write(unit,'(1x,a)') ' Yield:'
    write(unit,'(1x,4(a20))') 'Orbital', 'ionize', 'singlet', 'triplet'
    do io = 1, self%number
       write(unit,'(1x,15x,i5,3(8x,f12.4))') &
            io, self%yield_ionize(io), self%yield_singlet(io), self%yield_triplet(io)
    end do
    if (self%number > 1) then
       write(unit,'(1x,a20,3(8x,f12.4))') &
            'sum', sum(self%yield_ionize), sum(self%yield_singlet), sum(self%yield_triplet)
    end if

    write(unit,*)
    write(unit,'(1x,a)') ' G-value:'
    write(unit,'(1x,4(a20))') 'Orbital', 'ionize', 'singlet', 'triplet'
    do io = 1, self%number
       write(unit,'(1x,15x,i5,3(8x,f12.4))') &
            io, self%gvalue_ionize(io), self%gvalue_singlet(io), self%gvalue_triplet(io)
    end do
    if (self%number > 1) then
       write(unit,'(1x,a20,3(8x,f12.4))') &
            'sum', sum(self%gvalue_ionize), sum(self%gvalue_singlet), sum(self%gvalue_triplet)
    end if

    call print_results_degradation(self)

  end subroutine print_results

  subroutine print_results_degradation(self)
    use mod_file_utils, only: get_unused_unit, unit_stdout
    use mod_grid, only: egrid
    class(orbital) :: self
    integer :: ie, igen
    integer :: unit
    character (len=100) :: file

    file = 'results_degradation.dat'

    call get_unused_unit(unit)
    open(unit, file=trim(file))

    write(unit,'("#",2(1x,a20))', advance='no') 'Energy', 'Degradation (sum)'
    do igen = 1, ngen_degradation
       write(unit,'(1x, a16, 1x, i3)', advance='no') 'Degradation ', igen
    end do
    write(unit,*)

    do ie = 1, egrid%number
       write(unit,'(1x,2(1x,e20.12))', advance='no') egrid%val(ie), sum(self%degradation_gen(:, ie))
       do igen = 1, ngen_degradation
          write(unit,'(1x,e20.12)', advance='no') self%degradation_gen(igen, ie)
       end do
       write(unit,*)
    end do
    close(unit)

  end subroutine print_results_degradation

end module class_orbital

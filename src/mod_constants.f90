module mod_constants
  implicit none
  private
  public :: bb
  public :: bb_compat ! just for backward compatibility

  ! Pi to quad precision
  double precision, parameter :: dpi = &
       3.14159265358979323846264338327950288419716939938

  ! exact physical constants [SI unit]
  real, parameter :: light_speed_si = 2.99792458e8
  real, parameter :: elementary_charge_si = 1.602176634e-19

  ! elementary change in esu
  real, parameter :: elementary_charge_esu = elementary_charge_si * 10. * light_speed_si

  ! pi*e^4
  real, parameter :: bb = real(dpi) * (elementary_charge_esu**2 * 1.e-7/elementary_charge_si)**2
  real, parameter :: bb_compat = 6.51348e-14
  
end module mod_constants

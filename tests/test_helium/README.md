# test_helium

This test is to check calculation of degradation.

For simplicity, we use the parameters of Helium gas and consider
4 generations of electrons. (Further generations are obtained using
the same algorithm.)

The numerically obtained data are compared with those from 
WolfranEngine (WE). The WE data are also numerical using NIntegrate.
Since it is a heavy operation to calculate NIntegrate of NIntegrate-d
functions, we use Interpolation for the degradation of the previous
generation (e.g., for y3, we use interpolated y2 obained from NIntegrate
of y1).

To compare to data, we also use interpolation just in case if the dimensions
of two data are different. Only data in the valid energy range are compared,
otherwise the interpolation fails.
(This treatment leads to unnecessary complication and can be avoided.)

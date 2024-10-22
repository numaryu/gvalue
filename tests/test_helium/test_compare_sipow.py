import pytest
import numpy as np

eps = 1.e-6

def test_stoppow():
    test_helium = np.loadtxt("test_helium.dat")
    egrid_result = test_helium[:,0]
    stoppow_result = test_helium[:,1]
    print(stoppow_result[:10])
    
    stoppow = np.loadtxt("stoppow.dat")
    egrid_soln = stoppow[:,1]
    stoppow_soln = stoppow[:,3]
    print(stoppow_soln[:10])

    rms = np.sqrt(np.average( (egrid_result - egrid_soln)**2) )
    rms += np.sqrt(np.average( (stoppow_result - stoppow_soln)**2) )

    print(rms)
    assert rms < eps

def test_degradation():
    test_helium = np.loadtxt("test_helium_degradation.dat")
    egrid_result = test_helium[:,0]
    ytot_result = test_helium[:,1]
    y1_result = test_helium[:,2]
    y2_result = test_helium[:,3]
    y3_result = test_helium[:,4]
    y4_result = test_helium[:,5]
    y5_result = test_helium[:,6]
    y6_result = test_helium[:,7]
    print(y1_result[:10])
    print(y2_result[100:110])

    degradation = np.loadtxt("degradation320.dat")
    egrid_soln = degradation[:,1]
    y1_soln = degradation[:,2]
    y2_soln = degradation[:,3]

    y1_spl = interpolate.interp1d(egrid_soln, y1_soln, kind="cubic")
    y1_spl = interpolate.interp1d(egrid_soln, y1_soln, kind="cubic")
    print(y1_soln[:10])
    print(y2_soln[100:110])

    rms1 = np.sqrt(np.average( (y1_result - y1_soln)**2) )
    rms2 = np.sqrt(np.average( (y2_result - y2_soln)**2) )

    print(rms1)
    print(rms2)
    assert rms1 < eps
    assert rms2 < eps

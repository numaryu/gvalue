import pytest
import numpy as np

eps = 1.e-6

def test_stoppowstoppow():
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

import pytest
import numpy as np
import scipy.interpolate as interpolate

eps_spectra = 1.e-6
eps_gvalue = 1.e-4

def test_spectra():
    test_neon = np.loadtxt("test_neon.dat")
    spectra_result = np.zeros( (test_neon.shape[0], 8) )
    spectra_result[:, 0] = np.flipud(test_neon[:, 0])
    spectra_result[:, 1:7] = np.flipud(test_neon[:, 3:9])
    spectra_result[:, 7] = np.flipud(test_neon[:, 1])

    # data from WolframEngine
    spectra_data = np.loadtxt("spectra_neon.dat")

    rms = np.sqrt(np.mean((spectra_result-spectra_data)**2))

    # print(rms)
    assert rms < eps_spectra

def test_gvalue():
    gvalue_result = np.loadtxt("test_neon_gvalue.dat")
    gvalue_data = np.loadtxt("gvalue_neon.dat")

    rms = np.sqrt(np.mean((gvalue_result-gvalue_data)**2))

    # print(rms)
    assert rms < eps_gvalue

import pytest
import numpy as np

eps_spectra = 1.e-6
eps_gvalue = 1.e-4

def test_spectra():
    spectra_result = np.loadtxt("test_mixture.dat")
    spectra_ref = np.loadtxt("test_mixture_ref.dat")

    rms = np.sqrt(np.mean((spectra_result-spectra_ref)**2))

    # print(rms)
    assert rms < eps_spectra

def test_gvalue():
    gvalue_result = np.loadtxt("test_mixture_gvalue.dat")
    gvalue_ref = np.loadtxt("test_mixture_gvalue_ref.dat")

    rms = np.sqrt(np.mean((gvalue_result-gvalue_ref)**2))

    # print(rms)
    assert rms < eps_gvalue

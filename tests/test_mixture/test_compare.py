import pytest
import numpy as np

eps_spectra = 1.e-6
eps_gvalue = 1.e-4

def test_spectra():
    for i in range(6):
        spectra_result = np.loadtxt(f"test_mixture{i:d}.dat")
        spectra_ref = np.loadtxt(f"test_mixture{i:d}_ref.dat")

        rms = np.sqrt(np.mean((spectra_result-spectra_ref)**2))

        # print(rms)
        assert rms < eps_spectra

def test_gvalue():
    for i in range(6):
        gvalue_result = np.loadtxt(f"test_mixture{i:d}_gvalue.dat")
        gvalue_ref = np.loadtxt(f"test_mixture{i:d}_gvalue_ref.dat")

        rms = np.sqrt(np.mean((gvalue_result-gvalue_ref)**2))

        # print(rms)
        assert rms < eps_gvalue

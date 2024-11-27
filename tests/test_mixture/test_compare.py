import pytest
import numpy as np

eps_gvalue = 1.e-4

def test_gvalue():
    gvalue_result = np.loadtxt("test_mixture_gvalue.dat")
    gvalue_ref = np.loadtxt("test_mixture_gvalue_ref.dat")

    rms = np.sqrt(np.mean((gvalue_result-gvalue_ref)**2))

    # print(rms)
    assert rms < eps_gvalue

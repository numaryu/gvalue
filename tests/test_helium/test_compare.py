import pytest
import numpy as np
import scipy.interpolate as interpolate

eps_stoppow = 1.e-6
eps_degradation = 1.e-3

def test_stoppow():
    test_helium = np.loadtxt("test_helium.dat")
    egrid_result = np.flipud(test_helium[:, 0])
    stoppow_result = np.flipud(test_helium[:, 1])
    stoppow_result_interp = interpolate.CubicSpline(egrid_result, stoppow_result)

    # minimum energy value that gives finite stoppow
    energy_min_result = egrid_result[np.min(np.where(stoppow_result != 0))]

    # data from WolframEngine
    stoppow = np.loadtxt("stoppow_helium.dat")
    egrid_data = stoppow[:, 0]
    stoppow_data = stoppow[:, 1]
    stoppow_data_interp = interpolate.CubicSpline(egrid_data, stoppow_data)

    # minimum energy value that gives finite stoppow
    energy_min_data = egrid_data[np.min(np.where(stoppow_data != 0))]

    # minimum valid energy
    energy_min = max(energy_min_result, energy_min_data)

    if egrid_result.size <= egrid_data.size:
        egrid_check = egrid_result
    else:
        egrid_check = egrid_data

    # mask if egrid is less than valid energy
    madiff = np.ma.array((stoppow_result_interp(egrid_check) - \
                          stoppow_data_interp(egrid_check))**2,
                         mask = egrid_check < energy_min)
    rms = np.sqrt(np.ma.mean(madiff))

    # print(rms)
    assert rms < eps_stoppow

def test_degradation():
    ngen = 4

    test_helium = np.loadtxt("test_helium_degradation.dat")
    egrid_result = np.flipud(test_helium[:, 0])
    deggen_result = np.flipud(test_helium[:, 2:2+ngen])

    # minimum energy value that gives finite stoppow
    energy_min_result = np.zeros((ngen))
    energy_max_result = np.zeros((ngen))
    for i in range(ngen):
        energy_min_result[i] = egrid_result[np.min(np.where(deggen_result[:, i] != 0))]
        energy_max_result[i] = egrid_result[np.max(np.where(deggen_result[:, i] != 0))]

    # list of functions
    deggen_result_interp = []
    for i in range(ngen):
        deggen_result_interp.append(interpolate.CubicSpline(egrid_result, deggen_result[:, i]))

    # data from WolframEngine
    degradation = np.loadtxt("degradation_helium.dat")
    egrid_data = degradation[:, 0]
    deggen_data = degradation[:, 1:1+ngen]

    # minimum energy value that gives finite stoppow
    energy_min_data = np.zeros((ngen))
    energy_max_data = np.zeros((ngen))
    for i in range(ngen):
        energy_min_data[i] = egrid_data[np.min(np.where(deggen_data[:, i] != 0))]
        energy_max_data[i] = egrid_data[np.max(np.where(deggen_data[:, i] != 0))]

    # list of functions
    deggen_data_interp = []
    for i in range(ngen):
        deggen_data_interp.append(interpolate.CubicSpline(egrid_data, deggen_data[:, i]))

    # minimum valid energy
    energy_min = np.max(np.stack([energy_min_result, energy_min_data]), axis = 0)
    energy_max = np.min(np.stack([energy_max_result, energy_max_data]), axis = 0)

    if egrid_result.size <= egrid_data.size:
        egrid_check = egrid_result
    else:
        egrid_check = egrid_data

    rms = []
    for i in range(ngen):
        # mask if egrid is not in valid energy range
        madiff = np.ma.array((deggen_result_interp[i](egrid_check) - \
                              deggen_data_interp[i](egrid_check))**2,
                             mask = np.any([egrid_check < energy_min[i],
                                            egrid_check > energy_max[i]], axis=0))
        rms.append(np.sqrt(np.ma.mean(madiff)))

    # print(rms)
    for i in range(ngen):
        if i == 0:
            # for this case, degradation is inverse of stoppow
            assert rms[i] < eps_stoppow
        else:
            assert rms[i] < eps_degradation

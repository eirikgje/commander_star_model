"""Unzips spectrum files on the fly and calculates the response for the given spectrum per DIRBE band

NOTE THAT THE NORMALIZATION OF THE SPECTRUM IS NOT DONE HERE! No guarantees are given that the final amplitude should match anything.

"""

import gzip
import itertools
import numpy as np
import dirbe_utils
import astropy.units
import astropy.constants
import h5py
from scipy.interpolate import CubicSpline

splines = []

#for bandpass_num in range(11):
#    for wavelengths, response in dirbe_utils.get_bandpass(bandpass_num):
#        splines.append((CubicSpline(wavelengths, response), min(wavelengths), max(wavelengths)))

#outfile ='/home/eirik/data/commander_star_model/dirbe_integrated_bands.dat'
outfile ='/home/eirik/data/commander_star_model/commander_star_model.h5'
logg_vals = np.array([3.0, 3.5, 4.0, 4.5, 5.0])
teff_vals = np.array(list([i for i in range(3000, 4000, 50)]) +
                     list([i for i in range(4000, 12000, 100)]))
mh_vals = np.array([-2.0, -1.5, -0.5, 0.0, 0.5, 1.0])

def find_minmax(minval, maxval, target_array):
    minarg = np.searchsorted(target_array, minval)
    maxarg = np.searchsorted(target_array, maxval)
    return minarg, maxarg

param_vals = []
reported_vals = []

with h5py.File(outfile, 'w') as of:
    param_idx = 1
    for teff, logg, mh in itertools.product(teff_vals, logg_vals, mh_vals):
        notfound = False
        if mh < 0:
            fname = f'lte{teff:05d}-{logg:.2f}{mh}.7.gz'
        elif mh == 0:
            fname = f'lte{teff:05d}-{logg:.2f}-{mh}.7.gz'
        else:
            fname = f'lte{teff:05d}-{logg:.2f}+{mh}.7.gz'
        curr_spectrum = []
        curr_wavelengths = []
        try:
            with gzip.open(fname, 'rb') as f_in:
                for line in f_in:
                    currline = line.decode().split()
                    curr_wavelengths.append(float(currline[0].replace('D', 'E')))
                    curr_spectrum.append(10**float(currline[1].replace('D', 'E')))
        except FileNotFoundError:
            continue
        curr_spectrum = np.array(curr_spectrum) * (
                astropy.units.erg / (astropy.units.cm)**2 / astropy.units.s / astropy.units.cm)
#        curr_wavelengths = np.array(curr_wavelengths) * astropy.units.Angstrom
        curr_wavelengths = np.array(curr_wavelengths) * astropy.units.Angstrom
        break_idx = 239253
        np.savetxt('actual_wavelength.txt', np.array([curr_wavelengths[:break_idx], curr_spectrum[:break_idx]]))
#        1/0
        curr_spectrum *= curr_wavelengths **2 / astropy.constants.c
#        curr_spectrum *= (curr_wavelengths ** 2 / astropy.constants.c) * astropy.units.Jy
        curr_freqs = curr_wavelengths.to(astropy.units.Hz,
                                         astropy.units.spectral())[::-1]

        curr_spectrum = curr_spectrum.to(astropy.units.Jy)
#        print(curr_spectrum)
#        print(np.min(curr_freqs))
        break_idx = np.where(curr_freqs[1:] - curr_freqs[0:-1] <= 0)[0][0] + 1
#        break_idx = np.where(curr_wavelengths[1:] - curr_wavelengths[0:-1] <= 0)[0][0] + 1
#        print(break_idx)
#        1/0
#        print(curr_freqs)
#        print(curr_freqs[1:] - curr_freqs[0:-1])

#        print(np.where((curr_freqs[1:] - curr_freqs[0:-1]) <= 0))
#        print((curr_freqs)[239252])
#        print((curr_freqs)[239251])
#        print((curr_freqs)[239253])

        spline = CubicSpline(
                curr_freqs[:break_idx],
                curr_spectrum[:break_idx])
        np.savetxt('spline.txt', np.array([curr_freqs[:break_idx], spline(curr_freqs[:break_idx])]))
        np.savetxt('actual.txt', np.array([curr_freqs[:break_idx], curr_spectrum[:break_idx]]))
#        1/0

        curr_reported_vals = []

        for i in range(1, 11):
            wavelength_ref = dirbe_utils.BAND_TO_WAVELENGTH[i] * astropy.units.micron
#            print(wavelength_ref)
            freq_ref = wavelength_ref.to(astropy.units.Hz,
                                         astropy.units.spectral())
#            print(freq_ref)
            if spline(freq_ref) <= 0:
               curr_reported_vals.append(0)
               continue
            wavelength, response = dirbe_utils.get_bandpass(i)
            freqs = wavelength.to(astropy.units.Hz,
                                  astropy.units.spectral())[::-1]
            spectrum = spline(freqs)
            K_factor = (np.trapz(response * spectrum, freqs) /
                        np.trapz(response / freqs, freqs)) / (freq_ref * spline(freq_ref))
#            print(i)
#            print(spline(freq_ref))
            curr_reported_vals.append(K_factor * spline(freq_ref))
        param_vals.append((teff, logg, mh))
        reported_vals.append(curr_reported_vals)
        print(f"Calculated for {fname}")
        
    reported_vals = np.array(reported_vals)
    param_vals = np.array(param_vals)
    band_mapping = [f'{band_int:02d}' for band_int in range(1, 11)]
    of.create_dataset('reported_values', data=reported_vals)
    of.create_dataset('parameter_values', data=param_vals)
    of.create_dataset('band_column_mapping', data=band_mapping)

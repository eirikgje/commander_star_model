import numpy as np
import gzip
import requests
from pathlib import Path
import astropy.units
import astropy.constants
from astropy.io import fits as pf
import scipy.integrate
import dirbe_utils
import h5py
from scipy.interpolate import CubicSpline
from multiprocessing import Pool


"""
This assumes an already created cross-match file between GAIA and AllWISE
(or whatever), and it will then first find appropriate PHOENIX spectra, then
interpolate these to the actual parameter values. The way it does this is that
it first finds the closest 'cube' of values that surrounds the parameter point
in space, and then weights those eight spectra appropriately.

NB! This version of the script assumes you have downloaded all the PHOENIX
spectra already, and it's multiprocessed.

It also requires dirbe_utils to be in your path.
"""

new_mh_vals = np.array([-2.0, -1.5, -0.5, 0.0, 0.5, 1.0])
new_logg_vals = np.array([3.0, 3.5, 4.0, 4.5, 5.0])
new_teff_vals = np.array(list([i for i in range(3000, 4050, 50)]) +
                         list([i for i in range(4000, 12100, 100)]))

old_mh_vals = np.array([-4.0, -3.0, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0])
old_logg_vals = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0])
old_teff_vals = np.array(list([i for i in range(2300, 7100, 100)]) +
                         list([i for i in range(7000, 12200, 200)]))
mh_norm = 3.0
logg_norm = 2.0
teff_norm = 9000
new_points = []
old_points = []

old_wavelength_array = None

#starfile = '/home/eirik/data/commander_star_model/crossmatch_commander_unique.csv'
starfile = '/mn/stornext/u3/eirikgje/data/commander_star_model/crossmatch_commander_unique.csv'
#spectrum_dir = '/home/eirik/data/commander_star_model/phoenix_spectra/'
#outfile = '/home/eirik/data/commander_star_model/commander_star_model.h5'

new_spectrum_dir = '/mn/stornext/u3/eirikgje/data/commander_star_model/phoenix_spectra/'
old_spectrum_dir = '/mn/stornext/u3/eirikgje/data/commander_star_model/old_phoenix_spectra/'
outfile = '/mn/stornext/u3/eirikgje/data/commander_star_model/commander_star_model_including_old_spectra_with_scaling_40kstrongest.h5'
#outfile = '/mn/stornext/u3/eirikgje/data/commander_star_model/commander_star_model_test.h5'
coordinates_in_radians = True  # Otherwise, will be in degrees


def populate_grid():
#    mh_vals = list(new_mh_vals)
#    mh_vals.extend([el for el in old_mh_vals if el not in new_mh_vals])
#    logg_vals = list(new_logg_vals)
#    logg_vals.extend([el for el in old_logg_vals if el not in new_logg_vals])
#    teff_vals = list(new_teff_vals)
#    teff_vals.extend([el for el in old_teff_vals if el not in new_teff_vals])
    new_grid = {}
    old_grid = {}
    for teff in new_teff_vals:
        new_grid[teff] = {}
        for logg in new_logg_vals:
            new_grid[teff][logg] = {}
            for mh in new_mh_vals:
                new_grid[teff][logg][mh] = get_filename(teff, logg, mh, use_new=True)
    for teff in old_teff_vals:
        old_grid[teff] = {}
        for logg in old_logg_vals:
            old_grid[teff][logg] = {}
            for mh in old_mh_vals:
                old_grid[teff][logg][mh] = get_filename(teff, logg, mh, use_new=False)

    return new_grid, old_grid


def get_filename(teff, logg, mh, new_spectrum_dir=new_spectrum_dir,
                 old_spectrum_dir=old_spectrum_dir, use_new=True):
    teff = int(teff)
    if use_new:
        if mh < 0:
            fname = f'lte{teff:05d}-{logg:.2f}{mh}.dat'
        elif mh == 0:
            fname = f'lte{teff:05d}-{logg:.2f}-{mh}.dat'
        else:
            fname = f'lte{teff:05d}-{logg:.2f}+{mh}.dat'
#        if mh < 0:
#            fname = f'lte{teff:05d}-{logg:.2f}{mh}.7.gz'
#        elif mh == 0:
#            fname = f'lte{teff:05d}-{logg:.2f}-{mh}.7.gz'
#        else:
#            fname = f'lte{teff:05d}-{logg:.2f}+{mh}.7.gz'
        if Path(new_spectrum_dir + fname).exists():
            return Path(new_spectrum_dir + fname)
        return None

    if mh < 0:
        fname = f'lte{teff:05d}-{logg:.2f}{mh}.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
    elif mh == 0:
        fname = f'lte{teff:05d}-{logg:.2f}-{mh}.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
    else:
        fname = f'lte{teff:05d}-{logg:.2f}+{mh}.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'

    if Path(old_spectrum_dir + fname).exists():
        return Path(old_spectrum_dir + fname)
    return None

def get_source_params(line):
    "Parses one line from the source file to get the parameters"
    sp_line = line.split(',')
    if sp_line[5] == '':
        return None
    params = {}
    params['teff'] = float(sp_line[5])
    params['logg'] = float(sp_line[8])
    params['mh'] = float(sp_line[11])
    params['l'] = float(sp_line[3])
    params['b'] = float(sp_line[4])
    params['gaia_id'] = str(sp_line[1])
    if coordinates_in_radians:
        params['l'] = params['l'] * np.pi / 180
        params['b'] = params['b'] * np.pi / 180
    return params


def get_interpolated_spectrum(teff, logg, mh, param_points,
                              fname_grid, use_new=True):

    if use_new:
        cube, distances = find_nearest_cube(teff, logg, mh, param_points,
                                            fname_grid, new_teff_vals,
                                            new_logg_vals, new_mh_vals,
                                            use_new=True)
    else:
        cube, distances = find_nearest_cube(teff, logg, mh, param_points,
                                            fname_grid, old_teff_vals,
                                            old_logg_vals, old_mh_vals,
                                            use_new=False)
    cube_spectra = []
    norm = 0
    prev_wavelength = None
    for point, distance in zip(cube, distances):
        norm += 1/distance
        if use_new:
            wavelength, spectrum = get_spectrum_new(*point, fname_grid)
        else:
            wavelength, spectrum = get_spectrum_old(*point, fname_grid)
        if prev_wavelength is not None:
            if not np.all(wavelength == prev_wavelength):
                raise ValueError
        cube_spectra.append([spectrum, 1/distance])
        prev_wavelength = wavelength
    final_spectrum = 0
    for spectrum, weight in cube_spectra:
        final_spectrum += weight * spectrum
    final_spectrum /= norm
    return wavelength, final_spectrum


# This was before I re-wrote the new data
#def get_spectrum_new(teff, logg, mh, fname_grid):
#    fname = fname_grid[teff][logg][mh]
#    datafile = np.loadtxt(fname, converters=lambda s: s.replace('D', 'E'), usecols=(0, 1), encoding=None)
#
#    sorted_data = datafile[np.argsort(datafile[:, 0])]
#    return sorted_data[:, 0], 10**sorted_data[:, 1]


def get_spectrum_new(teff, logg, mh, fname_grid):
    fname = fname_grid[teff][logg][mh]
    data = np.loadtxt(fname)
    return data[:, 0], data[:, 1]


def get_spectrum_old(teff, logg, mh, fname_grid):
    fname = fname_grid[teff][logg][mh]

    spectrum = pf.open(fname)[0].data

    global old_wavelength_array
    if old_wavelength_array is None:
        global old_spectrum_dir
        fname = 'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits'
        old_wavelength_array = pf.open(f'{old_spectrum_dir}/{fname}')[0].data
    return old_wavelength_array, spectrum


def find_distances(target_point, points, normalization):
    point_arr = np.array(points)
    distances = (target_point - point_arr) / normalization
    return np.linalg.norm(distances, axis=1)


#def verify_param_point(param_point, spectrum_dir, use_new):
#
#    teff, logg, mh = param_point
#    teff = int(teff)
#    if use_new:
#        if mh < 0:
#            fname = f'lte{teff:05d}-{logg:.2f}{mh}.7.gz'
#        elif mh == 0:
#            fname = f'lte{teff:05d}-{logg:.2f}-{mh}.7.gz'
#        else:
#            fname = f'lte{teff:05d}-{logg:.2f}+{mh}.7.gz'
#    else:
#        if mh < 0:
#            fname = f'lte{teff:05d}-{logg:.2f}{mh}.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
#        elif mh == 0:
#            fname = f'lte{teff:05d}-{logg:.2f}-{mh}.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
#        else:
#            fname = f'lte{teff:05d}-{logg:.2f}+{mh}.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
#
#    # Does it exist in the folder?
#    if Path(spectrum_dir + fname).exists():
#        return True
#    return False


def find_nearest_cube(teff, logg, mh, possible_points,
                      fname_grid, teff_vals, logg_vals, mh_vals,
                      teff_norm=teff_norm, logg_norm=logg_norm,
                      mh_norm=mh_norm, use_new=True):
    dist_arr = find_distances(
        np.array([teff, logg, mh]), possible_points,
        np.array([teff_norm, logg_norm, mh_norm]))
    sort_args = dist_arr.argsort()
    sorted_points = np.array(possible_points)[sort_args]
    sorted_dists = dist_arr[sort_args]
    num_found_points = 0
    morethan = np.array([0, 0, 0])
    lessthan = np.array([0, 0, 0])
    found_points = []
    found_distances = []

    k = 0
    target_found_points = 8
    if teff < teff_vals[0] or teff > teff_vals[-1]:
        target_found_points /= 2
    if logg < logg_vals[0] or logg > logg_vals[-1]:
        target_found_points /= 2
    if mh < mh_vals[0] or mh > mh_vals[-1]:
        target_found_points /= 2
    target_found_points = int(target_found_points)

    while num_found_points < target_found_points:
        curr_trial = sorted_points[k]
        curr_dist = sorted_dists[k]
        k += 1
        teff_curr, logg_curr, mh_curr = curr_trial
        morethan_temp = np.array([0, 0, 0])
        lessthan_temp = np.array([0, 0, 0])
        if teff_curr > teff:
            morethan_temp[0] += 1
        else:
            lessthan_temp[0] += 1
        if logg_curr > logg:
            morethan_temp[1] += 1
        else:
            lessthan_temp[1] += 1
        if mh_curr > mh:
            morethan_temp[2] += 1
        else:
            lessthan_temp[2] += 1

        if not (any((morethan_temp + morethan) > 4) or
                any((lessthan_temp + lessthan > 4))):
            validated = fname_grid[teff_curr][logg_curr][mh_curr] is not None
            if validated:
                morethan += morethan_temp
                lessthan += lessthan_temp
                found_points.append(curr_trial)
                found_distances.append(curr_dist)
                num_found_points += 1
        if (k == len(sorted_points) and
                num_found_points < target_found_points):
            print("Warning: Could not find a cube for this source:")
            print(teff, logg, mh)
            print("Number of points expected:")
            print(target_found_points)
            print("Number of points found:")
            print(num_found_points)
            num_found_points = target_found_points # Hack to bypass the while loop - this means there are no more spectra available and we still want to find some approximation

    return found_points, found_distances


def get_dirbe_reported_values(wavelength, spectrum, dirbe_bandpasses):
    wavelength *= astropy.units.Angstrom
    spectrum *= (astropy.units.erg / (astropy.units.cm)**2 / astropy.units.s /
                 astropy.units.cm)
    spectrum *= wavelength ** 2 / astropy.constants.c
    freqs = wavelength.to(astropy.units.Hz, astropy.units.spectral())[::-1]
    spectrum = spectrum[::-1]

    spectrum = spectrum.to(astropy.units.MJy)
    sed_spline = CubicSpline(freqs, spectrum)
    reported_vals = []

    for i in range(1, 11):
        wavelength_ref = (dirbe_utils.BAND_TO_WAVELENGTH[i] *
                          astropy.units.micron)
        freq_ref = wavelength_ref.to(astropy.units.Hz,
                                     astropy.units.spectral())
        if freq_ref < freqs[0]:
            reported_vals.append(0)
            continue
#        dirbe_wavelength, response = dirbe_utils.get_bandpass(i)
        dirbe_wavelength, response = dirbe_bandpasses[i-1]
        dirbe_freqs = dirbe_wavelength.to(
                astropy.units.Hz, astropy.units.spectral())[::-1]
        response = response[::-1]
        K_factor = ((np.trapz(response * sed_spline(dirbe_freqs), dirbe_freqs) /
                    np.trapz(response / dirbe_freqs, dirbe_freqs)) /
                    (freq_ref * sed_spline(freq_ref)))
        reported_vals.append(K_factor.value * sed_spline(freq_ref))
    return reported_vals


def is_within_new_parameter_bounds(teff, logg, mh):
    if teff < new_teff_vals[0] or teff > new_teff_vals[-1]:
        return False
    if logg < new_logg_vals[0] or logg > new_logg_vals[-1]:
        return False
    if mh < new_mh_vals[0] or mh > new_mh_vals[-1]:
        return False
    return True


def scale_reported_vals(new_reported_vals, old_reported_vals):
#    threetofour_new = new_reported_vals[2] - new_reported_vals[3]
#    threetofour_old = old_reported_vals[2] - old_reported_vals[3]
    new_val = (new_reported_vals[2] + new_reported_vals[3]) * 0.5
    old_val = (old_reported_vals[2] + old_reported_vals[3]) * 0.5
    scaling = old_val / new_val
    scaled_reported_vals = old_reported_vals
    scaled_reported_vals[4] = new_reported_vals[4] * scaling
    scaled_reported_vals[5] = new_reported_vals[5] * scaling
    return scaled_reported_vals



for teff_val in new_teff_vals:
    for logg_val in new_logg_vals:
        for mh_val in new_mh_vals:
            new_points.append((teff_val, logg_val, mh_val))


for teff_val in old_teff_vals:
    for logg_val in old_logg_vals:
        for mh_val in old_mh_vals:
            old_points.append((teff_val, logg_val, mh_val))

dirbe_bandpasses= []
for i in range(1, 11):
    dirbe_bandpasses.append(dirbe_utils.get_bandpass(i))



#num_sources_to_process = 80
#num_sources_to_process = 2000
num_sources_to_process = None
with open(starfile, 'r') as source_file:
    lines = source_file.readlines()

edited_lines = []
for line in lines:
    sp_line = line.split(',')
    if sp_line[5] == '':
        continue
    edited_lines.append(line)

if num_sources_to_process is not None:
    edited_lines = edited_lines[:num_sources_to_process]

coordinates = []
param_vals = []
reported_vals = []
gaia_ids = []
#num_processes = None
num_processes = 128
#num_processes = 8
#num_processes = 100

new_fname_grid, old_fname_grid = populate_grid()

def multiprocess_core_func(source_line, new_fname_grid=new_fname_grid,
                           old_fname_grid=old_fname_grid,
                           new_points=new_points,
                           old_points=old_points,
                           dirbe_bandpasses=dirbe_bandpasses):
    #    print(f"Processing source {source_line}")
    params = get_source_params(source_line)
    can_use_new_spectra = is_within_new_parameter_bounds(params['teff'],
                                                         params['logg'],
                                                         params['mh'])
    new_wavelength, new_spectrum = get_interpolated_spectrum(params['teff'],
                                                             params['logg'],
                                                             params['mh'],
                                                             new_points,
                                                             new_fname_grid,
                                                             use_new=True)
    if not can_use_new_spectra:
        old_wavelength, old_spectrum = get_interpolated_spectrum(params['teff'],
                                                                 params['logg'],
                                                                 params['mh'],
                                                                 old_points,
                                                                 old_fname_grid,
                                                                 use_new=False)

    new_reported_vals = get_dirbe_reported_values(new_wavelength,
                                                  new_spectrum,
                                                  dirbe_bandpasses)

    if not can_use_new_spectra:
        old_reported_vals = get_dirbe_reported_values(old_wavelength,
                                                      old_spectrum,
                                                      dirbe_bandpasses)
        
        scaled_reported_vals = scale_reported_vals(new_reported_vals, old_reported_vals)


        return params, scaled_reported_vals
    return params, new_reported_vals


if num_processes is None:
    for line in edited_lines:
        params, curr_reported_vals = multiprocess_core_func(line)
        print(params)
        print(curr_reported_vals)
else:
    with Pool(processes=num_processes) as pool:
        chunksize = int(len(edited_lines) / num_processes)
        print(chunksize)
        for result in pool.map(multiprocess_core_func, edited_lines, chunksize):
        #        print(result)
            params = result[0]
            curr_reported_vals = result[1]
            param_vals.append((params['teff'], params['logg'], params['mh']))
            coordinates.append((params['l'], params['b']))
            gaia_ids.append(params['gaia_id'])
            reported_vals.append(curr_reported_vals)
    #    for source_line in source_file:
    #        params = get_source_params(source_line)
    #        if params is None:
    #            continue
    #        print(params)
    #        wavelength, spectrum = (
    #                get_interpolated_spectrum(params['teff'], params['logg'],
    #                                          params['mh'], points))
    #        coordinates.append((params['l'], params['b']))
    #        param_vals.append((params['teff'], params['logg'], params['mh']))
    #        curr_reported_vals = get_dirbe_reported_values(wavelength, spectrum)
    #        print(curr_reported_vals)
    #        print("Source processed")
    #        reported_vals.append(curr_reported_vals)
    #        counter += 1
    #        if counter == 15:
    #            break

with h5py.File(outfile, 'w') as of:
    coordinates = np.array(coordinates)
    param_vals = np.array(param_vals)
    reported_vals = np.array(reported_vals)
    band_mapping = [f'{band_int:02d}' for band_int in range(1, 11)]
    of.create_dataset('reported_values', data=reported_vals)
    of.create_dataset('parameter_values', data=param_vals)
    of.create_dataset('band_column_mapping', data=band_mapping)
    of.create_dataset('coordinates', data=coordinates)
    of.create_dataset('gaia_ids', data=gaia_ids)

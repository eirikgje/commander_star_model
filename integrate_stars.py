import numpy as np
import gzip
import requests
from pathlib import Path
import astropy.units
import astropy.constants
import scipy.integrate
import dirbe_utils
import h5py
from scipy.interpolate import CubicSpline


"""
This assumes an already created cross-match file between GAIA and AllWISE
(or whatever), and it will then first find appropriate PHOENIX spectra, then
interpolate these to the actual parameter values. The way it does this is that
it first finds the closest 'cube' of values that surrounds the parameter point
in space, and then weights those eight spectra appropriately. The script will
make several requests to the PHOENIX server unless all the needed spectra are
already downloaded.

It also requires dirbe_utils to be in your path.
"""

mh_vals = np.array([-2.0, -1.5, -0.5, 0.0, 0.5, 1.0])
logg_vals = np.array([3.0, 3.5, 4.0, 4.5, 5.0])
teff_vals = np.array(list([i for i in range(3000, 4000, 50)]) +
                     list([i for i in range(4000, 12000, 100)]))

mh_norm = 3.0
logg_norm = 2.0
teff_norm = 9000
points = []

starfile = '/home/eirik/data/commander_star_model/crossmatch_commander_unique.csv'
url = 'https://www.astro.uni-jena.de/Users/theory/for2285-phoenix/'
spectrum_dir = '/home/eirik/data/commander_star_model/phoenix_spectra/'
outfile = '/home/eirik/data/commander_star_model/commander_star_model_debug.h5'

coordinates_in_radians = True  # Otherwise, will be in degrees


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
    if coordinates_in_radians:
        params['l'] = params['l'] * np.pi / 180
        params['b'] = params['b'] * np.pi / 180
    return params


def get_interpolated_spectrum(teff, logg, mh, param_points,
                              verified_param_points,
                              spectrum_dir=spectrum_dir):

    cube, distances, param_points, verified_param_points = find_nearest_cube(
            teff, logg, mh, param_points, verified_param_points, spectrum_dir)
    cube_spectra = []
    norm = 0
    prev_wavelength = None
    for point, distance in zip(cube, distances):
        norm += 1/distance
        wavelength, spectrum = get_spectrum(*point)
        if prev_wavelength is not None:
            if not np.all(wavelength == prev_wavelength):
                raise ValueError
        cube_spectra.append([spectrum, 1/distance])
        prev_wavelength = wavelength
    final_spectrum = 0
    for spectrum, weight in cube_spectra:
        final_spectrum += weight * spectrum
    final_spectrum /= norm
#    np.savetxt('wavelength.dat', wavelength)
#    np.savetxt('spectrum.dat', final_spectrum)
    return wavelength, final_spectrum, param_points, verified_param_points


def get_spectrum(teff, logg, mh, spectrum_dir=spectrum_dir):
    teff = int(teff)
    if mh < 0:
        fname = f'lte{teff:05d}-{logg:.2f}{mh}.7.gz'
    elif mh == 0:
        fname = f'lte{teff:05d}-{logg:.2f}-{mh}.7.gz'
    else:
        fname = f'lte{teff:05d}-{logg:.2f}+{mh}.7.gz'

    curr_spectrum = []
    curr_wavelengths = []
    prev_wavelength = 0
    with gzip.open(spectrum_dir + fname, 'rb') as f_in:
        for line in f_in:
            currline = line.decode().split()
            curr_wavelength = float(currline[0].replace('D', 'E'))
            if curr_wavelength <= prev_wavelength:
                break
            curr_wavelengths.append(curr_wavelength)
            prev_wavelength = curr_wavelength
            curr_spectrum.append(10**float(currline[1].replace('D', 'E')))
    return np.array(curr_wavelengths), np.array(curr_spectrum)


def find_distances(target_point, points, normalization):
    point_arr = np.array(points)
    distances = (target_point - point_arr) / normalization
    return np.linalg.norm(distances, axis=1)


def verify_param_point(param_point, spectrum_dir=spectrum_dir):

    teff, logg, mh = param_point
    teff = int(teff)
    if mh < 0:
        fname = f'lte{teff:05d}-{logg:.2f}{mh}.7.gz'
    elif mh == 0:
        fname = f'lte{teff:05d}-{logg:.2f}-{mh}.7.gz'
    else:
        fname = f'lte{teff:05d}-{logg:.2f}+{mh}.7.gz'

    # Does it exist in the folder?
    if Path(spectrum_dir + fname).exists():
        return True

    r = requests.get(url+fname)
    if r.ok:
        open(spectrum_dir+fname, 'wb').write(r.content)
        return True
    return False


def find_nearest_cube(teff, logg, mh, possible_points, verified_param_points,
                      spectrum_dir=spectrum_dir, teff_norm=teff_norm,
                      logg_norm=logg_norm, mh_norm=mh_norm):
    dist_arr = find_distances(
        np.array([teff, logg, mh]), possible_points,
        np.array([teff_norm, logg_norm, mh_norm]))
    sort_args = dist_arr.argsort()
    sorted_points = np.array(possible_points)[sort_args]
    sorted_dists = dist_arr[sort_args]
    num_found_points = 0
    morethan = np.array([0, 0, 0])
    lessthan = np.array([0, 0, 0])
    points_to_be_deleted = []
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
        accepted = True
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

        if any((morethan_temp + morethan) > 4):
            accepted = False
            if k != len(sorted_points):
                continue
        if any((lessthan_temp + lessthan) > 4):
            accepted = False
            if k != len(sorted_points):
                continue

        if accepted and (tuple(curr_trial) not in verified_param_points):
            validated = verify_param_point(curr_trial)
            if not validated:
                points_to_be_deleted.append(tuple(curr_trial))
                continue
            else:
                verified_param_points.append(tuple(curr_trial))
        if accepted:
            morethan += morethan_temp
            lessthan += lessthan_temp
            found_points.append(curr_trial)
            found_distances.append(curr_dist)
            num_found_points += 1
        elif (k == len(sorted_points) and
              num_found_points < target_found_points):
            print("Warning: Could not find a cube for this source:")
            print(teff, logg, mh)
            print("Number of points expected:")
            print(target_found_points)
            print("Number of points found:")
            print(num_found_points)
            num_found_points = target_found_points # Hack to bypass the while loop - this means there are no more spectra available and we still want to find some approximation
        else:
            raise ValueError("This is not supposed to happen. Call your local"
                             " Python hotline.")

    for deletable_point in points_to_be_deleted:
        possible_points.remove(deletable_point)
    return (found_points, found_distances, possible_points,
            verified_param_points)


def get_dirbe_reported_values(wavelength, spectrum):
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
        dirbe_wavelength, response = dirbe_utils.get_bandpass(i)
        dirbe_freqs = dirbe_wavelength.to(
                astropy.units.Hz, astropy.units.spectral())[::-1]
        response = response[::-1]
        K_factor = ((np.trapz(response * sed_spline(dirbe_freqs), dirbe_freqs) /
                    np.trapz(response / dirbe_freqs, dirbe_freqs)) /
                    (freq_ref * sed_spline(freq_ref)))
        reported_vals.append(K_factor.value * sed_spline(freq_ref))
    return reported_vals


for teff_val in teff_vals:
    for logg_val in logg_vals:
        for mh_val in mh_vals:
            points.append((teff_val, logg_val, mh_val))

verified_points = []
coordinates = []
param_vals = []
reported_vals = []

counter = 0
with open(starfile, 'r') as source_file:
    source_file.readline()  # Skip the first line
    for source_line in source_file:
        params = get_source_params(source_line)
        if params is None:
            continue
        print(params)
        wavelength, spectrum, points, verified_points = (
                get_interpolated_spectrum(params['teff'], params['logg'],
                                          params['mh'], points,
                                          verified_points))
        coordinates.append((params['l'], params['b']))
        param_vals.append((params['teff'], params['logg'], params['mh']))
        curr_reported_vals = get_dirbe_reported_values(wavelength, spectrum)
        print(curr_reported_vals)
        print("Source processed")
        reported_vals.append(curr_reported_vals)
        counter += 1
        if counter == 20:
            break

with h5py.File(outfile, 'w') as of:
    coordinates = np.array(coordinates)
    param_vals = np.array(param_vals)
    reported_vals = np.array(reported_vals)
    band_mapping = [f'{band_int:02d}' for band_int in range(1, 11)]
    of.create_dataset('reported_values', data=reported_vals)
    of.create_dataset('parameter_values', data=param_vals)
    of.create_dataset('band_column_mapping', data=band_mapping)
    of.create_dataset('coordinates', data=coordinates)

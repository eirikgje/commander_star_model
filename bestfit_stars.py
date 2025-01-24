import numpy as np
import requests
import shutil

mh_vals = np.array([-2.0, -1.5, -0.5, 0.0, 0.5, 1.0])
logg_vals = np.array([3.0, 3.5, 4.0, 4.5, 5.0])
teff_vals = np.array(list([i for i in range(3000, 4000, 50)]) +
                     list([i for i in range(4000, 12000, 100)]))

mh_norm = 3.0
logg_norm = 2.0
teff_norm = 9000
points = []
for teff_val in teff_vals:
    for logg_val in logg_vals:
        for mh_val in mh_vals:
#             points.append((teff_val/teff_norm, logg_val/logg_norm, mh_val/mh_norm))
             points.append((teff_val, logg_val, mh_val))


#starfile = '/home/eirik/data/commander_star_model/crossmatch_commander_unique.csv'
starfile = '/mn/stornext/u3/eirikgje/data/commander_star_model/crossmatch_commander_unique.csv'
output_dir = '/mn/stornext/u3/eirikgje/data/commander_star_model/'

param_vals = {}


def find_nearest(val, array, excluded_indices=None, teff_norm=teff_norm,
        logg_norm=logg_norm, mh_norm=mh_norm):
    min_dist = 1000000000000
    min_idx = -1
    norm_val = np.array([val[0]/teff_norm, val[1]/logg_norm, val[2]/mh_norm])
    for idx, point in enumerate(array):
        if excluded_indices is not None and idx in excluded_indices:
            continue
        norm_point = np.array([point[0]/teff_norm, point[1]/logg_norm,
                               point[2]/mh_norm])
        if np.linalg.norm(norm_val - norm_point) < min_dist:
            min_dist = np.linalg.norm(norm_val - norm_point)
            min_idx = idx
    return min_idx, array[min_idx]


#k = 0
with open(starfile, 'r') as source_file:
    source_file.readline()
    for line in source_file:
#        if k == 100:
#            break
        if len(param_vals) % 50 == 0:
            print(len(param_vals))
        sp_line = line.split(',')
        if sp_line[5] == '':
            continue
        comm_source_id = int(sp_line[0])
        ang_sep = float(sp_line[2])
        l = float(sp_line[3])
        b = float(sp_line[4])
        teff = float(sp_line[5])
        teff_low = float(sp_line[6])
        teff_high = float(sp_line[7])
        logg = float(sp_line[8])
        logg_low = float(sp_line[9])
        logg_high = float(sp_line[10])
        mh = float(sp_line[11])
        mh_low = float(sp_line[12])
        mh_high = float(sp_line[13])
        idx, param_tuple = find_nearest((teff, logg, mh), points)
#        param_tuple = (listed_teff, listed_logg, listed_mh)
        if (idx, param_tuple) in param_vals.keys():
            param_vals[(idx, param_tuple)].append((comm_source_id,
                                                  (teff, logg, mh)))
        else:
            param_vals[(idx, param_tuple)] = [(comm_source_id,
                                              (teff, logg, mh))]
#        k += 1

print(param_vals)

url = 'https://www.astro.uni-jena.de/Users/theory/for2285-phoenix/'
excluded_indices = []
downloaded_indices = []

while True:
    any_not_found = False
    for key in param_vals.keys():
        idx = key[0]
        if idx in downloaded_indices:
            continue
        param_tuple = key[1]
        teff = param_tuple[0]
        logg = param_tuple[1]
        mh = param_tuple[2]
        print(f"{teff}, {logg}, {mh}")
        if mh < 0:
            fname = f'lte{teff:05d}-{logg:.2f}{mh}.7.gz'
        elif mh == 0:
            fname = f'lte{teff:05d}-{logg:.2f}-{mh}.7.gz'
        else:
            fname = f'lte{teff:05d}-{logg:.2f}+{mh}.7.gz'
        r = requests.get(url+fname)
        if r.ok:
            print("Okay")
            open(fname, 'wb').write(r.content)
            downloaded_indices.append(idx)
            shutil.move(fname, output_dir+fname)
        else:
            print("Not okay")
            excluded_indices.append(idx)
            comm_source_list = param_vals[key]
            print(comm_source_list)
            del param_vals[key]
            for source in comm_source_list:
                comm_source_id = source[0]
                teff, logg, mh = source[1]
                idx, param_tuple = find_nearest((teff, logg, mh),
                                                points,
                                                excluded_indices=excluded_indices)
                if (idx, param_tuple) in param_vals.keys():
                    param_vals[(idx, param_tuple)].append((comm_source_id,
                                                          (teff, logg, mh)))
                else:
                    param_vals[(idx, param_tuple)] = [(comm_source_id,
                                                      (teff, logg, mh))]
            any_not_found = True
            break
    if not any_not_found:
        break

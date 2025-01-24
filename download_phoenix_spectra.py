import numpy as np
from pathlib import Path
import requests
import itertools


output_dir = '/home/eirik/data/commander_star_model/phoenix_spectra/'
url = 'https://www.astro.uni-jena.de/Users/theory/for2285-phoenix/'

mh_vals = np.array([-2.0, -1.5, -0.5, 0.0, 0.5, 1.0])
logg_vals = np.array([3.0, 3.5, 4.0, 4.5, 5.0])
teff_vals = np.array(list([i for i in range(3000, 4000, 50)]) +
                     list([i for i in range(4000, 12000, 100)]))

for teff, logg, mh in itertools.product(teff_vals, logg_vals, mh_vals):
    if mh < 0:
        fname = f'lte{teff:05d}-{logg:.2f}{mh}.7.gz'
    elif mh == 0:
        fname = f'lte{teff:05d}-{logg:.2f}-{mh}.7.gz'
    else:
        fname = f'lte{teff:05d}-{logg:.2f}+{mh}.7.gz'
    currpath = Path(output_dir + fname)
    if currpath.exists():
        continue

    print(f"Downloading {fname}")

    r = requests.get(url+fname)
    if r.ok:
        open(output_dir + fname, 'wb').write(r.content)
    else:
        print(f"{fname} not downloaded")

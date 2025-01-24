#from astroquery.gaia import Gaia
#import astropy.units
#import astropy.coordinates

#Gaia.MAIN_GAIA_TABLE = 'gaiadr3.gaia_source'
#commander_star_fname = '/home/eirik/data/commander_star_model/AllWISE_catalog_W1_mag8_v1.dat'
#
#elements = []
#not_found = []
#
#with open(commander_star_fname, 'r') as starfile:
#    k = 0
#    for currstar in starfile:
#        #        if k >= 10:
#        #            break
#        star_params = currstar.split()
#        l, b = (star_params[0], star_params[1])
#        coords = astropy.coordinates.SkyCoord(l, b, unit='deg',
#                                              frame='galactic')
#        margin = 10 * astropy.units.arcsec
#        res = Gaia.query_object(coords, radius=margin)
#        k += 1
#        if not len(res):
#            print("No result")
#            print(f"{k}, {l}, {b}")
#            not_found.append((k, l, b))
#            continue
#        elements.append((res[0]['source_id'], star_params[-1]))

wise_gaia_mapping = {}

with open('/home/eirik/data/commander_star_model/crossref_wise_gaia.csv', 'r') as wise_crossref:
    for line in wise_crossref:
        spline = line.split(',')
        if spline[0] == 'source_id':
            continue
        if spline[1] not in wise_gaia_mapping.keys():
            wise_gaia_mapping[spline[1]] = spline[0]

star_ids = []
not_found = []
with open('/home/eirik/data/commander_star_model/AllWISE_catalog_W1_mag8_v1.dat', 'r') as comm_starfile:
    for line in comm_starfile:
        spline = line.split()
        try:
            star_ids.append(wise_gaia_mapping[spline[-1]])
        except KeyError:
            not_found.append(spline[-1])



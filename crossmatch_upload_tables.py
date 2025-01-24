from astroquery.gaia import Gaia
from astropy.table import Table
import time

username = 'egjerl01'
directory = '/home/eirik/data/commander_star_model/'
#fname = directory + 'AllWISE_catalog_W1_mag8_v1_reduced.dat'
fname = directory + 'AllWISE_catalog_W1_mag8_v1.dat'
arcsec = 5
#arcsec = 10
ls = []
bs = []

login = True
create_initial_table = True
create_crossmatch_table = True

if login:
    Gaia.login()

if create_initial_table:
    Gaia.delete_user_table('commander_sources')
    with open(fname, 'r') as starfile:
        for line in starfile:
            sp_line = line.split()
            ls.append(float(sp_line[0]))
            bs.append(float(sp_line[1]))

    table = Table([ls, bs], names=["l_galactic", "b_galactic"])
    Gaia.upload_table(upload_resource=table, table_name='commander_sources')

if create_crossmatch_table:
    Gaia.delete_user_table('commander_crossmatch')
    full_qualified_table_name = f"user_{username}.commander_sources"
    query = f'SELECT test.commander_sources_oid, gaia.source_id, DISTANCE(POINT(test.l_galactic, test.b_galactic),POINT(gaia.l, gaia.b)) AS ang_sep FROM {full_qualified_table_name} AS test JOIN gaiadr3.gaia_source AS gaia ON DISTANCE(test.l_galactic, test.b_galactic, gaia.l, gaia.b) < {arcsec}./3600'
    job = Gaia.launch_job_async(query)
    result = job.get_results()
    Gaia.upload_table(upload_resource=result, table_name='commander_crossmatch')

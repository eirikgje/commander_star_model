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
create_initial_table = False
create_crossmatch_table = False

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
#    time.sleep(60)


full_qualified_table_name = f"user_{username}.commander_crossmatch"

query = f'SELECT crossmatch.commander_sources_oid, crossmatch.source_id, crossmatch.ang_sep, gaia.l, gaia.b, gaia.teff_gspphot, gaia.teff_gspphot_lower, gaia.teff_gspphot_upper, gaia.logg_gspphot, gaia.logg_gspphot_lower, gaia.logg_gspphot_upper, gaia.mh_gspphot, gaia.mh_gspphot_lower, gaia.mh_gspphot_upper, gaia.non_single_star FROM {full_qualified_table_name} AS crossmatch JOIN gaiadr3.gaia_source AS gaia ON crossmatch.source_id = gaia.source_id'
job = Gaia.launch_job_async(query)
#table_1 = job.get_results()
result = job.get_results()
#
#query = f'SELECT crossmatch.commander_sources_oid, crossmatch.source_id, crossmatch.ang_sep, params.teff_msc1, params.teff_msc1_lower, params.teff_msc1_upper, params.teff_msc2, params.teff_msc2_lower, params.teff_msc2_upper, params.logg_msc1, params.logg_msc1_lower, params.logg_msc1_upper, params.logg_msc2, params.logg_msc2_lower, params.logg_msc2_upper, params.mh_msc, params.mh_msc_lower, params.mh_msc_upper FROM {full_qualified_table_name} AS crossmatch JOIN gaiadr3.astrophysical_parameters AS params ON crossmatch.source_id = params.source_id'
#
#job = Gaia.launch_job_async(query)
#table_2 = job.get_results()

#query = f'SELECT crossmatch.commander_sources_oid, crossmatch.source_id, crossmatch.ang_sep, gaia.l, gaia.b, gaia.teff_gspphot, gaia.teff_gspphot_lower, gaia.teff_gspphot_upper, gaia.logg_gspphot, gaia.logg_gspphot_lower, gaia.logg_gspphot_upper, gaia.mh_gspphot, gaia.mh_gspphot_lower, gaia.mh_gspphot_upper, gaia.non_single_star, params.teff_msc1, params.teff_msc1_lower, params.teff_msc1_upper, params.teff_msc2, params.teff_msc2_lower, params.teff_msc2_upper, params.logg_msc1, params.logg_msc1_lower, params.logg_msc1_upper, params.logg_msc2, params.logg_msc2_lower, params.logg_msc2_upper, params.mh_msc, params.mh_msc_lower, params.mh_msc_upper FROM {full_qualified_table_name} AS crossmatch JOIN gaiadr3.astrophysical_parameters AS params ON crossmatch.source_id = params.source_id JOIN gaiadr3.gaia_source AS gaia ON crossmatch.source_id = gaia.source_id'
job = Gaia.launch_job_async(query)
result = job.get_results()

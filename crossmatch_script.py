from astroquery.gaia import Gaia
from astropy.table import Table

username = 'egjerl01'
#directory = '/home/eirik/data/commander_star_model/'
directory = '/mn/stornext/u3/eirikgje/data/commander_star_model/'
#fname = directory + 'AllWISE_catalog_W1_mag8_v1_reduced.dat'
fname = directory + 'AllWISE_catalog_W1_mag8_v1.dat'
arcsec = 5
#arcsec = 10
ls = []
bs = []

login = True
create_table = True

if login:
    Gaia.login()

if create_table:
    Gaia.delete_user_table('commander_sources')
    with open(fname, 'r') as starfile:
        for line in starfile:
            sp_line = line.split()
            ls.append(float(sp_line[0]))
            bs.append(float(sp_line[1]))

    table = Table([ls, bs], names=["l_galactic", "b_galactic"])
    Gaia.upload_table(upload_resource=table, table_name='commander_sources')

full_qualified_table_name = f"user_{username}.commander_sources"

#query = f'SELECT test.*, gaia.*,DISTANCE(POINT(test.l_galactic, test.b_galactic),POINT(gaia.l, gaia.b)) AS ang_sep FROM {full_qualified_table_name} AS test JOIN gaiadr3.gaia_source AS gaia ON DISTANCE(test.l_galactic, test.b_galactic, gaia.l, gaia.b) < {arcsec}./3600'
query = f'SELECT test.commander_sources_oid, test.l_galactic, test.b_galactic,gaia.solution_id, gaia.DESIGNATION, gaia.source_id, gaia.l, gaia.b, gaia.teff_gspphot, gaia.teff_gspphot_lower, gaia.teff_gspphot_upper, gaia.logg_gspphot, gaia.logg_gspphot_lower, gaia.logg_gspphot_upper, gaia.mh_gspphot, gaia.mh_gspphot_lower, gaia.mh_gspphot_upper, DISTANCE(POINT(test.l_galactic, test.b_galactic),POINT(gaia.l, gaia.b)) AS ang_sep FROM {full_qualified_table_name} AS test JOIN gaiadr3.gaia_source AS gaia ON DISTANCE(test.l_galactic, test.b_galactic, gaia.l, gaia.b) < {arcsec}./3600'

job = Gaia.launch_job_async(query)
result = job.get_results()

#Gaia.update_user_table(table_name=full_qualified_table_name,
#                       list_of_changes=[['"Ra"', "flags", "Ra"],
#                                        ['"Dec"', "flags", "Dec"]])
#
#crossmatch_table = 'crossmatch_table'
#
#Gaia.cross_match(full_qualified_table_name_a=full_qualified_table_name,
#                 full_qualified_table_name_b='gaiadr3.gaia_source',
#                 results_table_name=crossmatch_table, radius=1.0)

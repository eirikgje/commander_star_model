import numpy as np

crossmatch_file = '/home/eirik/data/commander_star_model/crossmatch_commander.csv'
output_file = '/home/eirik/data/commander_star_model/crossmatch_commander_unique.csv'

data = {}

with open(crossmatch_file, 'r') as source_file:
    source_file.readline()
    for line in source_file:
        sp_line = line.split(',')
        comm_source_id = int(sp_line[0])
        ang_sep = float(sp_line[2])
#        l = float(sp_line[3])
#        b = float(sp_line[4])
#        teff = float(sp_line[5])
#        teff_low = float(sp_line[6])
#        teff_high = float(sp_line[7])
#        logg = float(sp_line[8])
#        logg_low = float(sp_line[9])
#        logg_high = float(sp_line[10])
#        mh = float(sp_line[11])
#        mh_low = float(sp_line[12])
#        mh_high = float(sp_line[13])
        insert = False
        if comm_source_id in data.keys():
            if ang_sep < data[comm_source_id][-1]:
                insert = True
        else:
            insert = True
        if insert:
            gaia_source_id = int(sp_line[1])
#            data[comm_source_id] = (gaia_source_id, l, b, ang_sep, teff, teff_low, teff_high, logg, logg_low, logg_high, mh, mh_low, mh_high)
            data[comm_source_id] = (line, ang_sep)

sorted_keys = np.sort(list(data.keys()))
with open(output_file, 'w') as outfile:
    for key in sorted_keys:
        outfile.write(data[key][0])
#        vals = data[key]
#        outline = f'{key}, {vals[0]}, {vals[1]}, {vals[2]}, {vals[3]}, {vals[4]}, {vals[5]}, {vals[6]}, {vals[7]}, {vals[8]}, {vals[9]}, {vals[10]}, {vals[11]}, {vals[12]}\n'
#        outfile.write(outline)

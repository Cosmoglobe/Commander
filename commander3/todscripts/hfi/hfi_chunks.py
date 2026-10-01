import sqlite3
import numpy as np

pid_database = '/mn/stornext/d23/cmbco/globe/orig/planck/aux/hfi_raw_rings_v3.db'

conn = sqlite3.connect(pid_database)
c = conn.cursor()

dbs = c.execute("SELECT * FROM ring_times_hfi WHERE pointID_unique LIKE '%S%'")

print(c.description)


scannum = 1
scanlist = []

for ring in dbs:
    scanlist.append(ring[1])

print(scanlist, len(scanlist))

np.savetxt("ring_mapping_hfi.txt", scanlist, fmt="%s")

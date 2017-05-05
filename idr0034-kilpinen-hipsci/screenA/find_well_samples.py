#!/usr/bin/env python

from os import environ
from sys import stderr
from pandas import read_csv
from psycopg2 import connect


pw = environ.get("PASSWORD", "secret")
name = "idr0034-kilpinen-hipsci/screenA"
df = read_csv("idr0034-screenA-wellsToExclude.txt", sep="\t", header=1)
conn = connect("dbname=idr user=omero password=%s host=192.168.21.5" % pw)
cur = conn.cursor()

plates = dict()
rows = list()
for idx, entry in df.iterrows():
    # Plate	Plate Column	Plate Row	Well	Well Number	Comments
    plate, col, row, well, num, comment = entry
    rows.append((plate, col, row, well, num, comment))
    if plate not in plates:
	    cur.execute(("select p.id from plate p, screen s, screenplatelink spl "
                         "where p.id = spl.child and spl.parent = s.id "
                         "  and s.name = %s and p.name = %s"),
			(name, plate))
            rv = cur.fetchall()
            if len(rv) != 1:
                print >>stderr, "Skipping", plate, rv
            else:
                plates[plate] = rv[0]

COLS = ["id", "permissions", "creation_id", "group_id", "owner_id", "update_id", "image", "plateacquisition", "well", "well_index"]
print "\t".join(COLS)
select = ",".join(["ws.%s" % c for c in COLS])

for row in rows:
    plate, col, row, well, num, comment = row
    if plate not in plates:
        continue
    pid = plates[plate]
    cur.execute(("select %s "
                 "from wellsample ws, plate p, well w "
                 "where ws.well = w.id and w.plate = p.id and p.id = %%s "
                 "  and w.column = %%s and w.row = %%s") % select,
                (pid, col-1, row-1))
    for rec in cur:
        print "\t".join([str(x) for x in rec])

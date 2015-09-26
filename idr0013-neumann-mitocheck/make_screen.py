#!/usr/bin/env python

from argparse import ArgumentParser
from os.path import basename

ROWS = 16
COLUMNS = 24
PLATE = """
[Plate]
Name = Mitocheck
Rows = %s
Columns = %s
Fields = 1
""" % (ROWS, COLUMNS)

WELL = """
[Well %(well)s]
Row = %(row)s 
Column = %(col)s
Field_0 = %(file)s
"""

parser = ArgumentParser()
parser.add_argument("file", nargs="+")
ns = parser.parse_args()

well_map = {}
def make_well_map_broken():
	for x in range(1, 385):
	    row = (x / COLUMNS) + 1
	    col = (x % COLUMNS) + 1
	    col = "?ABCDEFGHIJKLMNOPQRSTUVWXYZ"[col]
	    well_map[x] = "%s%s" % (col, row)
	return well_wap

def make_well_map2():
	with open("map.txt") as f:
		txt = f.read()
		for line in txt.split("\n"):
			parts = line.split()
			if parts:
				well_map[int(parts[0])] = parts[1]
make_well_map2()


print PLATE
for well in range(1, 385):
        f = "test/%05d_01.ch5" % well
	if f not in ns.file:
		print "# missing: %s" % f
	else:
		h5 = basename(f)
	cr = well_map[well]
	row = "ABCDEFGHIJKLMNOPQRSTUVWXYZ".find(cr[0])
	col = int(cr[1:])-1
	print WELL % {
		"well": well-1,
		"col": col,
		"row": row,
		"file": f,
	}

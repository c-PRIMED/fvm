#!/usr/bin/env python

"""
Run lammps and compare output with golden
"""

import sys, os

err = 0

if len(sys.argv) != 5:
    print "USAGE: %s lammps_executable mpi_args input golden_output" % sys.argv[0]
    sys.exit(-1)

# Read the next line.
# skip empty lines and comments
# and lines with "time" in them
def nextline(f):
    while True:
        line = f.readline()
        if line == '':
            break
        line = line.strip()
        if line == '' or line[0] == '#':
            continue
        if line.find("time") < 0:
            break
    return line

try:
    gf = open(sys.argv[4])
except:
    print "Failed to open gloden file %s" % sys.argv[4]
    sys.exit(-1)

try: os.remove("log.lammps")
except: pass

ret = os.system("mpirun %s %s < %s > /dev/null" % (sys.argv[2], sys.argv[1], sys.argv[3]))
if ret:
    print "ERROR: execution of mpirun failed: %s" % ret
    sys.exit(-1)

try:
    f = open("log.lammps")
except:
    print "Failed to open output file log.lammps"
    sys.exit(-1)

while True:
    line = nextline(f)
    gline = nextline(gf)
    if line != gline:
        print "ERROR: Got\n%s\nExpected\n%s\n" % (line, gline)
        err = 1
    if line == '' or gline == '':
        break

if err:
    print "FAILED"
else:
    print "OK"
sys.exit(err)

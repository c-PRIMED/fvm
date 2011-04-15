#!/usr/bin/env python

"""
Given script and destination directory,
this run the script in the destinationdirectory and compared
results in destination directory/GOLDEN and return 0 in success

Usage: ptest.py [options]
Options:
  --outdir       Directory where test data is written. Default is current dir.
  --np           Number of processors. Default is 1.
  --in           Input file (required).
  --golden       Directory where golden results are kept.
  --script       Script to run (required).
"""

import sys, os
from optparse import OptionParser
import numfile_compare

def usage():
    print __doc__
    sys.exit(-1)

def cleanup(ec):
    sys.exit(ec)

  
def check_tip_displacement(options):

    check_file     = os.path.join(options.outdir, 'tipDisplacement-se.dat')
    os.system("more "+check_file)
    reference_file = os.path.join(options.golden, 'tipDisplacement-se.dat')
    if numfile_compare.check_files(check_file, reference_file, 1.0e-8):
       print "convergence history files didn't match GOLDEN/* !!!!!!!!!!!!!"
       return  -1
      # break;
    print "convergence history files are ok"
    return 0

def main():
    funcs = {
        'mainCantilever2D_solid1_elec1.py':[check_tip_displacement],
        'mainCantilever2D_solid1_elec2.py':[check_tip_displacement],
        'mainCantilever2D_solid1_elec3.py':[check_tip_displacement],
        'mainCantilever2D_solid1_elec4.py':[check_tip_displacement],
        'mainCantilever2D_solid1_elec8.py':[check_tip_displacement],
        'mainCantilever2D_solid1_elec16.py':[check_tip_displacement],
        'mainCantilever2D_solid1_elec32.py':[check_tip_displacement]
        }
    parser = OptionParser()
    parser.set_defaults(np=1,outdir='',type='tri')
    parser.add_option("--in", dest='infile', help="Input file (required).")
    parser.add_option("--golden", help="Directory where golden files are stored.")
    parser.add_option("--np", type=int, help="Number of Processors.")
    parser.add_option("--outdir", help="Output directory.")
    parser.add_option("--script", help="Script to run.")
    parser.add_option("--type", help="'tri'[default], 'quad', 'hexa', or 'tetra'")
    parser.add_option("--xdmf", action='store_true', help="Dump data in xdmf")
    (options, args) = parser.parse_args()

    # convert path to absolute because we may cd to the test directory
    options.script   = os.path.abspath(options.script)
    options.outdir   = os.path.abspath(options.outdir)
    options.golden   = os.path.abspath(os.path.join(options.golden,'GOLDEN'))


    if options.outdir:
        if not os.path.isdir(options.outdir):
	    try:
		        os.makedirs(options.outdir)
            except:
                print "error creating directory " + options.outdir
                sys.exit(1)
        os.chdir(options.outdir)
    print options.outdir
    mpi = 'mpirun -np %s %s' % (options.np, options.script);
    print mpi

    if os.system( mpi ):
        cleanup(-1)

    for f in funcs[os.path.basename(options.script)]:
        err = f(options)
        if err:
            cleanup(err)

if __name__ == "__main__":
    main()


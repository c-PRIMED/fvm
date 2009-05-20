#!/usr/bin/env python

"""
Given script and destination directory,
this run the script in the destinationdirectory and compared
results in destination directory/GOLDEN and return 0 in success

Usage: test.py [options]
Options:
  --np           Number of processors. Default is 1.
  --mesh         Fluent case file.
  --nsweep       Number of sweep.
  --golden       Directory where golden files are stored.
  --outdir       Output directory
"""

import sys, os
from optparse import OptionParser
import numfile_compare

def usage():
    print __doc__
    sys.exit(-1)

def cleanup(ec):
    sys.exit(ec)

def check_particles(options):
    check_file = os.path.join(options.outdir, "cavity_nsweep%s.dat" % options.nsweep)
    reference_file = os.path.join(options.golden, "cavity_nsweep%s.dat" % options.nsweep)
    if numfile_compare.check_files(check_file, reference_file, 1.0e-5):
         print "File didn't match GOLDEN/* !!!!!!!!!!!!!"
         return  -1
    print "Particle file is ok"
    return 0


def main():
    parser = OptionParser()
    parser.set_defaults(np=1)
    parser.add_option("--np", type=int, help="Number of Processors.")
    parser.add_option("--golden", help="Directory where golden files are stored.")
    parser.add_option("--outdir", help="Output directory.")
    parser.add_option("--nsweep", type=int, help="number of sweep")
    parser.add_option("--mesh", help="fluent case file")
    (options, args) = parser.parse_args()

    # convert path to absolute because we may cd to the test directory
    options.outdir   = os.path.abspath(options.outdir)
    options.golden   = os.path.abspath(os.path.join(options.golden,'GOLDEN'))
    options.mesh   = os.path.abspath(options.mesh)
    script = os.path.abspath('./FVMParticleMarking.py')

    if options.outdir:
        if not os.path.isdir(options.outdir):
            try:
                os.makedirs(options.outdir)
            except:
                fatal("error creating directory " + options.outdir)
        os.chdir(options.outdir)

    mpi = 'mpirun -np %s %s %s %s' % (options.np, script, options.mesh, options.nsweep)    
    if os.system( mpi ):
        cleanup(-1)
    
    # OK.  Now check results.
    err = check_particles( options )
    cleanup(err)

if __name__ == "__main__":
    main()


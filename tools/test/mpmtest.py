#!/usr/bin/env python

"""
Given a directory containing MPM inputs and expected
outputs, run MPM and compare the outputs. Returns 0 on success.

Usage: mpmtest.py [options] input_dirname
Options:
  --datadir      Directory where test data is stored. Default is current dir.
  --np           Number of processors. Default in 1. 
"""

import sys, os, shutil
from optparse import OptionParser

def usage():
    print __doc__
    sys.exit(-1)

def cleanup(ec):
    sys.exit(ec)

def checknc(file, indir, mpi):

    # maybe this should be an option?
    maxerr = .001

    if mpi:
        use_part = "-p"
    else:
        use_part = ""

    count = 0
    exe = os.popen("ncdiff -e %s %s %s %s" % (maxerr, use_part, file, os.path.join(indir,file)))
    for line in exe:
        if count < 100:
            print line.rstrip()
            count += 1
        else:
            break
    if count == 100: print "[truncated]"
    if exe.close():
        return -1
    return count

def main():
    parser = OptionParser()
    parser.set_defaults(np=1)
    parser.add_option("--datadir",
                      action="store", dest="datadir", help="Data directory name.")
    parser.add_option("--np",
                      action="store", dest="np", type=int, help="Number of Processors.")
    (options, args) = parser.parse_args()

    if len(args) != 1:
        usage()

    dirname = os.path.abspath(args[0])
    if not os.path.isdir(dirname):
        print "%s is not a valid input directory." % dirname
        cleanup(-1)

    if options.np > 1:
        mpi = "mpirun -np %s " % options.np
    else:
        mpi = ""

    if options.datadir:
        if not os.path.isdir(options.datadir):
            try:
                os.makedirs(options.datadir)
            except:
                fatal("error creating directory " + options.datadir)
        os.chdir(options.datadir)

    # copy input data files to test directory
    try:
        shutil.copy(os.path.join(dirname,'pm2geometry'), '.')
        shutil.copy(os.path.join(dirname,'pm2input'), '.')
    except:
        print "Failed to copy pm2geometry and pm2input from %s" % dirname
        cleanup(-1)
            
    # run PM2-Pre and send stdout and stderr to a log file.
    if os.system("/bin/bash -c '%sPM2-Pre > PM2-Pre.out 2>&1'" % mpi):
        cleanup(-1)

    # run PM2 and send stdout and stderr to a log file.
    if os.system("/bin/bash -c '%sPM2 > PM2.out 2>&1'" % mpi):
        cleanup(-1)

    # OK.  MPM finished. Now check results.
    err = checknc('pmgrid.nc', dirname, False)
    if err:
        cleanup(err)

    err = checknc('pmpart.nc', dirname, mpi)
    cleanup(err)

if __name__ == "__main__":
    main()

#!/usr/bin/env python

"""
Given script and destination directory,
this run the script in the destinationdirectory and compared
results in destination directory/GOLDEN and return 0 in success

Usage: ptest.py [options] input_dirname
Options:
  --datadir      Directory where test data is stored. Default is current dir.
  --np           Number of processors. Default in 1. 
"""

import sys, os, shutil
from optparse import OptionParser
import filecmp


def usage():
    print __doc__
    sys.exit(-1)

def cleanup(ec):
    sys.exit(ec)

def check_particles(options):
    
    check_file = []
    reference_file = []
    check_file.append ( str(options.datadir) + "/cavity_nsweep" + str(options.nsweep) + ".dat" )
    reference_file.append( str(options.datadir) + "/GOLDEN/cavity_nsweep" + str(options.nsweep) + ".dat"  )
    print check_file
    print reference_file
    if  not(filecmp.cmp( check_file[0], reference_file[0], 0)) :
       print "mesh files didn't match GOLDEN/* !!!!!!!!!!!!!"
       return  -1
    print "mesh files are ok"
    return 0


def main():
    parser = OptionParser()
    parser.set_defaults(np=1)  
    parser.set_defaults(test="mesh")
    parser.set_defaults(runmode="run")
    parser.add_option("--script",
                      action="store", type="string", dest="scriptname", help="Name of script ")
    parser.add_option("--np",
                      action="store", dest="np", type=int, help="Number of Processors.")
    parser.add_option("--directory",
                      action="store", dest="datadir", help="working directory location")
    parser.add_option("--runmode", action="store", dest="mode", type="string", help="decide to which mode (run or  compare)")
    parser.add_option("--nsweep", action="store", dest="nsweep", type=int, help="number of sweep")
    parser.add_option("--mesh", action="store", dest="mesh", type="string", help="fluent case file")



    (options, args) = parser.parse_args()
    cwd = os.getcwd()
    options.scriptname = cwd + "/" + options.scriptname
    options.datadir    = cwd + "/" + options.datadir
    print options.datadir
    print options.scriptname

    
    mpi = "mpirun -np %s  " % options.np + " python " + options.scriptname + " " + options.mesh + " " + str(options.nsweep)

    if options.datadir:
        if not os.path.isdir(options.datadir):
            try:
                os.makedirs(options.datadir)
            except:
                fatal("error creating directory " + options.datadir)
        os.chdir(options.datadir)

    
    if os.system( mpi ):
        cleanup(-1)
    

    # OK.  Now check results.


    err = check_particles( options )
    cleanup(err)

 


if __name__ == "__main__":
    main()


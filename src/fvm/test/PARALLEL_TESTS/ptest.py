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

def check_mesh(options):
    
    check_file = []
    reference_file = []
    for np in range(0,options.np):
       check_file.append ( str(options.datadir) + "/mesh_proc" + str(np) + ".dat" )
       reference_file.append( str(options.datadir) +"/GOLDEN/" + "mesh_proc" + str(np) + ".dat" )
       if  not(filecmp.cmp( check_file[np], reference_file[np], 0)) :
         print "mesh files didn't match GOLDEN/* !!!!!!!!!!!!!"
         return  -1
         break;
    print "mesh files are ok"
    return 0

def check_parmetis(options):
    
    check_file = []
    reference_file = []
    for np in range(0,options.np):
       check_file.append ( str(options.datadir) + "/proc" + str(np) +"_debug_print" + ".dat" )
       reference_file.append( str(options.datadir) +"/GOLDEN/" + "proc" + str(np) +"_debug_print" + ".dat" )
       print check_file[np]
       print reference_file[np]
       if  not(filecmp.cmp( check_file[np], reference_file[np], 0)) :
         print "parmetis files didn't match GOLDEN/* !!!!!!!!!!!!!"
         return  -1
         break;
    print "parmetis files are ok"
    return 0

def check_mapping(options):
    
    check_file = []
    reference_file = []
    for np in range(0,options.np):
       check_file.append ( str(options.datadir) + "/mesh_proc" + str(np) +"_info" + ".dat" )
       reference_file.append( str(options.datadir) +"/GOLDEN/" + "mesh_proc" + str(np) + "_info" + ".dat" )
       print check_file[np]
       print reference_file[np]
       if  not(filecmp.cmp( check_file[np], reference_file[np], 0)) :
         print "mapping files didn't match GOLDEN/* !!!!!!!!!!!!!"
         return  -1
         break;
    print "mapping files are ok"
    return 0

def check_thermal_solver(options):
    
    check_file = []
    reference_file = []
    for np in range(0,options.np):
       check_file.append ( str(options.datadir) + "/temp_proc" + str(np) + ".dat" )
       reference_file.append( str(options.datadir) +"/GOLDEN/" + "temp_proc" + str(np) + ".dat" )
       if  not(filecmp.cmp( check_file[np], reference_file[np], 0)) :
         print "temp files didn't match GOLDEN/* !!!!!!!!!!!!!"
         return  -1
         break;
    print "temp files are ok"
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
    parser.add_option("--test", action="store", dest="testType", type="string", help="test function to be called")


    (options, args) = parser.parse_args()
    cwd = os.getcwd()
    options.scriptname = cwd + "/" + options.scriptname
    options.datadir    = cwd + "/" + options.datadir
    print options.datadir
    print options.scriptname

    
    mpi = "mpirun -np %s  " % options.np + " python " + options.scriptname

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
    print options.testType
    if options.testType == "parmetis" :
       err = check_parmetis( options )
       cleanup(err)

    if options.testType == "mesh" :
       err = check_mesh( options )
       cleanup(err)

    if options.testType == "mapping" :
       err = check_mapping( options )
       cleanup(err)
 
    if options.testType == "thermal_solver":
       err = check_thermal_solver( options )
       cleanup(err)



if __name__ == "__main__":
    main()


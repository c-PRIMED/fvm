#!/usr/bin/env python

"""
Given a directory containing MPM inputs and expected
outputs, run MPM and compare the outputs.

Usage: mpmtest.py [options] input_dirname
Options:
  --datafile     Directory where test data is stored. Default is current dir.
  --np           Number of processors. Default in 1. 
"""

import sys, os
from optparse import OptionParser

def usage():
    print __doc__
    sys.exit(-1)

def cleanup(ec, datafile, dl):
    sys.exit(ec)

def main():
    parser = OptionParser()
    parser.set_defaults(np=1,datafile=os.getcwd())
    parser.add_option("--datafile",
                      action="store", dest="datafile", help="Data directory name.")
    parser.add_option("--np",
                      action="store", dest="np", type=int, help="Number of Processors.")
    (options, args) = parser.parse_args()

    if len(args) != 1:
        usage()

    dirname = args[0]

    # run PM2 and send stdout and stderr to a log file.
    cmd = ' '.join(args) + " >" + datafile + " 2>&1"
    print "cmd=",cmd
    if os.system("/bin/bash -c '%s'" % cmd):
        cleanup(-1, datafile, delete_datafile)
        count = 0
        exe = os.popen("%s %s %s" % (diff, datafile, options.golden))
        for line in exe:
            if count < 100:
                print line.rstrip()
                count += 1
            else:
                break
        if count == 100: print "[truncated]"
        if exe.close():
            ret = 1




#        cmd = ' '.join(args)
#        ret = os.system(cmd)
#        if ret:
#            ret >>= 8

    cleanup(ret, datafile, delete_datafile)

if __name__ == "__main__":
    main()

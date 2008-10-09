#!/usr/bin/env python

"""
Run an executable with the supplied args and
compare the output with the expected output.
Expected output can be a string or a golden file, 
but not both.

Usage: mtest.py [options] exe [args]
Options:
  --expect       Expected string.
  --golden       Golden file
  --compare      Comparison executable. Default is "diff".
The options below are used by the test framework.
  --datafile     File where test data is stored.
"""

import sys, os, tempfile, re
from optparse import OptionParser

def usage():
    print __doc__
    sys.exit(-1)

def cleanup(ec, datafile, dl):
    if dl:
        os.remove(datafile)
    sys.exit(ec)

def main():
    parser = OptionParser()
    parser.add_option("--expect",
                      action="store", dest="expect", help="Expected output")
    parser.add_option("--golden",
                      action="store", dest="golden", help="Golden filename")
    parser.add_option("--datafile",
                      action="store", dest="datafile", help="Data filename")
    (options, args) = parser.parse_args()

    if options.golden and options.expect:
        print "Cannot have both --golden and --expect."
        usage()

    diff = 'diff'
    delete_datafile = 0
    datafile = options.datafile
    if not datafile:
        datafile = tempfile.mkstemp()[1]
        delete_datafile = 1


    ret = 0
    if options.golden:
        # run the test executable and send stdout and stderr to the data file.
        cmd = ' '.join(args) + " >>" + datafile + " 2>&1"
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
    elif options.expect:
        exp = re.compile(options.expect)
        exe = os.popen(' '.join(args))
        output = exe.read()
        match = exp.match(output)
        if not match:
            print "Expected '%s'" % options.expect
            print "Got '%s'" % output
            ret = 1
        if exe.close():
            ret = -1

    cleanup(ret, datafile, delete_datafile)

if __name__ == "__main__":
    main()

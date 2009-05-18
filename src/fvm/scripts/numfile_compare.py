#!/usr/bin/env python

"""
Very simple comparison of two files containing numbers.
Usage: numfile_compare.py file1 file2
"""

import sys, os

def usage():
    print __doc__
    sys.exit(-1)

# return 0 if two numbers are equal
def compare_nums(a,b,maxerr):
    global MaxErr
    if a == b:
        return 0
    try:
        fa = float(a)
        fb = float(b)
        if abs(fa - fb) <= maxerr:
            return 0
        if abs(fa-fb) > MaxErr:
            MaxErr = abs(fa-fb)
    except:
        pass
    return 1

def check_files(file1, file2, maxerr):
    try:
        f = open(file1, 'r')
    except:
        print "Cannot open '%s'" % file1
        return 1

    try:
        g = open(file2, 'r')
    except:
        print "Cannot open '%s'" % file2
        return 1

    lineno = 0
    while True:
        line = f.readline()
        gline = g.readline()
        lineno += 1
        if not line:
            if gline:
                print gline
                return 1
            break
        if gline == line:
            continue
        for a,b in zip(line.split(),gline.split()):
            res = compare_nums(a,b,maxerr)
            if res:
                print "Comparison of %s and %s failed on line %s" % (file1, file2, lineno)
                print line, gline
                return 1
    return 0

def main():
    global MaxErr
    if len(sys.argv) != 3:
        usage()
        
    MaxErr = 0.0
    if check_files(sys.argv[1], sys.argv[2], 1.e-5):
        print "Maximum Error =",MaxErr
        sys.exit(1)
    sys.exit(0)

if __name__ == "__main__":
    main()

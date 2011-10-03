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

def check_mesh(options):
    for np in range(0,options.np):
       check_file = os.path.join(options.outdir, 'mesh_proc%s.dat' % np)
       reference_file = os.path.join(options.golden, 'mesh_proc%s.dat' % np)
       if numfile_compare.check_files(check_file, reference_file, 1.0e-5):
         print "mesh files didn't match GOLDEN/* !!!!!!!!!!!!!"
         return  -1
         break;
    print "mesh files are ok"
    return 0

def check_parmetis(options):
    for np in range(0,options.np):
       check_file = os.path.join(options.outdir, 'proc%s_debug_print.dat' % np)
       reference_file = os.path.join(options.golden, 'proc%s_debug_print.dat' % np)
       if numfile_compare.check_files(check_file, reference_file, 1.0e-5):
         print "parmetis files didn't match GOLDEN/* !!!!!!!!!!!!!"
         return  -1
    print "parmetis files are ok"
    return 0

def check_mapping(options):
    for np in range(0,options.np):
       check_file = os.path.join(options.outdir, 'mesh_proc%s_info.dat' % np)
       reference_file = os.path.join(options.golden, 'mesh_proc%s_info.dat' % np)
       if numfile_compare.check_files(check_file, reference_file, 1.0e-5):
         print "mapping files didn't match GOLDEN/* !!!!!!!!!!!!!"
         return  -1
         break;
    print "mapping files are ok"
    return 0

def check_thermal_solver(options):
    for np in range(0,options.np):
       check_file = os.path.join(options.outdir, 'temp_proc%s.dat' % np)
       reference_file = os.path.join(options.golden, 'temp_proc%s.dat' % np)
       if numfile_compare.check_files(check_file, reference_file, 1.0e-5):
         print "temp files didn't match GOLDEN/* !!!!!!!!!!!!!"
         return  -1
         break;
    print "temp files are ok"
    return 0

def check_storage_site_merger(options):
    for np in range(0,options.np):
       check_file = os.path.join(options.outdir, 'proc%s_storage_site_merger.dat' % np)
       reference_file = os.path.join(options.golden, 'proc%s_storage_site_merger.dat' % np)
       if numfile_compare.check_files(check_file, reference_file, 1.0e-5):
         print "site merger files didn't match GOLDEN/* !!!!!!!!!!!!!"
         return  -1
         break;
    print "site merger files are ok"
    return 0
    
    
def check_convergence(options):
    check_file     = os.path.join(options.outdir, 'convergence.dat')
    os.system("more "+check_file)
    reference_file = os.path.join(options.golden, 'convergence.dat')
    if numfile_compare.check_files(check_file, reference_file, 1.0e-8):
       print "convergence history files didn't match GOLDEN/* !!!!!!!!!!!!!"
       return  -1
      # break;
    print "convergence history files are ok"
    return 0
  

def main():
    funcs = {
        'testPartMesh.py' : [check_parmetis, check_mesh, check_mapping],
        'testThermalParallel.py' : [check_thermal_solver], 
        'testMerger.py'                  : [check_storage_site_merger],
        'testThermalParallelJacobi.py':[check_convergence],
	'beamTest.py':[check_convergence],
	'beamTest3D.py':[check_convergence],
	'testFlowParallel.py':[check_convergence],
	'testPlate.py':[check_convergence],
	'testPlateTransient.py':[check_convergence],
	'testStructureModelTransient.py':[check_convergence]
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
    print options.infile
    # convert path to absolute because we may cd to the test directory
    options.script   = os.path.abspath(options.script)
    options.outdir   = os.path.abspath(options.outdir)
    options.golden   = os.path.abspath(os.path.join(options.golden,'GOLDEN'))
    options.infile   = os.path.abspath(options.infile)

    if options.outdir:
        if not os.path.isdir(options.outdir):
            try:
                os.makedirs(options.outdir)
            except:
                print "error creating directory " + options.outdir
                sys.exit(1)
        os.chdir(options.outdir)
    print "cwd = ", os.getcwd()
    print options.infile
    mpi = 'mpirun -np %s %s %s --type=%s' % (options.np, options.script, options.infile, options.type)
    print "MPI = ", mpi
    if options.xdmf:
        mpi += ' --xdmf'

    os.system(mpi)
    if os.system( mpi ):
        cleanup(-1)

    for f in funcs[os.path.basename(options.script)]:
        err = f(options)
        if err:
            cleanup(err)

if __name__ == "__main__":
    main()


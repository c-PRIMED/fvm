#! /usr/bin/env python
# encoding: utf-8
#
# Martin Hunt <mmh@purdue.edu>

"""
Functions for dealing with MOAB
"""

import os, config, time, subprocess

def start(bp, cname):
    if not config.config('Testing', 'MOAB'):
        return 0
    
    in_batch = 0
    try:
        if os.environ["SLURM_JOBID"]:
            in_batch = 1
    except:
        pass

    if in_batch:
        return 0

    msub(bp, cname)
    return 1


# write moab batch file and submit it
def msub(bp, cname):

    # first remove any old results
    os.system('/bin/rm -f run_tests.mob.*')
    
    f = open('run_tests.mob', 'w')
    f.write('#!/bin/bash -l\n')

    qname = config.config('Testing', 'queue')
    f.write('#MSUB -q %s\n' % qname)
    nodes = config.config('Testing', 'nodes')
    f.write('#MSUB -l nodes=%s\n' % nodes)
    f.write('#MSUB -l walltime=%s\n' % config.config('Testing', 'walltime'))
    f.write('#MSUB -l partition=%s\n' % config.config('Testing', 'partition'))

    f.write('source env-%s.sh\n' % cname)
    
    # Load Modules
    cmds = config.config('ALL', 'before')
    for c in cmds:
        f.write('%s\n' % c)            

    make = os.path.join(bp.topdir,'make.py')
    f.write('%s --test %s\n' % (make, cname))
    f.close()

    job = os.popen("msub run_tests.mob").read()
    job = job.split()[0]
    
    print "Job %s submitted to queue %s." % (job, qname)

    state = ''
    while True:
        try:
            time.sleep(15)
        except KeyboardInterrupt:
            print "\nInterrupted. Job will continue in the background."
            print "It you want to cancel it, do \"qdel %s\"\n" % job
            return

        so = file('/dev/null')
        st = subprocess.Popen('checkjob -A %s' % job, shell=True, stdout=subprocess.PIPE, stderr=so).stdout
        so.close()
        str = st.read().rstrip().rstrip(';')
        d = dict(x.split('=') for x in str.split(';'))
        newstate = d['STATE']
        
        try:
            if newstate != state:
                if newstate == 'Idle':
                    print "Waiting for job to start."
                elif newstate == 'Running':
                    print 'Job is running.'
                elif newstate == 'Completed':
                    print 'Job is finished.'
                    break
                else:
                    print "Job state is '%s'." % newstate
                state = newstate
        except:
            pass
        st.close()


    for filename in ['run_tests.mob.e%s' % job, 'run_tests.mob.o%s' % job]:
        try:
            st = os.stat(filename)
            if st.st_size:
                os.system('cat %s' % filename)
        except:
            pass

        

    

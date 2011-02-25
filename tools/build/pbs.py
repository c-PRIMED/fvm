#! /usr/bin/env python
# encoding: utf-8
#
# Martin Hunt <mmh@purdue.edu>

"""
Functions for dealing with PBS
"""

import re, os, config, time
import subprocess

def start(bp, cname):
    if not config.config('Testing', 'PBS'):
        return 0

    in_batch = 0
    try:
        if os.environ["PBS_ENVIRONMENT"] == 'PBS_BATCH':
            in_batch = 1
    except:
        pass

    if in_batch:
        return 0

    qsub(bp, cname)
    return 1


# write pbs batch file and submit it
def qsub(bp, cname):

    # first remove any old results
    os.system('/bin/rm -f run_tests.pbs.*')

    f = open('run_tests.pbs', 'w')
    f.write('#!/bin/bash\n')

    qname = config.config('Testing', 'queue')
    f.write('#PBS -q %s\n' % qname)

    cpus = config.config('Testing', 'cpus')
    nodes = config.config('Testing', 'nodes')
    f.write('#PBS -l select=%s:ncpus=%s:mpiprocs=%s\n' % (nodes,cpus,cpus))
    f.write('#PBS -l walltime=%s\n' % config.config('Testing', 'walltime'))

    f.write('cd  $PBS_O_WORKDIR\n')
    f.write('source env.sh\n')

    # Load Modules
    cmds = config.config('ALL', 'before')
    for c in cmds:
        f.write('%s\n' % c)

    make = os.path.join(bp.topdir,'make.py')
    f.write('%s --test %s\n' % (make, cname))
    f.close()

    job = os.popen("qsub run_tests.pbs").readline()
    job = job.split('.')[0]

    print "Job %s submitted to queue %s." % (job, qname)
    cmd = "qstat -f %s" % job
    print "cmd=%s" % cmd
    state = ''
    while True:
        try:
            time.sleep(15)
        except KeyboardInterrupt:
            print "\nInterrupted. Job will continue in the background."
            print "It you want to cancel it, do \"qdel %s\"\n" % job
            return

        st = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        saved_line = ''
        savelist = ['comment', # Gives time run and host name
                    'job_state', # 'F', 'Q' , 'R' , 'X' for finished with errors
                    'Exit_status', # 0 = normal, 271 for resources exceeded
                    # resources used
                    'resources_used.mem',  # string '280864kb'
                    'resources_used.vmem', # string '280864kb'
                    'resources_used.walltime', # string HH:MM:SS
                    # resources requested
                    'queue', # queue name
                    'Resource_List.walltime', # string HH:MM:SS
                    'Submit_arguments', # PBS file name. Useful for resubmitting
                    ]

        errline = st.stderr.read()
        if errline.startswith('qstat: Unknown Job Id'):
            print "Job no longer running"
            return

        d = {}
        for line in st.stdout:
            if line.startswith('Job Id:'):
                continue

            if line.startswith('\t'):
                saved_line += line[1:-1]
                continue

            if saved_line:
                k,v = re.findall(r'(\S+)\s*=\s*(.*)', saved_line)[0]
                if k in savelist:
                    d[k] = v

            if line.startswith('    '):
                line = line[4:-1]

            saved_line = line

        try:
            if d and d['Exit_status'] != '0':
                d['job_state'] = 'X'
        except:
            pass

        try:
            if d['job_state'] != state:
                state = d['job_state']
                print "State: %s" % state
            if d['job_state'] == 'F' or d['job_state'] == 'X':
                break
        except KeyError:
            pass

    for filename in ['run_tests.pbs.e%s' % job, 'run_tests.pbs.o%s' % job]:
        try:
            st = os.stat(filename)
            if st.st_size:
                os.system('cat %s' % filename)
        except:
            pass





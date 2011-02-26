#! /usr/bin/env python
# encoding: utf-8
#
# Memosa testing script
# Martin Hunt <mmh@purdue.edu>

"""
Functions for running tests
"""

import re, shutil
from build_packages import *
from build_utils import *
from subprocess import Popen, PIPE, STDOUT
from xml.sax.saxutils import escape
import xml.dom.minidom

class JTest:
    def __init__(self):
        self.ok = 0
        self.errs = 0

    def print_results(self):
        print '='*40
        print "PASSED  %d" % (self.ok)
        print "FAILED  %d" % (self.errs)

    def test_exec(self, cmd):
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT)
        res = p.stdout.read()
        rc = p.wait()
        return rc, res

    def fix_nose_xml(self, fname, pname):
        dom = xml.dom.minidom.parse(fname)
        ts = dom.getElementsByTagName("testsuite")[0]
        cases = dom.getElementsByTagName("testcase")
        for case in cases:
            cls = case.getAttribute("classname")
            case.setAttribute("classname", "%s.%s" % (pname,cls))
        fp = open(fname, 'w')
        fp.write(dom.toxml())
        fp.close

    def get_dom(self, pname, logdir):
        xmlname = os.path.join(logdir, "%s-test.xml" % pname)
        try:
            dom = xml.dom.minidom.parse(xmlname)
        except IOError:
            dom = xml.dom.minidom.parseString('<testsuite />')
        return dom

    #run all the tests in a file. return number of errors
    def run_tests(self, pname, fname, logdir):
        if not os.path.isdir(logdir):
            try:
                os.makedirs(logdir)
            except:
                fatal("error creating directory " + logdir)
        errs = ok = 0
        for line in open(fname):
            if line[0] == '\n' or line[0] == '#': continue
            tname = line.split()[0]
            pdir=''
            if line.find('TESTDIR') >= 0 or line.find('TESTOUT') >= 0:
                pdir = os.path.join(logdir, tname)
                if not os.path.isdir(pdir):
                    try:
                        os.makedirs(pdir)
                    except:
                        fatal("error creating directory " + pdir)
            line = line.replace('TESTDIR', pdir)
            line = line.replace('TESTOUT', os.path.join(pdir, tname + '.dat'))
            line = line.split()

            # force nosetests --with-xunit
            if line[1] == 'nosetests':
                nose = True
                xfile = os.path.join(logdir, '%s-test.xml' % tname)
                line[1] = 'nosetests --nologcapture --with-xunit --xunit-file=%s' % xfile
            else:
                nose = False
            cmd = ' '.join(line[1:])
            verbose(1, "Test %s: %s" % (tname, cmd))
            t = time.time()
            err, result_text = self.test_exec(cmd)
            t = time.time() - t
            if err:
                verbose_warn("%s Failed" % tname)
                errs += 1
                self.errs += 1
            else:
                ok += 1
                self.ok += 1
            if nose:
                self.fix_nose_xml(xfile, pname)
            else:
                # add testcase to the dom
                ts = self.dom.getElementsByTagName("testsuite")[0]
                tc = self.dom.createElement('testcase')
                tc.setAttribute('classname', pname)
                tc.setAttribute('name', tname)
                tc.setAttribute('time', str(t))
                if err:
                    f = self.dom.createElement('failure')
                    f.setAttribute('type', 'unknown')
                    f.setAttribute('message', result_text)
                    tc.appendChild(f)
                ts.appendChild(tc)
        return ok, errs


# generator that finds all the test files under a directory.
def find_tests(startdir):
    for root, dirs, files in os.walk(startdir):
        for x in files:
            if x == 'TESTS':
                yield root

# run all the tests in a directory.  return number of errors
def do_tests(tst, pname, startdir, logdir):
    shutil.rmtree(logdir, True)
    tst.dom = tst.get_dom(pname, logdir)
    errs = ok = 0
    for root in find_tests(startdir):
        os.chdir(root)
        _ok, _errs = tst.run_tests(pname, os.path.join(root, 'TESTS'), logdir)
        ok += _ok
        errs += _errs
    if ok or errs:
        ts = tst.dom.getElementsByTagName("testsuite")[0]
        ts.setAttribute('tests', str(ok+errs))
        ts.setAttribute('errors', '0')
        ts.setAttribute('failures', str(errs))
        ts.setAttribute('skip', '0')
        xmlname = os.path.join(logdir, "%s-test.xml" % pname)
        fp = open(xmlname, 'w')
        fp.write(tst.dom.toxml())
        fp.close
    return ok, errs

def run_all_tests(bld):
    testdir = os.path.join(bld.logdir, 'testing')
    try:
        os.mkdir(testdir)
    except OSError:
        pass

    # run any before commands
    run_commands('Testing', 'before')

    # add some dirs to our path
    tooldir = os.path.join(bld.topdir, "tools", "test")
    fix_path('PATH', tooldir, 1, 0)
    fix_path('PATH', bld.bindir, 1, 0)
    fix_path('LD_LIBRARY_PATH', bld.libdir, 1, 0)
    try:
        oldpypath = os.environ['PYTHONPATH']
    except:
        oldpypath = ''
    set_python_path(bld.blddir)

    # Create Test Object to hold and display results
    tst = JTest()

    # test each package
    for p in bld.packages:
        if not p.is_pkg:
            p.test(tst)

    tst.print_results()

    # restore paths
    if oldpypath:
        os.environ['PYTHONPATH'] = oldpypath
    else:
        del os.environ['PYTHONPATH']
    fix_path('LD_LIBRARY_PATH', bld.libdir, 1, 1)
    fix_path('PATH', bld.bindir, 1, 1)
    fix_path('PATH', tooldir, 1, 1)

if __name__ == "__main__":
    for testfile in find_tests(sys.argv[1]):
        print "testfile=", testfile


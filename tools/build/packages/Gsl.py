from build_packages import *
import cgi

# Not used anymore. Keeping as an example of another way to do tests.
class Gsl(BuildPkg):
    def _configure(self):
        return self.sys_log("%s/configure --prefix=%s" % (self.sdir, self.blddir))
    def _build(self):
        return self.sys_log("make -j%s" % jobs(self.name))
    def _install(self):
        return self.sys_log("make install")
    def _test(self):
        ok = errs = 0
        os.chdir(self.bdir)
        logfile = self.logfile.replace('xml','txt')
        os.system("make check > %s 2>&1" % logfile)
        for line in open(logfile):
            if line.find('PASS') == 0:
                ok += 1
            elif line.find('FAIL') == 0:
                errs += 1
        if errs:
            ostr = "<Test Status=\"failed\">\n"
        else:
            ostr = "<Test Status=\"passed\">\n"
        ostr += "\t<Name>gsl</Name>\n"
        ostr += "\t<Path>%s</Path>\n" % self.sdir
        ostr += "\t<FullName>gsl</FullName>\n"
        ostr += "\t<FullCommandLine>make check</FullCommandLine>\n"
        ostr += "\t<Results><Measurement><Value>"
        ostr += cgi.escape(open(logfile).read())
        ostr += "</Value></Measurement></Results></Test>\n"
        f = open(self.logfile,'w')
        f.write(ostr)
        f.close()
        return ok, errs

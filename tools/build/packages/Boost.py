from build_packages import *

class Boost(BuildPkg):
    def _installed(self):
        oname = 'boost_version'
        ename = "./" + oname
        cname = oname + '.cpp'
        src="#include <boost/version.hpp>\n#include <iostream>\nint main()\n{ std::cout << BOOST_LIB_VERSION; }"
        fp = open(cname, 'w')
        fp.write(src)
        fp.close()
        res = os.system('g++ -o %s %s' % (oname, cname))
        if res:
            try:
                os.remove(cname)            
            except:
                pass
            return False
        ver = subprocess.Popen(ename, shell=True, stdout=subprocess.PIPE).stdout.read()
        try:
            os.remove(cname)            
            os.remove(ename)
        except:
            pass
        a,b = re.findall(r'([^_]*)_([^\n]*)',ver)[0]
        #print "ver=%s.%s" % (a,b)        
        if int(a) > 1 or int(b) > 45:
            return True
        return False
    
    def _configure(self):
        return self.sys_log("./bootstrap.sh --prefix=%s" % (self.blddir))
    def _install(self):
        ret = self.sys_log("./b2 install")
        return ret

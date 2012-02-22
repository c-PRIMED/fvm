from build_packages import *
import build_utils

class Boost(BuildPkg):
    def _installed(self):
        if os.environ.has_key('BOOST_HOME'):
            binc = os.path.join(os.environ['BOOST_HOME'], 'include')
            blib = os.path.join(os.environ['BOOST_HOME'], 'lib')
            build_utils.fix_path('CPLUS_INCLUDE_PATH', binc, 0, 0)
            build_utils.fix_path('C_INCLUDE_PATH', binc, 0, 0)            
            build_utils.fix_path('LD_LIBRARY_PATH', blib, 0, 0)
            
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
        #print "ver=%s" % ver        
        a,b = re.findall(r'([^_]*)_([^_\n]*)',ver)[0]
        #print "ver=%s.%s" % (a,b)  
        # PCL needs 1.45 but nothing else does      
        if int(a) > 1 or int(b) > 40:
            return True
        return False
    
    def _configure(self):
        return self.sys_log("./bootstrap.sh --prefix=%s" % (self.blddir))
    def _build(self):
        ret = self.sys_log("./b2 install")
        return ret

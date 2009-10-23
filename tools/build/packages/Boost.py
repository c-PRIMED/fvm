from build_packages import *

# We don't build boost, just use the headers
class Boost(BuildPkg):
    def unpack_srcs(self):
        dst = os.path.join(self.blddir, 'include', 'boost')
        src = self.sdir
        suffix = src.split('.')[-1]
                
        compress = None
        if suffix == 'tgz' or suffix == 'gz':
            compress = 'z'
        elif suffix == 'bz2':
            compress = 'j'
        elif suffix == 'tar':
            compress = ''

        if compress == None:
            fatal('Cannot unpack %s' % src)

        # need to keep the build system happy
        self.sys_log("/bin/mkdir -p %s" % self.bdir)        
        
        if not os.path.isdir(dst):
            self.sys_log("/bin/mkdir -p %s" % dst)
        
        os.system('tar -C %s -%sxf %s' % (dst, compress, src))
        



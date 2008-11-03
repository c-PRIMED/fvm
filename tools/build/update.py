"""
Do subversion update and write out xml file.
"""

import cgi, time, datetime, os, socket, config
import sys, httplib, urlparse, build_utils

def svn_up(f):
   
   pass

def update(bp, cname, nightly):
   if nightly:
      bname = "NIGHTLY"
   else:
      bname = "EXPERIMENTAL"
      
   fname = os.path.join(bp.logdir, "Update.xml")
   f = open(fname, 'w')
   f.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
   f.write("<Update mode=\"Client\" Generator=\"make.py\">\n")
   
   Site = socket.gethostname()
   BuildName = cname
   BuildStamp = "%s-%s-%s" % (time.strftime('%Y%m%d'), os.getpid(), bname)
   StartDateTime = time.strftime("%b %d %H:%M %Z", time.localtime())
   StartTime = time.time()
   UpdateCommand = "svn status"
   UpdateType = "SVN"
   
   for v in ("Site", "BuildName", "BuildStamp", "StartDateTime", "StartTime", \
                "UpdateCommand", "UpdateType"): 
      f.write("\t<%s>%s</%s>\n" % (v, eval(v), v))
      
   #svn_up(f)
   f.write("</Update>\n")
   f.close()

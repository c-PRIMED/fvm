"""
Do subversion update and write out xml file.
"""

import cgi, time, datetime, os, socket

Dirname=None

def write_update(f, fname='', status=''):
   global Dirname
   dname = os.path.dirname(fname)
   bname = os.path.basename(fname)
   if fname == '' and Dirname:
      f.write("</Directory>\n")
      return
   if Dirname and Dirname != dname:
      f.write("</Directory>\n")
      Dirname = None
   if not Dirname:
      Dirname = dname
      f.write("<Directory>\n\t<Name>%s</Name>\n" % cgi.escape(dname))
        
   print "%s: %s" % (fname, status)
   f.write("<%s>\n\t<File Directory=\"%s\">%s</File>\n" % (status, dname, bname))

   Author = CheckinDate = Revision = PriorRevision = 'unknown'
   if status == 'Deleted':
      Log = "Deleted file"
   else:
      Log = ''
      # Get more info on the file
      for line in os.popen("svn info %s" % fname):
         if line.startswith('Revision:'):
            Revision = line.split(':')[1].strip()
         elif line.startswith('Last Changed Author:'):
            Author = line.split(':')[1].strip()
         elif line.startswith('Last Changed Rev:'):
            PriorRevision = line.split(':')[1].strip()
         elif line.startswith('Last Changed Date:'):
            CheckinDate = line[len('Last Chanegd Date:'):].strip()
            
      if status == "Modified":
         Log = "Locally modified file"
      else:
         for line in os.popen("svn log --xml -r %s" % Revision):
            start = line.find("<msg>")
            if start >= 0:
               start += len("<msg>")
               stop = line.find("</msg>")
               if stop >= 0:
                  Log = line[start:stop]
                  break

   Directory = dname
   FullName = fname

   for v in ("Directory", "FullName", "CheckinDate","Author","Log","Revision","PriorRevision"):
      f.write("\t<%s>%s</%s>\n" % (v, eval(v), v))
   f.write("</%s>\n" % status)

def svn_up(f):

   exe = os.popen("svn update")
   for line in exe:
      if (line[0] == 'U' or line[0] == 'A') and line[1] == ' ':
         write_update(f, line.split()[1], "Updated")
      elif line[0] == 'D' and line[1] == ' ':
         write_update(f, line.split()[1], "Deleted")
   ret = exe.close()
   if ret:
      debug("\"svn update\" Failed. Return code %s" % ret >> 8)

   exe = os.popen("svn status")
   for line in exe:
      if line[0] == 'M':
         write_update(f, line.split()[1], "Modified")
      elif line[0] == 'C':
         write_update(f, line.split()[1], "Conflicting")
   ret = exe.close()
   if ret:
      debug("\"svn status\" Failed. Return code %s" % ret >> 8)
   write_update(f)

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
      
   svn_up(f)
   f.write("</Update>\n")
   f.close()

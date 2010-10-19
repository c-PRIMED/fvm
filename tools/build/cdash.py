"""
Functions to submit build and test information to CDash
"""

import time, os, socket, re
import sys, httplib, urlparse, build_utils
from xml.sax.saxutils import escape

try:
    set
except NameError:
    # my python is old
    from sets import Set as set

def barf(msg):
    print >> sys.stderr, "Error! %s" % msg
    sys.exit(1)


def parseuri(uri):
    """Parse URI, return (host, port, path) tuple.

    >>> parseuri('http://example.org/testing?somequery#frag')
    ('example.org', 80, '/testing?somequery')
    >>> parseuri('http://example.net:8080/test.html')
    ('example.net', 8080, '/test.html')
    """

    scheme, netplace, path, query, fragid = urlparse.urlsplit(uri)

    if ':' in netplace:
        host, port = netplace.split(':', 2)
        port = int(port)
    else: host, port = netplace, 80

    if query: path += '?' + query

    return host, port, path

def putfile(f, uri, username=None, password=None):
    """HTTP PUT the file f to uri, with optional auth data."""
    host, port, path = parseuri(uri)

    redirect = set([301, 302, 307])
    authenticate = set([401])
    okay = set([200, 201, 204])

    authorized = False
    authorization = None
    tries = 0

    while True:
        # Attempt to HTTP PUT the data
        h = httplib.HTTPConnection(host, port)

        h.putrequest('PUT', path)

        h.putheader('User-Agent', 'put.py/1.0')
        h.putheader('Connection', 'keep-alive')
        h.putheader('Transfer-Encoding', 'chunked')
        h.putheader('Expect', '100-continue')
        h.putheader('Accept', '*/*')
        if authorization:
            h.putheader('Authorization', authorization)
        h.endheaders()

        # Chunked transfer encoding
        # Cf. 'All HTTP/1.1 applications MUST be able to receive and
        # decode the "chunked" transfer-coding'
        # - http://www.w3.org/Protocols/rfc2616/rfc2616-sec3.html
        while True:
            bytes = f.read(2048)
            if not bytes: break
            length = len(bytes)
            h.send('%X\r\n' % length)
            h.send(bytes + '\r\n')
        h.send('0\r\n\r\n')
        f.seek(0);

        resp = h.getresponse()
        status = resp.status # an int

        # Got a response, now decide how to act upon it
        if status in redirect:
            location = resp.getheader('Location')
            uri = urlparse.urljoin(uri, location)
            host, port, path = parseuri(uri)

            # We may have to authenticate again
            if authorization:
                authorization = None

        elif status in authenticate:
            # If we've done this already, break
            if authorization:
                # barf("Going around in authentication circles")
                barf("Authentication failed")

            if not (username and password):
                barf("Need a username and password to authenticate with")

            # Get the scheme: Basic or Digest?
            wwwauth = resp.msg['www-authenticate'] # We may need this again
            wauth = wwwauth.lstrip(' \t') # Hence use wauth not wwwauth here
            wauth = wwwauth.replace('\t', ' ')
            i = wauth.index(' ')
            scheme = wauth[:i].lower()

            if scheme in set(['basic', 'digest', 'Basic', 'Digest']):
                msg = "Performing %s Authentication..." % scheme.capitalize()
            else: 
                barf("Unknown authentication scheme: %s" % scheme)

            if scheme == 'basic' or scheme == 'Basic':
                import base64
                userpass = username + ':' + password
                userpass = base64.encodestring(userpass).strip()
                authorized, authorization = True, 'Basic ' + userpass

            elif scheme == 'digest':
                msg = "uses fragile, undocumented features in urllib2"
                build_utils.debug("Warning! Digest Auth %s" % msg)

                import urllib2 # See warning above

                passwd = type('Password', (object,), 
                              {
                               'find_user_password': lambda self, *args: (username, password),
                               'add_password': lambda self, *args: None
                               })()

                req = type('Request', (object,), 
                           {
                            'get_full_url': lambda self: uri,
                            'has_data': lambda self: None,
                            'get_method': lambda self: 'PUT',
                            'get_selector': lambda self: path
                            })()

                # Cf. urllib2.AbstractDigestAuthHandler.retry_http_digest_auth
                auth = urllib2.AbstractDigestAuthHandler(passwd)
                token, challenge = wwwauth.split(' ', 1)
                chal = urllib2.parse_keqv_list(urllib2.parse_http_list(challenge))
                userpass = auth.get_authorization(req, chal)
                authorized, authorization = True, 'Digest ' + userpass

        elif status in okay:
            if (username and password) and (not authorized):
                msg = "Warning! The supplied username and password went unused"
                print msg

            if build_utils.opt_debug:
                resultLine = "Success! Resource %s"
                statuses = {200: 'modified', 201: 'created', 204: 'modified'}
                build_utils.debug(resultLine % statuses[status])

                statusLine = "Response-Status: %s %s"
                build_utils.debug(statusLine % (status, resp.reason))

                body = resp.read(58)
                body = body.rstrip('\r\n')
                body = body.encode('string_escape')

                if len(body) >= 58:
                    body = body[:57] + '[...]'

                bodyLine = 'Response-Body: "%s"'
                build_utils.debug(bodyLine % body)
            break

        # @@ raise PutError, do the catching in main?
        else: 
            barf('Got "%s %s"' % (status, resp.reason))

        tries += 1
        if tries >= 50:
            barf("Too many redirects")

    return status, resp

def putname(fn, uri, username=None, password=None):
    """HTTP PUT the file with filename fn to uri, with optional auth data."""
    auth = {'username': username, 'password': password}

    if fn != '-':
        f = open(fn, 'rb')
        status, resp = putfile(f, uri, **auth)
        f.close()
    else: 
        status, resp = putfile(sys.stdin, uri, **auth)

    return status, resp

def site_str(BuildName, buildtype):
    ss = ""
    BuildStamp = "%s-%s-%s" % \
        (time.strftime('%Y%m%d'), os.getpid(), buildtype)
    (OSName, Name, OSRelease, OSVersion, OSPlatform) = os.uname()
    Hostname = socket.gethostname()
    Generator = "make.py"
    Is64Bits = 0
    if OSPlatform == 'x86_64': Is64Bits = 1
    NumberOfPhysicalCPU = 0

    for line in open('/proc/cpuinfo', 'r'):
        if line.startswith('processor'):
            NumberOfPhysicalCPU += 1
        elif line.startswith('vendor_id'):
            VendorID = line.split(':')[1].strip()
        elif line.startswith('model name'):
            VendorString = line.split(':')[1].strip()
        elif line.startswith('cpu family'):
            FamilyID = line.split(':')[1].strip()
        elif line.startswith('model'):
            ModelID = line.split(':')[1].strip()
        elif line.startswith('cache size'):
            ProcessorCacheSize = line.split(':')[1].strip().split(' ')[0]
        elif line.startswith('cpu MHz'):
            ProcessorClockFrequency = int(float(line.split(':')[1].strip()))

    TotalPhysicalMemory = 0
    for line in open('/proc/meminfo', 'r'):
        if line.startswith('MemTotal:'):
            TotalPhysicalMemory = int(line.split(':')[1].strip().split()[0])/1024

    for v in ("BuildName", "BuildStamp", "Name", "Generator", "OSName", \
                  "Hostname", "OSRelease", "OSVersion", "OSPlatform", \
                  "VendorString", "VendorID", "FamilyID", "ModelID", \
                  "ProcessorCacheSize", "NumberOfPhysicalCPU", "Is64Bits", \
                  "ProcessorClockFrequency", "TotalPhysicalMemory") :
        ss += "\t%s=\"%s\"\n" % (v, eval(v))
    return ss


saved_line = saved_file = ''
saved_lnum = 0
# parse a line of compiler output messages
def parse_line(line):
    global saved_line, saved_file, saved_lnum
    try:
        file, lnum, errstr = re.compile(r'(\w+.\w+):(\d+).*(warning|error):').findall(line)[0]
        saved_line = ''           
    except:
        try:            
            saved_file, saved_lnum  = re.compile(r'(\w+.\w+):(\d+)').findall(line)[0]
            saved_line = line
        except:
            if saved_line:
                saved_line += line
                if line.startswith('Warning:'):
                    line = saved_line
                    saved_line = ''
                    return saved_file, saved_lnum, 0, line
                if line.startswith('Error:'):
                    line = saved_line                    
                    saved_line = ''
                    return saved_file, saved_lnum, 1, line                
        return '', 0, 0, ''

    if errstr == 'error':
        err = 1
    else:
        err = 0
                
    return file, lnum, err, line

# given a partial path, find it in the sources and return
# the complete path relative to the toplevel sources directory.
def find_file(bp, fpath):
    where = len(bp.bld.topdir) + 1
    if os.path.isabs(fpath):
        return fpath
    fpath = os.path.normpath('/' + fpath)
    for root, dirs, files in os.walk(bp.sdir):
        for x in files:
            f = os.path.join(root, x)
            if f.endswith(fpath):
                if f.startswith(bp.bld.topdir):
                    f = f[where:]
                return f
    return fpath

def parse_builds(bp, filename):
    ostr = ''
    for line in open(filename):
        fname, lnum, err, line = parse_line(line)
        if fname:
            if err:
                ostr += "<Error>\n"
            else:
                ostr += "<Warning>\n"
            ostr += "\t<Text>%s</Text>\n" % escape(line)
            fname = find_file(bp, fname)
            ostr += "\t<SourceFile>%s</SourceFile>\n" % fname
            ostr += "\t<SourceLineNumber>%s</SourceLineNumber>\n" % lnum
            ostr += "<PreContext></PreContext><PostContext></PostContext><RepeatCount>0</RepeatCount><BuildLogLine>0</BuildLogLine>\n"
            if err:
                ostr += "</Error>\n"
            else:
                ostr += "</Warning>\n"
    return ostr

def build_str(bp, bstart, bstop, BuildCommand):
    StartDateTime = time.strftime("%b %d %H:%M %Z", time.localtime(bstart))
    StartBuildTime = bstart
    EndDateTime = time.strftime("%b %d %H:%M %Z", time.localtime(bstop))
    EndBuildTime = bstop
    ElapsedMinutes = (bstop - bstart) / 60.0
    bs = ""
    for v in ("StartDateTime", "StartBuildTime", "BuildCommand"):
        bs += "\t<%s>%s</%s>\n" % (v, eval(v), v)

    # warnings and errors go here
    for p in bp.packages:
        if not p.is_pkg:
            bs += parse_builds(p, os.path.join(bp.logdir, "%s-build.log" % p.name))
    bs += "\t<Log Encoding=\"base64\" Compression=\"/bin/gzip\"></Log>\n"

    for v in ("EndDateTime", "EndBuildTime", "ElapsedMinutes"):
        bs += "\t<%s>%s</%s>\n" % (v, eval(v), v)
    return bs

def parse_tests(fname):
    tl = []
    try:
        for line in open(fname):
            line = line.strip()
            if line.startswith("<Name>"):
                where = line.find("</Name>")
                line = line[6:where].strip()
                tl.append(line)
    except:
        pass
    return tl

def test_str(bp, tstart, tstop):
    testlist = []
    for p in bp.packages:
        testlist += parse_tests(os.path.join(bp.logdir, "%s-test.xml" % p.name))
    if os.path.isfile(os.path.join(bp.logdir, "MEMOSA-test.xml")):
        testlist += parse_tests(os.path.join(bp.logdir, "MEMOSA-test.xml"))

    StartDateTime = time.strftime("%b %d %H:%M %Z", time.localtime(tstart))
    StartTestTime = tstart
    ts = ""
    for v in ("StartDateTime", "StartTestTime"):
        ts += "<%s>%s</%s>\n" % (v, eval(v), v)

    ts += "<TestList>\n"
    for t in testlist:
        ts += "\t<Test>%s</Test>\n" % escape(t)
    ts += "</TestList>\n"

    for p in bp.packages:
        try:
            ts += open(os.path.join(bp.logdir, "%s-test.xml" % p.name)).read()
        except:
            pass

    if os.path.isfile(os.path.join(bp.logdir, "MEMOSA-test.xml")):
        ts += open(os.path.join(bp.logdir, "MEMOSA-test.xml")).read()

    EndDateTime = time.strftime("%b %d %H:%M %Z", time.localtime(tstop))
    EndTestTime = tstop
    ElapsedMinutes = (tstop - tstart) / 60.0
    for v in ("EndDateTime", "EndTestTime", "ElapsedMinutes"):
        ts += "<%s>%s</%s>\n" % (v, eval(v), v)
    return ts

def submit(bp, cname, cmd, nightly):
    cmd = ' '.join(cmd)
    if nightly:
        bname = "NIGHTLY"
    else:
        bname = "EXPERIMENTAL"
    ss = site_str(cname, bname)

    bs = float(open(bp.logdir + '/StartBuildTime').read())
    be = float(open(bp.logdir + '/EndBuildTime').read())
    bs = build_str(bp, bs, be, cmd)

    tests = True
    try:
        ts = float(open(bp.logdir + '/StartTestTime').read())
        te = float(open(bp.logdir + '/EndTestTime').read())
        ts = test_str(bp, ts, te)
    except:
        tests = False

    fname = os.path.join(bp.logdir, "Build.xml")
    f = open(fname, 'w')
    f.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
    f.write("<Site %s\n>\n<Build>\n%s\n</Build>\n</Site>" % (ss, bs))
    f.close()
    
    ret = putname(fname, "http://dash.prism.nanohub.org/cdash/submit.php?project=MEMOSA")
    if ret[0] != 200:
        build_utils.warning("http put returned %s" % ret[0])

    if tests:
        fname = os.path.join(bp.logdir, "Test.xml")
        f = open(fname, 'w')
        f.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
        f.write("<Site %s\n>\n<Testing>\n%s\n</Testing>\n</Site>" % (ss, ts))
        f.close()
        ret = putname(fname, "http://dash.prism.nanohub.org/cdash/submit.php?project=MEMOSA")
        if ret[0] != 200:
            build_utils.warning("http put returned %s" % ret[0])

    fname = os.path.join(bp.logdir, "Update.xml")
    if os.path.isfile(fname):
        ret = putname(fname, "http://dash.prism.nanohub.org/cdash/submit.php?project=MEMOSA")
        if ret[0] != 200:
            build_utils.warning("http put returned %s" % ret[0])



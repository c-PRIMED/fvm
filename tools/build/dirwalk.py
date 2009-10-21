#!/usr/bin/env python
import os

class DirWalk(object):
    ignore = ['.', '..', '.svn', '.cvsignore', 'lnx86_64', 'lnx86', 'lnxia64']
    ignore_suffixes = ['.pyc', '.bak', '~']

    def __init__(self, dir):
        self.dir = dir

    def ok(self, file):
        if file in self.ignore:
            return False
        for s in self.ignore_suffixes:
            if file.endswith(s):
                return False
        return True

    def walk(self, dir=None):
        """ walks a directory, and yields each file """
        if not dir:
            dir = self.dir
        dir = os.path.abspath(dir)
        for file in [file for file in os.listdir(dir) if self.ok(file)]:
            nfile = os.path.join(dir, file)
            if os.path.isdir(nfile):
                for nfile in self.walk(nfile):
                    yield nfile
            else:
                yield nfile

    def find_newest(self):
        """ find the most recently modified file in a directory """
        most_recent_time = 0
        most_recent_file = ''
        dir = os.path.abspath(self.dir)
        if not os.path.isdir(dir):
                return dir, os.path.getmtime(dir)
        for file in [file for file in os.listdir(dir) if self.ok(file)]:
            nfile = os.path.join(dir, file)
            if os.path.isdir(nfile):
                for nfile in self.walk(nfile):
                    if os.path.getmtime(nfile) > most_recent_time:
                        most_recent_time = os.path.getmtime(nfile)
                        most_recent_file = nfile
            else:
                if os.path.getmtime(nfile) > most_recent_time:
                    most_recent_time = os.path.getmtime(nfile)
                    most_recent_file = nfile
        return most_recent_file, most_recent_time

if __name__ == "__main__":
    #file, t = DirWalk('.').find_newest()
    #print file, time.strftime("%b %d %Y %H:%M %Z", time.localtime(t))
    for file in DirWalk('.').walk():
        print file

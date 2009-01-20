#include this file at any point in any of the .py scripts to start the
#debugger on the current process. Needs some work for ntx86

import sys
import time
import os
os.system('ddd %s %d& ' % (sys.executable,os.getpid()))
time.sleep(2)

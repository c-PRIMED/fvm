import distutils.sysconfig
import os

def generate(env):
    env.Append(LIBS=['umfpack', 'amd', 'ufconfig'])




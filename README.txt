MEMOSA: MEMS Overall Simulation Administrator

MEMOSA is developed as part of the PRISM project.
http://www.purdue.edu/research/prism/index.shtml


BUILDING
========

The Makefile just calls make.py.  

Build all the software
> make
or
> ./make.py

Build a system-specific set of software:
./make.py configname

Valid configs are stored in the config subdirectory. Current there
are configs for mmh (a linux desktop) and steele.  You can easily add your own.


To see verbose output
./make.py -v

---
Not working yet:

Run tests
> make test


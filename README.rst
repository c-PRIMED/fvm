************
Introduction
************

:Version: 1.0
:Authors: see src/fvm/AUTHORS
:Web site: http://c-PRIMED.github.io/fvm/
:Documentation: http://c-PRIMED.github.io/fvm/
:Copyright: This document has been placed in the public domain.
:License: MIT License.

Purpose
=======

A finite volume solver for fluid, thermal and charge transport.

Features
========

* Includes a Python interface.


Dependencies
============

Python 2.6+.


Installation
============

To build,
 ./make.py configname

Configurations are stored in the 'configs' subdirectory.  To build under
recent Ubuntu distributions,

 > sudo apt-get install g++ libboost-dev swig python-tk libopenmpi-dev python-mpi4py libmpfr-dev libhdf5-dev python-nose libpython-dev cmake python-h5py m4 libblas-dev

 > ./make.py ubuntu

When finished, there will be env-ubuntu.sh and env-ubuntu.csh. To use FVM, you will need to source the correct one for your shell to add the proper settings to your environment.

 > source env-ubuntu.sh


TESTING
=======

 ./make.py --test configname

History
=======

FVM is based upon work supported by the Department of Energy [National Nuclear Security Administration]
under Award Number DE-FC52-08NA28617.”


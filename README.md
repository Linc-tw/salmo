README
======
Chieh-An Lin (IfA Edinburgh)  
Date: 2019-12-11  


Description
-----------

This code generates fast mocks of galaxy tracers with or without lensing signals by sampling from HEALPix maps,
with regard to custom masks and redshift distributions in a very flexible way.


Requirements
------------

Required softwares:
- [cmake](https://cmake.org/cmake/resources/software.html)
- [pkg-config](https://www.freedesktop.org/wiki/Software/pkg-config/)
- [gcc](https://gcc.gnu.org/)
- [gsl](https://www.gnu.org/software/gsl/)
- [cfitsio](https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html)
- [chealpix](https://healpix.jpl.nasa.gov/index.shtml)

Optional software:
- [healpix_cxx](https://healpix.jpl.nasa.gov/index.shtml)

During the compilation, cmake uses pkg-config to find optional softwares. If they are missing, the compilation still continues without providing all functionalitites.


Compilation
-----------

For Mac users, do the follows before compilation:
```Bash
$ export CC=gcc
$ export CXX=g++
```
or use `setenv` command in tcsh.

To compile the package:
```Bash
$ cd build
$ cmake ..
$ make
```

To get program instructions:
```Bash
$ ./mockFootprint



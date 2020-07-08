Salmo
=====

Speedy Acquisition for Lensing and Matter Observables  
Chieh-An Lin (IfA Edinburgh)  
Date: 2020-07-08  


Description
-----------

Salmo is a software to generate fast mocks for cosmological analysis with the following characteristics:
- galaxy mocks, with or without lensing signals;
- curved-sky;
- can generate multiple tracers coherently;
- each tracer can have its own mask and redshift distribution.

Salmo is designed to be used with the [Flask]() software. 
The map outputs from Flask or those with the same format are required to generate catalogues.
See **Usage** and [Wikis]() for details.


Installation
------------

### Requirements

Required softwares:
- [cmake](https://cmake.org/cmake/resources/software.html)
- [pkg-config](https://www.freedesktop.org/wiki/Software/pkg-config/)
- [gcc](https://gcc.gnu.org/)
- [gsl](https://www.gnu.org/software/gsl/)
- [cfitsio](https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html)
- [chealpix](https://healpix.jpl.nasa.gov/index.shtml)

Optional software:
- [healpix_cxx](https://healpix.jpl.nasa.gov/index.shtml)

During the compilation, `cmake` uses `pkg-config` to find optional softwares. If they are missing, the compilation still continues without providing all functionalitites.


### Compilation

For Mac users, do the following before compilation:
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

To see if the compilation is successful, please execute:
```Bash
$ ./mockFootprint
```
and users will see usage instructions.


Usage
-----

### External files to provide

Salmo requires some external files as inputs to work:
- density & lensing maps,
- masks,
- redshift distributions, and 
- values like galaxy number densities & shape noise.

### Quick Run
(To be completed)

### Instructions Reminder

An example parameter file is provided: `param/MFPParam.par`. 

To get program instructions, please execute:
```Bash
$ ./mockFootprint
```


License
-------

Salmo is distributed under the terms of the [GNU General Public License Version 3 (GPLv3)](https://www.gnu.org/licenses/).

It gives to users the option to distribute an application which uses Salmo under the terms of GNU GPLv3.


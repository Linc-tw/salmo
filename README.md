Salmo
=====

Speedy Acquisition for Lensing and Matter Observables  
Chieh-An Lin (IfA Edinburgh)  
Date: 2020-07-08  


Description
-----------

Salmo is a software to generate fast mocks for cosmological analysis with the following characteristics:
- galaxy mocks, with or without lensing signals;
- curved sky;
- can generate multiple tracers coherently;
- each tracer can have its own mask and redshift distribution.

Salmo is designed to be used with the [Flask](https://github.com/hsxavier/flask) software. 
The map outputs from Flask or those with the same format are required to generate catalogues.
See **Usage** and [Wiki](https://github.com/Linc-tw/salmo/wiki) for details.


Installation
------------

### Requirements

- [cmake](https://cmake.org/cmake/resources/software.html)
- [pkg-config](https://www.freedesktop.org/wiki/Software/pkg-config/)
- [gcc](https://gcc.gnu.org/)
- [gsl](https://www.gnu.org/software/gsl/)
- [cfitsio](https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html)
- [chealpix](https://healpix.jpl.nasa.gov/index.shtml)


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
$ ./salmo
```
and users will see usage instructions.


Usage
-----

### External files to provide

Salmo requires some external files as inputs to work:
- density & lensing maps,
- masks, and 
- redshift distributions.

Other values to specify are galaxy number densities & shape noise.

### Quick Run

With the provided parameter file `param/salmoParam.par` and input files in `demo/input`, 
users can execute:
```Bash
$ ./salmo default 3
```
for a quick sample run and find catalogues in `demo/output`.

See [Wiki](https://github.com/Linc-tw/salmo/wiki) for detailed tutorials.

### Instructions reminder

To get program instructions, please execute:
```Bash
$ ./salmo
```


License
-------

Salmo is distributed under the terms of the [GNU General Public License Version 3 (GPLv3)](https://www.gnu.org/licenses/).

It gives to users the option to distribute an application which uses Salmo under the terms of GNU GPLv3.


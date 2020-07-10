Salmo
=====

Speedy Acquisition for Lensing and Matter Observables  
Chieh-An Lin (IfA Edinburgh)  
Date: 2020-07-10  


Description
-----------

_Salmo_ is a C-based software to generate fast mocks for cosmological analysis with the following characteristics:
- galaxy mocks, with or without lensing signals;
- curved sky;
- can generate multiple tracers coherently;
- each tracer can have its own mask and redshift distribution.

_Salmo_ is designed to be used with the [_Flask_](https://github.com/hsxavier/flask) software. 
The map outputs from _Flask_ or those with the same format are required to generate catalogues.
See **Usage** and [Wiki](https://github.com/Linc-tw/salmo/wiki) for details.


Installation
------------

### Requirements

- [cmake](https://cmake.org/cmake/resources/software.html)
- [pkg-config](https://www.freedesktop.org/wiki/Software/pkg-config/)
- [gcc](https://gcc.gnu.org/)
- [gsl](https://www.gnu.org/software/gsl/)
- [cfitsio](https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html)


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

_Salmo_ requires some external files as inputs to work:
- density & lensing maps,
- masks, and 
- redshift distributions.

Other values to specify are galaxy number densities & shape noise.

### Quick run

With the provided parameter file `param/salmoParam.par` and input files in `demo/input`, 
users can execute:
```Bash
$ ./salmo default 3
```
for a quick example run and find catalogues in `demo/output`.

See [Wiki](https://github.com/Linc-tw/salmo/wiki) for detailed tutorials.

### Instructions reminder

To get program instructions, please execute:
```Bash
$ ./salmo
```


License
-------

_Salmo_ is released under [GNU General Public License Version 3 (GPLv3)](https://www.gnu.org/licenses/).

It contains 2 pieces of code ([this](https://github.com/Linc-tw/salmo/blob/master/source/chealpix.c) 
and [that](https://github.com/Linc-tw/salmo/blob/master/source/chealpix.h))
taken from [HEALPix](https://healpix.sourceforge.io/index.php) and is released under 
[GNU General Public License Version 2 (GPLv2)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html).

This gives users the option to distribute an application which uses _Salmo_ under the same GPLv3 license.


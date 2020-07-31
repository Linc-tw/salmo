Salmo
=====

Speedy Acquisition for Lensing and Matter Observables  
Chieh-An Lin (Institute for Astronomy, University of Edinburgh)  
Date: 2020-07-30  
![logo](doc/Salmo_400px.png)


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

Attribution
-----------

Please cite the following paper for referencing _Salmo_:
- [Joachimi, Lin, et al.](https://arxiv.org/abs/2007.01844). Submitted to A&A. _KiDS-1000 Methodology: Modelling and inference for joint weak gravitational lensing and spectroscopic galaxy clustering analysis_.

Here are some other works in which _Salmo_ is involved:
- [Giblin et al.](https://arxiv.org/abs/2007.01845). Submitted to A&A. _KiDS-1000 catalogue: weak gravitational lensing shear measurements_.
- [Heydenreich et al. (2020)](https://doi.org/10.1051/0004-6361/201936966). A&A, 634, A104. _The effects of varying depth in cosmic shear surveys_.


Contributing to _Salmo_
-----------------------

No upcoming development is planned by the owner. 
Only maintenance will be done.
However, developers are welcome to contribute to _Salmo_ in various ways:
- draft an issue,
- e-mail to the owner at calin(at)roe.ac.uk,
- propose a pull request from a fork of this repository, or
- request for being a collaborator.


License
-------

_Salmo_ is released under [GNU General Public License Version 3 (GPLv3)](https://www.gnu.org/licenses/).

It contains 2 pieces of code ([this](https://github.com/Linc-tw/salmo/blob/master/source/chealpix.c) 
and [that](https://github.com/Linc-tw/salmo/blob/master/source/chealpix.h))
taken from [HEALPix](https://healpix.sourceforge.io/index.php) and is released under 
[GNU General Public License Version 2 (GPLv2)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html).

This gives users the option to distribute an application which uses _Salmo_ under the same GPLv3 license.


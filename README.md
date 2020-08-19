# dsp

dsp: functions for DSP (digital signal processing)

Erik Edwards (erik.edwards4@gmail.com)

================================================

This is a set of C programs, and associated command-line tools in C++,
that implement many of the usual DSP functions found in Octave, etc.
However, the intention is to support the voice and neural work
with a unified interface and set of functions. This is not over-complete,
with a vast number of DSP functions for all purposes.

The command-line programs are written in C++ with a consistent style and interface.
The low-level functions themselves are written in C for fastest performance (e.g., openBLAS).

The C functions are meant for the developer; the C++ command-line tools are meant for the end-user.
The interface to each C function is BLAS-like, meaning that one specifies the input and/or output dimensions,
the matrix order as row-major or column-major, and so on.

The C++ command-line programs are written in a consistent style that was developed for command-line tools in general.
All of these command-line tools use argtable2 (http://argtable.sourceforge.net/) for parsing
inputs and option flags. All of them allow -h (--help) as a flag to give description and usage info.

Input/output is supported for NumPy tensors (https://numpy.org/)
and several C++ tensor formats: Armadillo (http://arma.sourceforge.net/),
ArrayFire (https://arrayfire.com/), and a minimal format for Eigen (http://eigen.tuxfamily.org/).


## Dependencies
Requires argtable2, openBLAS, LAPACKE, FFTW.
For Ubuntu, these are available by apt-get:
```
sudo apt-get install libargtable2-0 libblas3 libopenblas-base liblapack3 liblapacke fftw3
```


## Installation
```
cd /opt/codee
git clone https://github.com/erikedwards4/dsp
cd /opt/codee/dsp
make
```

Each C function can also be compiled separately; see c subdirectory Makefile for details.
To make an archive library, do:
```
cd /opt/codee/dsp/c
make libdsp.a CC=clang
```
This creates /opt/codee/math/lib/libdsp.a with all of the C object files.
This could be useful if trying to use the C functions in other applications.
Change clang to clang++ to compile for use with C++ applications.


## Usage
See each resulting command-line tool for help (use -h or --help option).
For example:
```
/opt/codee/dsp/bin/fir --help
```


## Contributing
This is currently only to view the project in progress.


## License
[BSD 3-Clause](https://choosealicense.com/licenses/bsd-3-clause/)


# QMCkl: Quantum Monte Carlo Kernel Library (GPU addon)

<img src="https://trex-coe.eu/sites/default/files/styles/responsive_no_crop/public/2022-01/QMCkl%20code.png?itok=UvOUClA5" width=200>

This repository is a work-in-progress GPU addon to the [QMCkl library](https://github.com/TREX-CoE/qmckl). It provides alternatives to the standard QMCkl functions to perform computations on GPU. The library is based on OpenMP offload and OpenACC, as to support a wide variety of compilers and targets combinations.

This readme contains configure, build and installation instructions, and redirects to dedicated files for usage instructions and the troubleshooting section.

# Installation

The project uses GNU Autotools :

```
bash autogen.sh
./configure [GPU functions-related arguments & other arguments]
make
make install
```

The only other requirement for **building** the library is a compiler toolchain that supports OpenMP or OpenACC offloading for your chosen target(s). Linking with QMCkl to build QMCkl GPU is not needed (but it is when building an executable that uses QMCkl GPU).

For more details, se the dedicated INSTALL.md file.

## Enabling either OpenMP or OpenACC

Enabling either OpenMP and OpenACC is done at configure time, where (exactly) one of the two options has to be specified :

```
./configure --enable-[omp|acc]
```

In either case, the library interface is going to be exactly the same, as all of the OpenMP/OpenACC specific syntaxes are wrapped inside  QMCkl GPU's functions.

**Note:** Using exactly one of those arguments is mandatory, as doing otherwise would result in an empty library. If none is specified, Autotools will throw a warning message and attempt to build the library with `--enable-omp` as a fallback solution.


# Usage

See the dedicated [USAGE.md](https://github.com/TREX-CoE/qmckl_gpu/blob/main/doc/USAGE.md) file.


# Troubleshooting

During the development and testing of the library, we encountered some compiler related issues. This section contains a list of known issues and fixes if any. 

See the dedicated [TROUBLESHOOTING.md](https://github.com/TREX-CoE/qmckl_gpu/blob/main/doc/TROUBLESHOOTING.md) file.


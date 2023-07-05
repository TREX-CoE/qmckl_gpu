# QMCkl: Quantum Monte Carlo Kernel Library (GPU addon)

<img src="https://trex-coe.eu/sites/default/files/styles/responsive_no_crop/public/2022-01/QMCkl%20code.png?itok=UvOUClA5" width=200>

This repository is a GPU port of the [QMCkl library](https://github.com/TREX-CoE/qmckl). It provides alternatives to the standard QMCkl functions to perform computations on GPU. The library is based on OpenMP offload and OpenACC, as to support a wide variety of compilers and targets combinations.

This document contains configure, build and installation instructions, and redirects to dedicated files for usage instructions and the troubleshooting section.


# Build & install instructions

The project uses GNU Autotools :

```
bash autogen.sh
./configure [arguments]
make
make install
```

The only other requirement for the library is a compiler toolchain that supports OpenMP or OpenACC offloading for your chosen target(s). The library is now completely standalone, and it can be linked without QMCkl CPU.

You can also check that offloading works on your system with the `make check` command.


## Enabling or disabling TREXIO

The [TREXIO](https://github.com/TREX-CoE/trexio) integration enables the initialization of a device context from a TREXIO file.

By default, the configure will try to enable TREXIO by autodetecting it, but will leave it disabled if unable to find the library. This will be specified in the configure summary.

- If TREXIO can not be found, but you want to use it : specify the path by using the `--enable-trexio=...` configure option. Specify the TREXIO install path that contains the `lib` and `include` subdirectories. If you specified an incorrect path, the configure will fail.
- If TREXIO is found on your machine but you don't want to use it : explicitely disable it by using the `--disable-trexio` configure option.


## Enabling either OpenMP or OpenACC

Building the library with either OpenMP or OpenACC is done at configure time, where (exactly) one of the two options has to be specified :

```
./configure --enable-[omp|acc]
```

In either case, the library interface is going to be exactly the same, as all of the OpenMP/OpenACC specific syntaxes are wrapped inside QMCkl GPU's functions.

**Note:** If none of those arguments is specified, the configure will throw a warning message and attempt to build the library with `--enable-omp` as a fallback solution.


## Compilers support and autoflags

We currently support nvc, gcc and clang out of the box. This means we have succesfully built and run the library with one of these compilers, on hardware from at least one vendor. You can specify which compiler to use by specifying the `CC=...` variable to the configure (gcc should be the default).

When specifying a supported compiler, the configure also automatically tries to set the required flags to enable OpenMP or OpenACC offloading. In case the proposed flags don't work on your system, you can disable them with the `--disable-autoflags` configure option. Then, simply specify correct compiler flags manually in the `CFLAGS` variable.

### For AMD GPU

If you use AMD GPU you need to specified other thing. Change your compiler (HIPCC or clang with AMD gpu configuration), and add some specific compiler flags to your hardware.

```
./configure --enable-omp CC=/opt/rocm/hip/bin/hipcc  CFLAGS="-target x86_64-pc-linux-gnu  -fopenmp=libomp -fopenmp-targets=amdgcn-amd-amdhsa -Xopenmp-target=amdgcn-amd-amdhsa -march=gfx1030"  --prefix=$PWD/_INSTALL

```
You can found your specific arch for -march with ```  /opt/rocm/bin/rocminfo | grep gfx ```


# Usage

See the dedicated [USAGE.md](https://github.com/TREX-CoE/qmckl_gpu/blob/main/doc/USAGE.md) file.


# Troubleshooting

During the development and testing of the library, we encountered some compiler related issues. This section contains a list of known issues and fixes (if any).

See the dedicated [TROUBLESHOOTING.md](https://github.com/TREX-CoE/qmckl_gpu/blob/main/doc/TROUBLESHOOTING.md) file.


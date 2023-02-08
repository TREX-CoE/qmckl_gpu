# QMCkl: Quantum Monte Carlo Kernel Library (GPU addon)

<img src="https://trex-coe.eu/sites/default/files/styles/responsive_no_crop/public/2022-01/QMCkl%20code.png?itok=UvOUClA5" width=200>

This repository is a work-in-progress GPU addon to the [QMCkl library](https://github.com/TREX-CoE/qmckl). It provides alternatives to the standard QMCkl functions to perform computations on GPU. The library is based on OpenMP offload and OpenACC, as to support a wide variety of compilers and targets combinations.


# Basic installation

The project uses GNU Autotools :

```
bash autogen.sh
./configure [GPU functions-related arguments & other arguments]
make
make install
```

The only other requirement for **building** the library is a compiler toolchain that supports OpenMP or OpenACC offloading for your chosen target(s). Linking with QMCkl to build QMCkl GPU is not needed (but it is when building an executable that uses QMCkl GPU).


## Enabling either OpenMP or OpenACC

Enabling either OpenMP and OpenACC is done at configure time, where (exactly) one of the two options has to be specified :

```
./configure --enable-[omp|acc]
```

In either case, the library interface is going to be exactly the same, as all of the OpenMP/OpenACC specific syntaxes are wrapped inside  QMCkl GPU's functions.

**Note:** Using exactly one of those arguments is mandatory, as doing otherwise would result in an empty library. If none is specified, Autotools will throw a warning message and attempt to build the library with `--enable-omp` as a fallback solution.


# Basic usage

This library is an addon meant to be used on top of the main [QMCkl library](https://github.com/TREX-CoE/qmckl). While building QMCkl GPU doesn't require any linking with QMCkl, executables calling QMCkl GPU should be linked with

```
-lqmckl -lqmckl_gpu
```

Both `qmckl.h` and `qmckl_gpu.h` should also be included in C codes.

For instance, when linking their program with `-lqmckl`, users will have access to the "standard"

```
qmckl_exit_code
qmckl_get_ao_basis_ao_vgl(qmckl_context context,
                          double* const ao_vgl,
                          const int64_t size_max);
```

QMCkl function. And additionally, linking with `-lqmckl_gpu` will also expose

```
qmckl_exit_code
qmckl_get_ao_basis_ao_vgl_device (qmckl_context_device context,
                                  double* const ao_vgl,
                                  const int64_t size_max);
```

which performs the same task on a GPU, using memory allocated there.

**Notes :** 
- QMCkl GPU should be linked with an HPC-enabled QMCkl (configure QMCkl with `--enable-hpc`, see the [QMCkl README](https://github.com/TREX-CoE/qmckl/blob/master/README.md)), as well as with [TREXIO](https://github.com/TREX-CoE/trexio).
- While it is possible to use different compilers to build QMCkl, QMCkl GPU and your final executable, doing so might result in some missing symbols when linking your final executable. In that case, you would have to manually tell your linker where to find those symbols (which might come from different libc or libm implementations being simultaneously, for instance). 

The next sections provides more advanced explanations on how this library should be used compared to the standard CPU QMCkl.


## Using the "device" context and functions

The library works on a variant of the `qmckl_context` type: `qmckl_context_device`. It has the same structure as the default context, but contains a few additional and GPU-related informations, such as the OpenMP/OpenACC ID of the device you are allocating memory on and doing your offload to. This type comes with its own variants of the classic QMCkl functions. These functions work directly and (almost) exclusively on the GPU memory, referenced by "device pointers". As a rule of thumb, you should always pass such device pointers as arguments of these functions. This way, memory transfers are greatly reduced compared to a "naive" approach where memory is transferred back and forth everytime we need to perform a computation on GPU. Typically the entry data set would need to be allocated/transferred once from CPU to GPU at the beginning of a QMC simulation, and the results would be transferred on demand only when they are needed elsewhere. The main dowside is that **this will likely require the user to perform a few GPU allocations themself**, which requires a bit more care than using the CPU library.

**TODO** More detailed explanations on device pointer functions use, dedicated .md documentation files ?


## Troubleshooting

During the development and testing of the library, we encountered some compiler related issues. This section contains the fixes we used in case you run into the same errors : 

### Link error with nvc

When building an executable with an `nvc` built QMCkl GPU (when doing a `make check` for instance), you might run into this issue, or similar : 

```
/usr/bin/ld: /opt/nvidia/hpc_sdk/Linux_x86_64/22.11/compilers/lib/libgomp.so.1: undefined reference to `__pgi_nvomp_error_cuda_noversion'
```

It seems that by itself, nvc doesn't link correctly all of its libraries, resulting in missing symbols errors. We have been able to solve the issue by manually specifying the flags needed to link those libraries  : 

```
./configure CC="nvc" LDFLAGS="-L/[path]/[to]/nvidia/hpc_sdk/Linux_x86_64/[version]/compilers/lib -lpgc -lnvf"
```

... should solve the issue when building the tests. You can do something similar when building your own app.


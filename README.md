# QMCkl: Quantum Monte Carlo Kernel Library (GPU addon)

<img src="https://trex-coe.eu/sites/default/files/styles/responsive_no_crop/public/2022-01/QMCkl%20code.png?itok=UvOUClA5" width=200>

This repository is a GPU addon to the [QMCkl library](https://github.com/TREX-CoE/qmckl). It provides alternatives to the standard QMCkl functions to offload computations to the GPU.



# Installation

The project uses Autotools :

```
bash autogen.sh
./configure [GPU arguments & other arguments]
```

Where the GPU arguments can be :

| Argument | Description |
| ----------- | ----------- |
| **--enable-openmp** | Enable OpenMP offloaded functions |
| --with-device-pointers | Optionally enable OpenMP functions using device pointers (requires **--enable-openmp**) |
| **--enable-openacc** | Enable OpenACC offloaded functions |

**Note:** Using at least one GPU argument is mandatory. Doing otherwise would result in an empty library. If no GPU argument is specified, OpenMP offload will be enabled by default.

Then :

```
make
make install
```

The only requirement is a compiler with offloading support of your chosen GPU programming paradigm(s) to your desired target(s). Linking with QMCkl to build QMCkl GPU is not needed.



# Usage

This library is an addon to be used on top of the main [QMCkl library](https://github.com/TREX-CoE/qmckl). While building QMCkl GPU doesn't require any linking with QMCkl, programs calling QMCkl GPU should be linked with

```
-lqmckl -lqmckl_gpu
```

Both `qmckl.h` and `qmckl_gpu.h` should also be included in C codes.

**Note :** QMCkl GPU should be linked with an HPC-enabled QMCkl (configure QMCkl with `--enable-hpc`, see the [QMCkl README](https://github.com/TREX-CoE/qmckl/blob/master/README.md))


## Alternative get functions

Enabling one kind of GPU functions in the configure step will expose variants of the usual `get` functions. For instance, in order to perform atomic orbitals computations, and by configuring QMCkl GPU with

```
./configure --enable-openmp
```

users will still have access to the standard

```
qmckl_exit_code
qmckl_get_ao_basis_ao_vgl(qmckl_context context,
                          double* const ao_vgl,
                          const int64_t size_max);
```

function when linking their program with `-lqmckl`. However, additonally linking with `-lqmckl_gpu` will also expose

```
qmckl_exit_code
qmckl_get_ao_basis_ao_vgl_omp_offload (qmckl_context context,
                                       double* const ao_vgl,
                                       const int64_t size_max);
```

As a general rule, each kind of GPU functions comes with its own suffix. In order to call such functions, one can simply append the corresponding suffix to the standard `get` function name :

| Function type | Suffix |
| ----------- | ----------- |
| Simple OpenMP offload | `_omp_offload` |
| OpenMP + device pointers | `_device` |
| Simple OpenACC offload | `_acc_offload` |

The only exception to this rule are device pointers functions, which require special attention before being called.


## Device pointer functions

**TODO** Explain how to call the device pointers functions once they are added to the repository, as they require special care.

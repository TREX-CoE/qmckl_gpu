# QMCkl: Quantum Monte Carlo Kernel Library (GPU addon)

<img src="https://trex-coe.eu/sites/default/files/styles/responsive_no_crop/public/2022-01/QMCkl%20code.png?itok=UvOUClA5" width=200>

This repository is a GPU addon to the [QMCkl library](https://github.com/TREX-CoE/qmckl). It provides alternatives to the standard QMCkl functions (runningon CPU) to offload computations to the GPU.



# Installation

The project uses Autotools :

```
bash autogen.sh
./configure [GPU functions-related arguments & other arguments]
make
make install
```

The only requirement is a compiler with offloading support for your chosen GPU programming paradigm(s) and target(s). Linking with QMCkl to build QMCkl GPU is not needed.



# Usage

This library is an addon to be used on top of the main [QMCkl library](https://github.com/TREX-CoE/qmckl). While building QMCkl GPU doesn't require any linking with QMCkl, programs calling QMCkl GPU should be linked with

```
-lqmckl -lqmckl_gpu
```

Both `qmckl.h` and `qmckl_gpu.h` should also be included in C codes.

**Note :** QMCkl GPU should be linked with an HPC-enabled QMCkl (configure QMCkl with `--enable-hpc`, see the [QMCkl README](https://github.com/TREX-CoE/qmckl/blob/master/README.md))


## Enabling the different GPU functions

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

As a general rule, each type of GPU functions comes with its own suffix. In order to call such functions, one can simply append the corresponding suffix to the standard `get` function name :

| Function type | Suffix |
| ----------- | ----------- |
| Simple OpenMP offload | `_omp_offload` |
| Simple OpenACC offload | `_acc_offload` |
| OpenMP + device pointers | `_device` |

The only exception to this rule are device pointers functions, which require special attention before being called. These 3 functions types can be enabled at configure time with the following options : 

| Function type | Suffix |
| ----------- | ----------- |
| Simple OpenMP offload | `--enable-openmp` |
| Simple OpenACC offload | `--enable-openacc` |
| OpenMP + device pointers | `--enable-device` |

**Note:** Using at least one of those arguments is madatory, as doing otherwise would result in an empty library. If none is specified, `--enable-openmp` will be toggled on by default.


## Device pointer functions

The difference between the basic offload functions and the device functions is that device versions work all the way with device pointers. 

Indeed, when calling an offload `get` function, everything gets executed as in the CPU version of QMCkl, at the exception of some computation kernels that are offloaded "on the fly". This means that memory is allocated and transferred to the GPU just before the computation, only to be transferred back to the CPU and freed just after. This makes them easy to use, as the user can call them just like the standard CPU functions, but at the cost of performance when calling these kernels repeatedly (such as in a QMC simulation), which creates some costly and avoidable data transfers.

This is why the use of the device functions is advised to get the full benefits of GPU offloading. These functions work directly on the GPU memory, referenced by "device pointers". This way, memory transfers can be greatly reduced. Typically the entry data set would need to be done once from CPU to GPU at the beginning of a QMC simulation, theb the result would be transferred once at the very end. The only issue is that it will typically require the user to perform a few GPU allocations themself, which requires a bit more care.

**TODO** More detailed explanations on device pointer functions use

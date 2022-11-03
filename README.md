# QMCkl: Quantum Monte Carlo Kernel Library (GPU addon)

<img src="https://trex-coe.eu/sites/default/files/styles/responsive_no_crop/public/2022-01/QMCkl%20code.png?itok=UvOUClA5" width=200>

This repository is a work-in-progress GPU addon to the [QMCkl library](https://github.com/TREX-CoE/qmckl). It provides alternatives to the standard QMCkl functions to perform computations on GPU. The library is based on OpenMP offload and OpenACC, as to support a wide variety of compilers and targets.



# Basic installation

The project uses Autotools :

```
bash autogen.sh
./configure [GPU functions-related arguments & other arguments]
make
make install
```

The only other requirement is a compiler toolchain that supports OpenMP or OpenACC offloading for your chosen target(s). Linking with QMCkl to build QMCkl GPU is not needed.


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

As a general rule, each type of GPU functions comes with its own suffix. In order to call such functions, one can simply append the corresponding suffix to the standard `get` function name. Here is a table summarizing the different function types, their suffixes, and the configure flag used to enable them: 

| Function type | Suffix | Configure flag |
| ----------- | ----------- | ----------- |
| Simple OpenMP offload | `_omp_offload` | `--enable-openmp` |
| Simple OpenACC offload | `_acc_offload` | `--enable-openacc` |
| OpenMP + device pointers | `_device` | `--enable-device` |

**Note:** Using at least one of those arguments is mandatory, as doing otherwise would result in an empty library. If none is specified, `--enable-openmp` will be toggled on by default.


# Basic usage

This library is an addon meant to be used on top of the main [QMCkl library](https://github.com/TREX-CoE/qmckl). While building QMCkl GPU doesn't require any linking with QMCkl, executables calling QMCkl GPU should be linked with

```
-lqmckl -lqmckl_gpu
```

Both `qmckl.h` and `qmckl_gpu.h` should also be included in C codes.

**Notes :** 
- QMCkl GPU should be linked with an HPC-enabled QMCkl (configure QMCkl with `--enable-hpc`, see the [QMCkl README](https://github.com/TREX-CoE/qmckl/blob/master/README.md))
- While it is possible to use different compilers to build QMCkl, QMCkl GPU and your final executable, doing so might result in some missing symbols when linking your final executable. In that case, you would have to manually tell your linker where to find those symbols (which might come from different libc or libm implementations being simultaneously, for instance).

The next sections provide explanations on how this library should be used in your code.


## Difference between "offload" and "device pointers" functions

The library provides to main function types : the **basic "offload" functions**, which handle all data transfers between the CPU and the GPU internally, as well as **the more advanced "device pointers" functions**, which work all the way with memory allocatd on the GPU. The remainder of this section provides more details on the inner working of these two functions types, and outlines their advantages and inconvenients.

When calling an **offload function**, everything gets executed as in the CPU version of QMCkl, at the exception of some computation kernels that are offloaded "on the fly". This means that memory is allocated and transferred to the GPU just before the computation, only to be transferred back to the CPU and freed just after. This makes them easy to use, as **they can be called just like the standard CPU functions**, but at the cost of performance when calling these kernels repeatedly (such as in a QMC simulation), because it **repeats some costly and avoidable data transfers**.

This is why the use of the **device functions** is advised to get the **full benefits of GPU offloading**. These functions work directly and (almost) exclusively on the GPU memory, referenced by "device pointers". This way, memory transfers can be greatly reduced. Typically the entry data set would need to be transferred once from CPU to GPU at the beginning of a QMC simulation, and the results would be transferred once at the very end. The main dowside is that **this will typically require the user to perform a few GPU allocations themself as well as explicitly using the device version of each QMCkl function they might need**, which requires a bit more care.

**TODO** More detailed explanations on device pointer functions use, dedicated .md documentation files ?

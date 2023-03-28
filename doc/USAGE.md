# Usage

## Basic usage

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

which performs the same task on a GPU, using memory allocated on it.

Note that those functions accept a `qmckl_context_device` as parameter, as opposed to the `qmckl_context` used in the CPU library. A `qmckl_context_device` is associated to a single GPU.

**Notes :** 
- If you plan on using QMCkl CPU and GPU together : while it is possible to use different compilers to build QMCkl, QMCkl GPU and your final executable, doing so might result in some missing symbols when linking your final executable. In that case, you would have to manually tell your linker where to find those symbols (which might come from different libc implementations, for instance). 

The next sections provides more advanced explanations on how this library should be used compared to the standard CPU QMCkl.


## Using the "device" context and functions

The library works on a variant of the `qmckl_context` type: `qmckl_context_device`. It has the same structure as the default context, but contains a few additional and GPU-related informations, such as the OpenMP/OpenACC ID of the device you are allocating memory on and doing your offload to. This type comes with its own variants of the classic QMCkl functions. These functions work directly and (almost) exclusively on GPU memory, referenced by "device pointers". As a rule of thumb, you should always pass such device pointers as arguments of these functions (except for functions handling transfers between CPU and GPU). 

With this approach, memory transfers are greatly reduced compared to a "naive" approach where memory is transferred back and forth everytime we need to perform a computation on GPU. Typically the entry data set would need to be allocated/transferred once from CPU to GPU at the beginning of a QMC simulation, and the results would be transferred on demand only when they are needed elsewhere. The main dowside is that **this will likely require the user to perform a few GPU allocations themself**, which requires a bit more care than using the CPU library.


## Manipulating a GPU context

**TODO**

## Transferring data between CPU contextxs and GPU contexts

**TODO**

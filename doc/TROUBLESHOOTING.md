# Troubleshooting

During the development and testing of the library, we encountered some compiler related issues. This section contains a list of known issues and fixes if any. 


## Shared libraries 

The use of shared libraries can cause multiple runtime issues. On LLVM in particular, shared libaries are not supported at all. Running an executable linked against an LLVM-built QMCkl GPU would result in a runtime issue similar to this one : 

```
Libomptarget error: Host ptr 0x00007fedd3aa7380 does not have a matching target pointer.
Libomptarget error: Consult https://openmp.llvm.org/design/Runtimes.html for debugging options.
qmckl_trexio_omp.c:827:1: Libomptarget fatal error 1: failure of target construct while offloading is mandatory
```

Using GCC, we encoutered issues such as : 

```
/usr/bin/ld: tests/.libs/test_qmckl_ao: hidden symbol `__offload_funcs_end' in /var/lib/spack/packages/linux-endeavourosrolling-skylake/gcc-11.1.0/gcc-12.1.0-5radiwu4xtwuhqyxjhs4ue4qttn4gueq/lib/gcc/x86_64-pc-linux-gnu/12.1.0/crtoffloadend.o is referenced by DSO
/usr/bin/ld: final link failed: bad value
collect2: error: ld returned 1 exit status
```

While getting shared libaries to work is possible with compilers other than LLVM, it has been a common cause of runtime issues that you should be aware of. The best workaround is simply to link executables against the static version of the library. 

Generation of the shared library can be disabled at configure time using the `--disable-shared` option. Also, generation of the static library can be forced (even though this should be the default behaviour) with the `--enable-static` option.


## [nvc] pgi undefined reference when building executable

When building an executable with an `nvc` built QMCkl GPU (when doing a `make check` for instance), you might run into this issue, or similar : 

```
/usr/bin/ld: /opt/nvidia/hpc_sdk/Linux_x86_64/22.11/compilers/lib/libgomp.so.1: undefined reference to `__pgi_nvomp_error_cuda_noversion'
```

It seems that by itself, nvc doesn't link correctly all of its libraries, resulting in missing symbols errors. We have been able to solve the issue by manually specifying the flags needed to link those libraries  : 

```
./configure CC="nvc" LDFLAGS="-L/[path]/[to]/nvidia/hpc_sdk/Linux_x86_64/[version]/compilers/lib -lpgc -lnvf"
```

... should solve the issue when building the tests. You can do something similar when building your own app.


## [nvc] libgomp: TODO

On some systems, building the lib with nvc can result in the following runtime error (this has been observed with both OpenMP and OpenACC on different systems) : 

```
libgomp: TODO
```

Apparently, this is due to nvc linking executables with GCC's libgomp instead of its own. One way to fix this is to put the Nvidia libgomp path in `LD_LIBRARY_PATH` : 

```
export LD_LIBRARY_PATH=/[path]/[to]/nvidia/hpc_sdk/Linux_x86_64/[version]/compilers/lib:$LD_LIBRARY_PATH
```

As this affects the libraries' search path at runtime only, you should not need to recompile anything.


## [gcc] [OpenACC] lto1: internal compiler error: Segmentation fault

This error has been observed at link time with gcc (Spack GCC) 12.1.0 and OpenACC. This is apparently related to link time optimization (LTO) : 

```
lto1: internal compiler error: Segmentation fault
0xa4f83f crash_signal
	../../gcc/toplev.cc:322
0xd2b109 nvptx_record_offload_symbol
	../../gcc/config/nvptx/nvptx.cc:5911
0xd2b109 nvptx_record_offload_symbol
	../../gcc/config/nvptx/nvptx.cc:5891
0x93c54e omp_finish_file()
	../../gcc/omp-offload.cc:452
```

This issue has been fixed by explicitely enabling LTO in the build configuration : 

```
./configure CC="gcc -flto"
```

# Troubleshooting

During the development and testing of the library, we encountered some compiler related issues. This section contains a list of known issues and fixes if any. 

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


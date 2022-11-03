# QMCkl GPU: devtools

This folder contains helper scripts for developers. It will likely not be needed by an user building the library from sources.

## clang_format_all.sh

Requires a clang_format/LLVM install on your machine. Automatically formats (inplace) all your source and header files.

**Usage:**

```
[path]/[to]/clang_format_all.sh [path to qmckl_gpu]
```

## copy_qmckl_headers.sh

In order to be built, QMCkl generates some private header files (from the org files) which are not meant to be installed. QMCkl GPU also needs some of these headers, which are currently kept in the source tree and might therefore need to be updated along with the main QMCkl repository. This script will automatically copy/overwrite all the required private headers from a `qmckl` source tree to a `qmckl_gpu` source tree. 

**Usage:**

```
[path]/[to]/copy_qmckl_headers.sh [path to qmckl_gpu] [path to qmckl]
```

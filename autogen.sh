#!/bin/sh
install_qmckl()
{
	cd qmckl
	./autogen.sh
	./configure --enable-hpc --disable-doc
	make source -j 

}

test -f qmckl/configure || \
 	(git submodule update --init && install_qmckl)
autoreconf -is

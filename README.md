# libfasta
An open source library for parsing FASTA files

## Features
 * Supports files with multiple FASTA sequences, i.e. multi-FASTA
 * Support for reading sequence data into memory only on demand
 * Capable of indexing the FASTA files for faster repeated processing
 * API for processing user-defined coding sequences
 * No external dependencies

## Compilation

To compile the sources inside a release tarball, use the following commands:

    $ ./configure
    $ make
    $ sudo make install

If you want to compile the sources in a cloned repository, you'll have to
generate the configure script and other files, which aren't part of the repository.
To do that, use the `autogen.sh` script:

    $ ./autogen.sh

After completion of this command, you should be able to run the above mentioned commands
to compile and install the library

## Using the library in your project

To access the C API defined in the installed header file, use the following include statement in
your C/C++ source file:

    #include <libfasta/fasta.h>

When you install the library from a release tarball or a distribution package, it should also
install a *pkg-config* file which makes it easier to link to the library at compile-time, e.g.:

    $ gcc -o program $(pkg-config libfasta --cflags --libs) program.c

If you are using an Autotools or CMake based project, you can use the [**PKG_CHECK_MODULES** macro](https://autotools.io/pkgconfig/pkg_check_modules.html)
or [**PKG_CHECK_MODULES** CMake command](http://www.cmake.org/cmake/help/cmake2.6docs.html#module:FindPkgConfig) respectively.

Otherwise, you'll have to explicitly specify that you want to link to the library like this:

    $ gcc -o program -lfasta program.c

If you installed the files into a non-standard location, set the include and library path:

    $ gcc -o program -I/custom/location/include -L/custom/location/lib -lfasta program.c

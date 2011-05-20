/*
 * Copyright 2011 Daniel Kopecek. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are
 * permitted provided that the following conditions are met:
 *
 *    1. Redistributions of source code must retain the above copyright notice, this list of
 *       conditions and the following disclaimer.
 *
 *    2. Redistributions in binary form must reproduce the above copyright notice, this list
 *       of conditions and the following disclaimer in the documentation and/or other materials
 *       provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY Daniel Kopecek ``AS IS'' AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
 * FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL Daniel Kopecek OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 * Author: Daniel Kopecek <xkopecek@fi.muni.cz>
 *
 */
#ifndef HELPERS_H
#define HELPERS_H

/**
 * Memory allocation helpers
 */
#include <stdlib.h>

/**
 * Allocate storage for a variable of type `T'.
 */
#define alloc_type(T)                           \
        ((T *) malloc(sizeof(T)))

/**
 * Allocate an array of variables of type `T'.
 */
#define alloc_array(T, count)                   \
        ((T *) malloc(sizeof(T) * (size_t)(count)))

#define realloc_array(ptr, T, count)            \
        ((T *) realloc(ptr, sizeof(T) * (size_t)(count)))

/**
 * Create a symbol name with the current line as part
 * of its name making it a unique name in the current
 * file (provided it's used only once per line in a
 * file).
 */
#define __lsym(s) s ## __LINE__

/**
 * Various useful functions, macros, ...
 */
#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/stat.h>

/*
 * Wrappers for several file operations used for
 * easier portability.
 */
uint64_t file_get_offset(FILE *fp);
int file_set_offset(FILE *fp, uint64_t offset);
int file_get_stat(FILE *fp, struct stat *st);

/**
 * Save errno, execute the block, restore errno.
 */
#define errno_protect \
        for(int __lsym(errno) = errno, __lsym(exec) = 1; __lsym(exec) == 1; __lsym(exec) = 0, errno = __lsym(errno))

/**
 * Output to stderr in debug mode, otherwise do nothing.
 */
#ifndef NDEBUG
# define dP(...) fprintf(stderr, __VA_ARGS__)
#else
# define dP(...) while(0)
#endif

#endif /* HELPERS_H */

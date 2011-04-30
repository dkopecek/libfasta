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

#define __lsym(s) s ## __LINE__

/**
 * Various useful functions, macros, ...
 */

#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/stat.h>

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

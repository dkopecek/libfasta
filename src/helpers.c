#define _XOPEN_SOURCE /* for fileno */
#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include <unistd.h>
#include "helpers.h"

uint64_t file_get_offset(FILE *fp)
{
        return ((uint64_t)lseek(fileno(fp), 0, SEEK_CUR));
}

int file_set_offset(FILE *fp, uint64_t offset)
{
        if (lseek(fileno(fp), offset, SEEK_SET) == (off_t)(-1))
                return (-1);
        else
                return (0);
}

int file_get_stat(FILE *fp, struct stat *st)
{
        if (fstat(fileno(fp), st) != 0)
                return (-1);
        else
                return (0);
}

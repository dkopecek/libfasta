#include <sys/types.h>
#include <stdint.h>
#include "fasta.h"

FASTA *fasta_open(const char *path, uint32_t options)
{
	return (NULL);
}

uint32_t fasta_count(FASTA *fa)
{
	return (0);
}

FASTAREC *fasta_read(FASTA *fa)
{
	return (NULL);
}

int fasta_write(FASTA *fa, FASTAREC *farec)
{
	return (-1);
}

int fasta_rewind(FASTA *fa)
{
	return (-1);
}

int fasta_seeko(FASTA *fa, off_t off, int whence)
{
	return (-1);
}

off_t fasta_tello(FASTA *fa)
{
	return (-1);
}

void *fasta_apply(FASTA *fa, (void *)(*func)(FASTAREC *, void *), uint32_t options, void *funcarg)
{
	return (-1);
}

void fasta_close(FASTA *fa)
{
	return;
}

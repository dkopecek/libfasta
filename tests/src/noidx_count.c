#include <stdio.h>
#include <fasta.h>
#include <libgen.h>

int main(int argc, char *argv[])
{
	FASTA       *fa;
	FASTA_rec_t *farec;

	if (argc != 2) {
		fprintf(stderr, "Usage: %s <fasta-file>\n", basename(argv[0]));
		return (1);
	}

	fa = fasta_open(argv[1], FASTA_READ|FASTA_ONDEMSEQ, NULL);

	if (fa != NULL) {
		printf("Total records: %u\n", fasta_count(fa));
		fasta_close(fa);
	} else {
		printf("fasta_open => NULL\n");
		return (2);
	}

	return (0);
}

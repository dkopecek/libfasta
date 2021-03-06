#include <config.h>
#include <stdio.h>
#include <fasta.h>
#include <libgen.h>

int main(int argc, char *argv[])
{
	FASTA *fa;

	if (argc != 2) {
		fprintf(stderr, "Usage: %s <fasta-file>\n", basename(argv[0]));
		return (1);
	}

	fa = fasta_open(argv[1],
			FASTA_READ|FASTA_ONDEMSEQ|
			FASTA_USEINDEX|FASTA_GENINDEX|FASTA_CHKINDEX, NULL);

	if (fa != NULL) {
		printf("Total records: %u\n", fasta_count(fa));
		fasta_close(fa);
	} else {
		printf("fasta_open => NULL\n");
		return (2);
	}

	return (0);
}

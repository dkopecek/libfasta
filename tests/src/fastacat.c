#include <stdio.h>
#include <stdlib.h>
#include <fasta.h>
#include <libgen.h>

int main(int argc, char *argv[])
{
	FASTA       *fa;
	FASTA_rec_t *farec;

	if (argc < 2) {
		fprintf(stderr, "Usage: %s <fasta-file> [<seq #>]\n", basename(argv[0]));
		return (1);
	}

	fa = fasta_open(argv[1], FASTA_READ|FASTA_ONDEMSEQ, NULL);

	if (fa == NULL) {
		fprintf(stderr, "fasta_open(%s) => NULL\n", argv[1]);
		return (2);
	}

	if (argc > 2) {
		if (fasta_seeko(fa, atoi(argv[2]), SEEK_SET) != 0) {
			fprintf(stderr, "fasta_seeko(%s) != 0\n", argv[2]);
			return (3);
		}

		farec = fasta_read(fa, NULL, FASTA_INMEMSEQ|FASTA_CSTRSEQ, NULL);

		if (farec == NULL) {
			fprintf(stderr, "fasta_read => NULL\n");
			return (4);
		}

		printf("%s", farec->seq_mem);
		fasta_rec_free(farec);
	} else {
		while ((farec = fasta_read(fa, NULL, FASTA_INMEMSEQ|FASTA_CSTRSEQ, NULL)) != NULL) {
			printf("%s", farec->seq_mem);
			fasta_rec_free(farec);
		}
	}

	fasta_close(fa);

	return(0);
}

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

	fa = fasta_open(argv[1],
			FASTA_READ|FASTA_ONDEMSEQ|
			FASTA_USEINDEX|FASTA_GENINDEX|FASTA_CHKINDEX, NULL);

	if (fa != NULL) {
		printf("Total records: %u\n", fasta_count(fa));
	} else {
		printf("fasta_open => NULL\n");
		return (2);
	}

	while ((farec = fasta_read(fa, NULL, FASTA_INMEMSEQ|FASTA_CSTRSEQ, NULL)) != NULL) {
		printf("=> sequence: %s\n", farec->seq_mem);
		printf("     length: %lu\n", farec->seq_len);
		printf("     rawlen: %lu\n", farec->seq_rawlen);
		printf("      lines: %u\n", farec->seq_lines);
		fasta_rec_free(farec);
	}

	fasta_close(fa);

	return (0);
}

#include <stdio.h>
#include <fasta.h>
#include <libgen.h>

int main(int argc, char *argv[])
{
	FASTA       *fa;
	FASTA_rec_t *farec;
	atrans_t    *tr;

	if (argc != 2) {
		fprintf(stderr, "Usage: %s <fasta-file>\n", basename(argv[0]));
		return (1);
	}

	tr = atrans_new(8, 8, 0, 0);

	tr->tr_letter_s2d['A'] = '1';
	tr->tr_letter_s2d['a'] = '1';
	tr->tr_letter_s2d['T'] = '2';
	tr->tr_letter_s2d['t'] = '2';
	tr->tr_letter_s2d['C'] = '3';
	tr->tr_letter_s2d['c'] = '3';
	tr->tr_letter_s2d['G'] = '4';
	tr->tr_letter_s2d['g'] = '4';

	fa = fasta_open(argv[1], FASTA_READ|FASTA_ONDEMSEQ, tr);

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
	atrans_free(tr);

	return (0);
}

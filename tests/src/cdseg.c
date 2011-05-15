#include <fasta.h>
#include <libgen.h>
#include <inttypes.h>

int main(int argc, char *argv[])
{
	FASTA       *fa;
	FASTA_rec_t *farec;
        FASTA_CDS_t *faseg;

	if (argc != 2) {
		fprintf(stderr, "Usage: %s <fasta-file>\n", basename(argv[0]));
		return (1);
	}

	fa = fasta_open(argv[1],
			FASTA_READ|FASTA_ONDEMSEQ|FASTA_NASEQ, NULL);
        //FASTA_USEINDEX|FASTA_GENINDEX|FASTA_CHKINDEX, NULL);

	if (fa == NULL) {
		printf("fasta_open(%s) => NULL\n", argv[1]);
		return (2);
	}

	while ((farec = fasta_read(fa, NULL, FASTA_INMEMSEQ|FASTA_CSTRSEQ|FASTA_MAPCDSEG, NULL)) != NULL) {
                register size_t i;
                /*
                 * print:
                 *   seq_len cdseg_count
                 *   cdseg[0].a cdseg[0].b
                 *   ...
                 *   ...
                 *   cdseg[n].a cdseg[0].b
                 */
                printf("%"PRIu64" %zu\n", farec->seq_len, farec->cdseg_count);

                for (i = 0; i < farec->cdseg_count; ++i)
                        printf("%"PRIu64" %"PRIu64" %"PRIu64"\n",
                               farec->cdseg[i].a, farec->cdseg[i].b, farec->cdseg[i].b - farec->cdseg[i].a + 1);

                printf("-\n");
		fasta_rec_free(farec);
	}

	fasta_close(fa);

	return (0);
}

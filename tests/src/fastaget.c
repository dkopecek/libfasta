#define _XOPEN_SOURCE 800
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fasta.h>
#include <libgen.h>

int main(int argc, char *argv[])
{
        FILE *fp;
	FASTA *fa;
	FASTA_rec_t *farec;

        char *line;
        size_t llen;

	if (argc != 3) {
		fprintf(stderr, "Usage: %s <fasta-file> <positions>\n", basename(argv[0]));
		return (1);
	}

        fp = fopen(argv[2], "r");
	fa = fasta_open(argv[1], FASTA_READ|FASTA_ONDEMSEQ|FASTA_USEINDEX, NULL);

	if (fa == NULL) {
		fprintf(stderr, "fasta_open(%s) => NULL\n", argv[1]);
		return (2);
	}

        /*if (fasta_seeko(fa, seqn, SEEK_SET) != 0) {
                fprintf(stderr, "fasta_seeko(%s) != 0\n", argv[2]);
                return (3);
                }*/

        farec = fasta_read(fa, NULL, FASTA_INMEMSEQ, NULL);

        if (farec == NULL) {
                fprintf(stderr, "fasta_read => NULL\n");
                return (4);
        }

        while (getline(&line, &llen, fp) != -1) {
                char *b, *e;

                b = strtok(line, " ");
                e = strtok(NULL, " ");

                if (b == NULL || e == NULL)
                        continue;

                printf("%.*s\n", (int)(atoi(e) - atoi(b) + 1), farec->seq_mem + atoi(b));
        }

        fasta_rec_free(farec);
	fasta_close(fa);

	return(0);
}

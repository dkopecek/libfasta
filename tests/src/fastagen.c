#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <libgen.h>
#include <unistd.h>
#include <time.h>

static int rnd(int min, int max)
{
	return (min + (random() % ((max - min) + 1)));
}

typedef struct {
	char  *chars;
	size_t nchar;
} alphabet;

const alphabet aarray[] = {
	{ "ATCG",     4 },
	{ "atcg",     4 },
	{ "ATCGatcg", 8 }
};

int main(int argc, char *argv[])
{
	int seqs;
	int seqn, i;
	int seql, l, s;
	int lenv;
	int line;
	int linv;
	int lins;

	uint8_t *seq, *seq_orig;

	if (argc != 8) {
		fprintf(stderr, "Usage: %s <seq-seed> <line-seed> <numseq> <seqlen> <seqlen-variability> <linew> <linew-variability>\n", basename(argv[0]));
		return (1);
	}

	seqs = atoi(argv[1]);
	lins = atoi(argv[2]);
	seqn = atoi(argv[3]);
	seql = atoi(argv[4]);
	lenv = atoi(argv[5]);
	line = atoi(argv[6]);
	linv = atoi(argv[7]);

	fprintf(stderr, "params: %u %u %u %u %u %u %u\n", seqs, lins, seqn, seql, lenv, line, linv);

	for (i = 0; i < seqn; ++i) {
		if (seqs == 0)
			srandom((unsigned long)(getpid() ^ time(NULL) ^ clock()) + i);
		else
			srandom(seqs + i);
		/*
		 * generate length of the current sequence from the sequence length (seql)
		 * and sequence length variability (lenv) parameters
		 */
		l = seql + (random() % 2 ? rnd(0, lenv) : -rnd(0, lenv));

		/*
		 * emit the sequence header
		 */
		printf(">SEQUENCE_%u; length=%u\n", i, l);

		/*
		 * allocate memory for the sequence and generate the sequence
		 */
		seq = seq_orig = malloc(sizeof(uint8_t) * l);

		for (s = 0; s < l; ++s)
			seq[s] = aarray[0].chars[rnd(0, aarray[0].nchar - 1)];

		if (lins == 0)
			srandom((unsigned long)(getpid() ^ time(NULL) ^ clock()) + i);
		else
			srandom(lins + i);

		if (linv == 0) {
			/*
			 * constant line width
			 */
			int c;

			for (c = l / line; c > 0; --c) {
				for (s = line; s > 0; --s)
					putc(*seq++, stdout);

				putc('\n', stdout);
			}

			for (s = l % line; s > 0; --s)
				putc(*seq++, stdout);

			putc('\n', stdout);
		} else {
			/*
			 * variable line width
			 */
			int c;
			int w;

			while (l > 0) {
				/*
				 * Cycle the minimum number of times (i.e. cycle exactly the number of
				 * times it would require to cycle, if the line width would be maximized)
				 */
				for (c = l / (line + linv); c > 0; --c) {
					w = line + (random() % 2 ? rnd(0, linv) : -rnd(0, linv));

					for (s = w; s > 0; --s)
						putc(*seq++, stdout);

					putc('\n', stdout);

					l -= w;
				}

				if (l <= (line + linv)) {
					for (s = l; s > 0; --s)
						putc(*seq++, stdout);

					putc('\n', stdout);
					l = 0;
				}
			}
		}

		free(seq_orig);
	}

	return (0);
}

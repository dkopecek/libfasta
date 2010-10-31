#include <ctype.h>
#include <stdbool.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdint.h>
#include <limits.h>
#include <errno.h>
#include "assume.h"
#include "debug.h"
#include "sm_alloc.h"
#include "fasta.h"
#include "trans.h"

/*
 * XXX: not efficient... rewrite to bitmasks, separate
 *      alphabets + union
 */
static bool issequence(int ch)
{
	switch(toupper(ch)) {
	case 'A':
	case 'C':
	case 'G':
	case 'T':
	case 'U':
	case 'R':
	case 'Y':
	case 'K':
	case 'M':
	case 'S':
	case 'W':
	case 'B':
	case 'D':
	case 'H':
	case 'V':
	case 'N':
	case 'X':
	case '-':
		/* nucleic acid */
		return (true);
	case 'E':
	case 'F':
	case 'I':
	case 'L':
	case 'O':
	case 'P':
	case 'Q':
	case 'Z':
	case '*':
		/* amino acid */
		return (true);
	}

	return (false);
}

static int __fasta_read2(FILE *fp, FASTA_rec_t *dst, atrans_t *atr)
{
	return (-1);
}

static int __fasta_read1(FILE *fp, FASTA_rec_t *dst, atrans_t *atr)
{
	size_t alloc_size;
	char  *buffer;
	size_t buflen;

	register uint32_t l, n;
	register uint64_t i;

	if (fseeko(fp, dst->seq_start, SEEK_SET) != 0) {
		_D("Failed to seek to position %zu in %p\n", dst->seq_start, fp);
		return (-1);
	}

	_D("atr=%p\n", atr);

	if (atr != NULL)
		alloc_size = atrans_s2d_size(atr, dst->seq_len);
	else
		alloc_size = dst->seq_len;

	if (dst->flags & FASTA_CSTRSEQ)
		++alloc_size;

	_D("alloc_size=%zu\n", alloc_size);

	i = 0;
	dst->seq_mem = sm_alloc(alloc_size);

	if (atr != NULL) {
		bzero(dst->seq_mem, alloc_size);

		buffer = sm_alloc(sizeof(uint8_t) * dst->seq_linew);
		buflen = dst->seq_linew;

		for (l = dst->seq_lines - 1; l > 0; --l) {
			/*
			 * Read line into buffer
			 */
			if (fread(buffer, 1, buflen, fp) != buflen) {
				/* fail */
				sm_free(buffer);
				sm_free(dst->seq_mem);
				dst->seq_mem = NULL;

				return (-1);
			}

			/*
			 * Copy/Translate the buffer into seq_mem
			 */
			for (n = 0; n < buflen; ++n)
				atrans_letter_s2d(atr, buffer[n], i++, (uint8_t *)dst->seq_mem);

			getc(fp); /* skip the new-line */
		}

		buflen = dst->seq_lastw > 0 ? dst->seq_lastw : dst->seq_linew;
		/* no need to reallocate the buffer since lastw < linew */

		if (fread(buffer, 1, buflen, fp) != buflen) {
			/* fail */
			sm_free(buffer);
			sm_free(dst->seq_mem);
			dst->seq_mem = NULL;

			return (-1);
		}

		for (n = 0; n < buflen; ++n)
			atrans_letter_s2d(atr, buffer[n], i++, (uint8_t *)dst->seq_mem);

		sm_free(buffer);
	} else {
		buflen = dst->seq_linew;

		for (l = dst->seq_lines - 1; l > 0; --l) {
			/*
			 * Read line into buffer
			 */
			if (fread(dst->seq_mem + i, 1, buflen, fp) != buflen) {
				/* fail */
				sm_free(dst->seq_mem);
				dst->seq_mem = NULL;

				return (-1);
			}

			i += buflen;

			getc(fp); /* skip the new-line */
		}

		buflen = dst->seq_lastw > 0 ? dst->seq_lastw : dst->seq_linew;

		if (fread(dst->seq_mem + i, 1, buflen, fp) != buflen) {
			/* fail */
			sm_free(dst->seq_mem);
			dst->seq_mem = NULL;

			return (-1);
		}

		i += buflen;
	}

	if (dst->flags & FASTA_CSTRSEQ)
		dst->seq_mem[alloc_size - 1] = '\0';

	return (0);
}

static int __fasta_read0(FILE *fp, FASTA_rec_t *dst, uint32_t options, atrans_t *atr)
{
	int      ch;
	char    *buffer;
	uint32_t buflen;
	char    *buftok;
	uint32_t toklen;

	uint32_t plinew; /* previous line width */
	uint32_t clinew; /* current line width */

	uint32_t i;

	assume_d(fp  != NULL, -1);
	assume_d(dst != NULL, -1);

	dst->flags   = 0;
	dst->hdr     = NULL;
	dst->hdr_mem = NULL;
	dst->seq_mem = NULL;

	while (!feof(fp)) {
		ch = getc(fp);

		if (ch != '>')
			return (-1);

		/*
		 * Read all headers
		 */
		dst->hdr_start = (uint64_t)ftello(fp);
		dst->hdr_len   = 0;
		dst->hdr_cnt   = 1;

		buflen = 1024;
		buffer = sm_alloc(sizeof(char) * buflen);

		do {
			if (dst->hdr_len >= buflen) {
				buflen += 1024;
				buffer  = sm_realloc(buffer, sizeof(char) * buflen);
			}

			buffer[dst->hdr_len++] = ch = getc(fp);

			if (ch == 0x01) /* ^A */
				++dst->hdr_cnt;

		} while(ch != '\n');

		buffer[--dst->hdr_len] = '\0';
		buffer = sm_realloc(buffer, sizeof(char) * (dst->hdr_len + 1));
		dst->hdr_mem = buffer;


		_D(" Read header: \"%s\"\n", buffer);
		_D("Header count: %u\n", dst->hdr_cnt);

		/*
		 * Parse headers
		 */
		dst->hdr    = sm_alloc(sizeof(FASTA_rechdr_t) * dst->hdr_cnt);
		dst->flags |= FASTA_REC_FREEHDR;

		i = 0;

		while ((buftok = strsep(&buffer, "\x01")) != NULL) {
			/*
			 * Sanity check
			 */
			if (i >= dst->hdr_cnt) {
				_D("Insane value(s): i=%u, dst->hdr_cnt=%u\n", i, dst->hdr_cnt);
				goto fail;
			}

			/*
			 * Parse the header using the SeqID parser
			 */
			if ((dst->hdr[i].seqid_fmt = SeqID_parse(buftok, strlen(buftok),
								 &dst->hdr[i].seqid)) == SEQID_ERROR)
			{
				_D("SeqID returned an error: h=\"%s\" l=%zu\n", buftok, toklen);
				goto fail;
			}

			++i;
		}

		/*
		 * Analyze/Read the sequence
		 */
		if (options & FASTA_INMEMSEQ) {
			/* in-memory */
		} else {
			/*
			 * Analyze the sequence without storing it into main memory.
			 */
			bool linew_diff = false, linew_update = true;

			plinew = 0;
			clinew = 0;

			dst->seq_start  = (uint64_t)ftello(fp);

			dst->seq_len    = 0;
			dst->seq_rawlen = 0;
			dst->seq_lines  = 0;

			dst->seq_linew  = 0;
			dst->seq_lastw  = 0;

			/*
			 * Read in the first line of the sequence.
			 */
			for (;;) {
				ch = getc(fp);

				if (feof(fp)) {
					if (dst->seq_len == 0) {
						_D("Unexpected EOF: got header, but no sequence data\n");
						goto fail;
					} else {
						goto finalize_seq;
					}
				}

				++dst->seq_rawlen;

				if (issequence(ch)) {
					++plinew;
					++dst->seq_len;
				} else {
					if (ch == '\n' && dst->seq_len > 0) {
						++dst->seq_lines;
						break;
					} else if (isspace(ch) && dst->seq_len == 0) {
						dst->seq_start += dst->seq_rawlen;
						dst->seq_rawlen = 0;
					}
				}
			}

			/*
			 * Read the rest of the sequence lines.
			 */
			for (;;) {
				ch = getc(fp);

				if (feof(fp)) {
					if (clinew > 0 && !linew_diff)
						++dst->seq_lines;
					break;
				}

				++dst->seq_rawlen;

				if (issequence(ch)) {
					if (linew_update) {
						if (linew_diff) {
							clinew = 0;
							plinew = 0;
							linew_update = false;
							linew_diff   = false;
						} else
							++clinew;
					}

					++dst->seq_len;
				} else {
					switch (ch) {
					case '\n':
						if (linew_update && clinew > 0) {
							if (clinew != plinew)
								linew_diff = true;
							else {
								plinew = clinew;
								clinew = 0;
							}
						}

						++dst->seq_lines;

						break;
					case  '>':
						if (clinew == 0 || linew_diff == true) {
							ungetc('>', fp);
							goto finalize_seq;
						} else {
							_D("Unexpected '>': allowed only at the beginning of a line\n");
							goto fail;
						}
						break;
					}
				}
			}
		}
	}

finalize_seq:
	dst->seq_linew = plinew;
	dst->seq_lastw = clinew;
	dst->flags    |= FASTA_REC_MAGICFL;

	/*
	 * ret=0 - complete
	 * ret>0 - EOF (ret=1)
	 * ret<0 - error
	 */
	return (feof(fp) ? 1 : 0);
fail:
	if (dst->hdr != NULL)
		sm_free(dst->hdr);
	if (dst->hdr_mem != NULL)
		sm_free(dst->hdr_mem);
	if (dst->seq_mem != NULL)
		sm_free(dst->seq_mem);

	return (-1);
}

FASTA *fasta_open(const char *path, uint32_t options, atrans_t *atr)
{
	char   idx_path[PATH_MAX + 1];
	FASTA *fa;
	struct stat st;

	assume_r(path != NULL, NULL);

	fa             = sm_talloc(FASTA);
	fa->fa_options = options;
	fa->fa_path    = strdup(path);
	fa->fa_seqFP   = NULL;
	fa->fa_idxFP   = NULL;
	fa->fa_record  = NULL;
	fa->fa_rindex  = 0;
	fa->fa_rcount  = 0;
	fa->fa_atr     = atr;

	fa->fa_seqFP = fopen(path, "r");

	if (fa->fa_seqFP == NULL) {
		_D("Can't open the sequence file: %s\n", path);
		goto fail;
	}

	if (fstat(fileno(fa->fa_seqFP), &st) != 0) {
		int e = errno;
		_D("Failed to get stat information: fp=%p (fd=%d, path=%s), errno=%u, %s\n",
		   fa->fa_seqFP, fileno(fa->fa_seqFP), path, e, strerror(e));
		goto fail;
	}

	if (options & FASTA_USEINDEX) {
		if (strlen(path) + strlen(FASTA_INDEX_EXT) < (sizeof idx_path/sizeof(char))) {
			strcpy(idx_path, path);
			strcat(idx_path, FASTA_INDEX_EXT);
		} else {
			_D("Index path exceeds systems limits: idx_path_len=%zu, limit=%zu\n",
			   strlen(path) + strlen(FASTA_INDEX_EXT), sizeof idx_path/sizeof(char));
			goto fail;
		}

		fa->fa_idxFP = fopen(idx_path,
				     (options & FASTA_GENINDEX ? "r" : "r+"));

#if 0
		/* generate headers from the index */

		if (options & FATA_CHKINDEX_SLOW) {
			/* slow check */
		} else if (options & FASTA_CHKINDEX_FAST) {
			/* fast check */
		}
#endif
		/* continue */
		abort ();
	} else {
		register uint32_t i;
		int r;

		/*
		 * generate headers from the sequence file
		 */
		i = 0;
		fa->fa_rcount = 0;

		do {
			if (i >= fa->fa_rcount) {
				_D("=> pre-alloc: fa_rcount=%u\n", fa->fa_rcount);

				if (i == 0)
					fa->fa_rcount = 8;
				else if (i < 65535)
					fa->fa_rcount <<= 1;
				else
					fa->fa_rcount += 1024;

				fa->fa_record = sm_realloc(fa->fa_record, sizeof(FASTA_rec_t) * fa->fa_rcount);

				_D("<= pre-alloc: fa_rcount=%u\n", fa->fa_rcount);
			}

			_D("Reading sequence #%u\n", i);
		} while ((r = __fasta_read0(fa->fa_seqFP, fa->fa_record + i++, options, fa->fa_atr)) == 0);

		if (r < 0) {
			_D("An error ocured while reading the file \"%s\"\n", fa->fa_path);
			goto fail;
		}

		fa->fa_rcount = i;
		fa->fa_record = sm_realloc(fa->fa_record, sizeof(FASTA_rec_t) * fa->fa_rcount);
	}

	return (fa);
fail:
	for (; fa->fa_rcount > 0; --fa->fa_rcount)
		fasta_rec_free(fa->fa_record + fa->fa_rcount - 1);

	sm_free(fa);

	return (NULL);
}

uint32_t fasta_count(FASTA *fa)
{
	return (fa->fa_rcount);
}

FASTA_rec_t *fasta_read(FASTA *fa, FASTA_rec_t *dst, uint32_t flags, atrans_t *atr)
{
	FASTA_rec_t *farec;

	assume_d(fa != NULL, NULL);

	if (fa->fa_rindex >= fa->fa_rcount)
		return (NULL);

	if (atr == NULL)
		atr = fa->fa_atr;

	if (dst == NULL) {
		if (flags & FASTA_RAWREC)
			farec = fa->fa_record + fa->fa_rindex;
		else {
			farec = sm_talloc(FASTA_rec_t);
			memcpy(farec, fa->fa_record + fa->fa_rindex, sizeof(FASTA_rec_t));
			farec->flags = FASTA_REC_MAGICFL | FASTA_REC_FREEREC;
		}
	} else {
		farec = dst;
		memcpy(farec, fa->fa_record + fa->fa_rindex, sizeof(FASTA_rec_t));
		farec->flags = FASTA_REC_MAGICFL;
	}

	if (flags & FASTA_CSTRSEQ)
		farec->flags |= FASTA_CSTRSEQ;

	if ((flags & FASTA_INMEMSEQ) && farec->seq_mem == NULL) {
		farec->flags |= FASTA_REC_FREESEQ;

		if (farec->seq_linew != 0) {
			/*
			 * all the lines that form the sequence are of equal length
			 */
			if (__fasta_read1(fa->fa_seqFP, farec, atr) != 0) {
				/* fail */
			}
		} else {
			/*
			 * the lines have variable length, we have to look for new-lines
			 */
			if (__fasta_read2(fa->fa_seqFP, farec, atr) != 0) {
				/* fail */
			}
		}
	}

	++fa->fa_rindex;

	return (farec);
}

int fasta_write(FASTA *fa, FASTA_rec_t *farec)
{
	return (-1);
}

void fasta_rec_free(FASTA_rec_t *farec)
{
	if (farec->flags & FASTA_REC_FREEHDR) {
		sm_free(farec->hdr_mem);
		sm_free(farec->hdr);
		farec->hdr     = NULL;
		farec->hdr_cnt = 0;
	}

	if (farec->flags & FASTA_REC_FREESEQ) {
		sm_free(farec->seq_mem);
		farec->seq_mem = NULL;
		farec->flags  &= ~(FASTA_REC_FREESEQ);
	}

	if (farec->flags & FASTA_REC_FREEREC) {
		farec->flags = 0;
		sm_free(farec);
	}

	return;
}

int fasta_rewind(FASTA *fa)
{
	fa->fa_rindex = 0;
	return (0);
}

int fasta_seeko(FASTA *fa, off_t off, int whence)
{
	return (-1);
}

off_t fasta_tello(FASTA *fa)
{
	return (fa->fa_rindex);
}

void *fasta_apply(FASTA *fa, void * (*func)(FASTA_rec_t *, void *), uint32_t options, void *funcarg)
{
	return (NULL);
}

void fasta_close(FASTA *fa)
{
	if (fa->fa_record != NULL) {
		for (; fa->fa_rcount > 0; --fa->fa_rcount)
			fasta_rec_free(fa->fa_record + fa->fa_rcount - 1);
		sm_free(fa->fa_record);
	}

	sm_free(fa->fa_path);

	if (fa->fa_seqFP != NULL)
		fclose (fa->fa_seqFP);
	if (fa->fa_idxFP != NULL)
		fclose (fa->fa_idxFP);

	sm_free(fa);
	return;
}

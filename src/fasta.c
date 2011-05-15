#include <stdio.h>
#include <ctype.h>
#include <stdbool.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdint.h>
#include <inttypes.h>
#include <limits.h>
#include <errno.h>
#include <assert.h>

#include "helpers.h"
#include "fasta.h"
#include "trans.h"
#include "crc32.h"

#ifndef PATH_MAX
# define PATH_MAX 4096
#endif

/*
 * Nucleic Acid letter bitmask
 */
const uint32_t __NA_mask[] = {
	0x00000000, 0x00002000, 0x03fc699e, 0x03fc699e,
	0x00000000, 0x00000000, 0x00000000, 0x00000000
};

/*
 * Amino Acid letter bitmask
 */
const uint32_t __AA_mask[] = {
	0x00000000, 0x00002000, 0x07fffbfe, 0x07fffbfe,
	0x00000000, 0x00000000, 0x00000000, 0x00000000
};

/*
 * Sequence letter bitmask
 */
const uint32_t __SQ_mask[] = {
	0x00000000, 0x00002000, 0x07fffbfe, 0x07fffbfe,
	0x00000000, 0x00000000, 0x00000000, 0x00000000
};

/*
 * default NA coding letter bitmask
 */
const uint32_t __NA_CD_mask[] = {
        0x00000000, 0x00000000, 0x0010008a, 0x0010008a,
        0x00000000, 0x00000000, 0x00000000, 0x00000000
};

/*
 * default AA coding letter bitmask
 */
const uint32_t __AA_CD_mask[] = {
        0x00000000, 0x00000000, 0x06fffbfe, 0x02fffbfe,
        0x00000000, 0x00000000, 0x00000000, 0x00000000
};

static bool issequence(int ch)
{
	if (ch >= 0 && ch < 256)
		return (__SQ_mask[ch / (sizeof __SQ_mask[0] * 8)]) & (1 << (ch % (sizeof __SQ_mask[0] * 8)));
	else
		return (false);
}

static bool iscodingseq(FASTA *fa, int ch)
{
	if (ch >= 0 && ch < 256)
		return (fa->fa_CDSmask[ch / (sizeof fa->fa_CDSmask[0] * 8)]) & (1 << (ch % (sizeof fa->fa_CDSmask[0] * 8)));
	else
		return (false);
}

static int __index_write(FASTA *fa, const char *idxpath)
{
	register uint32_t i;
	struct stat st;

	assert(fa != NULL);
	assert(idxpath != NULL);

	fa->fa_idxFP = fopen(idxpath, "w");

	if (fa->fa_idxFP == NULL) {
		dP("Unable to open \"%s\" for writing\n", idxpath);
		return (-1);
	}

	if (file_get_stat(fa->fa_seqFP, &st) != 0) {
		fclose(fa->fa_idxFP);
		fa->fa_idxFP = NULL;
		return (-1);
	}

	fprintf(fa->fa_idxFP,
		";filesize=%"PRIu64"\n"
		";chksum=0x%08x\n"
		";rcount=%u\n",
		st.st_size, 0x0 /* TODO */, fa->fa_rcount);

	for (i = 0; i < fa->fa_rcount; ++i) {
		fprintf(fa->fa_idxFP,
			"%"PRIu64" "
			"%"PRIu32" "
			"%"PRIu64" "
			"%"PRIu64" "
			"%"PRIu64" "
			"%"PRIu32" "
			"%"PRIu32" "
			"%"PRIu32"\n",
			fa->fa_record[i].hdr_start,
			fa->fa_record[i].hdr_len,
			fa->fa_record[i].seq_start,
			fa->fa_record[i].seq_rawlen,
			fa->fa_record[i].seq_len,
			fa->fa_record[i].seq_lines,
			fa->fa_record[i].seq_linew,
			fa->fa_record[i].seq_lastw);
	}

	if (!(fa->fa_options & FASTA_KEEPOPEN)) {
		fclose(fa->fa_idxFP);
		fa->fa_idxFP = NULL;
	}

	return (0);
}

static int __fahdr_read0(FILE *fp, FASTA_rec_t *dst)
{
	register int ch;
	register uint32_t i;

	char    *buffer;
	uint32_t buflen;
	char    *buftok;

	ch = getc_unlocked(fp);

	if (ch != '>')
		return (-1);

	/*
	 * Read all headers
	 */
	dst->hdr_start = file_get_offset(fp);
	dst->hdr_len   = 0;
	dst->hdr_cnt   = 1;

	buflen = 1024;
	buffer = alloc_array(char, buflen);

	do {
		if (dst->hdr_len >= buflen) {
			buflen += 1024;
			buffer  = realloc_array(buffer, char, buflen);
		}

		buffer[dst->hdr_len++] = ch = getc_unlocked(fp);

		if (ch == 0x01) /* ^A */
			++dst->hdr_cnt;

	} while(ch != '\n');

	buffer[dst->hdr_len - 1] = '\0';
	buffer = realloc_array(buffer, char, dst->hdr_len);
	dst->hdr_mem = buffer;

	dP(" Read header: \"%s\"\n", buffer);
	dP("Header count: %u\n", dst->hdr_cnt);

	/*
	 * Parse headers
	 */
	dst->hdr    = alloc_array(FASTA_rechdr_t, dst->hdr_cnt);
	dst->flags |= FASTA_REC_FREEHDR;

	i = 0;

	while ((buftok = strsep(&buffer, "\x01")) != NULL) {
		/*
		 * Sanity check
		 */
		if (i >= dst->hdr_cnt) {
			dP("Insane value(s): i=%u, dst->hdr_cnt=%u\n", i, dst->hdr_cnt);
			goto fail;
		}

		/*
		 * Parse the header using the SeqID parser
		 */
		if ((dst->hdr[i].seqid_fmt = SeqID_parse(buftok, strlen(buftok),
							 &dst->hdr[i].seqid)) == SEQID_ERROR)
		{
			dP("SeqID returned an error: h=\"%s\" l=%zu\n", buftok, strlen(buftok));
			goto fail;
		}

		++i;
	}

	if (dst->hdr_cnt > 0)
		dst->rec_id = dst->hdr[0].seqid.common.id;

	return (0);
fail:
	if (dst->hdr != NULL)
		free(dst->hdr);
	if (dst->hdr_mem != NULL)
		free(dst->hdr_mem);

	return (-1);
}

/**
 * Read key=value pairs from the index file, ignoring unknown
 * keys
 */
static int __idxhdr_read0(FASTA *fa, FASTA_idxhdr_t *ihdr)
{
	char   buffer[2048+1];
	char  *bufptr;
	char  *buftok;
	size_t buflen;

	register int ch;

	while (!feof(fa->fa_idxFP)) {
		ch = getc_unlocked(fa->fa_idxFP);

		dP("ch=%c\n", ch);

		if (ch != ';') {
			ungetc(ch, fa->fa_idxFP);
			return(0);
		}

		if (fgets(buffer,
			  sizeof buffer, fa->fa_idxFP) == NULL)
		{
			break;
		}

		buflen = strlen(buffer);

		if (buffer[buflen - 1] != '\n')
			return (FASTA_ENOBUF);

		bufptr = buffer;
		buftok = strsep(&bufptr, "=");

		if (buftok == NULL || bufptr == NULL)
			continue; /* skip this line */

		if (strcmp(buftok, "filesize") == 0) {
			ihdr->filesize = strtoull(bufptr, NULL, 10);

			dP("filesize=%"PRIu64"\n", ihdr->filesize);

			if (errno == ERANGE || errno == EINVAL)
				return (FASTA_EINVAL);
		} else if (strcmp(buftok, "chksum") == 0) {
			ihdr->chksum = strtol(bufptr, NULL, 16);

			dP("chksum=0x%08x\n", ihdr->chksum);

			if (errno == ERANGE || errno == EINVAL)
				return (FASTA_EINVAL);
		} else if (strcmp(buftok, "rcount") == 0) {
			ihdr->rcount = strtol(bufptr, NULL, 10);

			dP("rcount=0x%08x (%s)\n", ihdr->rcount, bufptr);

			if (errno == ERANGE || errno == EINVAL)
				return (FASTA_EINVAL);
		}
	}

	return (FASTA_EUNEXPEOF);
}

static int __index_read0(FILE *idxFP, FILE *seqFP, FASTA_rec_t *dst)
{
	int r;

	dst->flags   = 0;
	dst->chksum  = 0;
	dst->hdr     = NULL;
	dst->hdr_mem = NULL;
	dst->seq_mem = NULL;

	if (!feof_unlocked(idxFP)) {
		switch (r = fscanf(idxFP,
			       "%"PRIu64" "
			       "%"PRIu32" "
			       "%"PRIu64" "
			       "%"PRIu64" "
			       "%"PRIu64" "
			       "%"PRIu32" "
			       "%"PRIu32" "
			       "%"PRIu32"\n",
			       &dst->hdr_start,
			       &dst->hdr_len,
			       &dst->seq_start,
			       &dst->seq_rawlen,
			       &dst->seq_len,
			       &dst->seq_lines,
			       &dst->seq_linew,
			       &dst->seq_lastw))
		{
		case 8:
			dP("index record\n"
			   "=> %"PRIu64"\n"
			   "=> %"PRIu32"\n"
			   "=> %"PRIu64"\n"
			   "=> %"PRIu64"\n"
			   "=> %"PRIu64"\n"
			   "=> %"PRIu32"\n"
			   "=> %"PRIu32"\n"
			   "=> %"PRIu32"\n",
			   dst->hdr_start,
			   dst->hdr_len,
			   dst->seq_start,
			   dst->seq_rawlen,
			   dst->seq_len,
			   dst->seq_lines,
			   dst->seq_linew,
			   dst->seq_lastw);

			dst->flags   = 0;
			dst->seq_mem = NULL;

			/*
			 * Seek to hdr_start - 1, because __fahdr_read0 expects the '>'
			 */
                        if (file_set_offset(seqFP, dst->hdr_start - 1) != 0) {
				dP("Failed to seek to position %zu in %p\n", dst->hdr_start, seqFP);
				return (-1);
			}

			return __fahdr_read0(seqFP, dst);
		case EOF:
			return (1);
#ifndef NDEBUG
		case 0:
		{
			char buffer[2048+1];

			fgets(buffer, sizeof buffer, idxFP);
			dP("Input doesn't match format:\n"
			   "=> %s\n", buffer);
		}
#endif
		default:
			dP("r=%d\n", r);
			return (-1);
		}
	}

	return (1);
}

static inline void __fasta_cdseg_process(FASTA *fa, FASTA_rec_t *dst, uint8_t ch, bool *in_cds, uint64_t i)
{
        if (iscodingseq(fa, ch) && ch != 0) {
                if (!*in_cds) {
                        /*
                         * transition from non-coding to coding
                         * => create a new coding segment entry
                         */
                        dst->cdseg = realloc_array(dst->cdseg, FASTA_u64p, ++dst->cdseg_count);

                        dst->cdseg[dst->cdseg_count - 1].a = i;
                        dst->cdseg[dst->cdseg_count - 1].b = i;

                        *in_cds = true;
                }
        } else {
                if (*in_cds) {
                        /*
                         * transition from coding to non-coding
                         * => finalize the current coding segment
                         *    entry
                         */
                        assert(i > 0);
                        assert(dst->cdseg_count > 0);

                        dst->cdseg[dst->cdseg_count - 1].b = i - 1;
                        *in_cds = false;
                }
        }
}

/**
 * Read a variable line length sequence record into memory.
 */
static int __fasta_read2(FASTA *fa, FASTA_rec_t *dst, atrans_t *atr)
{
	size_t   alloc_size;
        uint8_t *buffer;
	size_t   buflen;

	register uint32_t n;
	register uint64_t i;
        bool in_cds = false;

        dP("read2\n");

	if (file_set_offset(fa->fa_seqFP, dst->seq_start) != 0) {
		dP("Failed to seek to position %zu in %p\n", dst->seq_start, fa->fa_seqFP);
		return (-1);
	}

	dP("atr=%p\n", atr);

	if (atr != NULL)
		alloc_size = atrans_s2d_size(atr, dst->seq_rawlen);
	else
		alloc_size = dst->seq_rawlen;

	if (dst->flags & FASTA_CSTRSEQ)
		++alloc_size;

	dP("alloc_size=%zu\n", alloc_size);

	i = 0;
	dst->seq_mem = malloc(alloc_size);
	bzero(dst->seq_mem, alloc_size);

	buffer = alloc_array(uint8_t, FASTA_LINEBUFFER_SIZE);
	buflen = FASTA_LINEBUFFER_SIZE;

        dst->cdseg   = NULL;
        dst->cdseg_count = 0;

	for (; dst->seq_lines > 0;) {
		/*
		 * Read line into buffer
		 */
		buflen = fread(buffer, 1, buflen, fa->fa_seqFP);

		if (buflen == 0) {
			if (feof(fa->fa_seqFP) && dst->seq_lines < 2)
				break;
			else {
				dP("An error occured during fread(%p, 1, %zu, %p): errno=%d, %s.\n",
				   buffer, buflen, fa->fa_seqFP, errno, strerror(errno));

				free(buffer);
				free(dst->seq_mem);
				dst->seq_mem = NULL;

				return (-1);
			}
		}

		/*
		 * Copy/Translate the buffer into seq_mem
		 */
		if (atr != NULL) {
			for (n = 0; n < buflen; ++n) {
				if (issequence(buffer[n])) {
                                        if (dst->flags & FASTA_MAPCDSEG)
                                                __fasta_cdseg_process(fa, dst, buffer[n], &in_cds, i);

					atrans_letter_s2d(atr, buffer[n], i++, (uint8_t *)dst->seq_mem);
				} else {
					if (buffer[n] == '\n') {
						assert(dst->seq_lines > 0);
						--dst->seq_lines;

						if (dst->seq_lines == 0)
							goto __A_finish;
					} else {
						switch (buffer[n]) {
						case ' ':
							/* ignore */
							break;
						default:
							dP("Unexpected character: %c (%u)\n", (char)buffer[n], buffer[n]);

							free(buffer);
							free(dst->seq_mem);
							dst->seq_mem = NULL;

							return (-1);
						}
					}
				}
			}
		} else {
			for (n = 0; n < buflen; ++n) {
				if (issequence(buffer[n])) {
                                        if (dst->flags & FASTA_MAPCDSEG)
                                                __fasta_cdseg_process(fa, dst, buffer[n], &in_cds, i);

					((uint8_t *)(dst->seq_mem))[i++] = buffer[n];
				} else {
					if (buffer[n] == '\n') {
						assert(dst->seq_lines > 0);
						--dst->seq_lines;

						if (dst->seq_lines == 0)
							goto __A_finish;
					} else {
						switch (buffer[n]) {
						case ' ':
							/* ignore */
							break;
						default:
							dP("Unexpected character: %c (%u)\n", (char)buffer[n], buffer[n]);

							free(buffer);
							free(dst->seq_mem);
							dst->seq_mem = NULL;

							return (-1);
						}
					}
				}
			}
		}
	}
	__A_finish:
	free(buffer);
	dst->seq_len = i;

        if (dst->flags & FASTA_MAPCDSEG)
                __fasta_cdseg_process(fa, dst, 0, &in_cds, i);

	if (dst->flags & FASTA_CSTRSEQ) {
		dst->seq_mem[i] = '\0';
		dst->seq_mem = realloc_array(dst->seq_mem, uint8_t, dst->seq_len + 1);
	} else
		dst->seq_mem = realloc_array(dst->seq_mem, uint8_t, dst->seq_len);

	return (0);
}

/**
 * Read a equal line length sequence record into memory.
 */
static int __fasta_read1(FASTA *fa, FASTA_rec_t *dst, atrans_t *atr)
{
	size_t   alloc_size;
	uint8_t *buffer;
	size_t   buflen;

	register uint32_t l, n;
	register uint64_t i;
        bool in_cds = false;

	if (file_set_offset(fa->fa_seqFP, dst->seq_start) != 0) {
		dP("Failed to seek to position %zu in %p\n", dst->seq_start, fa->fa_seqFP);
		return (-1);
	}

	dP("atr=%p\n", atr);

	if (atr != NULL)
		alloc_size = atrans_s2d_size(atr, dst->seq_len);
	else
		alloc_size = dst->seq_len;

	if (dst->flags & FASTA_CSTRSEQ)
		++alloc_size;

	dP("alloc_size=%zu\n", alloc_size);

	i = 0;
	dst->seq_mem = malloc(alloc_size);
        dst->cdseg   = NULL;
        dst->cdseg_count = 0;

	if (atr != NULL) {
		/*
		 * An alphabet translation table was defined
		 * => Read & Translate
		 */
		bzero(dst->seq_mem, alloc_size);

		buffer = alloc_array(uint8_t, dst->seq_linew);
		buflen = dst->seq_linew;

		for (l = dst->seq_lines - 1; l > 0; --l) {
			/*
			 * Read line into buffer
			 */
			if (fread(buffer, 1, buflen, fa->fa_seqFP) != buflen) {
				/* fail */
				free(buffer);
				free(dst->seq_mem);
				dst->seq_mem = NULL;

				return (-1);
			}

			/*
			 * Copy/Translate the buffer into seq_mem and/or map coding segments
			 */
			for (n = 0; n < buflen; ++n) {
                                if (dst->flags & FASTA_MAPCDSEG)
                                        __fasta_cdseg_process(fa, dst, buffer[n], &in_cds, i);

				atrans_letter_s2d(atr, buffer[n], i++, (uint8_t *)dst->seq_mem);
                        }

			getc_unlocked(fa->fa_seqFP); /* skip the new-line */
		}

		buflen = dst->seq_lastw > 0 ? dst->seq_lastw : dst->seq_linew;
		/* no need to reallocate the buffer since lastw < linew */

		if (fread(buffer, 1, buflen, fa->fa_seqFP) != buflen) {
			/* fail */
			free(buffer);
			free(dst->seq_mem);
			dst->seq_mem = NULL;

			return (-1);
		}

		for (n = 0; n < buflen; ++n) {
                        if (dst->flags & FASTA_MAPCDSEG)
                                __fasta_cdseg_process(fa, dst, buffer[n], &in_cds, i);

			atrans_letter_s2d(atr, buffer[n], i++, (uint8_t *)dst->seq_mem);
                }

                /*
                 * If there is an open coding segment, this call will finalize it.
                 */
                if (dst->flags & FASTA_MAPCDSEG)
                        __fasta_cdseg_process(fa, dst, 0, &in_cds, i);

		free(buffer);
	} else {
		/*
		 * No alphabet translation defined
		 * => Read
		 */
		buflen = dst->seq_linew;

		for (l = dst->seq_lines - 1; l > 0; --l) {
			/*
			 * Read line into buffer
			 */
                        dP("l = %u, i=%"PRIu64", buflen=%"PRIu64"\n", l, i, buflen);

			if (fread(dst->seq_mem + i, 1, buflen, fa->fa_seqFP) != buflen) {
				/* fail */
				free(dst->seq_mem);
				dst->seq_mem = NULL;

				return (-1);
			}

                        if (dst->flags & FASTA_MAPCDSEG) {
                                for (n = 0; n < buflen; ++n, ++i)
                                        __fasta_cdseg_process(fa, dst, dst->seq_mem[i], &in_cds, i);
                        } else
                                i += buflen;

			getc_unlocked(fa->fa_seqFP); /* skip the new-line */
		}

		buflen = dst->seq_lastw > 0 ? dst->seq_lastw : dst->seq_linew;

		if (fread(dst->seq_mem + i, 1, buflen, fa->fa_seqFP) != buflen) {
			/* fail */
			free(dst->seq_mem);
			dst->seq_mem = NULL;

			return (-1);
		}

                if (dst->flags & FASTA_MAPCDSEG) {
                        for (n = 0; n < buflen; ++n, ++i)
                                __fasta_cdseg_process(fa, dst, dst->seq_mem[i], &in_cds, i);
                } else
                        i += buflen;

                /*
                 * If there is an open coding segment, this call will finalize it.
                 */
                if (dst->flags & FASTA_MAPCDSEG)
                        __fasta_cdseg_process(fa, dst, 0, &in_cds, i);
	}

	if (dst->flags & FASTA_CSTRSEQ)
		dst->seq_mem[alloc_size - 1] = '\0';

	return (0);
}

/**
 * Analyze a sequence record.
 */
static int __fasta_read0(FILE *fp, FASTA_rec_t *dst, uint32_t options, atrans_t *atr)
{
	int      ch;
	uint32_t plinew; /* previous line width */
	uint32_t clinew; /* current line width */

        (void)atr;

        dP("read0\n");

	assert(fp  != NULL);
	assert(dst != NULL);

	dst->flags   = 0;
	dst->chksum  = 0;
	dst->hdr     = NULL;
	dst->hdr_mem = NULL;
	dst->seq_mem = NULL;
        dst->cdseg   = NULL;
        dst->cdseg_count = 0;
	plinew = 0;
	clinew = 0;

	while (!feof_unlocked(fp)) {
		/*
		 * Read & Parse FASTA header(s)
		 */
		if (__fahdr_read0(fp, dst) != 0)
			return (-1);

		/*
		 * Analyze/Read the sequence
		 */
		if (options & FASTA_INMEMSEQ) {
			/*
                         * Not valid in this function
                         */
                        return (FASTA_EINVAL);
		} else {
			/*
			 * Analyze the sequence without storing it into main memory.
			 */
			bool linew_diff = false, linew_update = true;

			plinew = 0;
			clinew = 0;

			dst->seq_start  = file_get_offset(fp);

			dst->seq_len    = 0;
			dst->seq_rawlen = 0;
			dst->seq_lines  = 0;

			dst->seq_linew  = 0;
			dst->seq_lastw  = 0;

			/*
			 * Read in the first line of the sequence.
			 */
			for (;;) {
				ch = getc_unlocked(fp);

				if (feof_unlocked(fp)) {
					if (dst->seq_len == 0) {
						dP("Unexpected EOF: got header, but no sequence data\n");
						goto fail;
					} else {
						++dst->seq_lines;
						goto finalize_seq;
					}
				}

				++dst->seq_rawlen;
				dst->chksum = crc32(dst->chksum, (unsigned char *)&ch, 1);

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
				ch = getc_unlocked(fp);

				if (feof_unlocked(fp)) {
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
						if (linew_update /* && clinew > 0 */) {
							if (clinew != plinew) {
								if (linew_diff == true || clinew == 0) {
									linew_update = false;
									plinew = 0;
                                                                        clinew = 0;
								}

                                                                linew_diff = true;
							} else {
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
							dP("Unexpected '>': allowed only at the beginning of a line\n");
							goto fail;
						}
						break;
					case ' ': /* ignore */
						if (linew_update) {
							clinew = 0;
							plinew = 0;
							linew_update = false;
							linew_diff   = false;
						}
						break;
					default:
						dP("Unexpected character '%c' (%u)\n", ch, ch);
						goto fail;
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
	return (feof_unlocked(fp) ? 1 : 0);
fail:
	if (dst->hdr != NULL)
		free(dst->hdr);
	if (dst->hdr_mem != NULL)
		free(dst->hdr_mem);
	if (dst->seq_mem != NULL)
		free(dst->seq_mem);

	return (-1);
}

FASTA *fasta_open(const char *path, uint32_t options, atrans_t *atr)
{
	char   idx_path[PATH_MAX + 1];
	FASTA *fa;
	struct stat st;

	assert(path != NULL);

	fa             = alloc_type(FASTA);
	fa->fa_options = options;
	fa->fa_path    = strdup(path);
	fa->fa_seqFP   = NULL;
	fa->fa_idxFP   = NULL;
	fa->fa_record  = NULL;
	fa->fa_rindex  = 0;
	fa->fa_rcount  = 0;
	fa->fa_atr     = atr;
        fa->fa_CDSmask = (uint32_t *)__SQ_mask;

        fasta_setCDS(fa, options);

	fa->fa_seqFP = fopen(path, "r");

	if (fa->fa_seqFP == NULL) {
		dP("Can't open the sequence file: %s\n", path);
		goto fail;
	}

        setbuf(fa->fa_seqFP, NULL);
	flockfile(fa->fa_seqFP);

	if (file_get_stat(fa->fa_seqFP, &st) != 0) {
#ifndef NDEBUG
		int e = errno;
		dP("Failed to get stat information: fp=%p (path=%s), errno=%u, %s\n",
		   fa->fa_seqFP, path, e, strerror(e));
#endif
		goto fail;
	}

	if (options & FASTA_USEINDEX) {
		FASTA_idxhdr_t idxhdr;

		memset(&idxhdr, 0, sizeof (FASTA_idxhdr_t));

		if (strlen(path) + strlen(FASTA_INDEX_EXT) < (sizeof idx_path/sizeof(char))) {
			strcpy(idx_path, path);
			strcat(idx_path, FASTA_INDEX_EXT);
		} else {
			dP("Index path exceeds systems limits: idx_path_len=%zu, limit=%zu\n",
			   strlen(path) + strlen(FASTA_INDEX_EXT), sizeof idx_path/sizeof(char));
			goto fail;
		}

		fa->fa_idxFP = fopen(idx_path,
				     (options & FASTA_GENINDEX ? "r" : "r+"));

		if (fa->fa_idxFP == NULL)
			goto regen;

		flockfile(fa->fa_seqFP);

		switch (__idxhdr_read0(fa, &idxhdr)) {
		case 0:
			break;
		default:
			if (fa->fa_options & FASTA_CHKINDEX_FAIL)
				goto fail;
			else {
				funlockfile(fa->fa_idxFP);
				fclose(fa->fa_idxFP);
				fa->fa_idxFP = NULL;
				goto regen;
			}
		}

		dP("Read index header\n"
		   "=> filesize: %"PRIu64"\n"
		   "=>   chksum: 0x%08x\n"
		   "=>   rcount: %u\n", idxhdr.filesize, idxhdr.chksum, idxhdr.rcount);

		if ((uint64_t)st.st_size != idxhdr.filesize) {
			dP("Recorded (%zu) and actual (%zu) filesizes differ!\n", idxhdr.filesize, st.st_size);
			funlockfile(fa->fa_idxFP);
			fclose(fa->fa_idxFP);
			fa->fa_idxFP = NULL;
			goto regen;
		}

		if (options & FASTA_CHKINDEX_SLOW) {
			/* slow check */
		} else if (options & FASTA_CHKINDEX_FAST) {
			register uint32_t i;
			int r;

			i = 0;
			fa->fa_rcount = 0;
			fa->fa_record = NULL;

			do {
				if (i >= fa->fa_rcount) {
					dP("=> pre-alloc: fa_rcount=%u\n", fa->fa_rcount);

					if (i == 0)
						fa->fa_rcount = 8;
					else if (i < 65535)
						fa->fa_rcount <<= 1;
					else
						fa->fa_rcount += 1024;

					fa->fa_record = realloc_array(fa->fa_record, FASTA_rec_t, fa->fa_rcount);

					dP("<= pre-alloc: fa_rcount=%u\n", fa->fa_rcount);
				}

				dP("Reading index record #%u\n", i);
			} while ((r = __index_read0(fa->fa_idxFP, fa->fa_seqFP, fa->fa_record + i++)) == 0);
			if (r < 0) {
				dP("An error ocured while reading the file \"%s\"\n", idx_path);
				goto prefail;
			}

			fa->fa_rcount = --i;
			fa->fa_record = realloc_array(fa->fa_record, FASTA_rec_t, fa->fa_rcount);

			if (fa->fa_rcount != idxhdr.rcount) {
				dP("fa->fa_rcount (%u) != idxhdr.rcount (%u)\n", fa->fa_rcount, idxhdr.rcount);
			prefail:
				if (fa->fa_options & FASTA_CHKINDEX_FAIL)
					goto fail;
				else {
					funlockfile(fa->fa_idxFP);
					fclose(fa->fa_idxFP);
					fa->fa_idxFP = NULL;

					if (fa->fa_rcount > 0) {
						for (i = 0; i < fa->fa_rcount; ++i)
							fasta_rec_free(fa->fa_record + i);

						free(fa->fa_record);

						fa->fa_record = NULL;
						fa->fa_rcount = 0;
					}

					goto regen;
				}
			}
		}
	} else {
		register uint32_t i;
		int r;

	regen:
		/*
		 * generate headers from the sequence file
		 */
		i = 0;
		fa->fa_rcount = 0;

		rewind(fa->fa_seqFP);

		do {
			if (i >= fa->fa_rcount) {
				dP("=> pre-alloc: fa_rcount=%u\n", fa->fa_rcount);

				if (i == 0)
					fa->fa_rcount = 8;
				else if (i < 65535)
					fa->fa_rcount <<= 1;
				else
					fa->fa_rcount += 1024;

				fa->fa_record = realloc_array(fa->fa_record, FASTA_rec_t, fa->fa_rcount);

				dP("<= pre-alloc: fa_rcount=%u\n", fa->fa_rcount);
			}

			dP("Reading sequence #%u\n", i);
		} while ((r = __fasta_read0(fa->fa_seqFP, fa->fa_record + i++, options, fa->fa_atr)) == 0);

		if (r < 0) {
			dP("An error ocured while reading the file \"%s\"\n", fa->fa_path);

			/*
			 * Decrease `i' to prevent double-free. __fasta_read0 ensures that in
			 * case of an error the requested fasta record will not be initialized
			 * and therefore will not need to be freed.
			 */
			fa->fa_rcount = --i;

			goto fail;
		}

		fa->fa_rcount = i;
		fa->fa_record = realloc_array(fa->fa_record, FASTA_rec_t, fa->fa_rcount);

		if ((fa->fa_options & FASTA_USEINDEX) && (fa->fa_options & FASTA_GENINDEX))
			__index_write(fa, idx_path);
	}

	funlockfile(fa->fa_seqFP);

	if (fa->fa_idxFP != NULL)
		funlockfile(fa->fa_idxFP);

	if (!(options & FASTA_KEEPOPEN)) {
		if (fa->fa_idxFP != NULL)
			fclose(fa->fa_idxFP);

		fclose(fa->fa_seqFP);

		fa->fa_seqFP = NULL;
		fa->fa_idxFP = NULL;
	}

	return (fa);
fail:
	for (; fa->fa_rcount > 0; --fa->fa_rcount)
		fasta_rec_free(fa->fa_record + fa->fa_rcount - 1);

	if (fa->fa_seqFP != NULL) {
		funlockfile(fa->fa_seqFP);
		fclose (fa->fa_seqFP);
	}

	if (fa->fa_idxFP != NULL) {
		funlockfile(fa->fa_idxFP);
		fclose(fa->fa_idxFP);
	}

	free(fa);

	return (NULL);
}

uint32_t fasta_count(FASTA *fa)
{
	return (fa->fa_rcount);
}

int fasta_setCDS(FASTA *fa, uint32_t cds_flags)
{
        assert(fa != NULL);

        if (cds_flags & FASTA_NASEQ)
                fa->fa_CDSmask = (uint32_t *)__NA_CD_mask;
        else if (cds_flags & FASTA_AASEQ)
                fa->fa_CDSmask = (uint32_t *)__AA_CD_mask;

        return (0);
}

int fasta_setCDS_string(FASTA *fa, const char *letters)
{
        size_t    l, i;
        char     *s;
        uint32_t  b, m, *mask = alloc_array(uint32_t, 8);

        if (mask == NULL)
                return (-1);

        l = strlen(letters);

	for (s = (char *)letters, i = 0; i < l; ++i) {
		b = s[i] / (sizeof mask[0] * 8);
		m = s[i] % (sizeof mask[0] * 8);
		mask[b] |= 1 << m;
	}

        if (fa->fa_CDSmask != NULL && (fa->fa_options & FASTA_CDSFREEMASK))
                free(fa->fa_CDSmask);

        fa->fa_CDSmask  = mask;
        fa->fa_options |= FASTA_CDSFREEMASK;

        return (0);
}


FASTA_rec_t *fasta_read(FASTA *fa, FASTA_rec_t *dst, uint32_t flags, atrans_t *atr)
{
	FASTA_rec_t *farec;

	assert(fa != NULL);

	if (fa->fa_rindex >= fa->fa_rcount)
		return (NULL);

	if (atr == NULL)
		atr = fa->fa_atr;

	if (fa->fa_seqFP == NULL) {
		fa->fa_seqFP = fopen(fa->fa_path, "r");

		if (fa->fa_seqFP == NULL) {
			dP("Can't re-open the sequence file: %s\n", fa->fa_path);
			return (NULL);
		}

		flockfile(fa->fa_seqFP);
	}

	if (dst == NULL) {
		if (flags & FASTA_RAWREC)
			farec = fa->fa_record + fa->fa_rindex;
		else {
			farec = alloc_type(FASTA_rec_t);
			memcpy(farec, fa->fa_record + fa->fa_rindex, sizeof(FASTA_rec_t));
			farec->flags = FASTA_REC_MAGICFL | FASTA_REC_FREEREC;
		}
	} else {
		farec = dst;
		memcpy(farec, fa->fa_record + fa->fa_rindex, sizeof(FASTA_rec_t));
		farec->flags = FASTA_REC_MAGICFL;
	}

	if (flags & FASTA_MAPCDSEG)
		farec->flags |= FASTA_MAPCDSEG;

	if (flags & FASTA_CSTRSEQ)
		farec->flags |= FASTA_CSTRSEQ;

	if ((flags & FASTA_INMEMSEQ) && farec->seq_mem == NULL) {
		farec->flags |= FASTA_REC_FREESEQ;

		if (farec->seq_linew != 0) {
			/*
			 * all the lines that form the sequence are of equal length
			 */
			if (__fasta_read1(fa, farec, atr) != 0) {
				/* fail */
				fasta_rec_free(farec);
				farec = NULL;
			} else
				++fa->fa_rindex;
		} else {
			/*
			 * the lines have variable length, we have to look for new-lines
			 */
			if (__fasta_read2(fa, farec, atr) != 0) {
				/* fail */
				fasta_rec_free(farec);
				farec = NULL;
			} else
				++fa->fa_rindex;
		}
	}

	funlockfile(fa->fa_seqFP);

	if (!((fa->fa_options | flags) & FASTA_KEEPOPEN)) {
		fclose(fa->fa_seqFP);
		fa->fa_seqFP = NULL;
	}

	return (farec);
}

int fasta_write(FASTA *fa, FASTA_rec_t *farec)
{
        (void)fa;
        (void)farec;

        /* TODO */

	return (-1);
}

void fasta_rec_free(FASTA_rec_t *farec)
{
	if (farec->flags & FASTA_REC_FREEHDR) {
		free(farec->hdr_mem);
		free(farec->hdr);
		farec->hdr     = NULL;
		farec->hdr_cnt = 0;
	}

	if (farec->flags & FASTA_REC_FREESEQ) {
		free(farec->seq_mem);
		farec->seq_mem = NULL;
		farec->flags  &= ~(FASTA_REC_FREESEQ);
	}

	if (farec->flags & FASTA_REC_FREEREC) {
		farec->flags = 0;
		free(farec);
	}

	return;
}

FASTA_CDS_t *fasta_read_CDS(FASTA *fa, FASTA_rec_t *farec, FASTA_CDS_t *dst, uint32_t flags)
{
        (void)fa;
        (void)flags;

        if (farec->cdseg_count == 0 || farec->cdseg_index == farec->cdseg_count)
                return (NULL);

        if (dst == NULL) {
                dst = alloc_type(FASTA_CDS_t);
                dst->flags = FASTA_CDSFREEMASK;
        }

        dst->farec   = farec;
        dst->seg_idx = farec->cdseg_index;
        dst->seg_len = farec->cdseg[dst->seg_idx].b - farec->cdseg[dst->seg_idx].a + 1;
        dst->seg_mem = farec->seq_mem + farec->cdseg[dst->seg_idx].a;

        return (dst);
}

int fasta_rewind_CDS(FASTA *fa, FASTA_rec_t *farec)
{
        (void)fa;
        farec->cdseg_index = 0;
        return (0);
}

int fasta_seeko_CDS(FASTA *fa, FASTA_rec_t *farec, off_t off, int whence)
{
	off_t newpos;

        (void)fa;
	assert(fa != NULL);

	switch (whence) {
	case SEEK_SET:
		newpos = off;
		break;
	case SEEK_CUR:
		newpos = off + farec->cdseg_index;
		break;
	case SEEK_END:
		newpos = farec->cdseg_count - off - 1;
		break;
	default:
		errno = EINVAL;
		return (-1);
	}

	if (newpos < 0 || newpos >= (off_t)farec->cdseg_count) {
		errno = ERANGE;
		return (-1);
	}

	farec->cdseg_index = (size_t) newpos;

	return (0);
}

off_t fasta_tello_CDS(FASTA *fa, FASTA_rec_t *farec)
{
        (void)fa;
        return (farec->cdseg_index);
}

int fasta_rewind(FASTA *fa)
{
	fa->fa_rindex = 0;
	return (0);
}

int fasta_seeko(FASTA *fa, off_t off, int whence)
{
	off_t newpos;

	assert(fa != NULL);

	switch (whence) {
	case SEEK_SET:
		newpos = off;
		break;
	case SEEK_CUR:
		newpos = off + fa->fa_rindex;
		break;
	case SEEK_END:
		newpos = fa->fa_rcount - off - 1;
		break;
	default:
		errno = EINVAL;
		return (-1);
	}

	if (newpos < 0 || newpos >= fa->fa_rcount) {
		errno = ERANGE;
		return (-1);
	}

	fa->fa_rindex = (uint32_t) newpos;

	return (0);
}

off_t fasta_tello(FASTA *fa)
{
	return (fa->fa_rindex);
}

void *fasta_apply(FASTA *fa, void * (*func)(FASTA_rec_t *, void *), uint32_t options, void *funcarg)
{
	void        **result;
	uint32_t      resnum, i;
	FASTA_rec_t  *farec;

        (void)options;

	resnum = fasta_count(fa);

	if (resnum == 0 || fasta_rewind(fa) != 0)
		return (NULL);

	result = (void **) alloc_array(void *, resnum);

	for (i = 0; i < resnum; ++i) {
		farec     = fasta_read(fa, NULL, FASTA_INMEMSEQ|FASTA_CSTRSEQ, NULL);
		result[i] = func(farec, funcarg);
		fasta_rec_free(farec);
	}

	return (result);
}

void fasta_close(FASTA *fa)
{
	if (fa->fa_record != NULL) {
		for (; fa->fa_rcount > 0; --fa->fa_rcount)
			fasta_rec_free(fa->fa_record + fa->fa_rcount - 1);
		free(fa->fa_record);
	}

	free(fa->fa_path);

	if (fa->fa_seqFP != NULL)
		fclose (fa->fa_seqFP);
	if (fa->fa_idxFP != NULL)
		fclose (fa->fa_idxFP);

	free(fa);
	return;
}

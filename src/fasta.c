#if defined(__linux__)
# define _XOPEN_SOURCE
#endif
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
#include "assume.h"
#include "debug.h"
#include "sm_alloc.h"
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

static bool issequence(int ch)
{
	if (ch >= 0 && ch < 256)
		return (__SQ_mask[ch / (sizeof __SQ_mask[0] * 8)]) & (1 << (ch % (sizeof __SQ_mask[0] * 8)));
	else
		return (false);
}

static int __index_write(FASTA *fa, const char *idxpath)
{
	register uint32_t i;
	struct stat st;

	assume_d(fa != NULL, -1);
	assume_d(idxpath != NULL, -1);

	fa->fa_idxFP = fopen(idxpath, "w");

	if (fa->fa_idxFP == NULL) {
		_D("Unable to open \"%s\" for writing\n", idxpath);
		return (-1);
	}

	if (fstat(fileno(fa->fa_seqFP), &st) != 0) {
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
	uint32_t toklen;

	ch = getc_unlocked(fp);

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

		buffer[dst->hdr_len++] = ch = getc_unlocked(fp);

		if (ch == 0x01) /* ^A */
			++dst->hdr_cnt;

	} while(ch != '\n');

	buffer[dst->hdr_len - 1] = '\0';
	buffer = sm_realloc(buffer, sizeof(char) * dst->hdr_len);
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
			_D("SeqID returned an error: h=\"%s\" l=%zu\n", buftok, strlen(buftok));
			goto fail;
		}

		++i;
	}

	if (dst->hdr_cnt > 0)
		dst->rec_id = dst->hdr[0].seqid.common.id;

	return (0);
fail:
	if (dst->hdr != NULL)
		sm_free(dst->hdr);
	if (dst->hdr_mem != NULL)
		sm_free(dst->hdr_mem);

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

		_D("ch=%c\n", ch);

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

			_D("filesize=%"PRIu64"\n", ihdr->filesize);

			if (errno == ERANGE || errno == EINVAL)
				return (FASTA_EINVAL);
		} else if (strcmp(buftok, "chksum") == 0) {
			ihdr->chksum = strtol(bufptr, NULL, 16);

			_D("chksum=0x%08x\n", ihdr->chksum);

			if (errno == ERANGE || errno == EINVAL)
				return (FASTA_EINVAL);
		} else if (strcmp(buftok, "rcount") == 0) {
			ihdr->rcount = strtol(bufptr, NULL, 10);

			_D("rcount=0x%08x (%s)\n", ihdr->rcount, bufptr);

			if (errno == ERANGE || errno == EINVAL)
				return (FASTA_EINVAL);
		}
	}

	return (FASTA_EUNEXPEOF);
}

static int __index_read0(FILE *idxFP, FILE *seqFP, FASTA_rec_t *dst)
{
	int r;

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
			_D("index record\n"
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
			if (fseeko(seqFP, dst->hdr_start - 1, SEEK_SET) != 0) {
				_D("Failed to seek to position %zu in %p\n", dst->hdr_start, seqFP);
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
			_D("Input doesn't match format:\n"
			   "=> %s\n", buffer);
		}
#endif
		default:
			_D("r=%d\n", r);
			return (-1);
		}
	}

	return (1);
}

static int __fasta_read2(FILE *fp, FASTA_rec_t *dst, atrans_t *atr)
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
		alloc_size = atrans_s2d_size(atr, dst->seq_rawlen);
	else
		alloc_size = dst->seq_rawlen;

	if (dst->flags & FASTA_CSTRSEQ)
		++alloc_size;

	_D("alloc_size=%zu\n", alloc_size);

	i = 0;
	dst->seq_mem = sm_alloc(alloc_size);
	bzero(dst->seq_mem, alloc_size);

	buffer = sm_alloc(sizeof(uint8_t) * FASTA_LINEBUFFER_SIZE);
	buflen = FASTA_LINEBUFFER_SIZE;

	for (; dst->seq_lines > 0;) {
		/*
		 * Read line into buffer
		 */
		buflen = fread(buffer, 1, buflen, fp);

		if (buflen == 0) {
			if (feof(fp) && dst->seq_lines < 2)
				break;
			else {
				_D("An error occured during fread(%p, 1, %zu, %p): errno=%d, %s.\n",
				   buffer, buflen, fp, errno, strerror(errno));

				sm_free(buffer);
				sm_free(dst->seq_mem);
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
					atrans_letter_s2d(atr, buffer[n], i++, (uint8_t *)dst->seq_mem);
				} else {
					if (buffer[n] == '\n') {
						assume_d(dst->seq_lines > 0, -1);
						--dst->seq_lines;

						if (dst->seq_lines == 0)
							goto __A_finish;
					} else {
						switch (buffer[n]) {
						case ' ':
							/* ignore */
							break;
						default:
							_D("Unexpected character: %c (%u)\n", (char)buffer[n], buffer[n]);

							sm_free(buffer);
							sm_free(dst->seq_mem);
							dst->seq_mem = NULL;

							return (-1);
						}
					}
				}
			}
		} else {
			for (n = 0; n < buflen; ++n) {
				if (issequence(buffer[n])) {
					((uint8_t *)(dst->seq_mem))[i++] = buffer[n];
				} else {
					if (buffer[n] == '\n') {
						assume_d(dst->seq_lines > 0, -1);
						--dst->seq_lines;

						if (dst->seq_lines == 0)
							goto __A_finish;
					} else {
						switch (buffer[n]) {
						case ' ':
							/* ignore */
							break;
						default:
							_D("Unexpected character: %c (%u)\n", (char)buffer[n], buffer[n]);

							sm_free(buffer);
							sm_free(dst->seq_mem);
							dst->seq_mem = NULL;

							return (-1);
						}
					}
				}
			}
		}
	}
	__A_finish:
	sm_free(buffer);
	dst->seq_len = i;

	if (dst->flags & FASTA_CSTRSEQ) {
		dst->seq_mem[i] = '\0';
		dst->seq_mem = sm_realloc(dst->seq_mem, sizeof(uint8_t) * (dst->seq_len + 1));
	} else
		dst->seq_mem = sm_realloc(dst->seq_mem, sizeof(uint8_t) * dst->seq_len);

	return (0);
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
		/*
		 * An alphabet translation table was defined
		 * => Read & Translate
		 */
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

			getc_unlocked(fp); /* skip the new-line */
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
		/*
		 * No alphabet translation defined
		 * => Read
		 */
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

			getc_unlocked(fp); /* skip the new-line */
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
	uint32_t plinew; /* previous line width */
	uint32_t clinew; /* current line width */

	assume_d(fp  != NULL, -1);
	assume_d(dst != NULL, -1);

	dst->flags   = 0;
	dst->chksum  = 0;
	dst->hdr     = NULL;
	dst->hdr_mem = NULL;
	dst->seq_mem = NULL;
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
				ch = getc_unlocked(fp);

				if (feof_unlocked(fp)) {
					if (dst->seq_len == 0) {
						_D("Unexpected EOF: got header, but no sequence data\n");
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
					case ' ': /* ignore */
						if (linew_update) {
							clinew = 0;
							plinew = 0;
							linew_update = false;
							linew_diff   = false;
						}
						break;
					default:
						_D("Unexpected character '%c' (%u)\n");
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

	flockfile(fa->fa_seqFP);

	if (fstat(fileno(fa->fa_seqFP), &st) != 0) {
		int e = errno;
		_D("Failed to get stat information: fp=%p (fd=%d, path=%s), errno=%u, %s\n",
		   fa->fa_seqFP, fileno(fa->fa_seqFP), path, e, strerror(e));
		goto fail;
	}

	if (options & FASTA_USEINDEX) {
		FASTA_idxhdr_t idxhdr = { 0 };

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

		_D("Read index header\n"
		   "=> filesize: %llu\n"
		   "=>   chksum: 0x%08x\n"
		   "=>   rcount: %u\n", idxhdr.filesize, idxhdr.chksum, idxhdr.rcount);

		if (st.st_size != idxhdr.filesize) {
			_D("Recorded (%zu) and actual (%zu) filesizes differ!\n", idxhdr.filesize, st.st_size);
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

				_D("Reading index record #%u\n", i);
			} while ((r = __index_read0(fa->fa_idxFP, fa->fa_seqFP, fa->fa_record + i++)) == 0);

			if (r < 0) {
				_D("An error ocured while reading the file \"%s\"\n", idx_path);
				goto fail;
			}

			fa->fa_rcount = --i;
			fa->fa_record = sm_realloc(fa->fa_record, sizeof(FASTA_rec_t) * fa->fa_rcount);

			if (fa->fa_rcount != idxhdr.rcount) {
				_D("fa->fa_rcount (%u) != idxhdr.rcount (%u)\n", fa->fa_rcount, idxhdr.rcount);

				if (fa->fa_options & FASTA_CHKINDEX_FAIL)
					goto fail;
				else {
					funlockfile(fa->fa_idxFP);
					fclose(fa->fa_idxFP);
					fa->fa_idxFP = NULL;

					if (fa->fa_rcount > 0) {
						for (i = 0; i < fa->fa_rcount; ++i)
							fasta_rec_free(fa->fa_record + i);

						sm_free(fa->fa_record);

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

	if (fa->fa_seqFP == NULL) {
		fa->fa_seqFP = fopen(fa->fa_path, "r");

		if (fa->fa_seqFP == NULL) {
			_D("Can't re-open the sequence file: %s\n", fa->fa_path);
			return (NULL);
		}

		flockfile(fa->fa_seqFP);
	}

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
				fasta_rec_free(farec);
				farec = NULL;
			} else
				++fa->fa_rindex;
		} else {
			/*
			 * the lines have variable length, we have to look for new-lines
			 */
			if (__fasta_read2(fa->fa_seqFP, farec, atr) != 0) {
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
	off_t newpos;

	assume_d(fa != NULL, -1);

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

	resnum = fasta_count(fa);

	if (resnum == 0 || fasta_rewind(fa) != 0)
		return (NULL);

	result = (void **)sm_alloc(sizeof(void *) * resnum);

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

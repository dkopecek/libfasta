#ifndef FASTA_H
#define FASTA_H

#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

#define FASTA_KEEPOPEN      0x00000001 /**< Keep the FASTA file/index open */
#define FASTA_USEINDEX      0x00000002 /**< Use index, if present */
#define FASTA_GENINDEX      0x00000004 /**< Create the index, if not present */
#define FASTA_CHKINDEX_FAST 0x00000008 /**< Perform some index integrity check that aren't time consuming */
#define FASTA_CHKINDEX_SLOW 0x00000010 /**< Perform a full index integrity check */
#define FASTA_CHKINDEX_NONE 0x00000020 /**< Fully trust the index */
#define FASTA_CHKINDEX      FASTA_CHKINDEX_FAST
#define FASTA_INMEMSEQ      0x00000040 /**< Load all sequence data into memory */
#define FASTA_ONDEMSEQ      0x00000080 /**< Read a sequence into memory on-demand */
#define FASTA_PARALLEL      0x00000100 /**< Perform some operation in parallel, e.g. the _apply operation */
#define FASTA_READ          0x00000200 /**< Open for reading */
#define FASTA_WRITE         0x00000400 /**< Open for writing */
#define FASTA_RAWREC        0x00000800 /**< Return a pointer to an internally allocated FASTA record */
#define FASTA_CSTRSEQ       0x00001000 /**< Return the sequence as a C string (i.e. NUL terminated) */
#define FASTA_CHKINDEX_FAIL 0x00002000 /**< Fail if the index is corrupted, regenerate otherwise */

#define FASTA_INDEX_EXT ".index" /**< filename.fa.index */

#define FASTA_EUNEXPEOF 1
#define FASTA_EINVAL    2
#define FASTA_ENOBUF    3

#define FASTA_LINEBUFFER_SIZE 4096

typedef struct {
	uint64_t filesize; /**< filesize of the sequence file */
	uint32_t chksum; /**< checksum (CRC-32) of the sequence file */
	uint32_t rcount; /**< expected count of FASTA records */
} FASTA_idxhdr_t;

#include "seqid.h"
#include "trans.h"

typedef struct {
	SeqID_fmt_t seqid_fmt;
	SeqID_t     seqid;
} FASTA_rechdr_t;

#define FASTA_REC_MAGICFL 0xf0fa0000
#define FASTA_REC_FREESEQ 0x00000001 /**< Allowed to free the sequence memory */
#define FASTA_REC_FREEHDR 0x00000002 /**< Allowed to free the headers */
#define FASTA_REC_FREEREC 0x00000004 /**< Allowed to free the memory holding the FASTA record (FASTA_rec_t *) */

typedef struct {
	uint32_t flags;
	uint32_t chksum;

	FASTA_rechdr_t *hdr;     /**< parsed headers */
	uint32_t        hdr_cnt; /**< number of headers */
	void           *hdr_mem; /**< memory where all the headers are stored */
	char           *rec_id;  /**< ID guessed from the header information */

	uint64_t  hdr_start; /**< header file offset */
	uint32_t  hdr_len;   /**< lenght of the header */

	uint64_t  seq_start; /**< sequence file offset */
	uint64_t  seq_rawlen;/**< raw sequence length (including '\n', whitespaces, ...) */
	uint64_t  seq_len;   /**< sequence length */
	uint32_t  seq_lines; /**< sequence data line count */

	uint32_t  seq_linew; /**< line width */
	uint32_t  seq_lastw; /**< last line width */

	char     *seq_mem;  /**< in-memory sequence data */
} FASTA_rec_t;

typedef struct {
	uint32_t fa_options;
	char    *fa_path;

	FILE    *fa_seqFP;
	FILE    *fa_idxFP;

	atrans_t *fa_atr;

	FASTA_rec_t *fa_record;
	uint32_t     fa_rindex;
	uint32_t     fa_rcount;
} FASTA;

FASTA *fasta_open(const char *path, uint32_t options, atrans_t *atr);
uint32_t fasta_count(FASTA *fa);
FASTA_rec_t *fasta_read(FASTA *fa, FASTA_rec_t *dst, uint32_t flags, atrans_t *atr);
int fasta_write(FASTA *fa, FASTA_rec_t *farec);
void fasta_rec_free(FASTA_rec_t *farec);
int fasta_rewind(FASTA *fa);
int fasta_seeko(FASTA *fa, off_t off, int whence);
off_t fasta_tello(FASTA *fa);
void *fasta_apply(FASTA *fa, void * (*func)(FASTA_rec_t *, void *), uint32_t options, void *funcarg);
void fasta_close(FASTA *fa);

#ifdef __cplusplus
}
#endif

#endif /* FASTA_H */

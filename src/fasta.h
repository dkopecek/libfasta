#ifndef FASTA_H
#define FASTA_H

#define FASTA_KEEPOPEN      0x00000001 /**< Keep the FASTA file/index open */
#define FASTA_USEINDEX      0x00000002 /**< Use index, if present; Create it if not present */
#define FASTA_CHKINDEX_FAST 0x00000004 /**< Perform some index integrity check that aren't time consuming */
#define FASTA_CHKINDEX_LONG 0x00000008 /**< Perform a full index integrity check */
#define FASTA_CHKINDEX_NONE 0x00000010 /**< Fully trust the index */
#define FASTA_CHKINDEX      FASTA_CHKINDEX_FAST
#define FASTA_INMEMSEQ      0x00000020 /**< Load all sequence data into memory */
#define FASTA_ONDEMSEQ      0x00000040 /**< Read a sequence into memory on-demand */
#define FASTA_PARALLEL      0x00000080 /**< Perform some operation in parallel, e.g. the _apply operation */

#define FASTA_INDEX_EXT "index" /**< filename.fa.index */

#include "seqid.h"

typedef struct {
	SeqID_fmt_t seqid_fmt;
	SeqID       seqid;
} FASTA_rechdr_t;

typedef struct {
	FASTA_rechdr_t *hdr;     /**< parsed headers */
	uint32_t        hdr_cnt; /**< number of headers */

	uint64_t  hdr_start; /**< header file offset */
	uint32_t  hdr_len;   /**< lenght of the header */

	uint64_t  seq_start; /**< sequence file offset */
	uint64_t  seq_rawlen;/**< raw sequence length (including '\n', whitespaces, ...) */
	uint64_t  seq_len;   /**< sequence length */
	uint32_t  seq_lines; /**< sequence data line count */

	uint32_t  seq_linew; /**< line width */
	uint32_t  seq_lastw; /**< last line width */

	char     *seq_inmem; /**< in-memory sequence data */
} FASTA_rec_t;

typedef struct {
	uint32_t fa_options;
	char    *fa_path;

	FILE    *fa_seqFP;
	FILE    *fa_idxFP;

	FASTA_rec_t *fa_record;
	uint32_t     fa_rcount;
} FASTA;

#endif /* FASTA_H */

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

#define FASTAREC_SEQID_UNKNOWN   0
#define FASTAREC_SEQID_GENBANK   1
#define FASTAREC_SEQID_EMBL      2
#define FASTAREC_SEQID_DDJB      3
#define FASTAREC_SEQID_NBRFPIR   4
#define FASTAREC_SEQID_PRF       5
#define FASTAREC_SEQID_SWISSPROT 6
#define FASTAREC_SEQID_PDB1      7
#define FASTAREC_SEQID_PDB2      8
#define FASTAREC_SEQID_PATENTS   9
#define FASTAREC_SEQID_BBS       10
#define FASTAREC_SEQID_GNL       11
#define FASTAREC_SEQID_NCBIREF   12
#define FASTAREC_SEQID_LOCAL     13

typedef struct {
	uint32_t farec_SeqIDfmt; /**< format of SeqID */

	union {
		char *unknown;
	} farec_SeqID;

	char *farec_sequence;
} FASTAREC;

typedef struct {
	uint32_t fa_options;
	char    *fa_path;
	FILE    *fa_seqFP;
	FILE    *fa_idxFP;
} FASTA;

#endif /* FASTA_H */

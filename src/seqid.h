#ifndef SEQID_H
#define SEQID_H

typedef enum {
	SEQID_UNKNOWN,
	SEQID_GENBANK,
	SEQID_EMBL,
	SEQID_DDJB,
	SEQID_NBRFPIR,
	SEQID_PRF,
	SEQID_SWISSPROT,
	SEQID_PDB1,
	SEQID_PDB2,
	SEQID_PATENTS,
	SEQID_BBS,
	SEQID_GNL,
	SEQID_NCBIREF,
	SEQID_LOCAL
} SeqID_fmt_t;

typedef struct {
} SeqID_genbank_t;

typedef struct {

} SeqID_embl_t;

typedef struct {
} SeqID_ddjb_t;

typedef struct {
} SeqID_nbrfpir_t;

typedef struct {
} SeqID_prf_t;

typedef struct {
} SeqID_swissprot_t;

typedef struct {
} SeqID_pdb1_t;

typedef struct {
} SeqID_pdb2_t;

typedef struct {
} SeqID_patents_t;

typedef struct {
} SeqID_bbs_t;

typedef struct {
} SeqID_gnl_t;

typedef struct {
} SeqID_ncbiref_t;

typedef struct {
} SeqID_local_t;

typedef union {
	char             *unknown;
	SeqID_genbank_t   genbank;
	SeqID_embl_t      embl;
	SeqID_ddjb_t      ddjb;
	SeqID_nbrfpir_t   nbrfpir;
	SeqID_prf_t       prf;
	SeqID_swissprot_t swissprot;
	SeqID_pdb1_t      pdb1;
	SeqID_pdb2_t      pdb2;
	SeqID_patents_t   patents;
	SeqID_bbs_t       bbs;
	SeqID_gnl_t       gnl;
	SeqID_ncbiref_t   ncbiref;
	SeqID_local_t     local;
} SeqID;

#endif /* SEQID_H */

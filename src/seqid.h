#ifndef SEQID_H
#define SEQID_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>

typedef enum {
	SEQID_ERROR,
	SEQID_EMPTY,
	SEQID_UNKNOWN,
	SEQID_GENBANK,
	SEQID_EMBL,
	SEQID_DDBJ,
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
	char *id;
	char *rest;
	char *gi_number;
	char *accession;
	char *locus;
} SeqID_genbank_t;

typedef struct {
	char *id;
	char *rest;
	char *gi_number;
	char *accession;
	char *locus;
} SeqID_embl_t;

typedef struct {
	char *id;
	char *rest;
	char *gi_number;
	char *accession;
	char *locus;
} SeqID_ddbj_t;

typedef struct {
	char *id;
	char *rest;
	char *entry;
} SeqID_nbrfpir_t;

typedef struct {
	char *id;
	char *rest;
	char *name;
} SeqID_prf_t;

typedef struct {
	char *id;
	char *rest;
	char *accession;
	char *name;
} SeqID_swissprot_t;

typedef struct {
	char *id;
	char *rest;
	char *entry;
	char *chain;
} SeqID_pdb1_t;

typedef struct {
	char *id;
	char *rest;
	char *entry;
	char *chain;
	char *pdbid;
	char *sequence;
} SeqID_pdb2_t;

typedef struct {
	char *id;
	char *rest;
	char *country;
	char *number;
} SeqID_patents_t;

typedef struct {
	char *id;
	char *rest;
	char *number;
} SeqID_bbs_t;

typedef struct {
	char *id;
	char *rest;
	char *database;
	char *identifier;
} SeqID_gnl_t;

typedef struct {
	char *id;
	char *rest;
	char *accession;
	char *locus;
} SeqID_ncbiref_t;

typedef struct {
	char *id;
	char *rest;
	char *identifier;
} SeqID_local_t;

typedef struct {
	char *id;
	char *rest;
} SeqID_unknown_t;

typedef struct {
	char *id;
	char *rest;
} SeqID_common_t;

typedef union {
	SeqID_unknown_t   unknown;
	SeqID_genbank_t   genbank;
	SeqID_embl_t      embl;
	SeqID_ddbj_t      ddbj;
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
	SeqID_common_t    common;
} SeqID_t;

SeqID_fmt_t SeqID_parse(char *buffer, size_t buflen, SeqID_t *dst);

#ifdef __cplusplus
}
#endif

#endif /* SEQID_H */

/*
 * Copyright 2011 Daniel Kopecek. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are
 * permitted provided that the following conditions are met:
 *
 *    1. Redistributions of source code must retain the above copyright notice, this list of
 *       conditions and the following disclaimer.
 *
 *    2. Redistributions in binary form must reproduce the above copyright notice, this list
 *       of conditions and the following disclaimer in the documentation and/or other materials
 *       provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY Daniel Kopecek ``AS IS'' AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
 * FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL Daniel Kopecek OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 * Author: Daniel Kopecek <xkopecek@fi.muni.cz>
 *
 */
#ifndef SEQID_H
#define SEQID_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>

        /**
         * Format indentifiers for commonly used SeqID strings.
         */
        typedef enum {
                SEQID_ERROR,     /**< special value for error states */
                SEQID_EMPTY,     /**< special value for uninitialized state */
                SEQID_UNKNOWN,   /**< special value for unknown SeqID formats */
                SEQID_GENBANK,   /**< GenBank */
                SEQID_EMBL,      /**< EMBL Data Library */
                SEQID_DDBJ,      /**< DDBJ, DNA Database of Japan */
                SEQID_NBRFPIR,   /**< NBRF PIR */
                SEQID_PRF,       /**< Protein Research Foundation */
                SEQID_SWISSPROT, /**< SWISS-PROT */
                SEQID_PDB1,      /**< Brookhaven Protein Data Bank */
                SEQID_PDB2,      /**< Brookhaben Protein Data Bank */
                SEQID_PATENTS,   /**< Patents */
                SEQID_BBS,       /**< GenInfo Backbone Id */
                SEQID_GNL,       /**< General database identifier */
                SEQID_NCBIREF,   /**< NCBI Reference Sequence */
                SEQID_LOCAL      /**< Local Sequence Identifier */
        } SeqID_fmt_t;

        /*
         * The following structures define fields of each supported SeqID format
         */
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

        /**
         * Common type used in the source code.
         * Format of the SeqID string determines which
         * field is used to access header fields.
         */
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

        /**
         * Try to parse a string containing a SeqID string.
         * Returns SEQID_ERROR in case of an error, SEQID_UNKNOWN in case
         * of an unknown/unsupported identifier. If a SeqID was recognized,
         * then the relevant format constant is returned.
         */
        SeqID_fmt_t SeqID_parse(char *buffer, size_t buflen, SeqID_t *dst);

#ifdef __cplusplus
}
#endif

#endif /* SEQID_H */

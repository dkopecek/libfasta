#include <stdbool.h>
#include <stddef.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include "helpers.h"
#include "seqid.h"

static bool wspacep(char *s)
{
	while(*s != '\0')
		if (!isspace(*s++))
			return (false);
	return (true);
}

static SeqID_fmt_t SeqID_gi_parse(char *buffer, size_t buflen, SeqID_t *dst)
{
	char *gi_number, *db, *accession, *locus, *rest;
	size_t gi_number_len, db_len, accession_len, locus_len;

        (void)buflen;

	gi_number = strsep(&buffer, "|");

	if (gi_number != NULL) {
		dP("gi_number=\"%s\"\n", gi_number);

		gi_number_len = strlen(gi_number);
		db = strsep(&buffer, "|");

		if (db != NULL) {
			dP("db=\"%s\"\n", db);

			db_len = strlen(db);
			accession = strsep(&buffer, "|");

			if (accession != NULL) {
				dP("accession=\"%s\"\n", accession);

				accession_len = strlen(accession);
				locus = strsep(&buffer, " ");

				if (locus != NULL) {
					dP("locus=\"%s\"\n", locus);

					locus_len = strlen(locus);
					rest      = buffer;

					dP("rest=\"%s\"\n", rest);

					if (strcmp("gb", db) == 0) {
						dst->genbank.id        = gi_number;
						dst->genbank.gi_number = gi_number;
						dst->genbank.accession = accession;
						dst->genbank.locus     = locus;
						dst->genbank.rest      = rest;

						return (SEQID_GENBANK);
					}

					if (strcmp("emb", db) == 0) {
						dst->embl.id        = gi_number;
						dst->embl.gi_number = gi_number;
						dst->embl.accession = accession;
						dst->embl.locus     = locus;
						dst->embl.rest      = rest;

						return (SEQID_EMBL);
					}

					if (strcmp("dbj", db) == 0) {
						dst->ddbj.id        = gi_number;
						dst->ddbj.gi_number = gi_number;
						dst->ddbj.accession = accession;
						dst->ddbj.locus     = locus;
						dst->ddbj.rest      = rest;

						return (SEQID_DDBJ);
					}

					/* XXX: what is |ref| ?
					if (strcmp("ref", db) == 0) {
						dst->ddbj.gi_number = gi_number;
						dst->ddbj.accession = accession;
						dst->ddbj.locus     = locus;
						dst->ddbj.rest      = rest;

						return (SEQID_DDBJ);
					}
					*/

					locus[locus_len] = ' ';
				}
				accession[accession_len] = '|';
			}
			db[db_len] = '|';
		}
		gi_number[gi_number_len] = '|';
	}

	return (SEQID_UNKNOWN);
}

static SeqID_fmt_t SeqID_pir_parse(char *buffer, size_t buflen, SeqID_t *dst)
{
	char *empty, *entry, *rest;
	size_t empty_len;

        (void)buflen;

	empty = strsep(&buffer, "|");

	if (empty != NULL) {
		dP("empty=\"%s\"\n", empty);

		empty_len = strlen(empty);
		entry = strsep(&buffer, " ");

		if (entry != NULL) {
			rest = buffer;

			dP("entry=\"%s\"\n", entry);
			dP("rest=\"%s\"\n", rest);

			dst->nbrfpir.id    = entry;
			dst->nbrfpir.entry = entry;
			dst->nbrfpir.rest  = rest;

			return (SEQID_NBRFPIR);
		}

		empty[empty_len] = '|';
	}

	return (SEQID_UNKNOWN);
}

static SeqID_fmt_t SeqID_prf_parse(char *buffer, size_t buflen, SeqID_t *dst)
{
	char *empty, *name, *rest;
	size_t empty_len;

        (void)buflen;

	empty = strsep(&buffer, "|");

	if (empty != NULL) {
		dP("empty=\"%s\"\n", empty);

		empty_len = strlen(empty);
		name = strsep(&buffer, " ");

		if (name != NULL) {
			rest = buffer;

			dP("name=\"%s\"\n", name);
			dP("rest=\"%s\"\n", rest);

			dst->prf.id   = name;
			dst->prf.name = name;
			dst->prf.rest = rest;

			return (SEQID_PRF);
		}

		empty[empty_len] = '|';
	}

	return (SEQID_UNKNOWN);
}

static SeqID_fmt_t SeqID_sp_parse(char *buffer, size_t buflen, SeqID_t *dst)
{
	char *accession, *name, *rest;
	size_t accession_len;

        (void)buflen;

	accession = strsep(&buffer, "|");

	if (accession != NULL) {
		dP("accession=\"%s\"\n", accession);

		accession_len = strlen(accession);
		name = strsep(&buffer, " ");

		if (name != NULL) {
			rest = buffer;

			dP("name=\"%s\"\n", name);
			dP("rest=\"%s\"\n", rest);

			dst->swissprot.id        = accession;
			dst->swissprot.accession = accession;
			dst->swissprot.name      = name;
			dst->swissprot.rest      = rest;

			return (SEQID_SWISSPROT);
		}

		accession[accession_len] = '|';
	}

	return (SEQID_UNKNOWN);
}

static SeqID_fmt_t SeqID_pdb1_parse(char *buffer, size_t buflen, SeqID_t *dst)
{
	char *entry, *chain, *rest;
	size_t entry_len;

        (void)buflen;

	entry = strsep(&buffer, "|");

	if (entry != NULL) {
		dP("entry=\"%s\"\n", entry);

		entry_len = strlen(entry);
		chain = strsep(&buffer, " ");

		if (chain != NULL) {
			rest = buffer;

			dP("chain=\"%s\"\n", chain);
			dP("rest=\"%s\"\n", rest);

			dst->pdb1.id    = entry;
			dst->pdb1.entry = entry;
			dst->pdb1.chain = chain;
			dst->pdb1.rest  = rest;

			return (SEQID_PDB1);
		}

		entry[entry_len] = '|';
	}

	return (SEQID_UNKNOWN);
}

static SeqID_fmt_t SeqID_pdb2_parse(char *buffer, size_t buflen, char *ftok, size_t tlen, SeqID_t *dst)
{
	char *entry, *pdbid, *chain, *sequence, *rest;
	size_t entry_len, pdbid_len, chain_len;

        (void)tlen;
        (void)buflen;

	entry = strsep(&ftok, ":");

	if (entry != NULL) {
		dP("entry=\"%s\"\n", entry);

		entry_len = strlen(entry);
		pdbid = strsep(&buffer, "|");

		if (pdbid != NULL) {
			dP("pdbid=\"%s\"\n", pdbid);

			pdbid_len = strlen(pdbid);
			chain = strsep(&buffer, "|");

			if (chain != NULL) {
				dP("chain=\"%s\"\n", chain);

				chain_len = strlen(chain);
				sequence = strsep(&buffer, " ");

				if (sequence != NULL) {
					rest = buffer;

					dP("sequence=\"%s\"\n", sequence);
					dP("rest=\"%s\"\n", rest);

					dst->pdb2.id       = pdbid;
					dst->pdb2.entry    = entry;
					dst->pdb2.chain    = chain;
					dst->pdb2.pdbid    = pdbid;
					dst->pdb2.sequence = sequence;
					dst->pdb2.rest     = rest;

					return (SEQID_PDB2);
				}

				chain[chain_len] = '|';
			}
			pdbid[pdbid_len] = '|';
		}
		entry[entry_len] = ':';
	}

	return (SEQID_UNKNOWN);
}

static SeqID_fmt_t SeqID_pat_parse(char *buffer, size_t buflen, SeqID_t *dst)
{
	char *country, *number, *rest;
	size_t country_len;

        (void)buflen;

	country = strsep(&buffer, "|");

	if (country != NULL) {
		dP("country=\"%s\"\n", country);

		country_len = strlen(country);
		number = strsep(&buffer, " ");

		if (number != NULL) {
			rest = buffer;

			dP("number=\"%s\"\n", number);
			dP("rest=\"%s\"\n", rest);

			dst->patents.id      = number;
			dst->patents.country = country;
			dst->patents.number  = number;
			dst->patents.rest    = rest;

			return (SEQID_PATENTS);
		}

		country[country_len] = '|';
	}

	return (SEQID_UNKNOWN);
}

static SeqID_fmt_t SeqID_bbs_parse(char *buffer, size_t buflen, SeqID_t *dst)
{
	char *number, *rest;

        (void)buflen;

	number = strsep(&buffer, " ");

	if (number != NULL) {
		rest = buffer;

		dP("number=\"%s\"\n", number);
		dP("rest=\"%s\"\n", rest);

		dst->bbs.id     = number;
		dst->bbs.number = number;
		dst->bbs.rest   = rest;

		return (SEQID_BBS);
	}

	return (SEQID_UNKNOWN);
}

static SeqID_fmt_t SeqID_gnl_parse(char *buffer, size_t buflen, SeqID_t *dst)
{
	char *database, *identifier, *rest;
	size_t database_len;

        (void)buflen;

	database = strsep(&buffer, "|");

	if (database != NULL) {
		dP("database=\"%s\"\n", database);

		database_len = strlen(database);
		identifier = strsep(&buffer, " ");

		if (identifier != NULL) {
			rest = buffer;

			dP("identifier=\"%s\"\n", identifier);
			dP("rest=\"%s\"\n", rest);

			dst->gnl.id         = identifier;
			dst->gnl.database   = database;
			dst->gnl.identifier = identifier;
			dst->gnl.rest       = rest;

			return (SEQID_GNL);
		}

		database[database_len] = '|';
	}

	return (SEQID_UNKNOWN);
}

static SeqID_fmt_t SeqID_ref_parse(char *buffer, size_t buflen, SeqID_t *dst)
{
	char *accession, *locus, *rest;
	size_t accession_len;

        (void)buflen;

	accession = strsep(&buffer, "|");

	if (accession != NULL) {
		dP("accession=\"%s\"\n", accession);

		accession_len = strlen(accession);
		locus = strsep(&buffer, " ");

		if (locus != NULL) {
			rest = buffer;

			dP("locus=\"%s\"\n", locus);
			dP("rest=\"%s\"\n", rest);

			dst->ncbiref.id        = accession;
			dst->ncbiref.accession = accession;
			dst->ncbiref.locus     = locus;
			dst->ncbiref.rest      = rest;

			return (SEQID_NCBIREF);
		}

		accession[accession_len] = '|';
	}

	return (SEQID_UNKNOWN);
}

static SeqID_fmt_t SeqID_lcl_parse(char *buffer, size_t buflen, SeqID_t *dst)
{
	char *identifier, *rest;

        (void)buflen;

	identifier = strsep(&buffer, " ");

	if (identifier != NULL) {
		rest = buffer;

		dP("identifier=\"%s\"\n", identifier);
		dP("rest=\"%s\"\n", rest);

		dst->local.id         = identifier;
		dst->local.identifier = identifier;
		dst->local.rest       = rest;

		return (SEQID_LOCAL);
	}

	return (SEQID_UNKNOWN);
}

SeqID_fmt_t SeqID_parse(char *buffer, size_t buflen, SeqID_t *dst)
{
	char  *btok;
	size_t tlen;
	SeqID_fmt_t fmt = SEQID_UNKNOWN;
	char *orig = buffer;

        if (buffer == NULL || dst == SEQID_ERROR)
                return (SEQID_ERROR);
        if (buflen == 0)
                return (SEQID_EMPTY);

	btok = strsep(&buffer, "|");

	if (btok != NULL) {
		dP("btok=\"%s\"\n", btok);

		tlen = strlen(btok);

		switch (btok[0]) {
		case 'b':
			/* bbs - GenInfo Backbone Id */
			if (tlen >= 3) {
				if (btok[1] == 'b' && btok[2] == 's' && wspacep(btok + 3))
					fmt = SeqID_bbs_parse(buffer, buflen - tlen, dst);
			}
			break;
		case 'g':
			switch (btok[1]) {
			case 'i':
				/* gi  - (GenBank|EMBL|DDJB) */
				if (wspacep(btok + 2))
					fmt = SeqID_gi_parse(buffer, buflen - tlen, dst);
				break;
			case 'n':
				/* gnl - General database identifier */
				if (tlen >= 3) {
					if (btok[2] == 'l' && wspacep(btok + 3))
						fmt = SeqID_gnl_parse(buffer, buflen - tlen, dst);
				}
				break;
			}
			break;
		case 'l':
			/* lcl - Local sequence identifier */
			if (tlen >= 3) {
				if (btok[1] == 'c' && btok[2] == 'l' && wspacep(btok + 3))
					fmt = SeqID_lcl_parse(buffer, buflen - tlen, dst);
			}
			break;
		case 'p':
			if (tlen >= 3) {
				switch (btok[1]) {
				case 'd':
					/* pdb - Brookhaven Protein Data bank */
					if (btok[2] == 'b' && wspacep(btok + 3))
						fmt = SeqID_pdb1_parse(buffer, buflen - tlen, dst);
					break;
				case 'a':
					/* pat - Patents */
					if (btok[2] == 't' && wspacep(btok + 3))
						fmt = SeqID_pat_parse(buffer, buflen - tlen, dst);
					break;
				case 'i':
					/* pir - NBRF PIR */
					if (btok[2] == 'r' && wspacep(btok + 3))
						fmt = SeqID_pir_parse(buffer, buflen - tlen, dst);
					break;
				case 'r':
					/* prf - Protein Research Foundation */
					if (btok[2] == 'f' && wspacep(btok + 3))
						fmt = SeqID_prf_parse(buffer, buflen - tlen, dst);
					break;
				}
			}
			break;
		case 'r':
			/* ref - NCBI reference sequence */
			if (tlen >= 3) {
				if (btok[1] == 'e' && btok[2] == 'f' && wspacep(btok + 3))
					fmt = SeqID_ref_parse(buffer, buflen - tlen, dst);
			}
			break;
		case 's':
			/* sp - SWISS-PROT */
			if (tlen >= 2) {
				if (btok[1] == 'p' && wspacep(btok + 2))
					fmt = SeqID_sp_parse(buffer, buflen - tlen, dst);
			}
			break;
		default:
			if (strchr(btok, ':') != NULL) {
				/* entry:chain - Brookhaven Protein Data Bank */
				fmt = SeqID_pdb2_parse(buffer, buflen - tlen, btok, tlen, dst);
			}

			/* unknown */
		}

		if (fmt != SEQID_UNKNOWN)
			return (fmt);

		if (buffer != NULL)
			btok[tlen] = '|'; /* undo strsep */
	}

	dst->unknown.id = orig;
	orig = strchr(orig, ' ');

	if (orig != NULL) {
		*orig = '\0';
		dst->unknown.rest = ++orig;
	} else
		dst->unknown.rest = NULL;

	return (SEQID_UNKNOWN);
}

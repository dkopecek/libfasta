#include <stddef.h>
#include <string.h>
#include <ctype.h>
#include <assume.h>
#include "debug.h"
#include "seqid.h"

static bool wspacep(char *s)
{
	while(*s != '\0')
		if (!isspace(*s))
			return (false);
	return (true);
}

static SeqID_fmt_t SeqID_gi_parse(char *buffer, size_t buflen, SeqID_t *dst)
{
	char *gi_number, *db, *accession, *locus, *rest;
	size_t gi_number_len, db_len, accession_len, locus_len, rest_len;

	gi_number = strsep(&buffer, "|");

	if (gi_number != NULL) {
		_D("gi_number=\"%s\"\n", gi_number);

		gi_number_len = strlen(gi_number);
		db = strsep(&buffer, "|");

		if (db != NULL) {
			_D("db=\"%s\"\n", db);

			db_len = strlen(db);
			accession = strsep(&buffer, "|");

			if (accession != NULL) {
				_D("accession=\"%s\"\n", accession);

				accession_len = strsep(&buffer, "|");
				locus = strsep(&buffer, " ");

				if (locus != NULL) {
					_D("locus=\"%s\"\n", locus);

					locus_len = strlen(locus);
					rest      = buffer;
					rest_len  = buffer != NULL ? strlen(buffer) : 0;

					_D("rest=\"%s\"\n", rest);

					if (strcmp("gb", db) == 0) {
						dst->genbank->gi_number = gi_number;
						dst->genbank->accession = accession;
						dst->genbank->locus     = locus;
						dst->genbank->rest      = rest;

						return (SEQID_GENBANK);
					}

					if (strcmp("emb", db) == 0) {
						dst->embl->gi_number = gi_number;
						dst->embl->accession = accession;
						dst->embl->locus     = locus;
						dst->embl->rest      = rest;

						return (SEQID_EMBL);
					}

					if (strcmp("dbj", db) == 0) {
						dst->ddbj->gi_number = gi_number;
						dst->ddbj->accession = accession;
						dst->ddbj->locus     = locus;
						dst->ddbj->rest      = rest;

						return (SEQID_DDBJ);
					}

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
	size_t empty_len, entry_len, rest_len;

	empty = strsep(&buffer, "|");

	if (empty != NULL) {
		_D("empty=\"%s\"\n", empty);

		empty_len = strlen(empty);
		entry = strsep(&buffer, " ");

		if (entry != NULL) {
			_D("entry=\"%s\"\n", entry);

			rest      = buffer;
			rest_len  = buffer != NULL ? strlen(buffer) : 0;

			_D("rest=\"%s\"\n", rest);

			dst->nbrfpir->entry = entry;
			dst->nbrfpir->rest  = rest;

			return (SEQID_NBRFPIR);
		}

		empty[empty_len] = '|';
	}

	return (SEQID_UNKNOWN);
}

static SeqID_fmt_t SeqID_prf_parse(char *buffer, size_t buflen, SeqID_t *dst)
{
	char *empty, *name, *rest;
	size_t empty_len, name_len, rest_len;

	empty = strsep(&buffer, "|");

	if (empty != NULL) {
		_D("empty=\"%s\"\n", empty);

		empty_len = strlen(empty);
		name = strsep(&buffer, " ");

		if (name != NULL) {
			_D("name=\"%s\"\n", name);

			rest      = buffer;
			rest_len  = buffer != NULL ? strlen(buffer) : 0;

			_D("rest=\"%s\"\n", rest);

			dst->prf->name = name;;
			dst->prf->rest = rest;

			return (SEQID_PRF);
		}

		empty[empty_len] = '|';
	}

	return (SEQID_UNKNOWN);
}

static SeqID_fmt_t SeqID_sp_parse(char *buffer, size_t buflen, SeqID_t *dst)
{
	char *accession, *name, *rest;
	size_t accession_len, name_len, rest_len;

	accession = strsep(&buffer, "|");

	if (accession != NULL) {
		_D("accession=\"%s\"\n", accession);

		accession_len = strlen(accession);
		name = strsep(&buffer, " ");

		if (name != NULL) {
			_D("name=\"%s\"\n", name);

			name_len = strlen(name);
			rest     = buffer;
			rest_len = buffer != NULL ? strlen(buffer) : 0;

			_D("rest=\"%s\"\n", rest);

			dst->swissprot->accession = accession;
			dst->swissprot->name      = name;
			dst->swissprot->rest      = rest;

			return (SEQID_SWISSPROT);
		}

		accession[accession_len] = '|';
	}

	return (SEQID_UNKNOWN);
}

static SeqID_fmt_t SeqID_pdb1_parse(char *buffer, size_t buflen, SeqID_t *dst)
{
	char *entry, *chain;
	size_t entry_len, chain_len;

	entry = strsep(&buffer, "|");

	if (entry != NULL) {
		_D("entry=\"%s\"\n", entry);

		entry_len = strlen(entry);
		chain = strsep(&buffer, " ");

		if (chain != NULL) {
			_D("chain=\"%s\"\n", chain);

			chain_len = strlen(chain);
			rest      = buffer;
			rest_len  = buffer != NULL ? strlen(buffer) : 0;

			_D("rest=\"%s\"\n", rest);

			dst->pdb1->entry = entry;
			dst->pdb1->chain = chain;
			dst->pdb1->rest  = rest;

			return (SEQID_PDB1);
		}

		entry[entry_len] = '|';
	}

	return (SEQID_UNKNOWN);
}

static SeqID_fmt_t SeqID_pdb2_parse(char *buffer, size_t buflen, char *ftok, size_t tlen, SeqID_t *dst)
{
	char *entry, *pdbid, *chain, *sequence, *rest;
	size_t entry_len, pdbid_len, chain_len, sequence_len, rest_len;

	entry = strsep(&ftok, ":");

	if (entry != NULL) {
		_D("entry=\"%s\"\n", entry);

		entry_len = strlen(entry);
		pdbid = strsep(&buffer, "|");

		if (pdbid != NULL) {
			_D("pdbid=\"%s\"\n", pdbid);

			pdbid_len = strlen(pdbid);
			chain = strsep(&buffer, "|");

			if (chain != NULL) {
				_D("chain=\"%s\"\n", chain);

				chain_len = strlen(chain);
				sequence = strsep(&buffer, " ");

				if (sequence != NULL) {
					_D("sequence=\"%s\"\n", sequence);

					sequence_len = strlen(sequence);
					rest = buffer;
					rest_len = buffer != NULL ? strlen(buffer) : 0;

					_D("rest=\"%s\"\n", rest);

					dst->pdb2->entry    = entry;
					dst->pdb2->chain    = chain;
					dst->pdb2->pdbid    = pdbid;
					dst->pdb2->sequence = sequence;
					dst->pdb2->rest     = rest;

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
	char *country, *number;
	size_t country_len, number_len;

	country = strsep(&buffer, "|");

	if (country != NULL) {
		_D("country=\"%s\"\n", country);

		country_len = strlen(country);
		number = strsep(&buffer, " ");

		if (number != NULL) {
			_D("number=\"%s\"\n", number);

			number_len = strlen(number);
			rest       = buffer;
			rest_len   = buffer != NULL ? strlen(buffer) : 0;

			_D("rest=\"%s\"\n", rest);

			dst->patents->country = country;
			dst->patents->number  = number;
			dst->patents->rest    = rest;

			return (SEQID_PATENTS);
		}

		country[country_len] = '|';
	}

	return (SEQID_UNKNOWN);
}

static SeqID_fmt_t SeqID_bbs_parse(char *buffer, size_t buflen, SeqID_t *dst)
{
	char *number, *rest;

	number = strsep(&buffer, " ");

	if (number != NULL) {
		_D("number=\"%s\"\n", number);

		rest = buffer;
		rest_len = buffer != NULL ? strlen(buffer) : 0;

		_D("rest=\"%s\"\n", rest);

		dst->bbs->number = number;
		dst->bbs->rest   = rest;

		return (SEQID_BBS);
	}

	return (SEQID_UNKNOWN);
}

static SeqID_fmt_t SeqID_gnl_parse(char *buffer, size_t buflen, SeqID_t *dst)
{
	char *database, *identifier;
	size_t database_len, identifier_len;

	database = strsep(&buffer, "|");

	if (database != NULL) {
		_D("database=\"%s\"\n", database);

		database_len = strlen(database);
		identifier = strsep(&buffer, " ");

		if (identifier != NULL) {
			_D("identifier=\"%s\"\n", identifier);

			identifier_len = strlen(identifier);
			rest           = buffer;
			rest_len       = buffer != NULL ? strlen(buffer) : 0;

			_D("rest=\"%s\"\n", rest);

			dst->gnl->database   = database;
			dst->gnl->identifier = identifier;
			dst->gnl->rest       = rest;

			return (SEQID_GNL);
		}

		database[database_len] = '|';
	}

	return (SEQID_UNKNOWN);
}

static SeqID_fmt_t SeqID_ref_parse(char *buffer, size_t buflen, SeqID_t *dst)
{
	char *accession, *locus;
	size_t accession_len, locus_len;

	accession = strsep(&buffer, "|");

	if (accession != NULL) {
		_D("accession=\"%s\"\n", accession);

		accession_len = strlen(accession);
		locus = strsep(&buffer, " ");

		if (locus != NULL) {
			_D("locus=\"%s\"\n", locus);

			locus_len = strlen(locus);
			rest       = buffer;
			rest_len   = buffer != NULL ? strlen(buffer) : 0;

			_D("rest=\"%s\"\n", rest);

			dst->ncbiref->accession = accession;
			dst->ncbiref->locus     = locus;
			dst->ncbiref->rest      = rest;

			return (SEQID_NCBIREF);
		}

		accession[accession_len] = '|';
	}

	return (SEQID_UNKNOWN);
}

static SeqID_fmt_t SeqID_lcl_parse(char *buffer, size_t buflen, SeqID_t *dst)
{
	char *identifier, *rest;

	identifier = strsep(&buffer, " ");

	if (identifier != NULL) {
		_D("identifier=\"%s\"\n", identifier);

		rest = buffer;
		rest_len = buffer != NULL ? strlen(buffer) : 0;

		_D("rest=\"%s\"\n", rest);

		dst->local->identifier = identifier;
		dst->local->rest       = rest;

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

	assume_r(buffer != NULL, SEQID_ERROR);
	assume_r(buflen  > 0,    SEQID_EMPTY);
	assume_r(dst != NULL,    SEQID_ERROR);

	btok = strsep(&buffer, "|");

	if (btok != NULL) {
		_D("btok=\"%s\"\n", btok);

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

		btok[tlen] = '|'; /* undo strsep */
	}

	dst->unknown = orig;

	return (SEQID_UNKNOWN);
}

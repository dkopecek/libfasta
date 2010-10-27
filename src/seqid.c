#include <assume.h>
#include "seqid.h"

static SeqID_fmt_t SeqID_gi_parse(char *buffer, size_t buflen, SeqID_t *dst)
{
	return (SEQID_UNKNOWN);
}

static SeqID_fmt_t SeqID_pir_parse(char *buffer, size_t buflen, SeqID_t *dst)
{
	return (SEQID_UNKNOWN);
}

static SeqID_fmt_t SeqID_prf_parse(char *buffer, size_t buflen, SeqID_t *dst)
{
	return (SEQID_UNKNOWN);
}

static SeqID_fmt_t SeqID_sp_parse(char *buffer, size_t buflen, SeqID_t *dst)
{
	return (SEQID_UNKNOWN);
}

static SeqID_fmt_t SeqID_pdb1_parse(char *buffer, size_t buflen, SeqID_t *dst)
{
	return (SEQID_UNKNOWN);
}

static SeqID_fmt_t SeqID_pdb2_parse(char *buffer, size_t buflen, char *ftok, size_t tlen, SeqID_t *dst)
{
	return (SEQID_UNKNOWN);
}

static SeqID_fmt_t SeqID_pat_parse(char *buffer, size_t buflen, SeqID_t *dst)
{
	return (SEQID_UNKNOWN);
}

static SeqID_fmt_t SeqID_bbs_parse(char *buffer, size_t buflen, SeqID_t *dst)
{
	return (SEQID_UNKNOWN);
}

static SeqID_fmt_t SeqID_gnl_parse(char *buffer, size_t buflen, SeqID_t *dst)
{
	return (SEQID_UNKNOWN);
}

static SeqID_fmt_t SeqID_ref_parse(char *buffer, size_t buflen, SeqID_t *dst)
{
	return (SEQID_UNKNOWN);
}

static SeqID_fmt_t SeqID_lcl_parse(char *buffer, size_t buflen, SeqID_t *dst)
{
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
		tlen = strlen(btok);

		switch (btok[0]) {
		case 'b':
			/* bbs - GenInfo Backbone Id */
			fmt = SeqID_bbs_parse(buffer, buflen - tlen, dst);
			break;
		case 'g':
			switch (btok[1]) {
			case 'i':
				/* gi  - (GenBank|EMBL|DDJB) */
				fmt = SeqID_gi_parse(buffer, buflen - tlen, dst);
				break;
			case 'n':
				/* gnl - General database identifier */
				fmt = SeqID_gnl_parse(buffer, buflen - tlen, dst);
				break;
			}
			break;
		case 'l':
			/* lcl - Local sequence identifier */
			fmt = SeqID_lcl_parse(buffer, buflen - tlen, dst);
			break;
		case 'p':
			switch (btok[1]) {
			case 'd':
				/* pdb - Brookhaven Protein Data bank */
				fmt = SeqID_pdb1_parse(buffer, buflen - tlen, dst);
				break;
			case 'a':
				/* pat - Patents */
				fmt = SeqID_pat_parse(buffer, buflen - tlen, dst);
				break;
			case 'i':
				/* pir - NBRF PIR */
				fmt = SeqID_pir_parse(buffer, buflen - tlen, dst);
				break;
			case 'r':
				/* prf - Protein Research Foundation */
				fmt = SeqID_prf_parse(buffer, buflen - tlen, dst);
				break;
			}
			break;
		case 'r':
			/* ref - NCBI reference sequence */
			fmt = SeqID_ref_parse(buffer, buflen - tlen, dst);
			break;
		case 's':
			/* sp - SWISS-PROT */
			fmt = SeqID_sp_parse(buffer, buflen - tlen, dst);
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

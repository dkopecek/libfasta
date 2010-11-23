#include <string.h>
#include <stdint.h>
#include <errno.h>
#include "sm_alloc.h"
#include "trans.h"

atrans_t *atrans_new(uint8_t src_width, uint8_t dst_width, uint8_t src_unknown, uint8_t dst_unknown)
{
	atrans_t *atr = NULL;

	if (src_width > 8 || dst_width > 8) {
		errno = ERANGE;
		return (NULL);
	}

	atr = sm_talloc(atrans_t);
	atr->src_width = src_width;
	atr->dst_width = dst_width;

	memset(atr->tr_letter_s2d, src_unknown, sizeof atr->tr_letter_s2d);
	memset(atr->tr_letter_d2s, dst_unknown, sizeof atr->tr_letter_d2s);

	return (atr);
}

void atrans_free(atrans_t *atr)
{
	sm_free(atr);
}

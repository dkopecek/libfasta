#ifndef TRANS_H
#define TRANS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include "assume.h"

typedef struct {
	uint8_t src_width;
	uint8_t dst_width;
	uint8_t tr_letter_s2d[256];
	uint8_t tr_letter_d2s[256];
} atrans_t;

atrans_t *atrans_new(uint8_t src_width, uint8_t dst_width, uint8_t src_unknown, uint8_t dst_unknown);
void      atrans_free(atrans_t *atr);

static inline size_t atrans_s2d_size(atrans_t *atr, size_t size)
{
	register size_t dst_bw;

	assume_d(atr != NULL, 0);

	dst_bw = (((size * 8)/atr->src_width) * atr->dst_width);

	return (dst_bw / 8) + (dst_bw % 8 != 0 ? 1 : 0);
}

static inline size_t atrans_d2s_size(atrans_t *atr, size_t size)
{
	register size_t src_bw;

	assume_d(atr != NULL, 0);

	src_bw = (((size * 8)/atr->dst_width) * atr->src_width);

	return (src_bw / 8) + (src_bw % 8 != 0 ? 1 : 0);
}

static inline void atrans_letter_s2d(atrans_t *atr, uint8_t letter, uint32_t i, uint8_t *dst)
{
	register uint32_t di;

	assume_d(atr != NULL, /* void */);
	assume_d(dst != NULL, /* void */);

	di = i * atr->dst_width;
	dst[di/8] |= (atr->tr_letter_s2d[letter] & ((1 << atr->dst_width) - 1)) << (di % 8);

	return;
}

static inline void atrans_letter_d2s(atrans_t *atr, uint8_t letter, uint32_t i, uint8_t *dst)
{
	register uint32_t di;

	assume_d(atr != NULL, /* void */);
	assume_d(dst != NULL, /* void */);

	di = i * atr->src_width;
	dst[di/8] |= (atr->tr_letter_d2s[letter] & ((1 << atr->src_width) - 1)) << (di % 8);

	return;
}

#ifdef __cplusplus
}
#endif

#endif /* TRANS_H */

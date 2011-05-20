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
#ifndef TRANS_H
#define TRANS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <assert.h>

        /**
         * This structure defines a alphabet translation table.
         * Basically this defines how a 1-8bit number translates
         * to an other 1-8 number using the two uint8_t arrays of
         * size 256 (2^8).
         */
        typedef struct {
                uint8_t src_width;
                uint8_t dst_width;
                uint8_t tr_letter_s2d[256];
                uint8_t tr_letter_d2s[256];
        } atrans_t;

        /**
         * Allocate a new alphabet translation table and initialize required fields.
         */
        atrans_t *atrans_new(uint8_t src_width, uint8_t dst_width, uint8_t src_unknown, uint8_t dst_unknown);

        /**
         * Free an alphabet translation table.
         */
        void atrans_free(atrans_t *atr);

        /**
         * Transform size from the source alphabet to the size in the destination alphabet.
         */
        static inline size_t atrans_s2d_size(atrans_t *atr, size_t size)
        {
                register size_t dst_bw;

                assert(atr != NULL);

                dst_bw = (((size * 8)/atr->src_width) * atr->dst_width);

                return (dst_bw / 8) + (dst_bw % 8 != 0 ? 1 : 0);
        }

        /**
         * Transform size from the destination alphabet to the size in the source alphabet.
         */
        static inline size_t atrans_d2s_size(atrans_t *atr, size_t size)
        {
                register size_t src_bw;

                assert(atr != NULL);

                src_bw = (((size * 8)/atr->dst_width) * atr->src_width);

                return (src_bw / 8) + (src_bw % 8 != 0 ? 1 : 0);
        }

        /**
         * Transform a letter from the source alphabet to the destination alphabet.
         */
        static inline void atrans_letter_s2d(atrans_t *atr, uint8_t letter, uint32_t i, uint8_t *dst)
        {
                register uint32_t di;

                assert(atr != NULL);
                assert(dst != NULL);

                di = i * atr->dst_width;
                dst[di/8] |= (atr->tr_letter_s2d[letter] & ((1 << atr->dst_width) - 1)) << (di % 8);

                return;
        }

        /**
         * Transform a letter from the destination alphabet to the source alphabet.
         */
        static inline void atrans_letter_d2s(atrans_t *atr, uint8_t letter, uint32_t i, uint8_t *dst)
        {
                register uint32_t di;

                assert(atr != NULL);
                assert(dst != NULL);

                di = i * atr->src_width;
                dst[di/8] |= (atr->tr_letter_d2s[letter] & ((1 << atr->src_width) - 1)) << (di % 8);

                return;
        }

#ifdef __cplusplus
}
#endif

#endif /* TRANS_H */

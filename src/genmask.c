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
#include <stdio.h>
#include <string.h>
#include <stdint.h>

/*
 * This is a helper program for generating bitmaps used
 * by the various issomething(char *c) functions.
 */
int main(int argc, char *argv[])
{
	char    *s;
	uint32_t i, b, m;

	uint32_t __mask[] = {
		0x00000000, 0x00000000,	0x00000000, 0x00000000,
		0x00000000, 0x00000000,	0x00000000, 0x00000000
	};

	if (argc != 2) {
		fprintf(stderr, "Usage: %s <str>\n", argv[0]);
		return (1);
	}

	for (s = argv[1], i = 0; i < strlen(s); ++i) {
                /*
                 * Set a bit at a position equal to the
                 * value of s[i].
                 */
		b = s[i] / (sizeof __mask[0] * 8);
		m = s[i] % (sizeof __mask[0] * 8);
		__mask[b] |= 1 << m;
	}

	for (i = 0; i < (sizeof __mask / sizeof __mask[0]); ++i) {
		if (!(i % 4))
			printf("\n");
		printf("0x%08x, ", __mask[i]);
	}

	printf("\n");

	return (0);
}

#include <stdio.h>
#include <string.h>
#include <stdint.h>

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

#ifndef CRC32_H
#define CRC32_H

#include <stddef.h>
#include <stdint.h>

uint32_t crc32 (uint32_t crc, unsigned char *buf, size_t len);

#endif /* CRC32_H */

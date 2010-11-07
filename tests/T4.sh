#!/bin/sh

for file in ${srcdir}/data/*.fa; do
    ./T4_idx_read "${file}"

    if [ $? -ne 0 ]; then
	return 1
    fi
done

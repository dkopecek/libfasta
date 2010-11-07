#!/bin/sh

for file in ${srcdir}/data/*.fa; do
    ./T5_idx_count "${file}"

    if [ $? -ne 0 ]; then
	return 1
    fi
done

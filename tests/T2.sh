#!/bin/sh

for file in ${srcdir}/data/*.fa; do
    ./T2_noidx_read "${file}"

    if [ $? -ne 0 ]; then
	exit 1
    fi
done

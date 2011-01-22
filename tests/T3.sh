#!/bin/sh

for file in ${srcdir}/data/*.fa; do
    ./T3_noidx_trans "${file}"

    if [ $? -ne 0 ]; then
	exit 1
    fi
done

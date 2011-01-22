#!/bin/sh

for file in ${srcdir}/data/*.fa; do
    ./T1_noidx_count "${file}"

    if [ $? -ne 0 ]; then
	exit 1
    fi
done

#!/bin/sh

for file in data/*.fa; do
    ./T2_noidx_read "${file}"

    if [ $? -ne 0 ]; then
	return 1
    fi
done

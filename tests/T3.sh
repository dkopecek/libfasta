#!/bin/sh

for file in data/*.fa; do
    ./T3_noidx_trans "${file}"

    if [ $? -ne 0 ]; then
	return 1
    fi
done

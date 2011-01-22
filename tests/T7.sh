#!/bin/sh

for file in ${srcdir}/data/*.fa; do
    localname="$(basename "T7-${file}")"
    cp "${file}"  "${localname}"
    ./T4_idx_read "${localname}"

    SIZE=$(wc -c "${localname}" | sed -n 's|^\([0-9]*\).*$|\1|p')

echo \
";filesize=$SIZE
;chksum=0x00000000
;rcount=1000000
32 131 12 7112 1049 1333 543 0" > "${localname}.index"

    ./T4_idx_read "${localname}"

    if [ $? -ne 0 ]; then
	exit 1
    fi

    rm -f "${localname}" "${localname}.index"
done

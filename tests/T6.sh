#!/usr/bin/env bash

function chksum() {
	openssl md5 | sed -n 's|^.*\([0-9a-fA-F]*\)$|\1|p'
}

SEQ_LEN_START=3200000
SEQ_LEN_STOP=3200001

LINE_WIDTH_START=1
LINE_WIDTH_STOP=81

SEED_START=1
SEED_STOP=2

SQ_A_path="sequence_A.fa"
SQ_B_path="sequence_B.fa"
FP_path="sequence_params.log"

echo -n "This may take a while..."

for((l=$SEQ_LEN_START; l<=$SEQ_LEN_STOP; l++)); do
    echo
    for((w=$LINE_WIDTH_START; w<=$LINE_WIDTH_STOP; w++)); do
	echo -n .
	for((s=$SEED_START; s<=$SEED_STOP; s++)); do
	    echo "$s 0 1 $l $((l/2)) $w $((w/2))"     > $FP_path
	    ./fastagen $s 0 1 $l $((l/2)) $w $((w/2)) > $SQ_A_path 2> /dev/null
	    ./fastagen $s 0 1 $l $((l/2)) $w $((w/2)) > $SQ_B_path 2> /dev/null

	    SUM1=$(./fastacat $SQ_A_path | chksum)
	    SUM2=$(./fastacat $SQ_B_path | chksum)

	    if [[ "$SUM1" != "$SUM2" ]]; then
		echo "Checksums do not match!"
		echo "Keeping generated files for debugging purposes"
		ls -l $FP_path $SQ_A_path $SQ_B_path
		exit 1
	    fi
	done
    done
done

rm -f $FP_path $SQ_A_path $SQ_B_path
echo
exit 0

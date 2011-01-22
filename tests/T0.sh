#!/bin/sh

find "${srcdir}/data" -name '*.index' -exec rm -f {} \;
find "${srcdir}/"     -name '*.log'   -exec rm -f {} \;
find "${srcdir}/"     -name '*.out'   -exec rm -f {} \;

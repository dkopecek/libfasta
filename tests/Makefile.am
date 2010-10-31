TESTS= T1.sh T2.sh T3.sh
check_PROGRAMS= T1_noidx_count T2_noidx_read T3_noidx_trans

AM_CPPFLAGS= -I$(top_srcdir)/src
LDADD= $(top_builddir)/src/libfasta.la

EXTRA_DIST= data/multi.fa data/multi2.fa data/simple.fa

T1_noidx_count_SOURCES= src/noidx_count.c
T2_noidx_read_SOURCES=  src/noidx_read.c
T3_noidx_trans_SOURCES= src/noidx_trans.c

DISTCLEANFILES= *.log *.out


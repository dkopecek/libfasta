TESTS= T1.sh T2.sh T3.sh T4.sh T5.sh
check_PROGRAMS= T1_noidx_count T2_noidx_read T3_noidx_trans T4_idx_read T5_idx_count fastacat

AM_CPPFLAGS= -I$(top_srcdir)/src
LDADD= $(top_builddir)/src/libfasta-1.0.la

EXTRA_DIST= T1.sh T2.sh T3.sh T4.sh T5.sh data/multi.fa data/multi2.fa data/simple.fa

T1_noidx_count_SOURCES= src/noidx_count.c
T2_noidx_read_SOURCES=  src/noidx_read.c
T3_noidx_trans_SOURCES= src/noidx_trans.c
T4_idx_read_SOURCES= src/idx_read.c
T5_idx_count_SOURCES= src/idx_count.c
fastacat_SOURCES= src/fastacat.c

DISTCLEANFILES= *.log *.out

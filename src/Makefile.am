lib_LTLIBRARIES= libfasta.la

libfasta_la_SOURCES=\
	fasta.c	\
	seqid.c	\
	trans.c	\
	crc32.c	\
	crc32.h	\
	helpers.c \
	helpers.h

libfasta_la_CFLAGS=
libfasta_la_LDFLAGS=\
	-no-undefined \
	-version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE) \
	-Wl,--version-script,$(top_srcdir)/src/symbols.ver

library_includedir=$(includedir)/libfasta

library_include_HEADERS= fasta.h \
			 seqid.h \
			 trans.h

EXTRA_DIST=\
	symbols.ver

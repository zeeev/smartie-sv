.PHONY: all

all: bin/blasr bin/printGaps

bin/blasr:
	cd src/blasr && $(MAKE) && cp alignment/bin/blasr ../../bin/blasr

bin/printGaps: src/htslib/libhts.a
	cd src/print_gaps && $(MAKE)

src/htslib/libhts.a:
	cd src/htslib && $(MAKE)
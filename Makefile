.PHONY: all

all: bin/blasr bin/printGaps

bin/blasr: bin
	cd src/blasr && $(MAKE) && cp alignment/bin/blasr ../../bin/blasr && cp alignment/bin/sawriter ../../bin/sawriter

bin/printGaps: src/htslib/libhts.a bin
	cd src/print_gaps && $(MAKE)

src/htslib/libhts.a:
	cd src/htslib && $(MAKE)

bin:
	mkdir -p bin

clean:
	rm -rf bin
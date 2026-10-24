# Builds the static library libmd.a and the shared library (libmd.so, or
# libmd.dylib on macOS):
#
#     make
#
# builds and runs the example:
#
#     make test
#
# and checks that the work counters (-DMD_COUNT) leave the orderings
# unchanged:
#
#     make count-test

CC ?= cc
CFLAGS ?= -O3
AR ?= ar

ifeq ($(shell uname -s),Darwin)
SHARED = libmd.dylib
SHFLAGS = -dynamiclib -install_name @rpath/libmd.dylib
else
SHARED = libmd.so
SHFLAGS = -shared
endif

all: libmd.a $(SHARED)

md.o: md.c md.h md_impl.h
	$(CC) $(CFLAGS) -fPIC -c md.c -o md.o

libmd.a: md.o
	$(AR) rcs libmd.a md.o

$(SHARED): md.o
	$(CC) $(SHFLAGS) -o $(SHARED) md.o -lm

example: example.c md.h libmd.a
	$(CC) $(CFLAGS) example.c libmd.a -lm -o example

test: example
	./example

count_test: count_test.c md.c md.h md_impl.h
	$(CC) $(CFLAGS) count_test.c md.c -lm -o count_test

count_test_count: count_test.c md.c md.h md_impl.h
	$(CC) $(CFLAGS) -DMD_COUNT count_test.c md.c -lm -o count_test_count

count-test: count_test count_test_count
	./count_test > count_test.out
	./count_test_count > count_test_count.out
	cmp count_test.out count_test_count.out
	@echo "the orderings are the same with and without MD_COUNT"

clean:
	rm -f md.o libmd.a libmd.so libmd.dylib example
	rm -f count_test count_test_count count_test.out count_test_count.out

.PHONY: all test count-test clean

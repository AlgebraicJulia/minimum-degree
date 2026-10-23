# Builds the static library libmd.a and the shared library (libmd.so, or
# libmd.dylib on macOS):
#
#     make
#
# and builds and runs the example:
#
#     make test

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

clean:
	rm -f md.o libmd.a libmd.so libmd.dylib example

.PHONY: all test clean

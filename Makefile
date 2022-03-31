objects = test_main.o cmtx.o
options = -O3

test: test_main.o cmtx.o
	cc -o test $(options) test_main.o cmtx.o

test_main.o: test_main.c cmtx.h
	cc -c $(options) test_main.c

cmtx.o: cmtx.c
	cc -c $(options) cmtx.c

.PHONY: clean

clean:
	rm -rf test $(objects)

objects = test_main.o cmtx.o
options = -O3 

test: test_main.o cmtx.o
	cc -o test test_main.o cmtx.o

mul: test_main.o cmtx.o
	cc -o mul $(options) test_main.o cmtx.o

test_main.o: test_main.c cmtx.h

cmtx.o: cmtx.c
	cc -c $(options) cmtx.c

.PHONY: clean

clean:
	rm -rf mul test $(objects)

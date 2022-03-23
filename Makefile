objects = test_main.o

test: test_main.o
	cc -o test $(objects)

test_main.o: test_main.c cmtx.h

.PHONY: clean

clean:
	rm -rf test $(objects)

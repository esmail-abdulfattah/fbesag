CC = g++
CFLAGS = -Wall -fpic -g -O

all: fbesag.so

fbesag.so: fbesag.o
	$(CC) -shared -o fbesag.so fbesag.o

fbesag.o: fbesag.c
	$(CC) $(CFLAGS) -c -o fbesag.o fbesag.c

clean:
	rm -f fbesag.o fbesag.so


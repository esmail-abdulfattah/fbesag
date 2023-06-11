CC = gcc
CFLAGS = -Wall -fpic -g -O

# Check the platform and set the file extension accordingly
ifeq ($(OS),Windows_NT)
	SHARED_LIB_EXT = dll
else
	SHARED_LIB_EXT = so
endif

all: fbesag.$(SHARED_LIB_EXT)

fbesag.$(SHARED_LIB_EXT): fbesag.o
	$(CC) -shared -o fbesag.$(SHARED_LIB_EXT) fbesag.o

fbesag.o: fbesag.c
	$(CC) $(CFLAGS) -c -o fbesag.o fbesag.c

clean:
	rm -f fbesag.o fbesag.$(SHARED_LIB_EXT)


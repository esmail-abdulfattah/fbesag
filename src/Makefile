CC = g++
CFLAGS = -Wall -fpic -g -O $(shell Rscript -e "Rcpp:::CxxFlags()")
LDFLAGS = -shared $(shell Rscript -e "Rcpp:::LdFlags()")

all: fbesag.so

fbesag.so: fbesag.o
	$(CC) $(LDFLAGS) -o fbesag.so fbesag.o

fbesag.o: fbesag.c
	$(CC) $(CFLAGS) -I"$(R_HOME)/include" -I"$(R_HOME)/library/Rcpp/include" -c -o fbesag.o fbesag.c

clean:
	rm -f fbesag.o fbesag.so


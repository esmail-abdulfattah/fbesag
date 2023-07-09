# How to install fbesag package:

## Option 1: Install directly using devtools:

    library("devtools")
    install_github("esmail-abdulfattah/fbesag")

If it didn't work, try the next option.

## Option 2: Compile it yourself:

The c files are provided in:

https://github.com/esmail-abdulfattah/fbesag/tree/master/inst

The names of the files are cgeneric.h and pbesag.c. You need to compile it using:

Linux:

- g++ -Wall -fpic -g -O -c -o pbesag.o pbesag.c
- g++ -shared -o pbesag.so pbesag.o

Windows:

- g++ -Wall -fpic -g -O -c -o pbesag.o pbesag.c
- g++ -shared -o pbesag.dll pbesag.o

Then run the wrapper.R file in:

https://github.com/esmail-abdulfattah/fbesag/blob/master/inst/wrapper.R

If you are using windows you should replace "pbesag.so" by "pbesag.dll" in the get_pbesag() function in wrapper.R 

## Option 3: Using only R code.




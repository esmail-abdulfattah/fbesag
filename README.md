# How to install fbesag package:

## Option 1: Install directly using devtools:

    library("devtools")
    install_github("esmail-abdulfattah/fbesag")

If it didn't work, try the next option.
_run the R file "run_fbesag_option_1.R" in https://github.com/esmail-abdulfattah/fbesag/tree/master/inst_

## Option 2: Compile it yourself:

The c files are provided in:

https://github.com/esmail-abdulfattah/fbesag/tree/master/inst

The names of the files are cgeneric.h and fbesag.c. You need to compile it using:

Linux:

- gcc -Wall -fpic -g -O -c -o fbesag.o fbesag.c
- gcc -shared -o fbesag.so fbesag.o

Windows:

- gcc -Wall -fpic -g -O -c -o fbesag.o fbesag.c
- gcc -shared -o fbesag.dll fbesag.o

Then run the wrapper.R file in:

https://github.com/esmail-abdulfattah/fbesag/blob/master/inst/wrapper.R

If you are using windows you should replace "fbesag.so" by "fbesag.dll" in the get_fbesag() function in wrapper.R 

_run the R file "run_fbesag_option_2.R" in https://github.com/esmail-abdulfattah/fbesag/tree/master/inst_

## Option 3: Using only R code.

This option could be a bit slower, since you are using fully R code for fbesag, but it is the easiest if the previous two options didn't work. 

_run the R file "run_fbesag_option_2.R" in https://github.com/esmail-abdulfattah/fbesag/tree/master/inst_

# Contact me:

Please if you have any question, email me at esmail.abdulfattah@kaust.edu.sa

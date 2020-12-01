#Ccode: sam_convertor.c findbreaks.c find_localsnp.c
#Cython: lpa_init.pyx co_var.pyx

# the compiler: gcc for C program; the cythonize for Cython
CC = gcc
CFLAGS= -g -Wall -std=c99 -lm


all: convertor breaks localsnp Cython

convertor: 
        $(CC) $(CFLAGS) -o ./bin/sam_convertor ./src/C/sam_convertor.c 
                                                                                                   
breaks:
        $(CC) $(CFLAGS) -o ./bin/findbreaks ./src/C/findbreaks.c   

localsnp:
        $(CC) $(CFLAGS) -o ./lib/find_localsnp ./src/C/find_localsnp.c   

Cython:
        python setup.py build_ext --build-lib ./bin/
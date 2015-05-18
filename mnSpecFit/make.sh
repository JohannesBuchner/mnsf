#!/bin/sh
gcc -shared -o libsynchrotron.so -fPIC synchrotron.c -L/usr/local/lib -lgsl -lgslcblas

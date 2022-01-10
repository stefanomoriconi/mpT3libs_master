#!/bin/bash

for i in *.o; do
	rm $i
done

rm *.so

# Set the gcc version (e.g. gcc-6 or gcc-10) correctly configured with OpenMP
for i in *.c; do
    gcc -c -Wall -Werror -fpic -fopenmp $i
done

gcc -fopenmp -shared -o libmpT3libs.so *.o
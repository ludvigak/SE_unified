#!/bin/sh
icc -Wall -std=c99 -O3 -DBEENAKKER -DVERBOSE -msse4.1 -static ../mex/stresslet_real_rc.c test_real_rc.c -o test_real_rc -vec-report -openmp -lm -g
#gcc   -Wall -std=c99 -O3 -DBEENAKKER -DVERBOSE -static ../mex/stresslet_real_rc.c test_real_rc.c -o test_real_rc -g -pg -lm -ffast-math -fopenmp  -ftree-vectorizer-verbose=1

# -g for using valgrind
# -pg for profiling

# intel: amplxe-gui

# valgrind --tool=callgrind
# callgrind_annotate --auto=yes callgrind.out.pid

# valgrind --tool=cachegrind
# cg_annotate --auto=yes cachegrind.out.pid
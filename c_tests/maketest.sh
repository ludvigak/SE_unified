#!/bin/sh
icc  -Wall -std=c99 -O3 -DBEENAKKER -DVERBOSE -msse2 -static ../mex/stresslet_real_rc.c test_real_rc.c -o test_real_rc -vec-report -g -pg -openmp
#gcc  -Wall -std=c99 -O3 -DBEENAKKER -DVERBOSE -msse4.1 -static ../mex/stresslet_real_rc.c test_real_rc.c -o test_real_rc -g -pg -lm -ffast-math

# -g for using valgrind
# -pg for profiling

# intel: amplxe-gui

# valgrind --tool=callgrind
# callgrind_annotate --auto=yes callgrind.out.pid

# valgrind --tool=cachegrind
# cg_annotate --auto=yes cachegrind.out.pid
CC=icc

CFLAGS=-Wall -Wno-missing-field-initializers -Wno-sign-compare -std=c99
OPTIMFLAGS=-O3 -openmp -std=c99 -fast -fomit-frame-pointer -funroll-all-loops

a.out: FMM_expint.c
	$(CC) -o $@ $< $(CFLAGS) $(OPTIMFLAGS) -lm
clean:
	rm a.out

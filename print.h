// print.h
#ifndef PRINT_H
#define PRINT_H
#include "gf.h"

char *print_binary(unsigned x, unsigned n);
void print_rs_syn(gf c[], int n, int t, gf s[]);
void print_poly(char *name, gf p[], int n);
#endif

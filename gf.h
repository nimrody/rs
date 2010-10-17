// gf.h

#ifndef GF_H
#define GF_H

/* GF over 2^MM
 * The field has 256 elements from 0 to NN
 */
#define MM  8
#define NN  ((1 << MM) - 1)

typedef unsigned gf;

void generate_gf256();

gf gf_add(gf a, gf b);
gf gf_mul(gf a, gf b);
gf gf_div(gf a, gf b);
gf gf_inv(gf a);
gf gf_pow(gf a, int n);
gf gf_poly(gf p[], int n, gf a);

void print_alpha_to();
void print_index_of();
void print_multiply_by(gf m);
void print_roots_table(int t);
void print_inv_table();

gf alpha_to[NN + 1];
gf index_of[NN + 1];

#endif

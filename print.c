#include <assert.h>
#include <stdio.h>
#include <errno.h>
#include "print.h"
#include "gf.h"


char *print_binary(unsigned x, unsigned n)
{
	static char s[4][64];		/* circular buffer */
	static int  si = 0;			/* pointer into next available string */

	int  old_si = si;
	char *p = &s[si][0];
	unsigned mask = 1 << (n-1);

	while (mask) {
		*p++ = x & mask ? '1' : '0';
		mask >>= 1;
	}
	*p = '\0';

	si = (si + 1) % 4;
	return &s[old_si][0];
}


void print_rs_syn(gf c[], int n, int t, gf s[])
{
	int i;
#if defined(PRINT_SYNDROMS_CALC)
	/* evaluate 2t syndroms s_i = c(alpha^i), i=0..2T-1 */
	for(i=0; i<2*t; i++) {

		gf multiply_by = alpha_to[i];
		printf("s[%d]: multiply by %02x\n", i, multiply_by);
		s[i] = 0;
		for(j=n-1; j>=0; j--){
			gf feedback = mul_gf(s[i], multiply_by);
			printf("input[%3d] = %02x, s[%d] = %02x ", j, c[j], i, s[i]);
			s[i] = add_gf(c[j], feedback);
			printf("=> %02x\n", s[i]);
		}
		printf("=== s[%d] = %02x\n", i, s[i]);
	}
#else
	for(i=0; i<2*t; i++) 
		printf("=== s[%d] = %02x\n", i, s[i]);
#endif
}

void print_poly(char *name, gf p[], int n)
{
	int i;

	for(i=0; i<n; i++)
		printf("%s[%3d] = %02x\n", name, i, p[i]);
}


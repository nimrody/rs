#include <assert.h>
#include <stdio.h>
#include <errno.h>
#include "print.h"
#include "gf.h"

/* index->polynomial form conversion table */
gf alpha_to[NN + 1];

/* Polynomial->index form conversion table */
gf index_of[NN + 1];



/* Field generation and arithmetic operations *********************************/


void generate_gf256()
{
	int i;

	/* this function is specifically for GF(2^8) */
	assert(MM==8);

	/*
	 * Generate the standard base. Note that alpha = 0x02 
	 * 00000001 alpha^0
	 * 00000010 alpha^1
	 * 00000100 alpha^2
	 * 00001000 alpha^3
	 * 00010000 alpha^4
	 * 00100000 alpha^5
	 * 01000000 alpha^6
	 * 10000000 alpha^7
	 */
	
	for(i=0; i<MM; i++) 
		alpha_to[i] = 1 << i;

	/*
	 * Our generator polynomial is p(x) = 1+x^2+x^3+x^4+x^8.
	 * But in GF(256) the polynomial p(x) == 0 (mod p(x)), so
	 * we can write 0 = a^2+a^3+a^4+a^8
	 *
	 * If we substitute x=alpha (a for short) and regard the
	 * coefficients as binary digits (addition equals substruction)
	 * we have:
	 *
	 * a^8 = a^0 + a^2 + a^3 + a^4 or more generally, 
	 * a^n = a^(n-8) + a^(n-6) + a^(n-5) + a^(n-4),  n=8..254
	 *
	 * (a^0 .. a^254 are 255 element so a^255=1)
	 *
	 * Where '+' operations denote binary xor, of course.
	 */
	
	for(; i<NN; i++) {
		alpha_to[i] = alpha_to[i-8] ^ alpha_to[i-6] ^ alpha_to[i-5] ^
			          alpha_to[i-4];
	}


	/*
	 * we regard NN as infinity, and alpha^NN stands for the 
	 * zero element which actually has no alpha^n representation
	 */
	alpha_to[NN] = 0;

	/*
	 * now generate the inverse table: numeric representation
	 * to power of alpha (logarithm table)
	 */

	for(i=0; i<=NN; i++) {
		gf s = alpha_to[i];
		index_of[s] = i;
	}
	
}


gf gf_add(gf a, gf b)
{
	return a ^ b;
}

gf gf_mul(gf a, gf b)
{
	if (a == 0 || b == 0)
		return 0;
	else {
		int na = index_of[a];
		int nb = index_of[b];
		int nc = na + nb;

		/* note that this is % NN, since a^NN = 1 (this is
		 * unfortunate, as mod (NN+1) can be accomplished by
		 * a simple bitwise and with 0xFF */
		return alpha_to[nc % NN];
	}

}


gf gf_div(gf a, gf b)
{
	gf c = gf_inv(b);
	gf result = gf_mul(a, c);

	return result;
}

gf gf_inv(gf a)
{
	return gf_pow(a, -1);
}

gf gf_pow(gf a, int n)
{
	int ia;

	if (a==0)
		return 0;

	assert(n>=0 || (n<0 && a!=0));

	while(n<0)
		n += NN;

	ia = index_of[a];
	ia = (ia * n) % NN;

	return alpha_to[ia];	
}


/*
 * Compute p(x)
 *
 * p0 + p1*x + p2*x^2 + .. + p_n-1*x^n-1 =
 * (((p_n-1)*x + p_n-2)*x + .. )*x + p0;
 */

gf gf_poly(gf p[], int n, gf a)
{
	int i = n-1;
	gf acc = 0;

	while(i>=0)  {
		acc = gf_mul(acc, a);
		acc = gf_add(acc, p[i]);
		i--;
	}

	return acc;	
}

/**** Various output routines *******************************/

// just print the binary representation of the GF(2^MM) field
// for debugging
void print_alpha_to()
{
	int i;

	for (i=0; i<=NN; i++)
		printf("a^%3d = %s\n",i, print_binary(alpha_to[i], MM));
}

// The inverse table
void print_index_of()
{
	int i;

	for (i=0; i<=NN; i++)
		printf("%s = a^%d\n",print_binary(i, MM), index_of[i]);
}


void print_multiply_by(gf m)
{
	/* 
	 * to multiply by m we determine its logarithm n, (a^n = m in
	 * GF(256)), and add it to the table of powers.
	 */
	int i;
	for(i=0; i<=NN; i++) {
		printf("%02x x %02x = %02x\n",	i, m, gf_mul(i, m));
	}

}

/**** Output routines that generate Verilog code for the codec ************/


// Print the first 2*t roots of the RS generator polynomial.
// The output is in verilog format and is used by the 
// syndrom calculation unit, the chein search unit and 
// the forney algorithm unit.
void print_roots_table(int t)
{
	int i;

	for (i=0; i<=2*t-1; i++)
		printf("parameter alpha%d = 8'h%02x;\n",i, alpha_to[i]);
}

// Print the GF(2^MM) inverse table. This is in a verilog
// case format and is used to calculate the error magnitude
// (Forney formula)
void print_inv_table()
{
	int i;
	printf("case(x)\n");
	printf("  begin\n");

	for(i=1; i<=NN; i++) {
		int power = index_of[i];
		int inv_power   = (NN - power) % NN;
		gf  inv = alpha_to[inv_power];
		gf  mul = gf_mul(i, inv);
		assert(mul == 1);

		printf("    8'h%02x: y = 8'h%02x; // inv(a^%3d) = a^%d\n", i, inv, power, inv_power);
	}
	
	printf("  end\n");
	printf("endcase\n");
}

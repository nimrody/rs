#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "print.h"
#include "gf.h"



/*
 * Generate the coefficients of a Reed-Solomon generator
 * polynomial of degree n (n = 2t).
 *
 * In our case the roots are a^0 to a^n:
 *   p_2T(x) = (x+1)(x+a^1)(x+a^2)....(x+a^n)
 *
 * So we use the recursion:
 *   p_0(x)  = (x+1)
 *   p_1(x)  = (x+1)(x+a^1) = p_1(x)*x + p(x)*a^1
 *                           `--shift--'
 *   or ..
 *   p_k(x)  = p_(k-1)(x)*x + p_(k-1)(x)*a^k
 *
 *   and k goes from 0 to n-1 (k+1 is the degree of the
 *   resultant polynomial)
 */

void generate_rs(int t, gf g[])
{
	int n=2*t;
	gf  poly[n][n+1];
	int i, k;

	for(i=0; i<=n; i++)
		poly[0][i] = 0;

	poly[0][0] = 1;     /* the x^0 coefficient of x+1 is 1*/
	poly[0][1] = 1;     /* the x^1 coefficient of x+1 is 1*/

	for(k=1; k<n; k++) {
		/*
		 * first we generate the p_(k-1)(x)*x term by
		 * shifting the previous polynomial. we loop
		 * over all terms although many of them are 0
		 * just to make sure we generate the 0's for the
		 * new polynomial as well.
		 */
		for(i=n; i>0; i--)  {
			gf c = poly[k-1][i-1];
			poly[k][i] = c;
		}
		poly[k][0] = 0;

		/* 
		 * now add the p_(k-1)(x)*a^k term coefficient by
		 * coefficient
		 */
		for(i=0; i<=k; i++) {
			gf c  = poly[k-1][i];
			c  = gf_mul(c, alpha_to[k]);  /* c = c*a^k */
			poly[k][i] = gf_add(poly[k][i], c);
		}

	}

	for(i=0; i<n; i++) 
		g[i] = poly[n-1][i];

}

// print_rs_generator
// 
// Print the coefficients of the generator polynomial
// Of a t-error correcting RS code.
//
// The output is part of a Verilog case statement
// for inclusion in the encoder unit.

void print_rs_generator(int t)
{
	int n=2*t;
	gf 	g[n+1];
	int i;

	generate_rs(t, g);

	printf("%d:\n", t-1);
	printf("  begin\n");

	for(i=0;i<(32-2*t); i++) 
		 printf("    g%d <= %3d;\n", i,  0);

	for(i=0; i<n; i++) 
		 printf("    g%d <= %3d;\n", (32-2*t)+i,  g[i]); 

	printf("  end\n");

}


/* 
 * Encoding:
 *
 * m(x) is the data polynomial:
 *   m(x) = m[0] + m[1]x + ... + m[k-1]x^(k-1)
 *
 * where k sysmbols are to be encoded.
 *
 * The degree of the generator polynomial is b = n-k
 *   g(x) = g[0] + g[1]x + ... + g[n-k]x^(n-k)
 *
 * We could encode simply by multiplying m(x) by g(x):
 *   c(x) = m(x)g(x) = c[0] + c[1]x .. + c[n-1]x^(n-1)
 *
 * But this would result a non-systematic encoder (i.e., the output
 * codeword would not contain the input data symbols explicitly).
 *
 * Instead we represent x^b*m(x) as a divisor and reminder of 
 * a division by g(x) - the generator polynomial (this is always
 * possible, of course):
 *
 *   x^b * m(x) = f(x)g(x) + r(x)
 *
 * deg r(x) < deg g(x) = b, of course, so we have:
 *
 *   x^b * m(x) - r(x) = f(x)g(x)
 *
 * So we conclude that x^b*m(x)-r(x) is a codeword since it is
 * a multiple of g(x) and is of the right degree:
 *
 *   deg x^b * m(x) = n-1
 *   deg r(x) < b = n-k <= n-1, so we have
 *   deg x^b*m(x) - r(x) = n-1
 *
 *  Looking at the resulting codeword:
 *
 *  c(x) = x^b * m(x) - x^b*m(x) % g(x) 
 *  x^b*m(x) = m[0]x^b + m[1]x^(b+1) + ... + m[k-1]x^(n-1)
 * -x^b*m(x) % g(x) = r[0] + r[1]x + ... + r[b-1]x^(b-1)
 *  
 *  So if we start by transmitting m[k-1] m[k-2] .. m[0]
 *  we can calculate r(x) and then transmit r[b-1] .. r[0]
 *
 * Shortening:
 *   If we want to shorten the code (transmit less data with the
 *   same number of check bytes, b) we assume m[k-1] .. m[k-u] = 0
 *   and simply skip their transmission.A
 *
 * How do we calculate the reminder in the division x^b * m(x) / g(x)?
 * 
 *
 *     ,---------+----------+--------------+--------<--.
 *     |         |          |              |           |
 *     V         V          V              V           |
 *   .----.    .----.      ...          .-------.   .------.
 *   |x g0|    |x g1|                   |xg[b-1]|   | /g[b]|
 *   `----'    `----'                   `-------'   `------'
 *     |-        |-                        |-          ^
 *     |  .--.   V   .--.                  V   .----.  |
 *     `--|D0|--[+]->|D1|              -->[+]->|Db-1|-[+]--> reminder
 *        `--'       `--'                      `----'  ^     may be extracted
 *                                                     |     when m[0] has
 *           0 .. 0  m[0] ....  m[k] m[k-1]    --------'     been processed
 *
 *   Initially, the registers D0..Db are set to 0.
 *
 *   Notes:
 *     1. In our case g[b]=1 so 1/g[b] is simply a short circuit
 *     2. in GF(2^m) addition and substruction are the
 *        same, so all the '-' signs after the multiplications
 *        by g[0]..g[b-1] are not necessary.
 */

void encode_rs(gf m[], int k, gf g[], int b, gf c[])
{
	gf reg[b];
	int n=k+b;
	int i, j, p;

	/* initialize register to 0 */
	for(i=0; i<=b-1; i++)
		reg[i] = 0;

	/*
	 * p is the output pointer. Note that the output codeword
	 * is c(x) = c[n-1]x^n-1 ... c[0]x^0
	 * (i.e., the elements are arranged according to their
	 * power. That's why c[n-1] comes out first.
	 */
	p = n-1;  

	/*
	 * and now feed the data symbols one by one 
	 * (data is fed from highest power (m[k-1]) to lowest (m[0])
	 */
	for(i=k-1;i>=0; i--) {

		gf feedback = gf_add(reg[b-1], m[i]);

		/*
		 * systematic encoder - the first outputs are the 
		 * input data
		 */

		c[p--] = m[i];

		for(j=b-1; j>0; j--) {
			reg[j] = gf_mul(feedback, g[j]);
			reg[j] = gf_add(reg[j], reg[j-1]);
		}
		/* for the leftmost branch no addition is necessary: */
		reg[0] = gf_mul(feedback, g[0]);
	}

	/* finished with symbols. insert 0's and take out the reminder */
	for(i=b-1; i>=0 ;i--) {
		gf feedback = reg[i];
		c[p--] = feedback;
	}
	assert(p==-1);
}

/* Decoding:
 *
 * First we compute the syndroms. The syndroms are the result of evaluating
 * S[i] = c(alpha^i) where c(x) is the codeword polynomial and i = 0..2T-1 -
 * the roots of the code generator polynomial.
 *
 * If no errors occured, S[i]=0 for all i.
 *
 * For each syndrom calculation we use the following trick to evaluate
 * the polynomial:
 *
 * c(x) = c[n-1]x^n-1 + c[n-2]x^(n-1) + ... + c[0]x^0
 *
 * c_n      = 0
 * c_n-1(x) = c_n*x + c[n-1] = c[n-1]
 * c_n-2(x) = c_n-1(x)*x + c[n-2] = c[n-1]x + c[n-2]
 * c_n-3(x) = c_n-2(x)*x + c[n-3] = c[n-1]x^2 + c[n-2]*x + c[n-3]
 * ..
 * c(x)  = c_0(x) = c_1(x)*x + c[0]
 *
 * Or in hardware:
 *                   GF(256) multiplication by x        
 *                            .---.
 *                        .---| x |<--.
 *                        |   `---'   |
 *                        V   .----.  |
 * c[0]...c[n-2] c[n-1]->[+]->| Di |--'
 *                            `----'
 *
 * The register is initialized to 0 and then new symbols are
 * shifted in as they arrive.
 *
 * Adaptation to different RS(n,k):
 *   For n!=255
 *	   Same operation. Feed only the actual codeword
 *	   (no need to zero pad since the zeros would come
 *     *before* c[n] and will not affect the syndroms.
 *	
 *	 For t!=16
 *     Just make sure the extra syndroms are forced 
 *     to zero when transferred to the key equation
 *     solver.
 *
 */

void rs_syn(gf c[], int n, int t, gf s[])
{
	int i, j;

	/* evaluate 2t syndroms s_i = c(alpha^i), i=0..2T-1 */
	for(i=0; i<2*t; i++) {
		gf multiply_by = alpha_to[i];
		s[i] = 0;
		for(j=n-1; j>=0; j--){
			gf feedback = gf_mul(s[i], multiply_by);
			s[i] = gf_add(c[j], feedback);
		}
	}
	for(; i<2*16; i++)
		s[i]=0;
}

/* 
 *  Key equation solver (KES)

	 The PE1 processor: (i=0..2t-1)
                 ^wh_i
                 |
	             | .---.
	dhat_i[r]  <-+-| D |<---[+]----[X]<---+-- dhat_i+1[r]
                   `---'     ^      ^     |
                             |      |     |
	                         |      `-----)-- gamma[r]
                             |            |
	                        [X]<----------)-- delta[r]
                             |            |
					         |            |
					         `-----.      |
							       |      |
	                          .----)------'
	          .---.    ,----./ 1   |
	th_i[r] .-| D |<--< MUX |      |
			| `---'    `----'\_0   |
			|            |     |   |
			`------->----)-----+---'
				         |
						 ^
						 MC[r]
	

	  And the DC block:
                                                   ,---------.
	 ,------------------------------------------+->| control |
	 |                                          |  |         |
     |   wh0     wh1       wh_t-1               |  |         |
     |   ^       ^         ^                    |  `---------'
	 | ,----.   ,----.            ,----.        |    |  |
	 `-|PE1 |<--|PE1 |<--  ... <--|PE1 |<-- 0   |    |  V MC[r]
	   |0   |   |1   |            |2t-1|        |    | 
	   `----'   `----'            `----'        |    V gamma[r]
                                                V
	                                            delta[r]

		Initialization:												
		PE1_i = s_i, i=0..2t-1

		Output:
		wh0 .. wh_t-1 - error locator polynomial (after 2t iterations)


	 The PE0 processor: (i=t..0)
                                    ^ lambda_i
	                                |
                ,-------------------+-----.
				|                   |     |
	            |  .---.            V     |
	            `--| D |<---[+]----[X]    |
                   `---'     ^      ^     |
                             |      |     |
	                         |      `-----)-- gamma[r]
                             |            |
	                        [X]<----------)-- delta[r]
                             |            |
					         |            |
					         `-----.      V
							       |      |
	                          .----)------'
	           .---.   ,----./ 1   |
	B_i[r]  <--| D |<-< MUX |      |
			   `---'   `----'\ 0   |
			             |    `----+--------- B_i-1[r]
				         |
						 ^
						 MC[r]

	
		And the ELU block:

		 ^lambda_t  ^lambda_t-1      ^lambda1 ^lambda0
	   ,----.     ,----.            ,----.   ,----.
	   |PE0	|     |PE0 |            |PE0 |   |PE0 |
       |t   |<----|t-1 |<-- ...  <--|1   |<--|0   |<- 0
	   `----'     `----'            `----'   `----'

		Initialization:
		PE0_i = 0, i=1..t
		PE0_0 = 1

		Output:
		The error locator polynomial (after 2t iterations)
	

	Adaptation to N!=255, T!=16:

  	  No need to change anything. Just make sure
	  the extra syndroms were zeroed by the syndrom calculation
	  unit.

	  The number of iterations is of course 2*t (and not 2*16)
	
	  The extra lambda[] and omega[] symbols will be
	  zero and will not affect the forney/chien modules.
 *
 */

#define TT 16

void rs_kes(int n, int t, gf s[], gf l[TT+1], gf w[TT])
{
	int i;

   	/* r is the iteration count of the irBM algorithm */	
	int	r;

	/* PE1[i][r] state variables of PE1_i at iteration r */
	gf dhat[2*TT][2*t+1];
	gf th  [2*TT][2*t+1];

	/* PE0[i][r] state variables of PE0_i at iteration r */
	gf B     [TT+1][2*t+1];
	gf lambda[TT+1][2*t+1];

	/* control block signals (not all are registers!) */
	int kcnt [2*TT+1];
	gf  gamma[2*TT+1];
	gf  delta[2*TT+1];
	int MC   [2*TT+1];



	/* r is the iteration count. r = 0.. 2t-1 */
	r = 0; 

	/* initialize PE1 array (DC block) */
	for(i=0; i<2*TT; i++) {
		dhat[i][r] = s[i];
		th[i][r]   = s[i];
	}

	/* initialize PE0 array (ELU block */
	for(i=TT; i>0; i--)  {
		B[i][r]      = 0;
		lambda[i][r] = 0;
	}
	B[i][r] = 1;
	lambda[i][r] = 1;


	/* initialize the control block */
	kcnt[r] = 0;
	gamma[r] = 1;

	for(; r<2*t; r++) {

		/* control block */

		delta[r] = dhat[0][r];
		if (delta[r] != 0 && kcnt[r] >= 0)	
			MC[r] = 1;	
		else
			MC[r] = 0;

		if (MC[r] == 1) {
			gamma[r+1] = delta[r];
			kcnt[r+1] = -kcnt[r] - 1;
		} else {
			gamma[r+1] = gamma[r];
			kcnt[r+1] = kcnt[r] + 1;
		}

		//printf("%3d> delta=%02x gamma=%02x k=%d MC=%d\n", r, delta[r], gamma[r], kcnt[r], MC[r]);

		/* ELU block (PE0 processor) */

		for(i=0; i<=TT; i++) {
			gf m1, m2, prev;

			prev = (i > 0) ? B[i-1][r] : 0;

			if (MC[r])
				B[i][r+1] = lambda[i][r];
			else
				B[i][r+1] = prev;

			m1 = gf_mul( lambda[i][r], gamma[r] );
			m2 = gf_mul( prev, delta[r] );

			lambda[i][r+1] = gf_add(m1, m2);

		}

		/* DC block (PE1 processor) */
		

		for(i=0; i<2*TT; i++) {
			gf m1, m2, prev;

			prev = (i < 2*TT-1) ? dhat[i+1][r] : 0;

			m1 = gf_mul( prev, gamma[r] );
			m2 = gf_mul( th[i][r], delta[r] );

			dhat[i][r+1] = gf_add(m1, m2);

			if (MC[r])
				th[i][r+1] = prev;
			else
				th[i][r+1] = th[i][r];	
		}	

	}

	for(i=0; i<=TT; i++)
		l[i] = lambda[i][r];

	for(i=0; i<TT; i++)
		w[i] = dhat[i][r];

}

void decode_rs(gf c[], int n, int t, int *uncorrectable)
{
	int i, roots, deg_lambda;

	gf syndrom[2*TT];
	gf syndrom_check[2*TT];
	gf lambda[TT+1];
	gf omega[TT];
  
	gf lambda_odd[TT+1];

	// compute the syndroms: s[0]..s[2t-1] 
	rs_syn(c, n, t, syndrom);
	
	//print_poly("syn", syndrom, 2*16);

	// key equation solver
	rs_kes(n, t, syndrom, lambda, omega);

	printf("\n");
	//print_poly("lambda", lambda, TT+1);
	//print_poly("omega ", omega , TT);

	// chien search

	// Adaptaion for N!=255, T!=16
	// No need to change anything. However z^(2*t) is
	// computed with the real t of course.
	
	deg_lambda = 0;
	for(i=0; i<TT+1; i++) {
		lambda_odd[i] = (i%2 == 1) ? lambda[i] : 0;
		if (lambda[i] != 0) 
			deg_lambda = i;
	}

	/* 
	 * Calculate lambda[alpha^-i], i = n-1 .. 0
	 *
	 * lambda[a^-i] = l0 + l1*a^-i + l2*a^-i ... 
	 *
	 * If n=NN, i.e. n=2^m-1=255,
	 * alpha^-(n-1) = alpha^-(NN-1) = alpha^(NN-(NN-1)) =  alpha^1 
	 *
	 * If, however n<NN: n = NN-s
	 * alpha^-(n-1) = alpha^(NN-n+1) = alpha^(s+1)
	 *
	 * I.e., we must calculate alpha^s first
	 */
	roots = 0;
	for(i=n-1; i>=0; i--) {
		gf root = alpha_to[(255-i)%255];
		gf result;
		result = gf_poly(lambda, TT+1, root);

		if (result == 0) {
			gf den = gf_poly(lambda_odd, TT+1, root);
			gf num = gf_poly(omega, TT, root);
			gf z2t = gf_pow(root, 2*t);
			gf result = gf_mul(z2t, num);
			
			result = gf_div(result, den);

			//printf("lambda(alpha^%3d)=0 ", i);
			//printf("error=%02x\n", result);

			c[i] ^= result;

			roots++;
		}
	}

	*uncorrectable = (deg_lambda != roots);

	{
		int e=0;
		rs_syn(c, n, t, syndrom_check);
		for(i=0; i<2*t; i++)
			e += syndrom[i] != syndrom_check[i];

		if (!uncorrectable && e>0) {
			printf("----> could be detected (%d syndroms do not match\n", e);
			assert(0);
		}
	}
		
}

void test_encoder()
{
	int n = 255;
	int t = 16;
	int k = n-2*t;
	gf m[k];
	gf c[n];
	gf g[2*t + 1];
	int uncorrectable;
	
	FILE *fp;
	char line[100];
	int i;

	generate_gf256();
	generate_rs(t, g);

	fp = fopen("indata.txt","rt");
	assert(fp != 0);

	for(i=k-1; i>=0; i--)  {
		fgets(line, 100, fp);
		assert(!feof(fp) && errno == 0);
		sscanf(line, "%2x", &m[i]);
	}
	fclose(fp);

	encode_rs(m, k, g, n-k, c);
#if 1
	for(i=0; i<255; i++) {
		printf("i=%d\n", i);
		c[i] = 0;
		decode_rs(c, n, t, &uncorrectable);
//		c[i] = ~c[i];
	}
#else
	//printf("deleting c[0] = %02x\n", c[0]);
	//c[0] = 0;
	printf("deleting c[254] = %02x\n", c[254]);
	c[254] = 0;
	c[253] = 0;
	c[252] = 0;
	c[251] = 0;
	c[250] = 0;
	c[249] = 0;
	c[248] = 0;
	decode_rs(c, n, t, &uncorrectable);
#endif
	fp = fopen("outdatan.txt", "wt");
	for(i=n-1; i>=0; i--){
		fprintf(fp, "%02x\n", c[i]);
	}
	fclose(fp);
	
}

// add up to t errors
void add_errors(gf c[255], int n, int t)
{
	int i, j;

	for(i=0; i<t; i++) {
		j = rand() % n;
		c[j] = rand() & 0xFF;
	}
}


void test_codec()
{
	gf m[255];
	gf c[255];
	gf g[2*16 + 1];
	int n, k, t;
	int uncorrectable, errors;
	int i, j;

	generate_gf256();
	srand(200);

	for(i=0; i<10000; i++) {

		t = rand() % 16 + 1;
		k = rand() % (255-2*t) + 1;

		n = k + 2*t;
//		k=186; n=204; t=9;

		printf("RS(%d, %d, %d) ", n, k, t);

		for(j=k-1; j>=0; j--)
			m[j] = rand() & 0xFF;
		
		generate_rs(t, g);
		encode_rs(m, k, g, n-k, c);

		errors = rand() % (t+10) + 1;
		add_errors(c, n, t);
		decode_rs(c, n, t, &uncorrectable);

		errors = 0;
		for(j=0; j<k; j++)
			errors += m[j] != c[2*t+j];

		if (!uncorrectable && errors>0)
			printf(" %d errors undetected\n", errors);
		else
			printf(" \n");
	}
}

int main()
{
	test_codec();
//	print_inv_table();

//	print_roots_table(32);
	
//	test_decoder(16);
	
//	print_alpha_to();
//	print_index_of();
//	print_multiply_by(2);	

//	print_gf_multiplier();
//	print_roots_table(16);
	return 0;
}

#include "params.h"
#include "reduce.h"
#include "ntt.h"

const int16_t zetas[72] = {
	 -147,  1571,  -722, -1265,   867,   460,  1120,  -467,
	-1556,   556,   532,  -773,    58,  1200,  1722,   570,
	 -491,  1135,    36, -1443,  1664,  -952, -1020,  -599,
	  592,  1622,  1244, -1085,   566, -1369,   268,   781,
	   96,   211,   473,  -234, -1536,  -676,  -170,  -705,
	 -658, -1712,  1488,  1087,  1245,   791,  -875,   -70,
	  287,  -945,  -882,   206, -1421,   716,   838,  1035,
	  966,   -61,   404,   862, -1460,  1266, -1584,  -265,
	  905,  -619,   118,   286,  -194,  1229,  1608, -1669
};

/*************************************************
* Name:        fqmul
*
* Description: Multiplication followed by Montgomery reduction
*
* Arguments:   - int16_t a: first factor
*              - int16_t b: second factor
*
* Returns 16-bit integer congruent to a*b*R^{-1} mod q
**************************************************/
static int16_t fqmul(int16_t a, int16_t b)
{
	return montgomery_reduce((int32_t)a*b);
}

/*************************************************
* Name:        fqinv
*
* Description: Inversion
*
* Arguments:   - int16_t a: first factor a = x mod q
*
* Returns 16-bit integer congruent to x^{-1} * R^2 mod q
**************************************************/
static int16_t fqinv(int16_t a)
{
	int16_t t1,t2,t3;

	t1 = fqmul(a, a);    //10
	t2 = fqmul(t1, t1);  //100
	t2 = fqmul(t2, t2);  //1000
	t3 = fqmul(t2, t2);  //10000

	t1 = fqmul(t1, t2);  //1010

	t2 = fqmul(t1, t3);  //11010
	t2 = fqmul(t2, t2);  //110100
	t2 = fqmul(t2, a);   //110101

	t1 = fqmul(t1, t2);  //111111

	t2 = fqmul(t2, t2);  //1101010
	t2 = fqmul(t2, t2);  //11010100
	t2 = fqmul(t2, t2);  //110101000
	t2 = fqmul(t2, t2);  //1101010000
	t2 = fqmul(t2, t2);  //11010100000
	t2 = fqmul(t2, t2);  //110101000000
	t2 = fqmul(t2, t1);  //110101111111

	return t2;
}

/*************************************************
* Name:        ntt
*
* Description: number-theoretic transform (NTT) in Rq.
*
* Arguments:   - int16_t r[NTRUPLUS_N]: pointer to output vector of elements of Zq
*              - int16_t a[NTRUPLUS_N]: pointer to input vector of elements of Zq
**************************************************/
/*************************************************
* Name:        ntt
*
* Description: number-theoretic transform (NTT) in Rq.
*
* Arguments:   - int16_t r[NTRUPLUS_N]: pointer to output vector of elements of Zq
*              - int16_t a[NTRUPLUS_N]: pointer to input vector of elements of Zq
**************************************************/
void ntt(int16_t r[NTRUPLUS_N])
{
	int16_t t1,t2,t3;
	int16_t zeta1,zeta2;
	
	int k = 1;

	zeta1 = zetas[k++];

	for(int i = 0; i < NTRUPLUS_N/2; i++)
	{
		t1 = fqmul(zeta1, r[i + NTRUPLUS_N/2]);

		r[i + NTRUPLUS_N/2] = r[i] + r[i + NTRUPLUS_N/2] - t1;
		r[i               ] = r[i]                       + t1;
	}

	for(int step = NTRUPLUS_N/6; step >= 4; step = step/3)
	{
		for(int start = 0; start < NTRUPLUS_N; start += 3*step)
		{
			zeta1 = zetas[k++];
			zeta2 = zetas[k++];

			for(int i = start; i < start + step; i++)
			{
				t1 = fqmul(zeta1, r[i +   step]);
				t2 = fqmul(zeta2, r[i + 2*step]);
				t3 = fqmul(-886, t1 - t2);

				r[i + 2*step] = r[i] - t1 - t3;
				r[i +   step] = r[i] - t2 + t3;
				r[i         ] = r[i] + t1 + t2;
			}		
		}
	}

	for(int step = 2; step >= 1; step >>= 1)
	{
		for(int start = 0; start < NTRUPLUS_N; start += (step << 1))
		{
			zeta1 = zetas[k++];

			for(int i = start; i < start + step; i++)
			{
				t1 = fqmul(zeta1, r[i + step]);
				
				r[i + step] = (r[i] - t1);
				r[i       ] = (r[i] + t1);
			}
		}
	}
}

/*
void merge3_ntt(int16_t r[128])
{
	int16_t t;
	int16_t zeta[8];
	int16_t v[8];

	
    for (int i = 0; i < 8; i++)
    {
        zeta[i] = zetas[i+1];
    }
    
	for (int i = 0; i < 16; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			v[j] = r[16*j+i];
		}

		t = fqmul(zeta[0], v[4]);
		v[4] = (v[0] - t);
		v[0] = (v[0] + t);

		t = fqmul(zeta[0], v[5]);
		v[5] = (v[1] - t);
		v[1] = (v[1] + t);

		t = fqmul(zeta[0], v[6]);
		v[6] = (v[2] - t);
		v[2] = (v[2] + t);

		t = fqmul(zeta[0], v[7]);
		v[7] = (v[3] - t);
		v[3] = (v[3] + t);

		t = fqmul(zeta[1], v[2]);
		v[2] = (v[0] - t);
		v[0] = (v[0] + t);

		t = fqmul(zeta[1], v[3]);
		v[3] = (v[1] - t);
		v[1] = (v[1] + t);

		t = fqmul(zeta[2], v[6]);
		v[6] = (v[4] - t);
		v[4] = (v[4] + t);

		t = fqmul(zeta[2], v[7]);
		v[7] = (v[5] - t);
		v[5] = (v[5] + t);

		t = fqmul(zeta[3], v[1]);
		v[1] = (v[0] - t);
		v[0] = (v[0] + t);

		t = fqmul(zeta[4], v[3]);
		v[3] = (v[2] - t);
		v[2] = (v[2] + t);

		t = fqmul(zeta[5], v[5]);
		v[5] = (v[4] - t);
		v[4] = (v[4] + t);

		t = fqmul(zeta[6], v[7]);
		v[7] = (v[6] - t);
		v[6] = (v[6] + t);

		for (int j = 0; j < 8; j++)
		{
			r[16*j+i] = v[j];
		}
	}
}
*/

/*
void merge2_ntt(int16_t r[128])
{
	int16_t t;
	int32_t S0,S1;
	int32_t T0,T1,T2,T3,T4;
	int16_t zeta[3];
	int16_t v[4];

    for (int i = 0; i < 3; i++)
    {
        zeta[i] = zetas[i+1];
    }
    
	for (int i = 0; i < 32; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			v[j] = r[32*j+i];
		}

		S0 = R*v[0];
		S1 = zeta[0]*v[2];

		T0 = S0 + S1;
		T1 = S0 - S1;
		
		T2 = zeta[1]*v[1] + zeta[2]*v[3];
		T3 = zeta[2]*v[1] + zeta[1]*v[3];

		v[0] = montgomery_reduce(T0 + T2);
		v[1] = montgomery_reduce(T0 - T2);
		v[2] = montgomery_reduce(T1 + T3);
		v[3] = montgomery_reduce(T1 - T3);

		for (int j = 0; j < 4; j++)
		{
			r[32*j+i] = v[j];
		}
	}
}
*/
void merge6_ntt(int16_t r[72])
{
	int16_t t1,t2,t3;
	int16_t zeta;
	int16_t v[8];

	int16_t index;

	for (int j = 0; j < 12; j++)
	{
		for (int k = 0; k < 6; k++)
		{
			v[k] = r[12*k+j];
		}

		t1 = fqmul(zetas[1], v[3]);
		v[3] = (v[0] + v[3] - t1);
		v[0] = (v[0] + t1);

		t1 = fqmul(zetas[1], v[4]);
		v[4] = (v[1] + v[4] - t1);
		v[1] = (v[1] + t1);

		t1 = fqmul(zetas[1], v[5]);
		v[5] = (v[2] + v[5] - t1);
		v[2] = (v[2] + t1);

		t1 = fqmul(zetas[2], v[1]);
		t2 = fqmul(zetas[3], v[2]);
		t3 = fqmul(-886, t1 - t2);

		v[2] = v[0] - t1 - t3;
		v[1] = v[0] - t2 + t3;
		v[0] = v[0] + t1 + t2;

		t1 = fqmul(zetas[4], v[4]);
		t2 = fqmul(zetas[5], v[5]);
		t3 = fqmul(-886, t1 - t2);

		v[5] = v[3] - t1 - t3;
		v[4] = v[3] - t2 + t3;
		v[3] = v[3] + t1 + t2;

		for (int k = 0; k < 6; k++)
		{
			r[12*k+j] = v[k];
		}
	}

	index = 6;

	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			for (int k = 0; k < 6; k++)
			{
				v[k] = r[2*k+j+12*i];
			}

			t1 = fqmul(zetas[6+2*i], v[2]);
			t2 = fqmul(zetas[7+2*i], v[4]);
			t3 = fqmul(-886, t1 - t2);

			v[4] = v[0] - t1 - t3;
			v[2] = v[0] - t2 + t3;
			v[0] = v[0] + t1 + t2;

			t1 = fqmul(zetas[6+2*i], v[3]);
			t2 = fqmul(zetas[7+2*i], v[5]);
			t3 = fqmul(-886, t1 - t2);

			v[5] = v[1] - t1 - t3;
			v[3] = v[1] - t2 + t3;
			v[1] = v[1] + t1 + t2;			

			t1 = fqmul(zetas[18+3*i], v[1]);
			v[1] = (v[0] - t1);
			v[0] = (v[0] + t1);

			t1 = fqmul(zetas[19+3*i], v[3]);
			v[3] = (v[2] - t1);
			v[2] = (v[2] + t1);

			t1 = fqmul(zetas[20+3*i], v[5]);
			v[5] = (v[4] - t1);
			v[4] = (v[4] + t1);

			for (int k = 0; k < 6; k++)
			{
				r[2*k+j+12*i] = v[k];
			}
		}
	}

	index = 36;

	for(int start = 0; start < NTRUPLUS_N; start += 2)
	{
		zeta = zetas[index++];

		for(int i = start; i < start + 1; i++)
		{
			t1 = fqmul(zeta, r[i + 1]);
			
			r[i + 1] = (r[i] - t1);
			r[i    ] = (r[i] + t1);
		}
	}
}

/*************************************************
* Name:        invntt
*
* Description: inverse number-theoretic transform in Rq and
*              multiplication by Montgomery factor R = 2^16.
*
* Arguments:   - int16_t r[NTRUPLUS_N]: pointer to output vector of elements of Zq
*              - int16_t a[NTRUPLUS_N]: pointer to input vector of elements of Zq
**************************************************/
void invntt(int16_t r[NTRUPLUS_N], const int16_t a[NTRUPLUS_N])
{
	int16_t t1, t2, t3;
	int16_t zeta1, zeta2;
	int k = 71;

	for(int i = 0; i < NTRUPLUS_N; i++)
	{
		r[i] = a[i];
	}

	for(int step = 4; step <= 16; step <<= 1)
	{
		for(int start = 0; start < NTRUPLUS_N; start += (step << 1))
		{
			zeta1 = zetas[k--];

			for(int i = start; i < start + step; i++)
			{
				t1 = r[i + step];

				r[i + step] = fqmul(zeta1,  t1 - r[i]);
				r[i       ] = barrett_reduce(r[i] + t1);
			}
		}
	}

	for(int step = 32; step <= NTRUPLUS_N/6; step = 3*step)
	{
		for(int start = 0; start < NTRUPLUS_N; start += 3*step)
		{
			zeta2 = zetas[k--];
			zeta1 = zetas[k--];

			for(int i = start; i < start + step; i++)
			{
				t1 = fqmul(-886,  r[i +   step] - r[i]);
				t2 = fqmul(zeta1, r[i + 2*step] - r[i]        + t1);
				t3 = fqmul(zeta2, r[i + 2*step] - r[i + step] - t1);

				r[i         ] = barrett_reduce(r[i] + r[i + step] + r[i + 2*step]);
				r[i +   step] = t2;			
				r[i + 2*step] = t3;
			}
		}
	}

	for(int i = 0; i < NTRUPLUS_N/2; i++)
	{
		t1 = r[i] + r[i + NTRUPLUS_N/2];
		t2 = fqmul(-1665, r[i] - r[i + NTRUPLUS_N/2]);

		r[i               ] = fqmul(-132, t1 - t2);
		r[i + NTRUPLUS_N/2] = fqmul(-264, t2);
	}
}

/*************************************************
* Name:        baseinv
*
* Description: Inversion of polynomial in Zq[X]/(X^4-zeta)
*              used for inversion of element in Rq in NTT domain
*
* Arguments:   - int16_t r[4]: pointer to the output polynomial
*              - const int16_t a[4]: pointer to the input polynomial
*              - int16_t zeta: integer defining the reduction polynomial
**************************************************/
int baseinv(int16_t r[4], const int16_t a[4], int16_t zeta)
{
	int16_t t0, t1, t2;
	
	t0 = montgomery_reduce(a[2]*a[2] - (a[1]*a[3] << 1));
	t1 = montgomery_reduce(a[3]*a[3]);
	t0 = montgomery_reduce(a[0]*a[0] + t0*zeta);
	t1 = montgomery_reduce(((a[0]*a[2]) << 1) - a[1]*a[1] - t1*zeta);

	t2 = montgomery_reduce(t1*t1);
	t2 = montgomery_reduce(t0*t0 - t2*zeta);

	if(t2 == 0) return 1;

	t2 = fqinv(t2);
	t0 = fqmul(t0,t2);
	t1 = fqmul(t1,t2);
	t2 = fqmul(t1,zeta);
	
	r[0] = montgomery_reduce(a[0]*t0 - a[2]*t2);
	r[1] = montgomery_reduce(a[3]*t2 - a[1]*t0);
	r[2] = montgomery_reduce(a[2]*t0 - a[0]*t1);
	r[3] = montgomery_reduce(a[1]*t1 - a[3]*t0);

	return 0;
}

/*************************************************
* Name:        basemul
*
* Description: Multiplication of polynomials in Zq[X]/(X^4-zeta)
*              used for multiplication of elements in Rq in NTT domain
*
* Arguments:   - int16_t r[4]: pointer to the output polynomial
*              - const int16_t a[4]: pointer to the first factor
*              - const int16_t b[4]: pointer to the second factor
*              - int16_t zeta: integer defining the reduction polynomial
**************************************************/
void basemul(int16_t r[4], const int16_t a[4], const int16_t b[4], int16_t zeta)
{
	r[0] = montgomery_reduce(a[1]*b[3]+a[2]*b[2]+a[3]*b[1]);
	r[1] = montgomery_reduce(a[2]*b[3]+a[3]*b[2]);
	r[2] = montgomery_reduce(a[3]*b[3]);

	r[0] = montgomery_reduce(r[0]*zeta+a[0]*b[0]);
	r[1] = montgomery_reduce(r[1]*zeta+a[0]*b[1]+a[1]*b[0]);
	r[2] = montgomery_reduce(r[2]*zeta+a[0]*b[2]+a[1]*b[1]+a[2]*b[0]);
	r[3] = montgomery_reduce(a[0]*b[3]+a[1]*b[2]+a[2]*b[1]+a[3]*b[0]);
}

/*************************************************
* Name:        basemul_add
*
* Description: Multiplication then addition of polynomials in Zq[X]/(X^4-zeta)
*              used for multiplication of elements in Rq in NTT domain
*
* Arguments:   - int16_t c[4]: pointer to the output polynomial
*              - const int16_t a[4]: pointer to the first factor
*              - const int16_t b[4]: pointer to the second factor
*              - const int16_t c[4]: pointer to the third factor
*              - int16_t zeta: integer defining the reduction polynomial
**************************************************/
void basemul_add(int16_t r[4], const int16_t a[4], const int16_t b[4], const int16_t c[4], int16_t zeta)
{
	r[0] = montgomery_reduce(a[1]*b[3]+a[2]*b[2]+a[3]*b[1]);
	r[1] = montgomery_reduce(a[2]*b[3]+a[3]*b[2]);
	r[2] = montgomery_reduce(a[3]*b[3]);

	r[0] = montgomery_reduce(c[0]*(-147)+r[0]*zeta+a[0]*b[0]);
	r[1] = montgomery_reduce(c[1]*(-147)+r[1]*zeta+a[0]*b[1]+a[1]*b[0]);
	r[2] = montgomery_reduce(c[2]*(-147)+r[2]*zeta+a[0]*b[2]+a[1]*b[1]+a[2]*b[0]);
	r[3] = montgomery_reduce(c[3]*(-147)+a[0]*b[3]+a[1]*b[2]+a[2]*b[1]+a[3]*b[0]);
}

#include "params.h"
#include "reduce.h"
#include "ntt.h"

const int16_t zetas[128] = {
-1044, -758, -359, -1517, 1493, 1422, 287, 202, -171, 622, 1577, 182, 962, -1202, -1474, 1468, 573, -1325, 264, 383, -829, 1458, -1602, -130, -681, 1017, 732, 608, -1542, 411, -205, -1571, 1223, 652, -552, 1015, -1293, 1491, -282, -1544, 516, -8, -320, -666, -1618, -1162, 126, 1469, -853, -90, -271, 830, 107, -1421, -247, -951, -398, 961, -1508, -725, 448, -1065, 677, -1275, -1103, 430, 555, 843, -1251, 871, 1550, 105, 422, 587, 177, -235, -291, -460, 1574, 1653, -246, 778, 1159, -147, -777, 1483, -602, 1119, -1590, 644, -872, 349, 418, 329, -156, -75, 817, 1097, 603, 610, 1322, -1285, -1465, 384, -1215, -136, 1218, -1335, -874, 220, -1187, -1659, -1185, -1530, -1278, 794, -1510, -854, -870, 478, -108, -308, 996, 991, 958, -1460, 1522, 1628
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
* Name:        ntt
*
* Description: number-theoretic transform (NTT) in Rq.
*
* Arguments:   - int16_t r[NTRUPLUS_N]: pointer to output vector of elements of Zq
*              - int16_t a[NTRUPLUS_N]: pointer to input vector of elements of Zq
**************************************************/
void ntt(int16_t r[NTRUPLUS_N])
{
	int16_t t;
	int16_t zeta;
	
	int k = 1;

	for(int step = NTRUPLUS_N/2; step >= NTRUPLUS_N/4; step >>= 1)
	{
		for(int start = 0; start < NTRUPLUS_N; start += (step << 1))
		{
			zeta = zetas[k++];

			for(int i = start; i < start + step; i++)
			{
				t = fqmul(zeta, r[i + step]);
				
				r[i + step] = (r[i] - t);
				r[i       ] = (r[i] + t);
			}
		}
	}
}

void merge3_ntt(int16_t r[128])
{
	int16_t t;
	int16_t zeta[8];
	int16_t v[8];

	
    for (int i = 0; i < 8; i++)
    {
        zeta[i] = zetas[i+1];
    }
    
	for (int i = 0; i < 32; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			v[j] = r[32*j+i];
		}

		t = fqmul(zeta[0], v[2]);
		v[2] = (v[0] - t);
		v[0] = (v[0] + t);

		t = fqmul(zeta[0], v[3]);
		v[3] = (v[1] - t);
		v[1] = (v[1] + t);

		t = fqmul(zeta[1], v[1]);
		v[1] = (v[0] - t);
		v[0] = (v[0] + t);

		t = fqmul(zeta[2], v[3]);
		v[3] = (v[2] - t);
		v[2] = (v[2] + t);

		for (int j = 0; j < 4; j++)
		{
			r[32*j+i] = v[j];
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
		v[1] = barrett_reduce(v[0] - t);
		v[0] = barrett_reduce(v[0] + t);

		t = fqmul(zeta[4], v[3]);
		v[3] = barrett_reduce(v[2] - t);
		v[2] = barrett_reduce(v[2] + t);

		t = fqmul(zeta[5], v[5]);
		v[5] = barrett_reduce(v[4] - t);
		v[4] = barrett_reduce(v[4] + t);

		t = fqmul(zeta[6], v[7]);
		v[7] = barrett_reduce(v[6] - t);
		v[6] = barrett_reduce(v[6] + t);

		for (int j = 0; j < 8; j++)
		{
			r[16*j+i] = v[j];
		}
	}
}
*/
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
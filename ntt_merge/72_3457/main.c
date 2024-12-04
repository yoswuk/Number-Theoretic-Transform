#include <stdio.h>
#include <stdint.h>

#include "params.h"
#include "reduce.h"
#include "ntt.h"

#define TEST_LOOP 10000

static inline uint64_t cpucycles(void) 
{
	uint64_t result;
	
	__asm__ volatile ("rdtsc; shlq $32,%%rdx; orq %%rdx,%%rax"
	: "=a" (result) : : "%rdx");
	
	return result;
}

static void TEST_NTT()
{
    int16_t a[NTRUPLUS_N];
    int16_t b[NTRUPLUS_N];

    for (int i = 0; i < NTRUPLUS_N; i++)
    {
        a[i] = i;
    }
    
    ntt(a);

    printf("ntt\n");
    for (int i = 0; i < NTRUPLUS_N; i++)
    {
        printf("%d ", ((a[i] % NTRUPLUS_Q) + NTRUPLUS_Q) % NTRUPLUS_Q);
        if(i % 32 == 31) printf("\n");
    }
    printf("\n");
    
    for (int i = 0; i < NTRUPLUS_N; i++)
    {
        a[i] = i;
    }

    merge6_ntt(a);

    printf("merge6_ntt\n");
    for (int i = 0; i < NTRUPLUS_N; i++)
    {
        printf("%d ", ((a[i] % NTRUPLUS_Q) + NTRUPLUS_Q) % NTRUPLUS_Q);
        if(i % 32 == 31) printf("\n");
    }
    printf("\n");
}

static void TEST_REDUCE()
{
    for (int16_t i = -32768; i < 32767; i++)
    {
        int16_t r1 = ((((i % NTRUPLUS_Q) + NTRUPLUS_Q) + 1664) % NTRUPLUS_Q) - 1664;
        int16_t r2 = barrett_reduce(i);

        if(r1 != r2) printf("%d %d\n", r1, r2);
    }
}


static void TEST_SPEED_NTT()
{

	unsigned long long cycles;
	unsigned long long cycles1, cycles2;

    int16_t a[NTRUPLUS_N];
    int16_t b[NTRUPLUS_N];

    for (int i = 0; i < NTRUPLUS_N; i++)
    {
        a[i] = i;
    }
    
	printf("========= NTT SPEED TEST =========\n");
	
	cycles=0;
	for (int i = 0; i < TEST_LOOP; i++)
	{
		cycles1 = cpucycles();
        ntt(a);
		cycles2 = cpucycles();
		cycles += cycles2-cycles1;
	}

	printf("  NTT runs in ................. %8lld cycles", cycles/TEST_LOOP);
	printf("\n"); 


	printf("========= MERGE_NTT SPEED TEST =========\n");
	
	cycles=0;
	for (int i = 0; i < TEST_LOOP; i++)
	{
		cycles1 = cpucycles();
        merge6_ntt(a);
		cycles2 = cpucycles();
		cycles += cycles2-cycles1;
	}

	printf("  MERGE_NTT runs in ................. %8lld cycles", cycles/TEST_LOOP);
	printf("\n"); 
}

int main(void)
{
    //TEST_REDUCE();
    TEST_NTT();

    TEST_SPEED_NTT();

    return 0;
}
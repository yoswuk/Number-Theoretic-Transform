#include <stdio.h>
#include <stdint.h>

#define ALPHA 3
#define NTRUPLUS_Q 3457
#define NTRUPLUS_QPRIME 1951806081
#define NTRUPLUS_R 2590
#define NTRUPLUS_R2 1520


int add(int a, int b);

int32_t transform(int16_t a)
{
    return a*NTRUPLUS_R2;
}

int16_t plantard(int16_t a, int16_t b)
{
    int32_t r = a*b*NTRUPLUS_QPRIME;

    r = r >> 16;
    r = r + (1 << (ALPHA));
    r = (r * NTRUPLUS_Q) >> 16;

    return r;
}

int main(void)
{
    int16_t a = 1;
    int16_t b = 1;

    printf("%d\n", add(a,b));

    return 0;
}
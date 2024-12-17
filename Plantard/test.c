#include <stdint.h>

#define ALPHA 3
#define NTRUPLUS_Q 3457
#define NTRUPLUS_QPRIME 1951806081
#define NTRUPLUS_R 2590
#define NTRUPLUS_R2 1520

int16_t plantard(int16_t a, int16_t b)
{
    int32_t r = a*b*NTRUPLUS_QPRIME;

    r = r >> 16;
    r = r + (1 << (ALPHA));
    r = (r * NTRUPLUS_Q) >> 16;

    return r;
}
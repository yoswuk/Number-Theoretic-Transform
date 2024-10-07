#include <stdio.h>
#include <stdint.h>

#define NTRUPLUS1152

#define POW3_0 1
#define POW3_1 3
#define POW3_2 9
#define POW3_3 27

#define Q 3457
#define QINV 12929
#define Rmodq -147

#define CONCAT(a, b) a##b
#define POW3(k) CONCAT(POW3_, k)

#if defined(NTRUPLUS576)
    #define RADIX2 3
    #define RADIX3 2

    #define W 81
    #define W_ORDER 432
#endif

#if defined(NTRUPLUS768)
    #define RADIX2 5
    #define RADIX3 1

    #define W 22
    #define W_ORDER 576
#endif

#if defined(NTRUPLUS864) || defined(NTRUPLUS1152)
    #define RADIX2 4
    #define RADIX3 2

    #define W 9
    #define W_ORDER 864
#endif

#define TABLE_SIZE ((1 << (RADIX2+1))*POW3(RADIX3))

//basic
int16_t exp_table[W_ORDER] = {0};

//tree
int16_t tree[8][W_ORDER] = {0};

//ntt
int16_t ntt_tree[8][W_ORDER] = {0};
int16_t zetas_exp[TABLE_SIZE] = {0};

int center(int a)
{
    return (((((a % Q)+ Q) % Q) + (Q-1)/2) % Q) - (Q-1)/2;
}

void gen_exp()
{
    int a = W;

    exp_table[0] = (1 << 16) % Q;

    for (int i = 1; i < W_ORDER; i++)
    {
        exp_table[i] = center(exp_table[i-1] * a);
    }
}

void gen_tree()
{
    int t = 2;

    tree[0][0] = W_ORDER/6;
    tree[0][1] = 5*W_ORDER/6;

    for (int j = 0; j < RADIX3; j++)
    {
        for (int i = 0; i < t; i++)
        {
            tree[j+1][3*i+0] = tree[j][i]/3;
            tree[j+1][3*i+1] = tree[j+1][3*i+0] + W_ORDER/3;
            tree[j+1][3*i+2] = tree[j+1][3*i+1] + W_ORDER/3;
        }

        t = t*3;
    }

    for (int j = RADIX3; j <= RADIX2+RADIX3; j++)
    {
        for (int i = 0; i < t; i++)
        {
            tree[j+1][2*i+0] = tree[j][i]/2;
            tree[j+1][2*i+1] = tree[j+1][2*i+0] + W_ORDER/2;
        }

        t = t*2;
    }
}

void trans_tree()
{
    int t;

    ntt_tree[0][0] = exp_table[tree[0][0]];
    ntt_tree[0][1] = exp_table[tree[0][1]];

    t = 6;

    for (int j = 1; j < RADIX3+1; j++)
    {
        for (int i = 0; i < t; i++)
        {
            ntt_tree[j][2*i] = exp_table[tree[j][i]];
            ntt_tree[j][2*i+1] = exp_table[(tree[j][i] << 1) % W_ORDER];
        }

        t = t*3;
    }

    for (int j = RADIX3+1; j <= RADIX2+RADIX3; j++)
    {
        for (int i = 0; i < t; i++)
        {
            ntt_tree[j][i] = exp_table[tree[j][i]];
        }

        t = t*2;
    }

}

void ntt_encode()
{
    int k = 0;
    int t = 2;

    zetas_exp[k++] = Rmodq;
    zetas_exp[k++] = ntt_tree[0][0];

    for(int j = 1; j < RADIX3+1; j++)
    {
        for (int i = 0; i < t; i++)
        {                
            zetas_exp[k++] = ntt_tree[j][6*i];
            zetas_exp[k++] = ntt_tree[j][6*i+1];
        }

        t = t*3;
    }
    
    for(int j = RADIX3+1; j <= RADIX2+RADIX3; j++)
    {
        for (int i = 0; i < t; i++)
        {                
            zetas_exp[k++] = ntt_tree[j][2*i];
        }

        t = t*2;
    }

    //print
    printf("const int16_t zetas[%d] = {\n\t", k);
    for (int i = 0; i < k-1; i++)
    {
        printf("%5d, ", zetas_exp[i]);
        if(i%8 == 7) printf("\n\t");
    }
    printf("%5d\n", zetas_exp[k-1]);
    printf("};\n\n");
}

void init()
{
    int t;

    gen_exp();

    gen_tree();

    t = 2;
    for (int j = 0; j < RADIX3; j++)
    {
        printf("level %d\n", j);
        for (int i = 0; i < t; i++)
        {
            printf("%d ", tree[j][i]);
        }
        printf("\n\n");

        t = t*3;
    }

    for (int j = RADIX3; j <= RADIX2+RADIX3; j++)
    {
        printf("level %d\n", j);
        for (int i = 0; i < t; i++)
        {
            printf("%d ", tree[j][i]);
        }
        printf("\n\n");

        t = t*2;
    }

    trans_tree();
}

int main(void)
{
    init();

    ntt_encode();

    printf("w : %d\n", W);
    printf("(z - z^5)^-1 * 2^16 : %d\n", -1665);

#if defined(NTRUPLUS576)    
    printf("(3^RADIX3*2^(RADIX2+1))^-1*(2^32) : %d\n", -66);
    printf("(3^RADIX3*2^(RADIX2+1))^-1*(2^33) : %d\n", -132);
#endif
#if defined(NTRUPLUS768)
    printf("(3^RADIX3*2^(RADIX2+1))^-1*(2^32) : %d\n", 1679);
    printf("(3^RADIX3*2^(RADIX2+1))^-1*(2^33) : %d\n", -99);
#endif
#if defined(NTRUPLUS864) || defined(NTRUPLUS1152)
    printf("(3^RADIX3*2^(RADIX2+1))^-1*(2^32) : %d\n", -33);
    printf("(3^RADIX3*2^(RADIX2+1))^-1*(2^33) : %d\n", -66);
#endif
    return 0;
}
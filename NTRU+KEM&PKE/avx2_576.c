#include <stdio.h>
#include <stdint.h>

#define Q 3457
#define QINV 12929
#define W 81
#define W_ORDER 432

#define RADIX2 3
#define RADIX3 2

//basic
int16_t exp_table[W_ORDER] = {0};

//tree
int16_t tree[RADIX2+RADIX3+1][W_ORDER] = {0};

//ntt
int16_t ntt_tree[RADIX2+RADIX3+1][W_ORDER] = {0};
int16_t zetas_exp[1220] = {0};

//invntt
int16_t invntt_tree[RADIX2+RADIX3+1][W_ORDER] = {0};
int16_t zetas_inv_exp[1224] = {0};

//mul
int16_t zetas_mul_exp[1000] = {0};

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

void trans_tree_inv()
{
    int t;

    invntt_tree[0][0] = exp_table[W_ORDER - tree[0][0]];
    invntt_tree[0][1] = exp_table[W_ORDER - tree[0][1]];

    t = 6;

    for (int j = 1; j < RADIX3+1; j++)
    {
        for (int i = 0; i < t; i++)
        {
            invntt_tree[j][2*i] = exp_table[W_ORDER - tree[j][i]];
            invntt_tree[j][2*i+1] = exp_table[((W_ORDER - tree[j][i])*2) % W_ORDER];
        }

        t = t*3;
    }

    for (int j = RADIX3+1; j <= RADIX2+RADIX3; j++)
    {
        for (int i = 0; i < t; i++)
        {
            invntt_tree[j][i] = exp_table[W_ORDER - tree[j][i]];
        }

        t = t*2;
    }
}

void ntt_encode()
{
    int k = 0;
    int t = 0;

	//level0
    for (int i = 0; i < 1; i++)
    {
        zetas_exp[k++] = ntt_tree[0][0] * QINV;
        zetas_exp[k++] = ntt_tree[0][0] * QINV;
        zetas_exp[k++] = ntt_tree[0][0];
        zetas_exp[k++] = ntt_tree[0][0];
    }
    
	t = 2;

    //level1~2
    for(int j = 1; j < RADIX3 + 1; j++)
    {
        for (int i = 0; i < t; i++)
        {
            zetas_exp[k++] = ntt_tree[j][6*i] * QINV;
            zetas_exp[k++] = ntt_tree[j][6*i] * QINV;

            zetas_exp[k++] = ntt_tree[j][6*i];
            zetas_exp[k++] = ntt_tree[j][6*i];

            zetas_exp[k++] = ntt_tree[j][6*i+1] * QINV;
            zetas_exp[k++] = ntt_tree[j][6*i+1] * QINV;

            zetas_exp[k++] = ntt_tree[j][6*i+1];
            zetas_exp[k++] = ntt_tree[j][6*i+1];       
        }

        t = t*3;
    }

	//level3
    for (int i = 0; i < 9; i++)
    {
        zetas_exp[k++] = ntt_tree[3][(i << 2)    ] * QINV;
        zetas_exp[k++] = ntt_tree[3][(i << 2)    ] * QINV;
        zetas_exp[k++] = ntt_tree[3][(i << 2)    ] * QINV;
        zetas_exp[k++] = ntt_tree[3][(i << 2)    ] * QINV;
        zetas_exp[k++] = ntt_tree[3][(i << 2)    ] * QINV;
        zetas_exp[k++] = ntt_tree[3][(i << 2)    ] * QINV;
        zetas_exp[k++] = ntt_tree[3][(i << 2)    ] * QINV;
        zetas_exp[k++] = ntt_tree[3][(i << 2)    ] * QINV;
        zetas_exp[k++] = ntt_tree[3][(i << 2) + 2] * QINV;
        zetas_exp[k++] = ntt_tree[3][(i << 2) + 2] * QINV;
        zetas_exp[k++] = ntt_tree[3][(i << 2) + 2] * QINV;
        zetas_exp[k++] = ntt_tree[3][(i << 2) + 2] * QINV;
        zetas_exp[k++] = ntt_tree[3][(i << 2) + 2] * QINV;
        zetas_exp[k++] = ntt_tree[3][(i << 2) + 2] * QINV;
        zetas_exp[k++] = ntt_tree[3][(i << 2) + 2] * QINV;
        zetas_exp[k++] = ntt_tree[3][(i << 2) + 2] * QINV;

        zetas_exp[k++] = ntt_tree[3][(i << 2)    ];
        zetas_exp[k++] = ntt_tree[3][(i << 2)    ];
        zetas_exp[k++] = ntt_tree[3][(i << 2)    ];
        zetas_exp[k++] = ntt_tree[3][(i << 2)    ];
        zetas_exp[k++] = ntt_tree[3][(i << 2)    ];
        zetas_exp[k++] = ntt_tree[3][(i << 2)    ];
        zetas_exp[k++] = ntt_tree[3][(i << 2)    ];
        zetas_exp[k++] = ntt_tree[3][(i << 2)    ];
        zetas_exp[k++] = ntt_tree[3][(i << 2) + 2];
        zetas_exp[k++] = ntt_tree[3][(i << 2) + 2];
        zetas_exp[k++] = ntt_tree[3][(i << 2) + 2];
        zetas_exp[k++] = ntt_tree[3][(i << 2) + 2];
        zetas_exp[k++] = ntt_tree[3][(i << 2) + 2];
        zetas_exp[k++] = ntt_tree[3][(i << 2) + 2];
        zetas_exp[k++] = ntt_tree[3][(i << 2) + 2];
        zetas_exp[k++] = ntt_tree[3][(i << 2) + 2];
    }

	//level4
    for (int i = 0; i < 9; i++)
    {
        zetas_exp[k++] = ntt_tree[4][(i << 3)] * QINV;
        zetas_exp[k++] = ntt_tree[4][(i << 3)] * QINV;
        zetas_exp[k++] = ntt_tree[4][(i << 3)] * QINV;
        zetas_exp[k++] = ntt_tree[4][(i << 3)] * QINV;
        zetas_exp[k++] = ntt_tree[4][(i << 3) + 2] * QINV;
        zetas_exp[k++] = ntt_tree[4][(i << 3) + 2] * QINV;
        zetas_exp[k++] = ntt_tree[4][(i << 3) + 2] * QINV;
        zetas_exp[k++] = ntt_tree[4][(i << 3) + 2] * QINV;
        zetas_exp[k++] = ntt_tree[4][(i << 3) + 4] * QINV;
        zetas_exp[k++] = ntt_tree[4][(i << 3) + 4] * QINV;
        zetas_exp[k++] = ntt_tree[4][(i << 3) + 4] * QINV;
        zetas_exp[k++] = ntt_tree[4][(i << 3) + 4] * QINV;
        zetas_exp[k++] = ntt_tree[4][(i << 3) + 6] * QINV;
        zetas_exp[k++] = ntt_tree[4][(i << 3) + 6] * QINV;
        zetas_exp[k++] = ntt_tree[4][(i << 3) + 6] * QINV;
        zetas_exp[k++] = ntt_tree[4][(i << 3) + 6] * QINV;

        zetas_exp[k++] = ntt_tree[4][(i << 3)];
        zetas_exp[k++] = ntt_tree[4][(i << 3)];
        zetas_exp[k++] = ntt_tree[4][(i << 3)];
        zetas_exp[k++] = ntt_tree[4][(i << 3)];
        zetas_exp[k++] = ntt_tree[4][(i << 3) + 2];
        zetas_exp[k++] = ntt_tree[4][(i << 3) + 2];
        zetas_exp[k++] = ntt_tree[4][(i << 3) + 2];
        zetas_exp[k++] = ntt_tree[4][(i << 3) + 2];
        zetas_exp[k++] = ntt_tree[4][(i << 3) + 4];
        zetas_exp[k++] = ntt_tree[4][(i << 3) + 4];
        zetas_exp[k++] = ntt_tree[4][(i << 3) + 4];
        zetas_exp[k++] = ntt_tree[4][(i << 3) + 4];
        zetas_exp[k++] = ntt_tree[4][(i << 3) + 6];
        zetas_exp[k++] = ntt_tree[4][(i << 3) + 6];
        zetas_exp[k++] = ntt_tree[4][(i << 3) + 6];
        zetas_exp[k++] = ntt_tree[4][(i << 3) + 6];
    }
    


	//level5
    for (int i = 0; i < 9; i++)
    {
        zetas_exp[k++] = ntt_tree[5][(i << 4)] * QINV;
        zetas_exp[k++] = ntt_tree[5][(i << 4)] * QINV;
        zetas_exp[k++] = ntt_tree[5][(i << 4) + 2] * QINV;
        zetas_exp[k++] = ntt_tree[5][(i << 4) + 2] * QINV;
        zetas_exp[k++] = ntt_tree[5][(i << 4) + 4] * QINV;
        zetas_exp[k++] = ntt_tree[5][(i << 4) + 4] * QINV;
        zetas_exp[k++] = ntt_tree[5][(i << 4) + 6] * QINV;
        zetas_exp[k++] = ntt_tree[5][(i << 4) + 6] * QINV;
        zetas_exp[k++] = ntt_tree[5][(i << 4) + 8] * QINV;
        zetas_exp[k++] = ntt_tree[5][(i << 4) + 8] * QINV;
        zetas_exp[k++] = ntt_tree[5][(i << 4) + 10] * QINV;
        zetas_exp[k++] = ntt_tree[5][(i << 4) + 10] * QINV;
        zetas_exp[k++] = ntt_tree[5][(i << 4) + 12] * QINV;
        zetas_exp[k++] = ntt_tree[5][(i << 4) + 12] * QINV;
        zetas_exp[k++] = ntt_tree[5][(i << 4) + 14] * QINV;
        zetas_exp[k++] = ntt_tree[5][(i << 4) + 14] * QINV;

        zetas_exp[k++] = ntt_tree[5][(i << 4)];
        zetas_exp[k++] = ntt_tree[5][(i << 4)];
        zetas_exp[k++] = ntt_tree[5][(i << 4) + 2];
        zetas_exp[k++] = ntt_tree[5][(i << 4) + 2];
        zetas_exp[k++] = ntt_tree[5][(i << 4) + 4];
        zetas_exp[k++] = ntt_tree[5][(i << 4) + 4];
        zetas_exp[k++] = ntt_tree[5][(i << 4) + 6];
        zetas_exp[k++] = ntt_tree[5][(i << 4) + 6];
        zetas_exp[k++] = ntt_tree[5][(i << 4) + 8];
        zetas_exp[k++] = ntt_tree[5][(i << 4) + 8];
        zetas_exp[k++] = ntt_tree[5][(i << 4) + 10];
        zetas_exp[k++] = ntt_tree[5][(i << 4) + 10];
        zetas_exp[k++] = ntt_tree[5][(i << 4) + 12];
        zetas_exp[k++] = ntt_tree[5][(i << 4) + 12];
        zetas_exp[k++] = ntt_tree[5][(i << 4) + 14];
        zetas_exp[k++] = ntt_tree[5][(i << 4) + 14];        
    }

    printf("const int16_t zetas[%d] __attribute__((aligned(32))) = {\n\t", k);
    for (int i = 0; i < k-1; i++)
    {
        printf("%7d,", zetas_exp[i]);
        if(i%8 == 7) printf("\n\t");
    }
    printf("%7d\n", zetas_exp[k-1]);
    printf("};\n\n");
}

void invntt_encode()
{
    int t = 0;
    int k = 0;

	//level5
    for (int i = 0; i < 9; i++)
    {
        zetas_inv_exp[k++] = invntt_tree[5][(i << 4)] * QINV;
        zetas_inv_exp[k++] = invntt_tree[5][(i << 4)] * QINV;
        zetas_inv_exp[k++] = invntt_tree[5][(i << 4) + 2] * QINV;
        zetas_inv_exp[k++] = invntt_tree[5][(i << 4) + 2] * QINV;
        zetas_inv_exp[k++] = invntt_tree[5][(i << 4) + 4] * QINV;
        zetas_inv_exp[k++] = invntt_tree[5][(i << 4) + 4] * QINV;
        zetas_inv_exp[k++] = invntt_tree[5][(i << 4) + 6] * QINV;
        zetas_inv_exp[k++] = invntt_tree[5][(i << 4) + 6] * QINV;
        zetas_inv_exp[k++] = invntt_tree[5][(i << 4) + 8] * QINV;
        zetas_inv_exp[k++] = invntt_tree[5][(i << 4) + 8] * QINV;
        zetas_inv_exp[k++] = invntt_tree[5][(i << 4) + 10] * QINV;
        zetas_inv_exp[k++] = invntt_tree[5][(i << 4) + 10] * QINV;
        zetas_inv_exp[k++] = invntt_tree[5][(i << 4) + 12] * QINV;
        zetas_inv_exp[k++] = invntt_tree[5][(i << 4) + 12] * QINV;
        zetas_inv_exp[k++] = invntt_tree[5][(i << 4) + 14] * QINV;
        zetas_inv_exp[k++] = invntt_tree[5][(i << 4) + 14] * QINV;

        zetas_inv_exp[k++] = invntt_tree[5][(i << 4)];
        zetas_inv_exp[k++] = invntt_tree[5][(i << 4)];
        zetas_inv_exp[k++] = invntt_tree[5][(i << 4) + 2];
        zetas_inv_exp[k++] = invntt_tree[5][(i << 4) + 2];
        zetas_inv_exp[k++] = invntt_tree[5][(i << 4) + 4];
        zetas_inv_exp[k++] = invntt_tree[5][(i << 4) + 4];
        zetas_inv_exp[k++] = invntt_tree[5][(i << 4) + 6];
        zetas_inv_exp[k++] = invntt_tree[5][(i << 4) + 6];
        zetas_inv_exp[k++] = invntt_tree[5][(i << 4) + 8];
        zetas_inv_exp[k++] = invntt_tree[5][(i << 4) + 8];
        zetas_inv_exp[k++] = invntt_tree[5][(i << 4) + 10];
        zetas_inv_exp[k++] = invntt_tree[5][(i << 4) + 10];
        zetas_inv_exp[k++] = invntt_tree[5][(i << 4) + 12];
        zetas_inv_exp[k++] = invntt_tree[5][(i << 4) + 12];
        zetas_inv_exp[k++] = invntt_tree[5][(i << 4) + 14];
        zetas_inv_exp[k++] = invntt_tree[5][(i << 4) + 14];        
    }

	//level4
    for (int i = 0; i < 9; i++)
    {
        zetas_inv_exp[k++] = invntt_tree[4][(i << 3)] * QINV;
        zetas_inv_exp[k++] = invntt_tree[4][(i << 3)] * QINV;
        zetas_inv_exp[k++] = invntt_tree[4][(i << 3)] * QINV;
        zetas_inv_exp[k++] = invntt_tree[4][(i << 3)] * QINV;
        zetas_inv_exp[k++] = invntt_tree[4][(i << 3) + 2] * QINV;
        zetas_inv_exp[k++] = invntt_tree[4][(i << 3) + 2] * QINV;
        zetas_inv_exp[k++] = invntt_tree[4][(i << 3) + 2] * QINV;
        zetas_inv_exp[k++] = invntt_tree[4][(i << 3) + 2] * QINV;
        zetas_inv_exp[k++] = invntt_tree[4][(i << 3) + 4] * QINV;
        zetas_inv_exp[k++] = invntt_tree[4][(i << 3) + 4] * QINV;
        zetas_inv_exp[k++] = invntt_tree[4][(i << 3) + 4] * QINV;
        zetas_inv_exp[k++] = invntt_tree[4][(i << 3) + 4] * QINV;
        zetas_inv_exp[k++] = invntt_tree[4][(i << 3) + 6] * QINV;
        zetas_inv_exp[k++] = invntt_tree[4][(i << 3) + 6] * QINV;
        zetas_inv_exp[k++] = invntt_tree[4][(i << 3) + 6] * QINV;
        zetas_inv_exp[k++] = invntt_tree[4][(i << 3) + 6] * QINV;

        zetas_inv_exp[k++] = invntt_tree[4][(i << 3)];
        zetas_inv_exp[k++] = invntt_tree[4][(i << 3)];
        zetas_inv_exp[k++] = invntt_tree[4][(i << 3)];
        zetas_inv_exp[k++] = invntt_tree[4][(i << 3)];
        zetas_inv_exp[k++] = invntt_tree[4][(i << 3) + 2];
        zetas_inv_exp[k++] = invntt_tree[4][(i << 3) + 2];
        zetas_inv_exp[k++] = invntt_tree[4][(i << 3) + 2];
        zetas_inv_exp[k++] = invntt_tree[4][(i << 3) + 2];
        zetas_inv_exp[k++] = invntt_tree[4][(i << 3) + 4];
        zetas_inv_exp[k++] = invntt_tree[4][(i << 3) + 4];
        zetas_inv_exp[k++] = invntt_tree[4][(i << 3) + 4];
        zetas_inv_exp[k++] = invntt_tree[4][(i << 3) + 4];
        zetas_inv_exp[k++] = invntt_tree[4][(i << 3) + 6];
        zetas_inv_exp[k++] = invntt_tree[4][(i << 3) + 6];
        zetas_inv_exp[k++] = invntt_tree[4][(i << 3) + 6];
        zetas_inv_exp[k++] = invntt_tree[4][(i << 3) + 6];
    }
    


	//level3
    for (int i = 0; i < 9; i++)
    {
        zetas_inv_exp[k++] = invntt_tree[3][(i << 2)    ] * QINV;
        zetas_inv_exp[k++] = invntt_tree[3][(i << 2)    ] * QINV;
        zetas_inv_exp[k++] = invntt_tree[3][(i << 2)    ] * QINV;
        zetas_inv_exp[k++] = invntt_tree[3][(i << 2)    ] * QINV;
        zetas_inv_exp[k++] = invntt_tree[3][(i << 2)    ] * QINV;
        zetas_inv_exp[k++] = invntt_tree[3][(i << 2)    ] * QINV;
        zetas_inv_exp[k++] = invntt_tree[3][(i << 2)    ] * QINV;
        zetas_inv_exp[k++] = invntt_tree[3][(i << 2)    ] * QINV;
        zetas_inv_exp[k++] = invntt_tree[3][(i << 2) + 2] * QINV;
        zetas_inv_exp[k++] = invntt_tree[3][(i << 2) + 2] * QINV;
        zetas_inv_exp[k++] = invntt_tree[3][(i << 2) + 2] * QINV;
        zetas_inv_exp[k++] = invntt_tree[3][(i << 2) + 2] * QINV;
        zetas_inv_exp[k++] = invntt_tree[3][(i << 2) + 2] * QINV;
        zetas_inv_exp[k++] = invntt_tree[3][(i << 2) + 2] * QINV;
        zetas_inv_exp[k++] = invntt_tree[3][(i << 2) + 2] * QINV;
        zetas_inv_exp[k++] = invntt_tree[3][(i << 2) + 2] * QINV;

        zetas_inv_exp[k++] = invntt_tree[3][(i << 2)    ];
        zetas_inv_exp[k++] = invntt_tree[3][(i << 2)    ];
        zetas_inv_exp[k++] = invntt_tree[3][(i << 2)    ];
        zetas_inv_exp[k++] = invntt_tree[3][(i << 2)    ];
        zetas_inv_exp[k++] = invntt_tree[3][(i << 2)    ];
        zetas_inv_exp[k++] = invntt_tree[3][(i << 2)    ];
        zetas_inv_exp[k++] = invntt_tree[3][(i << 2)    ];
        zetas_inv_exp[k++] = invntt_tree[3][(i << 2)    ];
        zetas_inv_exp[k++] = invntt_tree[3][(i << 2) + 2];
        zetas_inv_exp[k++] = invntt_tree[3][(i << 2) + 2];
        zetas_inv_exp[k++] = invntt_tree[3][(i << 2) + 2];
        zetas_inv_exp[k++] = invntt_tree[3][(i << 2) + 2];
        zetas_inv_exp[k++] = invntt_tree[3][(i << 2) + 2];
        zetas_inv_exp[k++] = invntt_tree[3][(i << 2) + 2];
        zetas_inv_exp[k++] = invntt_tree[3][(i << 2) + 2];
        zetas_inv_exp[k++] = invntt_tree[3][(i << 2) + 2];
    }

    //level2~1
    t = 6;

    for(int j = RADIX3; j >= 1; j--)
    {
        for (int i = 0; i < t; i++)
        {
            zetas_inv_exp[k++] = invntt_tree[j][6*i] * QINV;
            zetas_inv_exp[k++] = invntt_tree[j][6*i] * QINV;

            zetas_inv_exp[k++] = invntt_tree[j][6*i];
            zetas_inv_exp[k++] = invntt_tree[j][6*i];

            zetas_inv_exp[k++] = invntt_tree[j][6*i+1] * QINV;
            zetas_inv_exp[k++] = invntt_tree[j][6*i+1] * QINV;

            zetas_inv_exp[k++] = invntt_tree[j][6*i+1];
            zetas_inv_exp[k++] = invntt_tree[j][6*i+1];       
        }

        t = t/3;
    }
		
	//level0
    for (int i = 0; i < 1; i++)
    {
        //z = w^(W_ORDER/6)
        //(z-z^5)^-1*(2^16)
        zetas_inv_exp[k++] = (int16_t) (-1665 * QINV);
        zetas_inv_exp[k++] = (int16_t) (-1665 * QINV);
        zetas_inv_exp[k++] = -1665;
        zetas_inv_exp[k++] = -1665;

        //(3^RADIX3*2^(RADIX2+1))^-1*(2^32)
        zetas_inv_exp[k++] = (int16_t)(-66 * QINV);
        zetas_inv_exp[k++] = (int16_t)(-66 * QINV);
        zetas_inv_exp[k++] = -66;
        zetas_inv_exp[k++] = -66;
    }

    printf("const int16_t zetas_inv[%d] __attribute__((aligned(32))) = {\n\t", k);
    for (int i = 0; i < k-1; i++)
    {   
        printf("%7d,", zetas_inv_exp[i]);   
        if(i%8==7) printf("\n\t");
    }
    printf("%7d\n", zetas_inv_exp[k-1]);
    printf("};\n\n");
}

void mul_encode()
{
    int k =0;

    for (int i = 0; i < 9; i++)
    {
        zetas_mul_exp[k++] = (ntt_tree[5][(i << 4)     ] * QINV);
        zetas_mul_exp[k++] = (ntt_tree[5][(i << 4) +  1] * QINV);
        zetas_mul_exp[k++] = (ntt_tree[5][(i << 4) +  2] * QINV);
        zetas_mul_exp[k++] = (ntt_tree[5][(i << 4) +  3] * QINV);
        zetas_mul_exp[k++] = (ntt_tree[5][(i << 4) +  4] * QINV);
        zetas_mul_exp[k++] = (ntt_tree[5][(i << 4) +  5] * QINV);
        zetas_mul_exp[k++] = (ntt_tree[5][(i << 4) +  6] * QINV);
        zetas_mul_exp[k++] = (ntt_tree[5][(i << 4) +  7] * QINV);
        zetas_mul_exp[k++] = (ntt_tree[5][(i << 4) +  8] * QINV);
        zetas_mul_exp[k++] = (ntt_tree[5][(i << 4) +  9] * QINV);
        zetas_mul_exp[k++] = (ntt_tree[5][(i << 4) + 10] * QINV);
        zetas_mul_exp[k++] = (ntt_tree[5][(i << 4) + 11] * QINV);
        zetas_mul_exp[k++] = (ntt_tree[5][(i << 4) + 12] * QINV);
        zetas_mul_exp[k++] = (ntt_tree[5][(i << 4) + 13] * QINV);
        zetas_mul_exp[k++] = (ntt_tree[5][(i << 4) + 14] * QINV);
        zetas_mul_exp[k++] = (ntt_tree[5][(i << 4) + 15] * QINV);

        zetas_mul_exp[k++] = (ntt_tree[5][(i << 4)     ]);
        zetas_mul_exp[k++] = (ntt_tree[5][(i << 4) +  1]);
        zetas_mul_exp[k++] = (ntt_tree[5][(i << 4) +  2]);
        zetas_mul_exp[k++] = (ntt_tree[5][(i << 4) +  3]);
        zetas_mul_exp[k++] = (ntt_tree[5][(i << 4) +  4]);
        zetas_mul_exp[k++] = (ntt_tree[5][(i << 4) +  5]);
        zetas_mul_exp[k++] = (ntt_tree[5][(i << 4) +  6]);
        zetas_mul_exp[k++] = (ntt_tree[5][(i << 4) +  7]);
        zetas_mul_exp[k++] = (ntt_tree[5][(i << 4) +  8]);
        zetas_mul_exp[k++] = (ntt_tree[5][(i << 4) +  9]);
        zetas_mul_exp[k++] = (ntt_tree[5][(i << 4) + 10]);
        zetas_mul_exp[k++] = (ntt_tree[5][(i << 4) + 11]);
        zetas_mul_exp[k++] = (ntt_tree[5][(i << 4) + 12]);
        zetas_mul_exp[k++] = (ntt_tree[5][(i << 4) + 13]);
        zetas_mul_exp[k++] = (ntt_tree[5][(i << 4) + 14]);
        zetas_mul_exp[k++] = (ntt_tree[5][(i << 4) + 15]);
    }

    printf("const int16_t zetas_mul[%d] __attribute__((aligned(32))) = {\n\t", k);
    for (int i = 0; i < k-1; i++)
    {
        printf("%7d,", zetas_mul_exp[i]);
        if(i%8 == 7) printf("\n\t");
    }
    printf("%7d\n", zetas_mul_exp[k-1]);
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
    trans_tree_inv();
}

int main(void)
{
    init();

    ntt_encode();
    invntt_encode();
    mul_encode();

    return 0;
}

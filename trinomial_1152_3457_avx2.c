#include <stdio.h>
#include <stdint.h>

#define Q 3457
#define QINV 12929
#define W 9
#define W_ORDER 864

//basic
int16_t exp_table[W_ORDER] = {0};

//tree
int16_t tree[8][288] = {0};
int16_t ntt_tree[8][288] = {0};
int16_t invntt_tree[8][288] = {0};

//invntt
int16_t zetas_exp[860] = {0};
int16_t zetas_inv_exp[864] = {0};

int center(int a)
{
    return (((((a % Q)+ Q) % Q) + (Q-1)/2) % Q) - (Q-1)/2;
}

void gen_exp()
{
    int a = W;

    exp_table[0] = (1 << 16) % Q;
    //exp_table[0] = 1;
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

    for (int j = 0; j < 2; j++)
    {
        for (int i = 0; i < t; i++)
        {
            tree[j+1][3*i+0] = tree[j][i]/3;
            tree[j+1][3*i+1] = tree[j+1][3*i+0] + W_ORDER/3;
            tree[j+1][3*i+2] = tree[j+1][3*i+1] + W_ORDER/3;
        }

        t = t*3;
    }

    for (int j = 2; j < 6; j++)
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

    for (int j = 1; j < 3; j++)
    {
        for (int i = 0; i < t; i++)
        {
            ntt_tree[j][2*i] = exp_table[tree[j][i]];
            ntt_tree[j][2*i+1] = exp_table[(tree[j][i] << 1) % W_ORDER];
        }

        t = t*3;
    }

    for (int j = 3; j < 7; j++)
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

    for (int j = 1; j < 3; j++)
    {
        for (int i = 0; i < t; i++)
        {
            invntt_tree[j][2*i] = exp_table[W_ORDER - tree[j][i]];
            invntt_tree[j][2*i+1] = exp_table[((W_ORDER - tree[j][i])*2) % W_ORDER];
        }

        t = t*3;
    }

    for (int j = 3; j < 7; j++)
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

    printf("ntt_encode\n");
//level0
    for (int i = 0; i < 1; i++)
    {
        zetas_exp[k++] = ntt_tree[0][0] * QINV;
        zetas_exp[k++] = ntt_tree[0][0] * QINV;
        zetas_exp[k++] = ntt_tree[0][0];
        zetas_exp[k++] = ntt_tree[0][0];
    }
    printf("k : %d\n", k);    
//level1
    for (int i = 0; i < 2; i++)
    {
        zetas_exp[k++] = ntt_tree[1][6*i] * QINV;
        zetas_exp[k++] = ntt_tree[1][6*i] * QINV;

        zetas_exp[k++] = ntt_tree[1][6*i];
        zetas_exp[k++] = ntt_tree[1][6*i];

        zetas_exp[k++] = ntt_tree[1][6*i+1] * QINV;
        zetas_exp[k++] = ntt_tree[1][6*i+1] * QINV;

        zetas_exp[k++] = ntt_tree[1][6*i+1];
        zetas_exp[k++] = ntt_tree[1][6*i+1];        
    }
    printf("k : %d\n", k);
//level2
    for (int i = 0; i < 6; i++)
    {
        zetas_exp[k++] = ntt_tree[2][6*i] * QINV;
        zetas_exp[k++] = ntt_tree[2][6*i] * QINV;

        zetas_exp[k++] = ntt_tree[2][6*i];
        zetas_exp[k++] = ntt_tree[2][6*i];

        zetas_exp[k++] = ntt_tree[2][6*i+1] * QINV;
        zetas_exp[k++] = ntt_tree[2][6*i+1] * QINV;

        zetas_exp[k++] = ntt_tree[2][6*i+1];
        zetas_exp[k++] = ntt_tree[2][6*i+1];        
    }
    printf("k : %d\n", k);

//level3
    for (int i = 0; i < 18; i++)
    {
        zetas_exp[k++] = ntt_tree[3][2*i] * QINV;
        zetas_exp[k++] = ntt_tree[3][2*i] * QINV;

        zetas_exp[k++] = ntt_tree[3][2*i];
        zetas_exp[k++] = ntt_tree[3][2*i];
    }
    printf("k : %d\n", k);

//level4
    for (int i = 0; i < 36; i++)
    {
        zetas_exp[k++] = ntt_tree[4][2*i] * QINV;
        zetas_exp[k++] = ntt_tree[4][2*i] * QINV;

        zetas_exp[k++] = ntt_tree[4][2*i];
        zetas_exp[k++] = ntt_tree[4][2*i];
    }
    printf("k : %d\n", k);

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

    printf("k : %d\n", k);

//level6
    for (int i = 0; i < 9; i++)
    {
        zetas_exp[k++] = ntt_tree[6][(i << 4)     ] * QINV;
        zetas_exp[k++] = ntt_tree[6][(i << 4) +  2] * QINV;
        zetas_exp[k++] = ntt_tree[6][(i << 4) +  4] * QINV;
        zetas_exp[k++] = ntt_tree[6][(i << 4) +  6] * QINV;
        zetas_exp[k++] = ntt_tree[6][(i << 4) +  8] * QINV;
        zetas_exp[k++] = ntt_tree[6][(i << 4) + 10] * QINV;
        zetas_exp[k++] = ntt_tree[6][(i << 4) + 12] * QINV;
        zetas_exp[k++] = ntt_tree[6][(i << 4) + 14] * QINV;
        zetas_exp[k++] = ntt_tree[6][(i << 4) + 16] * QINV;
        zetas_exp[k++] = ntt_tree[6][(i << 4) + 18] * QINV;
        zetas_exp[k++] = ntt_tree[6][(i << 4) + 20] * QINV;
        zetas_exp[k++] = ntt_tree[6][(i << 4) + 22] * QINV;
        zetas_exp[k++] = ntt_tree[6][(i << 4) + 24] * QINV;
        zetas_exp[k++] = ntt_tree[6][(i << 4) + 26] * QINV;
        zetas_exp[k++] = ntt_tree[6][(i << 4) + 28] * QINV;
        zetas_exp[k++] = ntt_tree[6][(i << 4) + 30] * QINV;

        zetas_exp[k++] = ntt_tree[6][(i << 4)     ];
        zetas_exp[k++] = ntt_tree[6][(i << 4) +  2];
        zetas_exp[k++] = ntt_tree[6][(i << 4) +  4];
        zetas_exp[k++] = ntt_tree[6][(i << 4) +  6];
        zetas_exp[k++] = ntt_tree[6][(i << 4) +  8];
        zetas_exp[k++] = ntt_tree[6][(i << 4) + 10];
        zetas_exp[k++] = ntt_tree[6][(i << 4) + 12];
        zetas_exp[k++] = ntt_tree[6][(i << 4) + 14];
        zetas_exp[k++] = ntt_tree[6][(i << 4) + 16];
        zetas_exp[k++] = ntt_tree[6][(i << 4) + 18];
        zetas_exp[k++] = ntt_tree[6][(i << 4) + 20];
        zetas_exp[k++] = ntt_tree[6][(i << 4) + 22];
        zetas_exp[k++] = ntt_tree[6][(i << 4) + 24];
        zetas_exp[k++] = ntt_tree[6][(i << 4) + 26];
        zetas_exp[k++] = ntt_tree[6][(i << 4) + 28];
        zetas_exp[k++] = ntt_tree[6][(i << 4) + 30];       
    }

    printf("k : %d\n\n", k);
}

void invntt_encode()
{
    int k = 0;

    printf("invntt_encode\n");
//level6
    for (int i = 0; i < 9; i++)
    {
        zetas_inv_exp[k++] = invntt_tree[6][(i << 4)     ] * QINV;
        zetas_inv_exp[k++] = invntt_tree[6][(i << 4) +  2] * QINV;
        zetas_inv_exp[k++] = invntt_tree[6][(i << 4) +  4] * QINV;
        zetas_inv_exp[k++] = invntt_tree[6][(i << 4) +  6] * QINV;
        zetas_inv_exp[k++] = invntt_tree[6][(i << 4) +  8] * QINV;
        zetas_inv_exp[k++] = invntt_tree[6][(i << 4) + 10] * QINV;
        zetas_inv_exp[k++] = invntt_tree[6][(i << 4) + 12] * QINV;
        zetas_inv_exp[k++] = invntt_tree[6][(i << 4) + 14] * QINV;
        zetas_inv_exp[k++] = invntt_tree[6][(i << 4) + 16] * QINV;
        zetas_inv_exp[k++] = invntt_tree[6][(i << 4) + 18] * QINV;
        zetas_inv_exp[k++] = invntt_tree[6][(i << 4) + 20] * QINV;
        zetas_inv_exp[k++] = invntt_tree[6][(i << 4) + 22] * QINV;
        zetas_inv_exp[k++] = invntt_tree[6][(i << 4) + 24] * QINV;
        zetas_inv_exp[k++] = invntt_tree[6][(i << 4) + 26] * QINV;
        zetas_inv_exp[k++] = invntt_tree[6][(i << 4) + 28] * QINV;
        zetas_inv_exp[k++] = invntt_tree[6][(i << 4) + 30] * QINV;

        zetas_inv_exp[k++] = invntt_tree[6][(i << 4)     ];
        zetas_inv_exp[k++] = invntt_tree[6][(i << 4) +  2];
        zetas_inv_exp[k++] = invntt_tree[6][(i << 4) +  4];
        zetas_inv_exp[k++] = invntt_tree[6][(i << 4) +  6];
        zetas_inv_exp[k++] = invntt_tree[6][(i << 4) +  8];
        zetas_inv_exp[k++] = invntt_tree[6][(i << 4) + 10];
        zetas_inv_exp[k++] = invntt_tree[6][(i << 4) + 12];
        zetas_inv_exp[k++] = invntt_tree[6][(i << 4) + 14];
        zetas_inv_exp[k++] = invntt_tree[6][(i << 4) + 16];
        zetas_inv_exp[k++] = invntt_tree[6][(i << 4) + 18];
        zetas_inv_exp[k++] = invntt_tree[6][(i << 4) + 20];
        zetas_inv_exp[k++] = invntt_tree[6][(i << 4) + 22];
        zetas_inv_exp[k++] = invntt_tree[6][(i << 4) + 24];
        zetas_inv_exp[k++] = invntt_tree[6][(i << 4) + 26];
        zetas_inv_exp[k++] = invntt_tree[6][(i << 4) + 28];
        zetas_inv_exp[k++] = invntt_tree[6][(i << 4) + 30];       
    }

    printf("k : %d\n", k);

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

    printf("k : %d\n", k);

//level4
    for (int i = 0; i < 36; i++)
    {
        zetas_inv_exp[k++] = invntt_tree[4][2*i] * QINV;
        zetas_inv_exp[k++] = invntt_tree[4][2*i] * QINV;

        zetas_inv_exp[k++] = invntt_tree[4][2*i];
        zetas_inv_exp[k++] = invntt_tree[4][2*i];
    }
    printf("k : %d\n", k);

//level3
    for (int i = 0; i < 18; i++)
    {
        zetas_inv_exp[k++] = invntt_tree[3][2*i] * QINV;
        zetas_inv_exp[k++] = invntt_tree[3][2*i] * QINV;

        zetas_inv_exp[k++] = invntt_tree[3][2*i];
        zetas_inv_exp[k++] = invntt_tree[3][2*i];
    }
    printf("k : %d\n", k);

//level2
    for (int i = 0; i < 6; i++)
    {
        zetas_inv_exp[k++] = invntt_tree[2][6*i] * QINV;
        zetas_inv_exp[k++] = invntt_tree[2][6*i] * QINV;

        zetas_inv_exp[k++] = invntt_tree[2][6*i];
        zetas_inv_exp[k++] = invntt_tree[2][6*i];

        zetas_inv_exp[k++] = invntt_tree[2][6*i+1] * QINV;
        zetas_inv_exp[k++] = invntt_tree[2][6*i+1] * QINV;

        zetas_inv_exp[k++] = invntt_tree[2][6*i+1];
        zetas_inv_exp[k++] = invntt_tree[2][6*i+1];

    }
    printf("k : %d\n", k);

//level1
    for (int i = 0; i < 2; i++)
    {
        zetas_inv_exp[k++] = invntt_tree[1][6*i] * QINV;
        zetas_inv_exp[k++] = invntt_tree[1][6*i] * QINV;

        zetas_inv_exp[k++] = invntt_tree[1][6*i];
        zetas_inv_exp[k++] = invntt_tree[1][6*i];

        zetas_inv_exp[k++] = invntt_tree[1][6*i+1] * QINV;
        zetas_inv_exp[k++] = invntt_tree[1][6*i+1] * QINV;

        zetas_inv_exp[k++] = invntt_tree[1][6*i+1];
        zetas_inv_exp[k++] = invntt_tree[1][6*i+1];
    }
    printf("k : %d\n", k);

//level0
    for (int i = 0; i < 1; i++)
    {
        //(z-z^5)^-1
        zetas_inv_exp[k++] = (int16_t) (-1665 * QINV);
        zetas_inv_exp[k++] = (int16_t) (-1665 * QINV);
        zetas_inv_exp[k++] = -1665;
        zetas_inv_exp[k++] = -1665;

        //기억이 안남.... 전체적인 기법이랑, 3*2^x의 역수 포함
        zetas_inv_exp[k++] = (int16_t)(-33 * QINV);
        zetas_inv_exp[k++] = (int16_t)(-33 * QINV);
        zetas_inv_exp[k++] = -33;
        zetas_inv_exp[k++] = -33; //&3391
    }
    printf("k : %d\n", k);
}

void init()
{
    int t;

    gen_exp();

    gen_tree();

    printf("tree\n");
    t = 2;
    for (int j = 0; j < 2; j++)
    {
        printf("level %d\n", j);
        for (int i = 0; i < t; i++)
        {
            printf("%d ", tree[j][i]);
            if(i % 12 == 11) printf("\n");
        }
        printf("\n\n");

        t = t*3;
    }

    for (int j = 2; j < 7; j++)
    {
        printf("level %d\n", j);
        for (int i = 0; i < t; i++)
        {
            printf("%3d ", tree[j][i]);
            if(i % 6 == 5) printf("\n");
        }
        printf("\n\n");

        t = t*2;
    }
}

void ntt()
{
    trans_tree();
    ntt_encode();
    printf("int16_t zetas[860] = {\n");
    for (int i = 0; i < 859; i++)
    {
        if(i%8==0) printf("\t");
        printf("%6d, ", zetas_exp[i]);
        if(i%8 == 7) printf("\n");
    }
    printf("%6d", zetas_exp[859]);
    printf("\n};");
    printf("\n\n");
}

void invntt()
{
    trans_tree_inv();
    invntt_encode();

    printf("int16_t zetas_inv[864] = {\n\t");
    for (int i = 0; i < 863; i++)
    {   
        printf("%6d, ", zetas_inv_exp[i]);   
        if(i%8==7) printf("\n\t");
    }
    printf("%6d", zetas_inv_exp[863]);
    printf("\n};");
    printf("\n\n");
}



int main(void)
{
    init();

    ntt();
    invntt();

    return 0;
}

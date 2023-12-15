#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <openacc.h>  
#include "string.h"


/*
生成A和B矩阵，Omega矩阵，2*Omega-A-B矩阵，A-B矩阵,q1向量,q2向量，xm向量
a1_1, a2_1, a3_1, a4_1, a5_1表示A矩阵从左到右元素，b_1是相对于a3_1的偏移量
a1_2, a2_2, a3_2, a4_2, a5_2表示B矩阵从左到右元素，b_2是相对于a3_2的偏移量
m表示阶数
*/
void createSparseMatrix(double a1_1, double a2_1, double a3_1, double a4_1, double a5_1, double b_1,
    double a1_2, double a2_2, double a3_2, double a4_2, double a5_2, double b_2,
    int* A_rows_1, int* A_cols_1, double* A_vals_1, int* A_dataLen_1, double* q_1,
    int* A_rows_2, int* A_cols_2, double* A_vals_2, int* A_dataLen_2, double* q_2,
    int* Omega_rows, int* Omega_cols, double* Omega_vals,
    int* A1SUBA2_rows, int* A1SUBA2_cols, double* A1SUBA2_vals, int* A1SUBA2_dataLen,
    int* OmegaA1A2_rows, int* OmegaA1A2_cols, double* OmegaA1A2_vals, int* OmegaA1A2_dataLen,
    double alpha_1, double alpha_2, double* xm, int m) {

    int n = m * m;
    double* z = (double*)malloc(sizeof(double) * (n + 4));
    double* w1 = (double*)malloc(sizeof(double) * (n + 3));
    double* w2 = (double*)malloc(sizeof(double) * (n + 3));

    memset(z, 0, sizeof(double) * (n + 4));
    memset(w1, 0, sizeof(double) * (n));
    memset(w2, 0, sizeof(double) * (n));
    int i;
    for (i = 0; i < n / 3; i++)
    {
        if (3 * i >= n)
        {
            break;
        }
        z[i * 3] = 1.0;
        w1[i * 3 + 1] = 1.0;
        w2[i * 3 + 2] = 1.0;
    }
    /*
    for (int i = 0; i < n; i++)
    {
        printf("%d, %lf\n",i, z[i]);
    }
    */

    int tag_1 = 0;
    int lenPerRow_1 = 0;
    A_dataLen_1[n] = 0;

    int tag_2 = 0;
    int lenPerRow_2 = 0;
    A_dataLen_2[n] = 0;

    int tag = 0;
    int lenPerRow = 0;
    A1SUBA2_dataLen[n] = 0;

    int tag_3 = 0;
    int lenPerRow_3 = 0;
    OmegaA1A2_dataLen[n] = 0;

    double q_tmp = 0;

    int rows_tag_tmp = 0;
    int cols_tag_tmp = 0;
    double A1SUBA2_tmp = 0;
    double OmegaA1A2_tmp = 0;

    for (i = 0; i < n; i++)
    {
        q_tmp = 0;
        if ((i - m) >= 0)
        {
            rows_tag_tmp = 0;
            cols_tag_tmp = 0;
            A1SUBA2_tmp = 0;
            OmegaA1A2_tmp = 0;
            //if (fabs(a1_1) > 1.0e-16)
           // {
            A_rows_1[tag_1] = i;
            A_cols_1[tag_1] = i - m;
            A_vals_1[tag_1] = a1_1;
            q_tmp = q_tmp + A_vals_1[tag_1] * z[i - m];
            rows_tag_tmp = i;
            cols_tag_tmp = i - m;
            A1SUBA2_tmp = A_vals_1[tag_1];
            OmegaA1A2_tmp = -A_vals_1[tag_1];
            tag_1++;
            lenPerRow_1++;
            //    }

                //if (fabs(a1_2) > 1.0e-16)
            //    {
            A_rows_2[tag_2] = i;
            A_cols_2[tag_2] = i - m;
            A_vals_2[tag_2] = a1_2;
            rows_tag_tmp = i;
            cols_tag_tmp = i - m;
            A1SUBA2_tmp = A1SUBA2_tmp - A_vals_2[tag_2];
            OmegaA1A2_tmp = OmegaA1A2_tmp - A_vals_2[tag_2];
            tag_2++;
            lenPerRow_2++;
            //    }

                //if (fabs(A1SUBA2_tmp) > 1.0e-16)
           //     {
            A1SUBA2_rows[tag] = rows_tag_tmp;
            A1SUBA2_cols[tag] = cols_tag_tmp;
            A1SUBA2_vals[tag] = A1SUBA2_tmp;
            tag++;
            lenPerRow++;
            //    }

                //if (fabs(OmegaA1A2_tmp) > 1.0e-16)
             //   {
            OmegaA1A2_rows[tag_3] = rows_tag_tmp;
            OmegaA1A2_cols[tag_3] = cols_tag_tmp;
            OmegaA1A2_vals[tag_3] = OmegaA1A2_tmp;
            tag_3++;
            lenPerRow_3++;
            //     }
        }
        if ((i % m) > 0)
        {

            rows_tag_tmp = 0;
            cols_tag_tmp = 0;
            A1SUBA2_tmp = 0;
            OmegaA1A2_tmp = 0;
            //if (fabs(a2_1) > 1.0e-16)
         //   {
            A_rows_1[tag_1] = i;
            A_cols_1[tag_1] = i - 1;
            A_vals_1[tag_1] = a2_1;
            q_tmp = q_tmp + A_vals_1[tag_1] * z[i - 1];
            rows_tag_tmp = i;
            cols_tag_tmp = i - 1;
            A1SUBA2_tmp = A_vals_1[tag_1];
            OmegaA1A2_tmp = -A_vals_1[tag_1];
            tag_1++;
            lenPerRow_1++;
            //  }

              //if (fabs(a2_2) > 1.0e-16)
           //   {
            A_rows_2[tag_2] = i;
            A_cols_2[tag_2] = i - 1;
            A_vals_2[tag_2] = a2_2;
            rows_tag_tmp = i;
            cols_tag_tmp = i - 1;
            A1SUBA2_tmp = A1SUBA2_tmp - A_vals_2[tag_2];
            OmegaA1A2_tmp = OmegaA1A2_tmp - A_vals_2[tag_2];
            tag_2++;
            lenPerRow_2++;
            //    }

                //if (fabs(A1SUBA2_tmp) > 1.0e-16)
             //   {
            A1SUBA2_rows[tag] = rows_tag_tmp;
            A1SUBA2_cols[tag] = cols_tag_tmp;
            A1SUBA2_vals[tag] = A1SUBA2_tmp;
            tag++;
            lenPerRow++;
            //   }

               //if (fabs(OmegaA1A2_tmp) > 1.0e-16)
            //   {
            OmegaA1A2_rows[tag_3] = rows_tag_tmp;
            OmegaA1A2_cols[tag_3] = cols_tag_tmp;
            OmegaA1A2_vals[tag_3] = OmegaA1A2_tmp;
            tag_3++;
            lenPerRow_3++;
            //    }
        }

        rows_tag_tmp = 0;
        cols_tag_tmp = 0;
        A1SUBA2_tmp = 0;
        OmegaA1A2_tmp = 0;

        //if (fabs(a3_1 + b_1) > 1.0e-16)
       // {
        A_rows_1[tag_1] = i;
        A_cols_1[tag_1] = i;
        A_vals_1[tag_1] = a3_1 + b_1;
        q_tmp = q_tmp + A_vals_1[tag_1] * z[i];
        rows_tag_tmp = i;
        cols_tag_tmp = i;
        A1SUBA2_tmp = A_vals_1[tag_1];
        OmegaA1A2_tmp = -A_vals_1[tag_1];
        tag_1++;
        lenPerRow_1++;
        //  }

          //if (fabs(a3_2 + b_2) > 1.0e-16)
         // {
        A_rows_2[tag_2] = i;
        A_cols_2[tag_2] = i;
        A_vals_2[tag_2] = a3_2 + b_2;
        rows_tag_tmp = i;
        cols_tag_tmp = i;
        A1SUBA2_tmp = A1SUBA2_tmp - A_vals_2[tag_2];
        OmegaA1A2_tmp = OmegaA1A2_tmp - A_vals_2[tag_2];
        tag_2++;
        lenPerRow_2++;
        //  }

        Omega_rows[i] = i;
        Omega_cols[i] = i;
        Omega_vals[i] = ((a3_1 + b_1) / alpha_1 + (a3_2 + b_2) / alpha_2) / 2;

        rows_tag_tmp = i;
        cols_tag_tmp = i;
        OmegaA1A2_tmp = OmegaA1A2_tmp + 2 * Omega_vals[i];

        //if (fabs(A1SUBA2_tmp) > 1.0e-16)
    //    {
        A1SUBA2_rows[tag] = rows_tag_tmp;
        A1SUBA2_cols[tag] = cols_tag_tmp;
        A1SUBA2_vals[tag] = A1SUBA2_tmp;
        tag++;
        lenPerRow++;
        //    }

            //if (fabs(OmegaA1A2_tmp) > 1.0e-16)
         //   {
        OmegaA1A2_rows[tag_3] = rows_tag_tmp;
        OmegaA1A2_cols[tag_3] = cols_tag_tmp;
        OmegaA1A2_vals[tag_3] = OmegaA1A2_tmp;
        tag_3++;
        lenPerRow_3++;
        //     }

        if ((i % m) < m - 1)
        {
            rows_tag_tmp = 0;
            cols_tag_tmp = 0;
            A1SUBA2_tmp = 0;
            OmegaA1A2_tmp = 0;

            //if (fabs(a4_1) > 1.0e-16)
       //     {
            A_rows_1[tag_1] = i;
            A_cols_1[tag_1] = i + 1;
            A_vals_1[tag_1] = a4_1;
            q_tmp = q_tmp + A_vals_1[tag_1] * z[i + 1];
            rows_tag_tmp = i;
            cols_tag_tmp = i + 1;
            A1SUBA2_tmp = A_vals_1[tag_1];
            OmegaA1A2_tmp = -A_vals_1[tag_1];
            tag_1++;
            lenPerRow_1++;
            //       }

                   //if (fabs(a4_2) > 1.0e-16)
               //    {
            A_rows_2[tag_2] = i;
            A_cols_2[tag_2] = i + 1;
            A_vals_2[tag_2] = a4_2;
            rows_tag_tmp = i;
            cols_tag_tmp = i + 1;
            A1SUBA2_tmp = A1SUBA2_tmp - A_vals_2[tag_2];
            OmegaA1A2_tmp = OmegaA1A2_tmp - A_vals_2[tag_2];
            tag_2++;
            lenPerRow_2++;
            //    }

                //if (fabs(A1SUBA2_tmp) > 1.0e-16)
              //  {
            A1SUBA2_rows[tag] = rows_tag_tmp;
            A1SUBA2_cols[tag] = cols_tag_tmp;
            A1SUBA2_vals[tag] = A1SUBA2_tmp;
            tag++;
            lenPerRow++;
            //   }

               //if (fabs(OmegaA1A2_tmp) > 1.0e-16)
            //   {
            OmegaA1A2_rows[tag_3] = rows_tag_tmp;
            OmegaA1A2_cols[tag_3] = cols_tag_tmp;
            OmegaA1A2_vals[tag_3] = OmegaA1A2_tmp;
            tag_3++;
            lenPerRow_3++;
            //   }
        }

        if ((i + m) < n)
        {
            rows_tag_tmp = 0;
            cols_tag_tmp = 0;
            A1SUBA2_tmp = 0;
            OmegaA1A2_tmp = 0;

            //if (fabs(a5_1) > 1.0e-16)
         //   {
            A_rows_1[tag_1] = i;
            A_cols_1[tag_1] = i + m;
            A_vals_1[tag_1] = a5_1;
            q_tmp = q_tmp + A_vals_1[tag_1] * z[i + m];
            rows_tag_tmp = i;
            cols_tag_tmp = i + m;
            A1SUBA2_tmp = A_vals_1[tag_1];
            OmegaA1A2_tmp = -A_vals_1[tag_1];
            tag_1++;
            lenPerRow_1++;
            //  }

              //if (fabs(a5_2) > 1.0e-16)
           //   {
            A_rows_2[tag_2] = i;
            A_cols_2[tag_2] = i + m;
            A_vals_2[tag_2] = a5_2;
            rows_tag_tmp = i;
            cols_tag_tmp = i + m;
            A1SUBA2_tmp = A1SUBA2_tmp - A_vals_2[tag_2];
            OmegaA1A2_tmp = OmegaA1A2_tmp - A_vals_2[tag_2];
            tag_2++;
            lenPerRow_2++;
            //   }

               //if (fabs(A1SUBA2_tmp) > 1.0e-16)
           //    {
            A1SUBA2_rows[tag] = rows_tag_tmp;
            A1SUBA2_cols[tag] = cols_tag_tmp;
            A1SUBA2_vals[tag] = A1SUBA2_tmp;
            tag++;
            lenPerRow++;
            //   }

               //if (fabs(OmegaA1A2_tmp) > 1.0e-16)
            //   {
            OmegaA1A2_rows[tag_3] = rows_tag_tmp;
            OmegaA1A2_cols[tag_3] = cols_tag_tmp;
            OmegaA1A2_vals[tag_3] = OmegaA1A2_tmp;
            tag_3++;
            lenPerRow_3++;
            //  }
        }

        A_dataLen_1[i] = lenPerRow_1;
        lenPerRow_1 = 0;

        A_dataLen_2[i] = lenPerRow_2;
        lenPerRow_2 = 0;

        A1SUBA2_dataLen[i] = lenPerRow;
        lenPerRow = 0;

        OmegaA1A2_dataLen[i] = lenPerRow_3;
        lenPerRow_3 = 0;

        A_dataLen_1[n] += A_dataLen_1[i];

        A_dataLen_2[n] += A_dataLen_2[i];

        A1SUBA2_dataLen[n] += A1SUBA2_dataLen[i];

        OmegaA1A2_dataLen[n] += OmegaA1A2_dataLen[i];

        q_1[i] = w1[i] - q_tmp;

        q_2[i] = w2[i] - q_tmp;

        //printf("%d\t%lf\t %lf\n" , i , q_1[i],q_2[i]);

        // xm[i] = -0.5;
        // if (i % 2 == 0)
        // {
        //     xm[i] = -1 * xm[i];
        // }
        xm[i] = 1.0;

    }

    free(z);
    free(w1);
    free(w2);
}

/*工具方法*/
int calSSum(int* s, int len)
{
    int sum = 0;
    int i;
    for ( i = 0; i < len; i++)
    {
        sum = sum + s[i];
    }
    return sum;
}



/*
该函数是返回算法迭代的步数

参数说明：
runAlgorithmTag 表示算法执行与否，2表示执行算法3 1表示执行算法2，0表示不执行算法
a1_1, a2_1, a3_1, a4_1, a5_1表示A矩阵从左到右元素，b_1是相对于a3_1的偏移量
a1_2, a2_2, a3_2, a4_2, a5_2表示B矩阵从左到右元素，b_2是相对于a3_2的偏移量

alpha_1, beta_1,gamma_1, alpha_2, beta_2, gamma_2,表示论文中提到的系数

res是一个数组，记录每次迭代的残差值
maxit表示最大迭代步数
numworkers表示开启的线程数
m表示阶数
tol表示收敛阈值
*/
int MSM(int runAlgorithmTag, double a1_1, double a2_1, double a3_1, double a4_1, double a5_1, double b_1,
    double a1_2, double a2_2, double a3_2, double a4_2, double a5_2, double b_2,
    double alpha_1, double beta_1, double gamma_1, double alpha_2, double beta_2, double gamma_2,
    double* res, int maxit, int numworkers, int m, double tol, double* const_time)
{
    int i, j = 0;

    int n = m * m;
    int* A_rows_1 = (int*)malloc(sizeof(int) * n * 5);
    int* A_cols_1 = (int*)malloc(sizeof(int) * n * 5);
    double* A_vals_1 = (double*)malloc(sizeof(double) * n * 5);
    int* A_dataLen_1 = (int*)malloc(sizeof(int) * (n + 1));
    double* q_1 = (double*)malloc(sizeof(double) * n);

    int* A_rows_2 = (int*)malloc(sizeof(int) * n * 5);
    int* A_cols_2 = (int*)malloc(sizeof(int) * n * 5);
    double* A_vals_2 = (double*)malloc(sizeof(double) * n * 5);
    int* A_dataLen_2 = (int*)malloc(sizeof(int) * (n + 1));
    double* q_2 = (double*)malloc(sizeof(double) * n);

    int* Omega_rows = (int*)malloc(sizeof(int) * n);
    int* Omega_cols = (int*)malloc(sizeof(int) * n);
    double* Omega_vals = (double*)malloc(sizeof(double) * n);

    int* A1SUBA2_rows = (int*)malloc(sizeof(int) * n * 5);
    int* A1SUBA2_cols = (int*)malloc(sizeof(int) * n * 5);
    double* A1SUBA2_vals = (double*)malloc(sizeof(double) * n * 5);
    int* A1SUBA2_dataLen = (int*)malloc(sizeof(int) * (n + 1));

    int* OmegaA1A2_rows = (int*)malloc(sizeof(int) * n * 5);
    int* OmegaA1A2_cols = (int*)malloc(sizeof(int) * n * 5);
    double* OmegaA1A2_vals = (double*)malloc(sizeof(double) * n * 5);
    int* OmegaA1A2_dataLen = (int*)malloc(sizeof(int) * (n + 1));

    double* xm = (double*)malloc(sizeof(double) * n);

    createSparseMatrix(a1_1, a2_1, a3_1, a4_1, a5_1, b_1, a1_2, a2_2, a3_2, a4_2, a5_2, b_2,
        A_rows_1, A_cols_1, A_vals_1, A_dataLen_1, q_1,
        A_rows_2, A_cols_2, A_vals_2, A_dataLen_2, q_2,
        Omega_rows, Omega_cols, Omega_vals,
        A1SUBA2_rows, A1SUBA2_cols, A1SUBA2_vals, A1SUBA2_dataLen,
        OmegaA1A2_rows, OmegaA1A2_cols, OmegaA1A2_vals, OmegaA1A2_dataLen,
        alpha_1, alpha_2,
        xm, m);

    /*

    for (size_t i = 0; i < n * 5; i++)
    {
        printf("%d,%d,%lf\n", A_rows_1[i], A_cols_1[i], A_vals_1[i]);
    }
    */


    /*
    for (size_t i = 0; i < n * 5; i++)
    {
        printf("%d,%d,%lf\n", A_rows_2[i], A_cols_2[i], A_vals_2[i]);

    }
    */

    //------------------------------初始化完成

    //------------------------------并行内存分配变量初始化
    int* s = (int*)malloc(sizeof(int) * numworkers);
    int* ind1 = (int*)malloc(sizeof(int) * numworkers);
    int* ind2 = (int*)malloc(sizeof(int) * numworkers);
    double* x = (double*)malloc(sizeof(double) * n);
    double* z = (double*)malloc(sizeof(double) * n);
    
    int maxRowLen = 0;
    memset(s, 0, sizeof(int) * numworkers);
    for (i = 0; i < numworkers; i++)
    {
        if ((n % numworkers) <= i)
        {
            s[i] = (int)floor(n / numworkers);
        }
        else
        {
            s[i] = (int)floor(n / numworkers) + 1;
        }
        if (s[i] > maxRowLen)
        {
            maxRowLen = s[i];
        }
        ind1[i] = 1 + calSSum(s, i);
        ind2[i] = calSSum(s, i) + s[i];
    }


    //------------------------------定义运算静态变量
    //------------------------------这里静态指不变的变量  

    int maxRowLen5numworkers = maxRowLen * 5 * numworkers;
    int maxRowLen1numworkers = (maxRowLen + 1) * numworkers;


    int* OAA_rows = (int*)malloc(sizeof(int) * maxRowLen5numworkers);
    int* OAA_cols = (int*)malloc(sizeof(int) * maxRowLen5numworkers);
    double* OAA_vals = (double*)malloc(sizeof(double) * maxRowLen5numworkers);
    int* OAA_dataLen = (int*)malloc(sizeof(int) * maxRowLen1numworkers);

    int* AA_rows = (int*)malloc(sizeof(int) * maxRowLen5numworkers);
    int* AA_cols = (int*)malloc(sizeof(int) * maxRowLen5numworkers);
    double* AA_vals = (double*)malloc(sizeof(double) * maxRowLen5numworkers);
    int* AA_dataLen = (int*)malloc(sizeof(int) * maxRowLen1numworkers);

    int* M_rows_1 = (int*)malloc(sizeof(int) * maxRowLen5numworkers);
    int* M_cols_1 = (int*)malloc(sizeof(int) * maxRowLen5numworkers);
    double* M_vals_1 = (double*)malloc(sizeof(double) * maxRowLen5numworkers);
    int* M_dataLen_1 = (int*)malloc(sizeof(int) * maxRowLen1numworkers);

    int* N_rows_1 = (int*)malloc(sizeof(int) * maxRowLen5numworkers);
    int* N_cols_1 = (int*)malloc(sizeof(int) * maxRowLen5numworkers);
    double* N_vals_1 = (double*)malloc(sizeof(double) * maxRowLen5numworkers);
    int* N_dataLen_1 = (int*)malloc(sizeof(int) * maxRowLen1numworkers);

    int* M_rows_2 = (int*)malloc(sizeof(int) * maxRowLen5numworkers);
    int* M_cols_2 = (int*)malloc(sizeof(int) * maxRowLen5numworkers);
    double* M_vals_2 = (double*)malloc(sizeof(double) * maxRowLen5numworkers);
    int* M_dataLen_2 = (int*)malloc(sizeof(int) * maxRowLen1numworkers);

    int* N_rows_2 = (int*)malloc(sizeof(int) * maxRowLen5numworkers);
    int* N_cols_2 = (int*)malloc(sizeof(int) * maxRowLen5numworkers);
    double* N_vals_2 = (double*)malloc(sizeof(double) * maxRowLen5numworkers);
    int* N_dataLen_2 = (int*)malloc(sizeof(int) * maxRowLen1numworkers);


    //额外开辟M1~M2~N1~N2~内存
    int* M_rows_1_ = (int*)malloc(sizeof(int) * maxRowLen5numworkers);
    int* M_cols_1_ = (int*)malloc(sizeof(int) * maxRowLen5numworkers);
    double* M_vals_1_ = (double*)malloc(sizeof(double) * maxRowLen5numworkers);
    int* M_dataLen_1_ = (int*)malloc(sizeof(int) * maxRowLen1numworkers);

    int* N_rows_1_ = (int*)malloc(sizeof(int) * maxRowLen5numworkers);
    int* N_cols_1_ = (int*)malloc(sizeof(int) * maxRowLen5numworkers);
    double* N_vals_1_ = (double*)malloc(sizeof(double) * maxRowLen5numworkers);
    int* N_dataLen_1_ = (int*)malloc(sizeof(int) * maxRowLen1numworkers);

    int* M_rows_2_ = (int*)malloc(sizeof(int) * maxRowLen5numworkers);
    int* M_cols_2_ = (int*)malloc(sizeof(int) * maxRowLen5numworkers);
    double* M_vals_2_ = (double*)malloc(sizeof(double) * maxRowLen5numworkers);
    int* M_dataLen_2_ = (int*)malloc(sizeof(int) * maxRowLen1numworkers);

    int* N_rows_2_ = (int*)malloc(sizeof(int) * maxRowLen5numworkers);
    int* N_cols_2_ = (int*)malloc(sizeof(int) * maxRowLen5numworkers);
    double* N_vals_2_ = (double*)malloc(sizeof(double) * maxRowLen5numworkers);
    int* N_dataLen_2_ = (int*)malloc(sizeof(int) * maxRowLen1numworkers);
    //end

    int* MM_rows = (int*)malloc(sizeof(int) * maxRowLen5numworkers);
    int* MM_cols = (int*)malloc(sizeof(int) * maxRowLen5numworkers);
    double* MM_vals = (double*)malloc(sizeof(double) * maxRowLen5numworkers);
    int* MM_dataLen = (int*)malloc(sizeof(int) * maxRowLen1numworkers);

    //额外产生的M1~-M2~
    int* MM_rows_ = (int*)malloc(sizeof(int) * maxRowLen5numworkers);
    int* MM_cols_ = (int*)malloc(sizeof(int) * maxRowLen5numworkers);
    double* MM_vals_ = (double*)malloc(sizeof(double) * maxRowLen5numworkers);
    int* MM_dataLen_ = (int*)malloc(sizeof(int) * maxRowLen1numworkers);


    //end
    double* T = (double*)malloc(sizeof(double) * maxRowLen * numworkers);
    double* v = (double*)malloc(sizeof(double) * maxRowLen * numworkers);

    memset(OAA_rows, 0, sizeof(int) * maxRowLen5numworkers);
    memset(OAA_cols, 0, sizeof(int) * maxRowLen5numworkers);
    memset(OAA_vals, 0, sizeof(double) * maxRowLen5numworkers);
    memset(OAA_dataLen, 0, sizeof(int) * maxRowLen1numworkers);
    memset(AA_rows, 0, sizeof(int) * maxRowLen5numworkers);
    memset(AA_cols, 0, sizeof(int) * maxRowLen5numworkers);
    memset(AA_vals, 0, sizeof(double) * maxRowLen5numworkers);
    memset(AA_dataLen, 0, sizeof(int) * maxRowLen1numworkers);

    memset(M_rows_1, 0, sizeof(int) * maxRowLen5numworkers);
    memset(M_cols_1, 0, sizeof(int) * maxRowLen5numworkers);
    memset(M_vals_1, 0, sizeof(double) * maxRowLen5numworkers);
    memset(M_dataLen_1, 0, sizeof(int) * maxRowLen1numworkers);
    memset(N_rows_1, 0, sizeof(int) * maxRowLen5numworkers);
    memset(N_cols_1, 0, sizeof(int) * maxRowLen5numworkers);
    memset(N_vals_1, 0, sizeof(double) * maxRowLen5numworkers);
    memset(N_dataLen_1, 0, sizeof(int) * maxRowLen1numworkers);
    memset(M_rows_2, 0, sizeof(int) * maxRowLen5numworkers);
    memset(M_cols_2, 0, sizeof(int) * maxRowLen5numworkers);
    memset(M_vals_2, 0, sizeof(double) * maxRowLen5numworkers);
    memset(M_dataLen_2, 0, sizeof(int) * maxRowLen1numworkers);
    memset(N_rows_2, 0, sizeof(int) * maxRowLen5numworkers);
    memset(N_cols_2, 0, sizeof(int) * maxRowLen5numworkers);
    memset(N_vals_2, 0, sizeof(double) * maxRowLen5numworkers);
    memset(N_dataLen_2, 0, sizeof(int) * maxRowLen1numworkers);

    memset(M_rows_1_, 0, sizeof(int) * maxRowLen5numworkers);
    memset(M_cols_1_, 0, sizeof(int) * maxRowLen5numworkers);
    memset(M_vals_1_, 0, sizeof(double) * maxRowLen5numworkers);
    memset(M_dataLen_1_, 0, sizeof(int) * maxRowLen1numworkers);
    memset(N_rows_1_, 0, sizeof(int) * maxRowLen5numworkers);
    memset(N_cols_1_, 0, sizeof(int) * maxRowLen5numworkers);
    memset(N_vals_1_, 0, sizeof(double) * maxRowLen5numworkers);
    memset(N_dataLen_1_, 0, sizeof(int) * maxRowLen1numworkers);
    memset(M_rows_2_, 0, sizeof(int) * maxRowLen5numworkers);
    memset(M_cols_2_, 0, sizeof(int) * maxRowLen5numworkers);
    memset(M_vals_2_, 0, sizeof(double) * maxRowLen5numworkers);
    memset(M_dataLen_2_, 0, sizeof(int) * maxRowLen1numworkers);
    memset(N_rows_2_, 0, sizeof(int) * maxRowLen5numworkers);
    memset(N_cols_2_, 0, sizeof(int) * maxRowLen5numworkers);
    memset(N_vals_2_, 0, sizeof(double) * maxRowLen5numworkers);
    memset(N_dataLen_2_, 0, sizeof(int) * maxRowLen1numworkers);

    memset(MM_rows, 0, sizeof(int) * maxRowLen5numworkers);
    memset(MM_cols, 0, sizeof(int) * maxRowLen5numworkers);
    memset(MM_vals, 0, sizeof(double) * maxRowLen5numworkers);
    memset(MM_dataLen, 0, sizeof(int) * maxRowLen1numworkers);

    memset(MM_rows_, 0, sizeof(int) * maxRowLen5numworkers);
    memset(MM_cols_, 0, sizeof(int) * maxRowLen5numworkers);
    memset(MM_vals_, 0, sizeof(double) * maxRowLen5numworkers);
    memset(MM_dataLen_, 0, sizeof(int) * maxRowLen1numworkers);


    memset(T, 0, sizeof(double) * maxRowLen * numworkers);
    memset(v, 0, sizeof(double) * maxRowLen * numworkers);


    //------------------------------初始化运算静态变量

    for (i = 0; i < numworkers; i++)
    {
        int M_rows = ind2[i] - ind1[i] + 1;

        int threadCur1 = i * (maxRowLen + 1);
        int threadCur2 = i * maxRowLen * 5;

        int OAA_tag = 0;
        int OAA_dataLenPerRow = 0;
        int OAA_curRow = ind1[i] - 1;
        int OAA_DataLenTag = 0;
        OAA_dataLen[threadCur1 + M_rows] = 0;

        int AA_tag = 0;
        int AA_dataLenPerRow = 0;
        int AA_curRow = ind1[i] - 1;
        int AA_DataLenTag = 0;
        AA_dataLen[threadCur1 + M_rows] = 0;

        for (j = 0; j < OmegaA1A2_dataLen[n]; j++)
        {
            if (OmegaA1A2_rows[j] >= (ind1[i] - 1) && OmegaA1A2_rows[j] < ind2[i])
            {
                if (OAA_curRow != OmegaA1A2_rows[j])
                {
                    OAA_dataLen[threadCur1 + OAA_DataLenTag] = OAA_dataLenPerRow;
                    OAA_DataLenTag++;
                    OAA_dataLen[threadCur1 + M_rows] = OAA_dataLen[threadCur1 + M_rows] + OAA_dataLenPerRow;
                    OAA_curRow = OmegaA1A2_rows[j];
                    OAA_dataLenPerRow = 0;
                }
                OAA_rows[threadCur2 + OAA_tag] = OmegaA1A2_rows[j] - (ind1[i] - 1);
                OAA_cols[threadCur2 + OAA_tag] = OmegaA1A2_cols[j];
                OAA_vals[threadCur2 + OAA_tag] = OmegaA1A2_vals[j];
                OAA_tag++;
                OAA_dataLenPerRow++;
            }
        }
        OAA_dataLen[threadCur1 + OAA_DataLenTag] = OAA_dataLenPerRow;
        OAA_dataLen[threadCur1 + M_rows] = OAA_dataLen[threadCur1 + M_rows] + OAA_dataLenPerRow;
        
        /*
        int cursorOAAal = 0;
        for (size_t j = 0; j < M_rows; j++)
        {
            for (l = 0; l < OAA_dataLen[j]; l++)
            {
                printf("%d ,%d ,%d,%lf\n", OAA_dataLen[j], OAA_rows[cursorOAAal + l], OAA_cols[cursorOAAal + l], OAA_vals[cursorOAAal + l]);
                //printf("%d ,%d ,%d,%lf\n", MM_dataLen_[j], threadCur1_t + MM_rows_[cursorMMVal + l], threadCur1_t + MM_cols_[cursorMMVal + l], MM_vals_[cursorMMVal + l]);
            }
            cursorOAAal = cursorOAAal + OAA_dataLen[j];
        }
        */

        for (j = 0; j < A1SUBA2_dataLen[n]; j++)
        {
            if (A1SUBA2_rows[j] >= (ind1[i] - 1) && A1SUBA2_rows[j] < ind2[i])
            {
                if (AA_curRow != A1SUBA2_rows[j])
                {
                    AA_dataLen[threadCur1 + AA_DataLenTag] = AA_dataLenPerRow;
                    AA_DataLenTag++;
                    AA_dataLen[threadCur1 + M_rows] = AA_dataLen[threadCur1 + M_rows] + AA_dataLenPerRow;
                    AA_curRow = A1SUBA2_rows[j];
                    AA_dataLenPerRow = 0;
                }
                AA_rows[threadCur2 + AA_tag] = A1SUBA2_rows[j] - (ind1[i] - 1);
                AA_cols[threadCur2 + AA_tag] = A1SUBA2_cols[j];
                AA_vals[threadCur2 + AA_tag] = A1SUBA2_vals[j];
                AA_tag++;
                AA_dataLenPerRow++;
            }
        }
        AA_dataLen[threadCur1 + AA_DataLenTag] = AA_dataLenPerRow;
        AA_dataLen[threadCur1 + M_rows] = AA_dataLen[threadCur1 + M_rows] + AA_dataLenPerRow;

        //初始化M1M2N1N2 M1-M2

        int Mtag_1 = 0;
        int MdataLenPerRow_1 = 0;
        int McurRow_1 = (ind1[i] - 1);
        int MDataLenTag_1 = 0;
        M_dataLen_1[threadCur1 + M_rows] = 0;

        int Ntag_1 = 0;
        int NdataLenPerRow_1 = 0;
        int NcurRow_1 = (ind1[i] - 1);
        int NDataLenTag_1 = 0;
        N_dataLen_1[threadCur1 + M_rows] = 0;

        double M_N_vals_1_tmp = 0;

        for (j = 0; j < A_dataLen_1[n]; j++)
        {

            if (A_rows_1[j] >= ind2[i])
            {
                break;
            }
            M_N_vals_1_tmp = 0;
            if (A_rows_1[j] >= (ind1[i] - 1) && A_rows_1[j] < ind2[i] )
            {

                if (A_rows_1[j] >= A_cols_1[j] && A_cols_1[j] >= (ind1[i] - 1))
                {
                    if (McurRow_1 != A_rows_1[j])
                    {
                        M_dataLen_1[threadCur1 + MDataLenTag_1] = MdataLenPerRow_1;
                        MDataLenTag_1++;
                        M_dataLen_1[threadCur1 + M_rows] = M_dataLen_1[threadCur1 + M_rows] + MdataLenPerRow_1;
                        McurRow_1 = A_rows_1[j];
                        MdataLenPerRow_1 = 0;
                    }

                    M_rows_1[threadCur2 + Mtag_1] = A_rows_1[j] - (ind1[i] - 1);
                    M_cols_1[threadCur2 + Mtag_1] = A_cols_1[j] - (ind1[i] - 1);
                    if (A_rows_1[j] == A_cols_1[j])
                    {
                        M_vals_1[threadCur2 + Mtag_1] = A_vals_1[j] / alpha_1;
                    }
                    else
                    {
                        M_vals_1[threadCur2 + Mtag_1] = beta_1 * A_vals_1[j] / alpha_1;
                    }

                    M_N_vals_1_tmp = M_vals_1[threadCur2 + Mtag_1];


                    Mtag_1++;
                    MdataLenPerRow_1++;
                }

                M_N_vals_1_tmp = M_N_vals_1_tmp - A_vals_1[j];

                //if (fabs(M_N_vals_1_tmp) > 1.0e-16)
            //    {
                if (NcurRow_1 != A_rows_1[j])
                {
                    N_dataLen_1[threadCur1 + NDataLenTag_1] = NdataLenPerRow_1;
                    NDataLenTag_1++;
                    N_dataLen_1[threadCur1 + M_rows] = N_dataLen_1[threadCur1 + M_rows] + NdataLenPerRow_1;
                    NcurRow_1 = A_rows_1[j];
                    NdataLenPerRow_1 = 0;
                }
                N_rows_1[threadCur2 + Ntag_1] = A_rows_1[j] - (ind1[i] - 1);
                N_cols_1[threadCur2 + Ntag_1] = A_cols_1[j];
                N_vals_1[threadCur2 + Ntag_1] = M_N_vals_1_tmp;
                Ntag_1++;
                NdataLenPerRow_1++;

                //  }
            }
        }

        M_dataLen_1[threadCur1 + MDataLenTag_1] = MdataLenPerRow_1;
        M_dataLen_1[threadCur1 + M_rows] = M_dataLen_1[threadCur1 + M_rows] + MdataLenPerRow_1;
        N_dataLen_1[threadCur1 + NDataLenTag_1] = NdataLenPerRow_1;
        N_dataLen_1[threadCur1 + M_rows] = N_dataLen_1[threadCur1 + M_rows] + NdataLenPerRow_1;

        int Mtag_2 = 0;
        int MdataLenPerRow_2 = 0;
        int McurRow_2 = (ind1[i] - 1);
        int MDataLenTag_2 = 0;
        M_dataLen_2[threadCur1 + M_rows] = 0;

        int Ntag_2 = 0;
        int NdataLenPerRow_2 = 0;
        int NcurRow_2 = (ind1[i] - 1);
        int NDataLenTag_2 = 0;
        N_dataLen_2[threadCur1 + M_rows] = 0;

        double M_N_vals_2_tmp = 0;

        for (j = 0; j < A_dataLen_2[n]; j++)
        {
            if (A_rows_2[j] >= ind2[i])
            {
                break;
            }
            M_N_vals_2_tmp = 0;
            if (A_rows_2[j] >= (ind1[i] - 1) && A_rows_2[j] < ind2[i])
            {
                if (A_rows_2[j] >= A_cols_2[j] && A_cols_2[j] >= (ind1[i] - 1))
                {
                    if (McurRow_2 != A_rows_2[j])
                    {
                        M_dataLen_2[threadCur1 + MDataLenTag_2] = MdataLenPerRow_2;
                        MDataLenTag_2++;
                        M_dataLen_2[threadCur1 + M_rows] = M_dataLen_2[threadCur1 + M_rows] + MdataLenPerRow_2;
                        McurRow_2 = A_rows_2[j];
                        MdataLenPerRow_2 = 0;
                    }

                    M_rows_2[threadCur2 + Mtag_2] = A_rows_2[j] - (ind1[i] - 1);
                    M_cols_2[threadCur2 + Mtag_2] = A_cols_2[j] - (ind1[i] - 1);
                    if (A_rows_2[j] == A_cols_2[j])
                    {
                        M_vals_2[threadCur2 + Mtag_2] = A_vals_2[j] / alpha_2;
                    }
                    else
                    {
                        M_vals_2[threadCur2 + Mtag_2] = beta_2 * A_vals_2[j] / alpha_2;
                    }

                    M_N_vals_2_tmp = M_vals_2[threadCur2 + Mtag_2];

                    Mtag_2++;
                    MdataLenPerRow_2++;
                }

                M_N_vals_2_tmp = M_N_vals_2_tmp - A_vals_2[j];

                //if (fabs(M_N_vals_2_tmp) > 1.0e-16)
            //    {
                if (NcurRow_2 != A_rows_2[j])
                {
                    N_dataLen_2[threadCur1 + NDataLenTag_2] = NdataLenPerRow_2;
                    NDataLenTag_2++;
                    N_dataLen_2[threadCur1 + M_rows] = N_dataLen_2[threadCur1 + M_rows] + NdataLenPerRow_2;
                    NcurRow_2 = A_rows_2[j];
                    NdataLenPerRow_2 = 0;
                }
                N_rows_2[threadCur2 + Ntag_2] = A_rows_2[j] - (ind1[i] - 1);
                N_cols_2[threadCur2 + Ntag_2] = A_cols_2[j];
                N_vals_2[threadCur2 + Ntag_2] = M_N_vals_2_tmp;
                Ntag_2++;
                NdataLenPerRow_2++;
                //   }
            }
        }

        M_dataLen_2[threadCur1 + MDataLenTag_2] = MdataLenPerRow_2;
        M_dataLen_2[threadCur1 + M_rows] = M_dataLen_2[threadCur1 + M_rows] + MdataLenPerRow_2;
        N_dataLen_2[threadCur1 + NDataLenTag_2] = NdataLenPerRow_2;
        N_dataLen_2[threadCur1 + M_rows] = N_dataLen_2[threadCur1 + M_rows] + NdataLenPerRow_2;

        int MMtag = 0;
        int MMdataLenPerRow = 0;
        int MMcurRow = 0;
        int MMDataLenTag = 0;
        MM_dataLen[threadCur1 + M_rows] = 0;

        double M_M_vals_tmp = 0;
        int M1_cursor = 0;
        int M2_cursor = 0;
        // int OmegaA1A2_cursor = 0;

        for (j = 0; j < M_rows; j++)
        {
            if (MMcurRow != j)
            {
                MM_dataLen[threadCur1 + MMDataLenTag] = MMdataLenPerRow;
                MMDataLenTag++;
                MM_dataLen[threadCur1 + M_rows] = MM_dataLen[threadCur1 + M_rows] + MMdataLenPerRow;
                MMcurRow = j;
                MMdataLenPerRow = 0;
            }
            if (j - m >= 0)
            {
                M_M_vals_tmp = 0;

                if (M_cols_1[threadCur2 + M1_cursor] == j - m)
                {
                    M_M_vals_tmp += M_vals_1[threadCur2 + M1_cursor];
                    M1_cursor += 1;
                }

                if (M_cols_2[threadCur2 + M2_cursor] == j - m)
                {
                    M_M_vals_tmp += M_vals_2[threadCur2 + M2_cursor];
                    M2_cursor++;
                }

                //if (fabs(M_M_vals_tmp) > 1.0e-16)
            //    {
                MM_rows[threadCur2 + MMtag] = j;
                MM_cols[threadCur2 + MMtag] = j - m;
                MM_vals[threadCur2 + MMtag] = M_M_vals_tmp;

                MMtag++;
                MMdataLenPerRow++;
                //    }
            }
            if (j - 1 >= 0)
            {
                M_M_vals_tmp = 0;

                if (M_cols_1[threadCur2 + M1_cursor] == j - 1)
                {
                    M_M_vals_tmp += M_vals_1[threadCur2 + M1_cursor];
                    M1_cursor++;
                }

                if (M_cols_2[threadCur2 + M2_cursor] == j - 1)
                {
                    M_M_vals_tmp += M_vals_2[threadCur2 + M2_cursor];
                    M2_cursor++;
                }

                //if (fabs(M_M_vals_tmp) > 1.0e-16)
            //    {
                MM_rows[threadCur2 + MMtag] = j;
                MM_cols[threadCur2 + MMtag] = j - 1;
                MM_vals[threadCur2 + MMtag] = M_M_vals_tmp;
                MMtag++;
                MMdataLenPerRow++;
                //   }
            }

            M_M_vals_tmp = 2 * Omega_vals[j + (ind1[i] - 1)];

            if (M_cols_1[threadCur2 + M1_cursor] == j)
            {
                M_M_vals_tmp += M_vals_1[threadCur2 + M1_cursor];
                M1_cursor++;
            }

            if (M_cols_2[threadCur2 + M2_cursor] == j)
            {
                M_M_vals_tmp += M_vals_2[threadCur2 + M2_cursor];
                M2_cursor++;
            }

            //if (fabs(M_M_vals_tmp) > 1.0e-16)
          //  {
            MM_rows[threadCur2 + MMtag] = j;
            MM_cols[threadCur2 + MMtag] = j;
            MM_vals[threadCur2 + MMtag] = M_M_vals_tmp;
            MMtag++;
            MMdataLenPerRow++;
            //    }
        }

        MM_dataLen[threadCur1 + MMDataLenTag] = MMdataLenPerRow;



        //初始化 M1~M2~N1~N2~ M1~ - M2~ 
        Mtag_1 = 0;
        MdataLenPerRow_1 = 0;
        McurRow_1 = (ind1[i] - 1);
        MDataLenTag_1 = 0;
        M_dataLen_1_[threadCur1 + M_rows] = 0;

        Ntag_1 = 0;
        NdataLenPerRow_1 = 0;
        NcurRow_1 = (ind1[i] - 1);
        NDataLenTag_1 = 0;
        N_dataLen_1_[threadCur1 + M_rows] = 0;

        M_N_vals_1_tmp = 0;

        for (j = 0; j < A_dataLen_1[n]; j++)
        {
            if (A_rows_1[j] >= ind2[i])
            {
                break;
            }
            M_N_vals_1_tmp = 0;
            if (A_cols_1[j] >= ind2[i])continue;
            if (A_rows_1[j] >= (ind1[i] - 1) && A_rows_1[j] < ind2[i] && A_cols_1[j] < ind2[i])
            {

                if (A_rows_1[j] <= A_cols_1[j] && A_cols_1[j] >= (ind1[i] - 1))
                {
                    if (McurRow_1 != A_rows_1[j])
                    {
                        M_dataLen_1_[threadCur1 + MDataLenTag_1] = MdataLenPerRow_1;
                        MDataLenTag_1++;
                        M_dataLen_1_[threadCur1 + M_rows] = M_dataLen_1_[threadCur1 + M_rows] + MdataLenPerRow_1;
                        McurRow_1 = A_rows_1[j];
                        MdataLenPerRow_1 = 0;
                    }

                    M_rows_1_[threadCur2 + Mtag_1] = A_rows_1[j] - (ind1[i] - 1);
                    M_cols_1_[threadCur2 + Mtag_1] = A_cols_1[j] - (ind1[i] - 1);
                    if (A_rows_1[j] == A_cols_1[j])
                    {
                        M_vals_1_[threadCur2 + Mtag_1] = A_vals_1[j] / alpha_1;
                    }
                    else
                    {
                        M_vals_1_[threadCur2 + Mtag_1] = beta_1 * A_vals_1[j] / alpha_1;
                    }

                    M_N_vals_1_tmp = M_vals_1_[threadCur2 + Mtag_1];
                    Mtag_1++;
                    MdataLenPerRow_1++;

                }

                M_N_vals_1_tmp = M_N_vals_1_tmp - A_vals_1[j];

                if (NcurRow_1 != A_rows_1[j])
                {
                    N_dataLen_1_[threadCur1 + NDataLenTag_1] = NdataLenPerRow_1;
                    NDataLenTag_1++;
                    N_dataLen_1_[threadCur1 + M_rows] = N_dataLen_1_[threadCur1 + M_rows] + NdataLenPerRow_1;
                    NcurRow_1 = A_rows_1[j];
                    NdataLenPerRow_1 = 0;
                }
                N_rows_1_[threadCur2 + Ntag_1] = A_rows_1[j] - (ind1[i] - 1);
                N_cols_1_[threadCur2 + Ntag_1] = A_cols_1[j];
                N_vals_1_[threadCur2 + Ntag_1] = M_N_vals_1_tmp;
                Ntag_1++;
                NdataLenPerRow_1++;
                //  }
            }
        }
        M_dataLen_1_[threadCur1 + MDataLenTag_1] = MdataLenPerRow_1;
        M_dataLen_1_[threadCur1 + M_rows] = M_dataLen_1[threadCur1 + M_rows] + MdataLenPerRow_1;
        N_dataLen_1_[threadCur1 + NDataLenTag_1] = NdataLenPerRow_1;
        N_dataLen_1_[threadCur1 + M_rows] = N_dataLen_1[threadCur1 + M_rows] + NdataLenPerRow_1;


        Mtag_2 = 0;
        MdataLenPerRow_2 = 0;
        McurRow_2 = (ind1[i] - 1);
        MDataLenTag_2 = 0;
        M_dataLen_2[threadCur1 + M_rows] = 0;

        Ntag_2 = 0;
        NdataLenPerRow_2 = 0;
        NcurRow_2 = (ind1[i] - 1);
        NDataLenTag_2 = 0;
        N_dataLen_2[threadCur1 + M_rows] = 0;
        M_N_vals_2_tmp = 0;

        for (j = 0; j < A_dataLen_2[n]; j++)
        {
            if (A_rows_2[j] < ind1[i] - 1)continue;
            if (A_rows_2[j] >= ind2[i])break;
            M_N_vals_2_tmp = 0;

            if (A_rows_2[j] <= A_cols_2[j] && A_cols_2[j] < ind2[i])
            {
                if (McurRow_2 != A_rows_2[j])
                {
                    M_dataLen_2_[threadCur1 + MDataLenTag_2] = MdataLenPerRow_2;
                    MDataLenTag_2++;
                    M_dataLen_2_[threadCur1 + M_rows] = M_dataLen_2_[threadCur1 + M_rows] + MdataLenPerRow_2;
                    McurRow_2 = A_rows_2[j];
                    MdataLenPerRow_2 = 0;
                }

                M_rows_2_[threadCur2 + Mtag_2] = A_rows_2[j] - (ind1[i] - 1);
                M_cols_2_[threadCur2 + Mtag_2] = A_cols_2[j] - (ind1[i] - 1);
                if (A_rows_2[j] == A_cols_2[j])
                {
                    M_vals_2_[threadCur2 + Mtag_2] = A_vals_2[j] / alpha_2;
                }
                else
                {
                    M_vals_2_[threadCur2 + Mtag_2] = beta_2 * A_vals_2[j] / alpha_2;
                }

                M_N_vals_2_tmp = M_vals_2_[threadCur2 + Mtag_2];

                Mtag_2++;
                MdataLenPerRow_2++;
            }

            M_N_vals_2_tmp = M_N_vals_2_tmp - A_vals_2[j];

            //if (fabs(M_N_vals_2_tmp) > 1.0e-16)
        //    {
            if (NcurRow_2 != A_rows_2[j])
            {
                N_dataLen_2_[threadCur1 + NDataLenTag_2] = NdataLenPerRow_2;
                NDataLenTag_2++;
                N_dataLen_2_[threadCur1 + M_rows] = N_dataLen_2[threadCur1 + M_rows] + NdataLenPerRow_2;
                NcurRow_2 = A_rows_2[j];
                NdataLenPerRow_2 = 0;
            }


            N_rows_2_[threadCur2 + Ntag_2] = A_rows_2[j] - (ind1[i] - 1);
            N_cols_2_[threadCur2 + Ntag_2] = A_cols_2[j];
            N_vals_2_[threadCur2 + Ntag_2] = M_N_vals_2_tmp;
            Ntag_2++;
            NdataLenPerRow_2++;
            //   }
        }
        M_dataLen_2_[threadCur1 + MDataLenTag_2] = MdataLenPerRow_2;
        M_dataLen_2_[threadCur1 + M_rows] = M_dataLen_2_[threadCur1 + M_rows] + MdataLenPerRow_2;
        N_dataLen_2_[threadCur1 + NDataLenTag_2] = NdataLenPerRow_2;
        N_dataLen_2_[threadCur1 + M_rows] = N_dataLen_2_[threadCur1 + M_rows] + NdataLenPerRow_2;


        MMtag = 0;
        MMdataLenPerRow = 0;
        MMcurRow = 0;
        MMDataLenTag = 0;
        MM_dataLen[threadCur1 + M_rows] = 0;

        M_M_vals_tmp = 0;
        M1_cursor = 0;
        M2_cursor = 0;
        // int OmegaA1A2_cursor = 0;

        for (j = 0; j < M_rows; j++)
        {
            if (j - 1 >= 0)
            {
                M_M_vals_tmp = 0;

                if (M_cols_1_[threadCur2 + M1_cursor] == j)
                {
                    M_M_vals_tmp += M_vals_1_[threadCur2 + M1_cursor];
                    M1_cursor++;
                }

                if (M_cols_2_[threadCur2 + M2_cursor] == j)
                {
                    M_M_vals_tmp += M_vals_2_[threadCur2 + M2_cursor];
                    M2_cursor++;
                }


                MM_rows_[threadCur2 + MMtag] = j - 1;
                MM_cols_[threadCur2 + MMtag] = j;
                MM_vals_[threadCur2 + MMtag] = M_M_vals_tmp;
                MMtag++;
                MMdataLenPerRow++;

                if (j + m - 1 < M_rows) {
                    M_M_vals_tmp = 0;

                    if (M_cols_1_[threadCur2 + M1_cursor] == j - 1 + m)
                    {
                        M_M_vals_tmp += M_vals_1_[threadCur2 + M1_cursor];
                        M1_cursor++;
                    }

                    if (M_cols_2_[threadCur2 + M2_cursor] == j - 1 + m)
                    {
                        M_M_vals_tmp += M_vals_2_[threadCur2 + M2_cursor];
                        M2_cursor++;
                    }


                    MM_rows_[threadCur2 + MMtag] = j - 1;
                    MM_cols_[threadCur2 + MMtag] = j - 1 + m;
                    MM_vals_[threadCur2 + MMtag] = M_M_vals_tmp;
                    MMtag++;
                    MMdataLenPerRow++;
                }
            }


            if (MMcurRow != j)
            {
                MM_dataLen_[threadCur1 + MMDataLenTag] = MM_dataLen_[threadCur1 + MMDataLenTag] + MMdataLenPerRow;
                MMDataLenTag++;
                MM_dataLen_[threadCur1 + M_rows] = MM_dataLen_[threadCur1 + M_rows] + MMdataLenPerRow;
                MMcurRow = j;
                MMdataLenPerRow = 0;
            }


            M_M_vals_tmp = 2 * Omega_vals[j + (ind1[i] - 1)];

            if (M_cols_1_[threadCur2 + M1_cursor] == j)
            {
                M_M_vals_tmp += M_vals_1_[threadCur2 + M1_cursor];
                M1_cursor++;
            }

            if (M_cols_2_[threadCur2 + M2_cursor] == j)
            {
                M_M_vals_tmp += M_vals_2_[threadCur2 + M2_cursor];
                M2_cursor++;
            }

            MM_rows_[threadCur2 + MMtag] = j;
            MM_cols_[threadCur2 + MMtag] = j;
            MM_vals_[threadCur2 + MMtag] = M_M_vals_tmp;
            MMtag++;
            MMdataLenPerRow++;

        }
        MM_dataLen_[threadCur1 + MMDataLenTag] = MMdataLenPerRow;
        MM_dataLen_[threadCur1 + M_rows] = MM_dataLen_[threadCur1 + M_rows] + MMdataLenPerRow;
        
        /*
        int cursorMMVal = 0;
        for (size_t j = 0; j < M_rows; j++)
        {
            for (l = 0; l < MM_dataLen_[j]; l++)
            {
                MM_dataLen_[j];
                MM_rows_[cursorMMVal + l];
                if (MM_cols_[cursorMMVal + l] > M_rows)
                    printf("cccc");
                MM_vals_[cursorMMVal + l];
            }
            cursorMMVal += MM_dataLen_[j];
        }
        */
    }


    //------------------------------迭代
    int* z_cursor1 = (int*)malloc(sizeof(int) * n);
    int* z_cursor2 = (int*)malloc(sizeof(int) * n);
    for (i = 0; i < numworkers; i++)
    {
        int j;
        z_cursor1[i] = 0;
        z_cursor2[i] = 0;
        for (j = 0; j < ind1[i] - 1; j++)
        {
            z_cursor1[i] = z_cursor1[i] + A_dataLen_1[j];
            z_cursor2[i] = z_cursor2[i] + A_dataLen_2[j];
        }
    }
    
    int k, l,M_rows,iters,threadCursor_1,threadCursor_2,threadCursor_3,cursorN1Val,cursorN2Val,cursorOAAVal,cursorMMVal,cursorAAVal = 0;
    double AAValTmp = 0.;
    *const_time = clock();
    
    #pragma acc data copyin(ind1 [0:numworkers], ind2 [0:numworkers], T [0:maxRowLen * numworkers], v [0:maxRowLen * numworkers], q_1 [0:n], q_2 [0:n],xm [0:n],\
            N_dataLen_1 [0:maxRowLen1numworkers], N_dataLen_1_ [0:maxRowLen1numworkers], N_cols_1 [0:maxRowLen5numworkers], \
            N_cols_1_ [0:maxRowLen5numworkers], N_vals_1 [0:maxRowLen5numworkers],N_vals_1_[0:maxRowLen5numworkers],                                                                                \
            N_dataLen_2 [0:maxRowLen1numworkers], N_dataLen_2_ [0:maxRowLen1numworkers], N_cols_2 [0:maxRowLen5numworkers], \
            N_cols_2_ [0:maxRowLen5numworkers], N_vals_2 [0:maxRowLen5numworkers],N_vals_2_[0:maxRowLen5numworkers],                                                                                                  \
            A_dataLen_1[0:n],A_vals_1[0:n*5-1],A_cols_1[0:n*5-1], \
            A_dataLen_2[0:n],A_vals_2[0:n*5-1],A_cols_2[0:n*5-1], \
            OAA_dataLen [0:maxRowLen1numworkers], OAA_vals [0:maxRowLen5numworkers], OAA_cols [0:maxRowLen5numworkers],                                                                                                             \
            AA_dataLen [0:maxRowLen1numworkers], AA_vals [0:maxRowLen5numworkers], AA_cols [0:maxRowLen5numworkers],                                                                                                                \
            MM_dataLen_ [0:maxRowLen1numworkers], MM_vals_ [0:maxRowLen5numworkers], MM_cols_ [0:maxRowLen5numworkers],                                                                                                              \
            MM_dataLen [0:maxRowLen1numworkers], MM_vals [0:maxRowLen5numworkers], MM_cols [0:maxRowLen5numworkers])                                                                                                                \
            copy(x [0:n], z [0:n])
    {
        for (iters = 0; iters < maxit; iters++){
            #pragma acc parallel loop num_gangs(numworkers) num_workers(1) vector_length(1)\
            private(i,j, k, l,M_rows,iters,threadCursor_1,threadCursor_2,threadCursor_3,cursorN1Val,\
            cursorN2Val,cursorOAAVal,cursorMMVal,cursorAAVal,AAValTmp)\
            firstprivate(gamma_1, gamma_2, maxRowLen)
            for (i = 0; i < numworkers; i++){
                if (runAlgorithmTag > 0){
                    M_rows = ind2[i] - ind1[i] + 1;

                    threadCursor_1 = i * maxRowLen;
                    threadCursor_2 = i * (maxRowLen + 1);
                    threadCursor_3 = i * maxRowLen * 5;

                    cursorN1Val = 0;
                    cursorN2Val = 0;
                    cursorOAAVal = 0;
                    cursorAAVal = 0;
                    AAValTmp = 0;
                    cursorMMVal = 0;

                    for (j = 0; j < M_rows; j++)
                    {
                        T[threadCursor_1 + j] = -gamma_1 * q_1[j + ind1[i] - 1] - gamma_2 * q_2[j + ind1[i] - 1];

                        for (k = 0; k < N_dataLen_1[threadCursor_2 + j]; k++)
                        {
                            T[threadCursor_1 + j] = T[threadCursor_1 + j] + N_vals_1[threadCursor_3 + cursorN1Val + k] * xm[N_cols_1[threadCursor_3 + cursorN1Val + k]];
                        }
                        cursorN1Val = cursorN1Val + N_dataLen_1[threadCursor_2 + j];
                        for (k = 0; k < N_dataLen_2[threadCursor_2 + j]; k++)
                        {
                            T[threadCursor_1 + j] = T[threadCursor_1 + j] + N_vals_2[threadCursor_3 + cursorN2Val + k] * xm[N_cols_2[threadCursor_3 + cursorN2Val + k]];
                        }
                        cursorN2Val = cursorN2Val + N_dataLen_2[threadCursor_2 + j];

                        for (k = 0; k < OAA_dataLen[threadCursor_2 + j]; k++)
                        {
                            T[threadCursor_1 + j] = T[threadCursor_1 + j] + OAA_vals[threadCursor_3 + cursorOAAVal + k] * fabs(xm[OAA_cols[threadCursor_3 + cursorOAAVal + k]]);
                        }
                        cursorOAAVal = cursorOAAVal + OAA_dataLen[threadCursor_2 + j];
                        AAValTmp = 0;
                        for (k = 0; k < AA_dataLen[threadCursor_2 + j]; k++)
                        {
                            AAValTmp = AAValTmp + AA_vals[threadCursor_3 + cursorAAVal + k] * (fabs(xm[AA_cols[threadCursor_3 + cursorAAVal + k]]) + xm[AA_cols[threadCursor_3 + cursorAAVal + k]]);
                        }
                        T[threadCursor_1 + j] = T[threadCursor_1 + j] + fabs(AAValTmp + gamma_1 * q_1[j + ind1[i] - 1] - gamma_2 * q_2[j + ind1[i] - 1]);
                        cursorAAVal = cursorAAVal + AA_dataLen[threadCursor_2 + j];
                        //printf("%d  \t %1.15f\n", threadCursor_1 + j, T[threadCursor_1 + j]);

                        v[threadCursor_1 + j] = T[threadCursor_1 + j];
                        for (l = 0; l < MM_dataLen[threadCursor_2 + j] - 1; l++)
                        {
                            v[threadCursor_1 + j] = v[threadCursor_1 + j] - MM_vals[threadCursor_3 + cursorMMVal + l] * v[threadCursor_1 + MM_cols[threadCursor_3 + cursorMMVal + l]];
                        }
                        v[threadCursor_1 + j] = v[threadCursor_1 + j] / MM_vals[threadCursor_3 + cursorMMVal + MM_dataLen[threadCursor_2 + j] - 1];

                        cursorMMVal = cursorMMVal + MM_dataLen[threadCursor_2 + j];

                        x[j + ind1[i] - 1] = v[threadCursor_1 + j];
                        //printf("%d  \t %1.15f\n", iters, x[j + ind1[i] - 1]);
                        z[j + ind1[i] - 1] = x[j + ind1[i] - 1] + fabs(x[j + ind1[i] - 1]);
                    }
                }
            }

            //printf("算法3\n");
            //算法三中这里的x将只是xi二分之一
            //-------------------后续将xi看做xi二分之一参与运算最后重新赋值------------------------------------------
            #pragma acc parallel loop num_gangs(numworkers) num_workers(1) vector_length(1)\
            private(i,j, k, l,M_rows,iters,threadCursor_1,threadCursor_2,threadCursor_3,cursorN1Val,\
            cursorN2Val,cursorOAAVal,cursorMMVal,cursorAAVal,AAValTmp)\
            firstprivate(gamma_1, gamma_2, maxRowLen)   
            for (i = 0; i < numworkers; i++){
                if (runAlgorithmTag > 1){
                    M_rows = ind2[i] - ind1[i] + 1;

                    threadCursor_1 = i * maxRowLen;
                    threadCursor_2 = i * (maxRowLen + 1);
                    threadCursor_3 = i * maxRowLen * 5;

                    cursorN1Val = 0;
                    cursorN2Val = 0;
                    cursorOAAVal = 0;
                    cursorAAVal = 0;
                    cursorMMVal = 0;

                        
                    for (j = 0; j < M_rows; j++){
                        
                        T[threadCursor_1 + j] = -gamma_1 * q_1[j + ind1[i] - 1] - gamma_2 * q_2[j + ind1[i] - 1];
                        
                        for (k = 0; k < N_dataLen_1_[threadCursor_2 + j]; k++)
                        {
                            T[threadCursor_1 + j] = T[threadCursor_1 + j] + N_vals_1_[threadCursor_3 + cursorN1Val + k] * x[N_cols_1_[threadCursor_3 + cursorN1Val + k]];
                        }
                        
                        cursorN1Val = cursorN1Val + N_dataLen_1_[threadCursor_2 + j];
                        for (k = 0; k < N_dataLen_2_[threadCursor_2 + j]; k++)
                        {
                            T[threadCursor_1 + j] = T[threadCursor_1 + j] + N_vals_2_[threadCursor_3 + cursorN2Val + k] * x[N_cols_2_[threadCursor_3 + cursorN2Val + k]];
                        }

                        cursorN2Val = cursorN2Val + N_dataLen_2_[threadCursor_2 + j];
                        for (k = 0; k < OAA_dataLen[threadCursor_2 + j]; k++)
                        {
                            T[threadCursor_1 + j] = T[threadCursor_1 + j] + OAA_vals[threadCursor_3 + cursorOAAVal + k] * fabs(x[OAA_cols[threadCursor_3 + cursorOAAVal + k]]);
                        }
                        cursorOAAVal = cursorOAAVal + OAA_dataLen[threadCursor_2 + j];

                        AAValTmp = 0;
                        for (k = 0; k < AA_dataLen[threadCursor_2 + j]; k++)
                        {
                            double kk = AA_vals[threadCursor_3 + cursorAAVal + k] * (fabs(x[AA_cols[threadCursor_3 + cursorAAVal + k]]) + x[AA_cols[threadCursor_3 + cursorAAVal + k]]);
                            AAValTmp = AAValTmp + AA_vals[threadCursor_3 + cursorAAVal + k] * (fabs(x[AA_cols[threadCursor_3 + cursorAAVal + k]]) + x[AA_cols[threadCursor_3 + cursorAAVal + k]]);
                        }
                        T[threadCursor_1 + j] = T[threadCursor_1 + j] + fabs(AAValTmp + gamma_1 * q_1[j + ind1[i] - 1] - gamma_2 * q_2[j + ind1[i] - 1]);
                        //printf("%lf\n", T[threadCursor_1 + j]);
                        cursorAAVal = cursorAAVal + AA_dataLen[threadCursor_2 + j];
                        
                    }
                    
                    cursorMMVal = MM_dataLen_[threadCursor_2 + M_rows] - 1;
                    for (j = M_rows - 1; j > -1; j--) {
                        v[threadCursor_1 + j] = T[threadCursor_1 + j];
                        //printf("%lf\t", T[threadCursor_1 + j]);
                        for (l = 0; l < MM_dataLen_[threadCursor_2 + j] - 1; l++)
                        {
                            v[threadCursor_1 + j] = v[threadCursor_1 + j] - MM_vals_[threadCursor_3 + cursorMMVal - l] * v[threadCursor_1 + MM_cols_[threadCursor_3 + cursorMMVal - l]];
                        }
                        //printf("%lf\t", v[threadCursor_1 + j]);
                        //printf("%d  \t %lf\n", threadCursor_1 + j, MM_vals_[threadCursor_3 + cursorMMVal - l]);
                        v[threadCursor_1 + j] = v[threadCursor_1 + j] / MM_vals_[threadCursor_3 + cursorMMVal - l];

                        cursorMMVal = cursorMMVal - MM_dataLen_[threadCursor_2 + j];

                        x[j + ind1[i] - 1] = v[threadCursor_1 + j];

                        //printf("%d  \t %lf\n", threadCursor_1 + j, x[j + ind1[i] - 1]);
                        z[j + ind1[i] - 1] = x[j + ind1[i] - 1] + fabs(x[j + ind1[i] - 1]);
                    }
                    
                }
            }


            //printf('%d' ,iters);
            //终止条件
            if (runAlgorithmTag > 0)
            {
                double* z_sum = (double*)malloc(sizeof(double) * numworkers);
                double zminval, ztmpval,zsum = 0;
                int cursor_1,cursor_2;
            
                #pragma acc parallel loop num_gangs(numworkers) num_workers(1) vector_length(1) \
                private(i,j, k, l,zminval, ztmpval,zsum,cursor_1,cursor_2)
                for (i = 0; i < numworkers; i++ )
                {

                    cursor_1 = z_cursor1[i], cursor_2 = z_cursor2[i];
                    zminval = 0, ztmpval = 0;
                    z_sum[i] = 0;

                    for (j = ind1[i] - 1; j < ind2[i]; j++)
                    {
                        xm[j] = x[j];
                        zminval = z[j];
                        ztmpval = 0;
                        
                        for (k = 0; k < A_dataLen_1[j]; k++)
                        {
                            ztmpval = ztmpval + A_vals_1[cursor_1 + k] * z[A_cols_1[cursor_1 + k]];
                        }
                        ztmpval += q_1[j];
                        cursor_1 = cursor_1 + A_dataLen_1[j];
                        if (ztmpval < zminval)
                        {
                            zminval = ztmpval;
                        }
                        ztmpval = 0;
                        for (k = 0; k < A_dataLen_2[j]; k++)
                        {
                            ztmpval = ztmpval + A_vals_2[cursor_2 + k] * z[A_cols_2[cursor_2 + k]];
                        }
                        ztmpval += q_2[j];
                        cursor_2 = cursor_2 + A_dataLen_2[j];
                        if (ztmpval < zminval)
                        {
                            zminval = ztmpval;
                        }
                        z_sum[i] = z_sum[i] + (zminval * zminval);
                    }

                }

                for (i = 0; i < numworkers; i++)
                {
                    zsum = zsum + z_sum[i];
                }
                
            
                res[iters] = sqrt(zsum);

                if (res[iters] < tol)
                {
                    iters += 1;
                    break;
                }
            }  
        }    
    }

    *const_time = clock() - *const_time;

    free(v);
    free(T);
    free(MM_dataLen_);
    free(MM_vals_);
    free(MM_cols_);
    free(MM_rows_);
    free(N_dataLen_2_);
    free(N_vals_2_);
    free(N_cols_2_);
    free(N_rows_2_);
    free(M_dataLen_2_);
    free(M_vals_2_);
    free(M_cols_2_);
    free(M_rows_2_);
    free(N_dataLen_1_);
    free(N_vals_1_);
    free(N_cols_1_);
    free(N_rows_1_);
    free(M_dataLen_1_);
    free(M_vals_1_);
    free(M_cols_1_);
    free(M_rows_1_);

    free(MM_dataLen);
    free(MM_vals);
    free(MM_cols);
    free(MM_rows);
    free(N_dataLen_2);
    free(N_vals_2);
    free(N_cols_2);
    free(N_rows_2);
    free(M_dataLen_2);
    free(M_vals_2);
    free(M_cols_2);
    free(M_rows_2);
    free(N_dataLen_1);
    free(N_vals_1);
    free(N_cols_1);
    free(N_rows_1);
    free(M_dataLen_1);
    free(M_vals_1);
    free(M_cols_1);
    free(M_rows_1);

    free(AA_dataLen);
    free(AA_vals);
    free(AA_cols);
    free(AA_rows);

    free(OAA_dataLen);
    free(OAA_vals);
    free(OAA_cols);
    free(OAA_rows);

    free(z);
    free(x);
    free(ind2);
    free(ind1);
    free(s);
    free(xm);

    free(OmegaA1A2_dataLen);
    free(OmegaA1A2_vals);
    free(OmegaA1A2_cols);
    free(OmegaA1A2_rows);
    free(A1SUBA2_dataLen);
    free(A1SUBA2_vals);
    free(A1SUBA2_cols);
    free(A1SUBA2_rows);
    free(Omega_vals);
    free(Omega_cols);
    free(Omega_rows);
    free(q_2);
    free(A_dataLen_2);
    free(A_vals_2);
    free(A_cols_2);
    free(A_rows_2);
    free(q_1);
    free(A_dataLen_1);
    free(A_vals_1);
    free(A_cols_1);
    free(A_rows_1);

    return iters;

}


int main(int argc,char *argv[])
{
    
    int numworkers = atoi(argv[1]);
    int m = atoi(argv[2]);
    int runAlgorithmTag = atoi(argv[3]);
    int gpu_id = atoi(argv[4]);
    acc_set_device_num(gpu_id,acc_device_nvidia);

    char outfile[200];
    sprintf(outfile,"%s%s%s%s%s%s%s","out_put/numworkers_",argv[1],"_m_",argv[2],"_runAlgorithmTag_",argv[3],".txt");
    
    int maxit = 1000;
    double tol = 1.0e-6;
    double* res = (double*)malloc(sizeof(double) * maxit);
    
    double const_time = 0.0;
    int iters = MSM(runAlgorithmTag, 0, -1, 5, -1, 0, 1, -1, -1, 5, -1, -1, 1, 1, 1, 1, 1, 1, 1, res, maxit, numworkers, m, tol, &const_time);
    //int iters = MSM(runAlgorithmTag, 0, -1.5, 4, -0.5, 0, 0, -1.5, -1.5, 4, -0.5, -0.5, 1,1, 1, 1, 1, 1, 1,res, maxit, numworkers, m, tol);

    // int iters = MSM(runAlgorithmTag, 0, -1, 4, -1, 0, 1, -1, -1, 4, -1, -1, 1, 1, 0, 1, 1, 0, 1, res, maxit, numworkers, m, tol);
//  int iters = MSM(runAlgorithmTag, 0, -1, 4, -1, 0, 1, -1, -1, 4, -1, -1, 1, 1, 1, 1, 1, 1, 1, res, maxit, numworkers, m, tol);
//  int iters = MSM(runAlgorithmTag, 0, -1, 4, -1, 0, 1, -1, -1, 4, -1, -1, 1, 1.1, 1.1, 1, 1.1, 1.1, 1, res, maxit, numworkers, m, tol);

 //int iters = MSM(runAlgorithmTag, 0, -1.5, 4, -0.5, 0, 1, -1.5, -1.5, 4, -0.5, -0.5, 1, 1, 0, 1, 1, 0, 1, res, maxit, numworkers, m, tol);
 //int iters = MSM(runAlgorithmTag, 0, -1.5, 4, -0.5, 0, 1, -1.5, -1.5, 4, -0.5, -0.5, 1, 1, 1, 1, 1, 1, 1, res, maxit, numworkers, m, tol);
//   int iters = MSM(runAlgorithmTag, 0, -1.5, 4, -0.5, 0, 1, -1.5, -1.5, 4, -0.5, -0.5, 1, 1.1, 1.1, 1, 1.1, 1.1, 1, res, maxit, numworkers, m, tol);
    FILE * pFile;
    pFile = fopen(outfile, "w");
    fprintf(pFile,"%lfs\n", const_time / CLOCKS_PER_SEC);
    fprintf(pFile,"%d\n", iters);
    
    for (int i = 0; i < iters; i++)
    {
        fprintf(pFile,"%d: %0.15f\n", i + 1, res[i]);
    }
    
    free(res);
    return 0;
}


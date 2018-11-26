#include <carvig.h>

#define M 6
#define N 5
#define LDA M

int main(int argc, char **argv)
{
    double U[M*M],W[M],V[N*N];
    double a[LDA*N] = {
            8.79,  6.11, -9.15,  9.57, -3.49,  9.84,
            9.93,  6.91, -7.93,  1.64,  4.02,  0.15,
            9.83,  5.04,  4.86,  8.83,  9.80, -8.99,
            5.45, -0.27,  4.85,  0.74, 10.00, -6.02,
            3.16,  7.98,  3.01,  5.80,  4.27, -5.31
    };
    /* svd */
    svd(a,M,N,U,W,V);

    tracemat(3,U,M,M,12,6);

    tracemat(3,V,N,N,12,6);

    tracemat(3,W,1,N,12,6);
}
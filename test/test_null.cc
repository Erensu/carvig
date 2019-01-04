#include <carvig.h>

int main(int argc, char **argv)
{
    double A[18]={
            1,3,3,
            2,4,7,
            5,1,8,
            8,2,9,
            12,7,11,
            1,4,2
    };
    double a[30] = {
            8.79,  6.11, -9.15,  9.57, -3.49,  9.84,
            9.93,  6.91, -7.93,  1.64,  4.02,  0.15,
            9.83,  5.04,  4.86,  8.83,  9.80, -8.99,
            5.45, -0.27,  4.85,  0.74, 10.00, -6.02,
            3.16,  7.98,  3.01,  5.80,  4.27, -5.31
    };
    double N[18]; int p,q;
    double Na[30]; int pa,qa;
    double at[30];

    matt(a,6,5,at);

    null(A,3,6,N,&p,&q);

    null(at,5,6,Na,&pa,&qa);

    tracemat(3,at,5,6,12,6);

    tracemat(3,N,p,q,12,6);
    tracemat(3,Na,pa,qa,12,6);

    return 0;
}
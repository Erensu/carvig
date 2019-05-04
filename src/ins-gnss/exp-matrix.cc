/* ----------------------------------------------------------------------------
 * Dubious Ways to Implement the Exponential of a Matrix
 * a small C++ library for computing the exponential of a matrix. It is inspired
 * by the paper "Nineteen Dubious Ways to Compute the Exponential of a Matrix,
 * Twenty-Five Years Later" by Cleve Moler and Charles Van Loan. Note that I have
 * no affiliation with either author.
 *
 * The library is distributed under the [Boost Software License]
 * (<http://www.boost.org/users/license.html>).
 * ---------------------------------------------------------------------------*/
#include "carvig.h"
#include "r8lib.h"

/* calculate the matrix exponential (e^A)-------------------------------------
 * args:  double *A  input matrix: A (col-major)
 *        int n      rows and cols of A
 * return: matrix of (e^A)
 * ---------------------------------------------------------------------------*/
extern double* expm(const double *A,int n)
{
    double *a2,a_norm,c,*d,*e,t,*x,*a=mat(n,n),*E=mat(n,n);
    int ee,k,s,p;
    const double one=1.0;
    const int q=6;

    matt(A,n,n,a);

    a2=r8mat_copy_new(n,n,a);
    a_norm=r8mat_norm_li(n,n,a2);
    ee=(int)(r8_log_2(a_norm))+1;
    s=i4_max(0,ee+1);
    t=1.0/pow(2.0,s);

    r8mat_scale(n,n,t,a2);
    x=r8mat_copy_new(n,n,a2);
    c=0.5;

    e=r8mat_identity_new(n);
    r8mat_add(n,n,one,e,c,a2,e);
    d=r8mat_identity_new(n);
    r8mat_add(n,n,one,d,-c,a2,d);
    p=1;

    for (k=2;k<=q;k++) {
        c=c*(double)(q-k+1)/(double)(k*(2*q-k+1));
        r8mat_mm(n,n,n,a2,x,x);

        r8mat_add(n,n,c,x,one,e,e);
        if (p) {
            r8mat_add(n,n,c,x,one,d,d);
        }
        else {
            r8mat_add(n,n,-c,x,one,d,d);
        }
        p=!p;
    }
    r8mat_minvm(n,n,d,e,e);
    for (k=1;k<=s;k++) {
        r8mat_mm(n,n,n,e,e,e);
    }
    matt(e,n,n,E);
    free(a2);
    free(a);
    free(d);
    free(x); free(e);
    return E;
}

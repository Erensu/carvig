#include <carvig.h>

/* compute camera project matrix----------------------------------------------*/
static int prjmatrix(const double *C,const double *t,const double *K,double *P)
{
    double T[16],M[12];
    int i,j;

    rt2tf(C,t,T);
    for (i=0;i<3;i++) for (j=0;j<4;j++) M[i+j*3]=T[i+j*4];
    matmul("NN",3,4,3,1.0,K,M,0.0,P);
    return 1;
}

static void test_1()
{
    double C[9],t[3],P1[12],P2[12];
    double dT[16]={
            0.999999000000000,7.80000000000000e-05,-0.00158500000000000,0,
            -7.80000000000000e-05,1,-0.000187000000000000,0,
            0.00158500000000000,0.000187000000000000,0.999999000000000,0,
            0.00147100000000000,0.00359500000000000,0.0923580000000000,1
    };
    double J[16],V[16];
    double K[9]={
            500,0,0,
            0,500,0,
            250,250,1
    };
    double uv1[2]={305,301};
    double uv2[2]={298.2413, 299.5223};

    double uv3[2]={150,346};
    double uv4[2]={ 147.9870,  351.0920};

    double uv5[2]={-248,451};
    double uv6[2]={-252.4763,450.5836};

    double pf[2];
    int j;

    matinv(dT,4); tf2rt(dT,C,t);

    tracemat(3,C,3,3,12,6);
    tracemat(3,t,1,3,12,3);

    /* camera project matrix */
    prjmatrix(C,t,K,P2); seteye(C,3); setzero(t,1,3);
    prjmatrix(C,t,K,P1);
    for (j=0;j<4;j++) {
        J[0+4*j]=P1[2+3*j]*uv1[0]-P1[0+3*j];
        J[1+4*j]=P1[2+3*j]*uv1[1]-P1[1+3*j];
        J[2+4*j]=P2[2+3*j]*uv2[0]-P2[0+3*j];
        J[3+4*j]=P2[2+3*j]*uv2[1]-P2[1+3*j];
    }
    if (!svd(J,4,4,NULL,NULL,V)) return;

    /* return false if this point is at infinity */
    if (fabs(V[3+3*4])<1E-10) {
        trace(2,"feature point is at infinity\n");
        return;
    }
    for (j=0;j<3;j++) pf[j]=V[j+3*4]/V[3+3*4];
    trace(3,"feature position: %8.4lf  %8.4lf  %8.4lf\n",pf[0],pf[1],pf[2]);
}

static int bearingvector(const double *uv,
                         const double *K,double *bv)
{
    double img_p[2],uimg_p[2];

    trace(2,"bearingvector:\n");

    /* unscale and center */
    img_p[0]=(uv[0]-K[6])/K[0];
    img_p[1]=(uv[1]-K[7])/K[4];

    /* project 1 into z direction and normalization*/
    bv[0]=img_p[0];
    bv[1]=img_p[1];
    bv[2]=1.0;
    return 1;
}

static void test_2()
{
    double R[9],t[3],P1[12],P2[12],b1[3],b2[3];
    double dT[16]={
            0.999999000000000,7.80000000000000e-05,-0.00158500000000000,0,
            -7.80000000000000e-05,1,-0.000187000000000000,0,
            0.00158500000000000,0.000187000000000000,0.999999000000000,0,
            0.00147100000000000,0.00359500000000000,0.0923580000000000,1
    };
    double K[9]={
            500,0,0,
            0,500,0,
            250,250,1
    };
    double uv1[2]={-251,452};
    double uv2[2]={ -250.4763, 449.5836};

    double pf[2],bt[3],b[3],lambda[3],A[4],xm[3],xn[3];
    int i;

    tf2rt(dT,R,t);

    /* bearing vector: cam1->uv1, cam2->uv2 */
    bearingvector(uv1,K,b1);
    bearingvector(uv2,K,b2);

    matmul("NN",3,1,3,1.0,R,b2,0.0,bt);

    b[0]=dot(t,b1,3);
    b[1]=dot(t,bt,3);

    A[0]= dot(b1,b1,3); A[1]=dot(b1,bt,3); A[2]=-A[1];
    A[3]=-dot(bt,bt,3);

    if (matinv(A,2)) return;
    matmul("NN",2,1,2,1.0,A,b,0.0,lambda);

    for (i=0;i<3;i++) {
        xm[i]=lambda[0]*b1[i]; xn[i]=t[i]+lambda[1]*bt[i];
        pf[i]=(xm[i]+xn[i])/2.0;
    }
    trace(3,"feature position: %8.4lf  %8.4lf  %8.4lf\n",pf[0],pf[1],pf[2]);
}

static void test_3()
{
    double R[9],t[3];
    double dT[16]={
            0.999999000000000,7.80000000000000e-05,-0.00158500000000000,0,
            -7.80000000000000e-05,1,-0.000187000000000000,0,
            0.00158500000000000,0.000187000000000000,0.999999000000000,0,
            0.00147100000000000,0.00359500000000000,1.5924,1
    };
    double K[9]={
            500,0,0,
            0,500,0,
            250,250,1
    };
    double uv[2]={150,350};
    double uvp[2];
    tf2rt(dT,R,t);

    predictfeat(R,t,K,uv,uvp);

    trace(3,"feature position: %8.4lf  %8.4lf \n",uvp[0],uvp[1]);
}

int main(int argc, char **argv)
{
    test_1();
    test_2();
    test_3();
    return 1;
}
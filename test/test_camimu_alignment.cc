#include <carvig.h>
int main(int argc, char **argv)
{
    double Cbc[9]={
            0,0,1,
            1,0,0,
            0,1,0
    };
    double rpy[3];

    dcm2rpy(Cbc,rpy);

    tracemat(3,rpy,1,3,12,6);

    rpy2dcm(rpy,Cbc);

    tracemat(3,Cbc,3,3,12,6);

    /*-------------------------------------------------------*/
    double rpy0[3]={51.0*D2R,3*D2R,7*D2R};
    double rpy1[3]={10.0*D2R,15*D2R,20*D2R};
    double C0[9],C1[9];
    double Cc0[9],Cc1[9];
    double dC[9],dCc[9],dCz[9],phic[3];
    double phi[3],phiz[3];

    rpy2dcm(rpy0,C0);
    rpy2dcm(rpy1,C1);
    matmul("NN",3,3,3,1.0,Cbc,C0,0.0,Cc0);
    matmul("NN",3,3,3,1.0,Cbc,C1,0.0,Cc1);

    matmul("TN",3,3,3,1.0,C0,C1,0.0,dC);
    matmul("TN",3,3,3,1.0,Cc0,Cc1,0.0,dCc);

    so3_log(dCc,phic,NULL);
    //matmul("TN",3,1,3,1.0,Cbc,phic,0.0,phiz);

    so3_exp(phic,dCz);

    so3_log(dC ,phi ,NULL);
    so3_log(dCz,phiz,NULL);

    double dCC[9];
    matmul("TN",3,3,3,1.0,dC,dCz,0.0,dCC);

    trace(3,"dCC =\n"); tracemat(3,dCC ,3,3,12,5);

    trace(3,"phi =\n"); tracemat(3,phi ,1,3,12,5);
    trace(3,"phiz=\n"); tracemat(3,phiz,1,3,12,5);

    return 0;
}
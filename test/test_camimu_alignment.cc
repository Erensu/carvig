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

    return 0;
}
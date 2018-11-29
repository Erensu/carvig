#include <carvig.h>
#include <iostream>
#include <png++/png.hpp>
using namespace std;

/* generates a uniform distributed random number between min and max --------*/
static double getuniform(double min,double max)
{
    return 1.0*rand()/RAND_MAX*(max-min)+min;
}
/* creates gaussian distributed random numbers (Box-Müller method)-----------*/
static double getgaussian(double std)
{
    if (std<0.0) std=-std;
    double x1,x2,w,y1;
    do {
        x1=getuniform(-1.0,1.0); x2=getuniform(-1.0,1.0); w=x1*x1+x2*x2;
    } while (w>=1.0);
    w=sqrt((-2.0*log(w))/w); y1=x1*w; return std*y1;
}
int main(int argc, char **argv)
{
    feature feats[6];
    double R[9]={
            -0.3536,0.9330,-0.0670,
            0.3536,0.0670,-0.9330,
            -0.8660,-0.3536,-0.3536
    };
    double t[3]={20,-30,100};
    double Re[9],te[3];

    double xp[3*6]={
            0,0,0,
            20,0,40,
            10,-10,0,
            15,-5,24,
            14,-2,21,
            16,-5,13
    };
    voopt_t opt={0};
    seteye(opt.cam.K,3);

    double mp[3*6];
    for (int i=0;i<6;i++) {
        matmul("NN",3,1,3,1.0,R,xp+3*i,0.0,mp+3*i);
        for (int j=0;j<3;j++) mp[3*i+j]+=t[j];

        feats[i].u=mp[3*i+0]/mp[3*i+2];
        feats[i].v=mp[3*i+1]/mp[3*i+2];
    }
    p3pthree(feats,6,xp,6,&opt,Re,te);

    tracemat(3,Re,3,3,12,5);
    tracemat(3,te,1,3,12,5);
}
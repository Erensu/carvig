#include <include/carvig.h>
#include "carvig.h"

static void test_C()
{
    double C1[9]={
            0.9999690 ,  -0.0008225  ,  0.0078361  ,
            0.0008216 ,   0.9999997  ,  0.0001192  ,
            -0.0078362,   -0.0001128 ,   0.9999693 ,
    };
    double C2[9]={
            0.999999   ,  0.000078   , -0.001585  ,
            -0.000078  ,   1.000000  ,  -0.000187 ,
            0.001585   ,  0.000187   ,  0.999999
    };
    double C1t[9],C2t[9];
    matt(C1,3,3,C1t);
    matt(C2,3,3,C2t);

    double rpy1[3],rpy2[3];

    dcm2rot(C1t,rpy1);
    dcm2rot(C2t,rpy2);

    rpy1[0]=-rpy1[0];
    rpy1[1]=-rpy1[1];
    rpy1[2]=-rpy1[2];

    rv2m(rpy1,C1t); tracemat(3,C1t,3,3,12,8);

    tracemat(3,rpy1,3,1,12,8);
    tracemat(3,rpy2,3,1,12,8);
}

static void test_kinematicsecef()
{
    double Cbe1[9]={
            0.0627179158608598,	0.769273836348678,	-0.635833490577099    ,
            0.998031293612615,	-0.0483424238536104,	0.0399568145112667 ,
            -1.14444047162599e-10,	-0.637087729253472,	-0.770791298105171 ,
    };
    double Cbe0[9]={
            0.0627179403692973,	0.769273835161543,	-0.635833489595888   ,
            0.998031292072465,	-0.0483424427445007	,0.0399568301252915   ,
            -1.14444047162599e-10,	-0.637087729253472,	-0.770791298105171
    };
    double ve1[3]={
            0.627132905021703,
            9.97957692406307 ,
            0                ,
    };
    double ve0[3]={
            0.627179402311078,
            9.98031292081149 ,
            0
    };
    double re[3]={
            4063553.78987791 ,
            -255360.453867144,
            4893080.18618536
    };
    double Cbe1t[9],Cbe0t[9],fb[3],omgb[3];
    matt(Cbe1,3,3,Cbe1t);
    matt(Cbe0,3,3,Cbe0t);

    kinematicsecef(0.01,Cbe1t,Cbe0t,ve1,ve0,re,fb,omgb);

    tracemat(3,fb,3,1,12,8);
    tracemat(3,omgb,3,1,12,8);
}

static void test_align0()
{
    double vn0[3]={
            0,10.0,0
    };
    double rn0[3]={
            50.425,-3.5958333,50
    };
}

static int test_file(char **file)
{
    fprintf(stderr,"%s\n",file[0]);
}

int main(int argc, char **argv)
{

    test_kinematicsecef();

    test_C();

    const char *datafile="/home/sujinglan/carvig/test/data/Profile_2.csv";
    const char *vofile="./vo-sim.out";
    const char *imufile="./imu-sim.out";
    const char *gpsfile="./gps-sim.out";

    imu_err_t err={0};
    insopt_t insopt={0};
    cam_t cam={0};

    cam.K[0]=500; cam.K[4]=500;
    cam.K[6]=250; cam.K[7]=250;
    cam.K[8]=1.0;

#if 1
    err.ba[0]=900*9.80665E-5;
    err.ba[1]=-1300*9.80665E-5;
    err.ba[2]=800*9.80665E-5;

    err.bg[0]=-900.0*D2R/3600;
    err.bg[1]=1300*D2R/3600;
    err.bg[2]=-800*D2R/3600;

    err.wgn[0]=0.01*D2R/60.0;
    err.wgn[1]=0.01*D2R/60.0;
    err.wgn[2]=0.01*D2R/60.0;

    err.wan[0]=100*9.80665E-6;
    err.wan[1]=100*9.80665E-6;
    err.wan[2]=100*9.80665E-6;

    /*
    err.Ma[0]=500*1E-6; err.Ma[3]=-300*1E-6; err.Ma[6]=200*1E-6;
    err.Ma[1]=-150*1E-6; err.Ma[4]=-600*1E-6; err.Ma[7]=250*1E-6;
    err.Ma[2]=-250*1E-6; err.Ma[5]=100*1E-6; err.Ma[8]=450*1E-6;

    err.Mg[0]=400*1E-6; err.Mg[3]=-300*1E-6; err.Mg[6]=250*1E-6;
    err.Mg[1]=0.0; err.Mg[4]=-300*1E-6; err.Mg[7]=-150*1E-6;
    err.Mg[2]=0.0; err.Mg[5]=0.0; err.Mg[8]=-350*1E-6;

    err.Gg[0]=0.9; err.Gg[3]=-1.1; err.Gg[6]=-0.6;
    err.Gg[1]=-0.5; err.Gg[4]=1.9; err.Gg[7]=-1.6;
    err.Gg[2]=0.3; err.Gg[5]=-1.1; err.Gg[8]=-1.3;

    for (int i=0;i<9;i++) err.Gg[i]*=(D2R/(3600*9.80665));
     */
#endif

    insopt.voopt.ebc[0]=1.57079632;
    insopt.voopt.ebc[1]=0.0;
    insopt.voopt.ebc[2]=1.57079632;

    insopt.voopt.lbc[0]=0.065;
    insopt.voopt.lbc[1]=-0.206;
    insopt.voopt.lbc[2]=0.002;

    insopt.lever[0]=-0.701;
    insopt.lever[1]=-0.247;
    insopt.lever[2]=0.030;

    generatepath(datafile,&cam,&err,&insopt,imufile,gpsfile,vofile);
    return 0;
}

/*--------------------------------------------------------------------------------
 * ins-gnss-vo-lc.cc : ins-gnss-vo coupled common functions
 *
 * reference :
 *    [01] P.D.Groves, Principles of GNSS, Intertial, and Multisensor Integrated
 *         Navigation System, Artech House, 2008
 *    [02] Tedaldi D, Pretto A, Menegatti E. A robust and easy to implement method
 *         for IMU calibration without external equipments,2014.
 *    [03] Skog I, Handel P. Effects of time synchronization errors in GNSS-aided
 *         INS 2008.
 *    [04] Shin E H. Accuracy Improvement of Low Cost INS/GPS for Land Applications
 *    [05] Li M, Mourikis A I. High-precision, consistent EKF-based visual–inertial
 *         odometry[J].International Journal of Robotics Research,2013,32(6):690-711.
 *    [06] Monocular Visual Inertial Odometry on a Mobile Device.
 *    [07] Mourikis A / , Roumeliotis S / . A Multi-State Constraint Kalman Filter
 *         for Vision-aided Inertial Navigation[C]// IEEE International Conference
 *         on Robotics & Automation. IEEE, 2007.
 *    [08] Forster C , Carlone L , Dellaert F , et al. On-Manifold Preintegration
 *         for Real-Time Visual-Inertial Odometry[J]. IEEE Transactions on Robotics,
 *         2015, 33(1):1-21.
 *    [09] Observability-constrained vision-aided inertial navigation.
 *    [10] Li M ,Mourikis A I. Improving the accuracy of EKF-based visual-inertial
 *         odometry[C]// IEEE International Conference on Robotics & Automation.
 *         IEEE, 2012.
 *    [11] Li M ,Mourikis A I. High-precision, consistent EKF-based visual-inertial
 *         odometry[M].Sage Publications, Inc. 2013.
 *    [12] Pizzoli M , Forster C , Scaramuzza D . REMODE: Probabilistic, Monocular
 *         Dense Reconstruction in Real Time[C]// IEEE International Conference on
 *         Robotics and Automation (ICRA), Hong Kong, 2014. IEEE, 2014.
 *    [13] Li M , Kim B H , Mourikis A I . Real-time motion tracking on a cellphone
 *         using inertial sensing and a rolling-shutter camera[C]// IEEE International
 *         Conference on Robotics & Automation. IEEE, 2013.
 *    [14] Monocular Visual-Inertial SLAM and Self Calibration for Long Term Autonomy
 *    [15] Lupton T , Sukkarieh S . Visual-Inertial-Aided Navigation for High-Dynamic
 *         Motion in Built Environments Without Initial Conditions[J]. IEEE Transactions
 *         on Robotics, 2012, 28(1):61-76.
 *
 * version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
 * history : 2018/12/02 1.0 new
 *-----------------------------------------------------------------------------*/
#include "carvig.h"

/* constants-------------------------------------------------------------------*/
#define MIN_TRACK_LEN           2             /* min length of tracking data */
#define MAX_TRACK_LEN           50            /* max length of tracking data */
#define MAXRES_POSE             10.0*D2R      /* max residual for pose filter (rad) */
#define MAXRES_POS              3.00          /* max residual for position filter (m) */
#define VARPOSE                 3.00*D2R      /* variance of pose measurement (rad) */
#define VARPOS                  0.50          /* variance of position measurement (m) */
#define SWAP(type,x,y)          {type tmp; tmp=x; x=y; y=tmp;}

typedef struct {                              /* filter workspace */
    double Cbe[2][9],re[2][3];                /* ins states in precious and current epoch */
    double Cce[2][9],rc[2][3];                /* camera states in precious and current epoch */
    double *Px;                               /* error states covariance matrix */
    double scale;                             /* mono visual odometry scale */
    int nx,n,flag;                            /* number of error states */
} filt_t;

/* global variables------------------------------------------------------------*/
static track_t  tracks={0};                   /* all tracking feature points data */
static match_t  matchs={0};                   /* match feature points data */
static filt_t   filts ={0};                   /* vo aid filter workspace */

/* initial vo-aid-------------------------------------------------------------*/
extern void initvoaidlc(insopt_t *opt)
{
    trace(3,"initvoaidlc:\n");

    opt->voopt.match.f =opt->voopt.calib.f;
    opt->voopt.match.fu=opt->voopt.calib.fu;
    opt->voopt.match.fv=opt->voopt.calib.fv;

    opt->voopt.match.cu=opt->voopt.calib.cu;
    opt->voopt.match.cv=opt->voopt.calib.cv;

    init_match(&matchs,&opt->voopt.match);
    inittrackimgbuf(&opt->voopt);
    return;
}
/* free vo-aid----------------------------------------------------------------*/
extern void freevoaidlc()
{
    trace(3,"freevoaidlc:\n");

    freetrackset(&tracks);
    free_match(&matchs);
    if (filts.Px) free(filts.Px); filts.Px=NULL;
    freetrackimgbuf();
}
/* initial filter workspace----------------------------------------------------*/
static int initfilt(const insstate_t *ins)
{
   if (filts.Px==NULL) {

       filts.nx=2*ins->nx+1;
       filts.Px=mat(filts.nx,filts.nx);
       filts.scale=1.0;
   }
}
/* set filter ins states-------------------------------------------------------*/
static void setfiltinsstat(const insstate_t *ins,int pos)
{
    matcpy(filts.Cbe[pos],ins->Cbe,3,3);
    matcpy(filts.re [pos],ins->re ,3,1);
}
/* set filter camera states----------------------------------------------------*/
static void setfiltcamstat(const vostate_t *vo,int pos)
{
    matcpy(filts.Cce[pos],vo->Cce,3,3);
    matcpy(filts.rc [pos],vo->rc ,3,1);
}
/* set filter states-----------------------------------------------------------*/
static void setfiltstate(const insstate_t *ins,const vostate_t *vo,int pos)
{
    setfiltinsstat(ins,pos%2);
    setfiltcamstat(vo ,pos%2);
    filts.n++;
}
/* set filter covariance matrix of ins states----------------------------------*/
static void setfiltcov(const insstate_t *ins)
{
    int i,j,nx=ins->nx;
    double *T;

    T=mat(nx,nx);
    matmul("NN",nx,nx,nx,1.0,ins->F,ins->P,0.0,T);

    for (i=0;i<nx;i++) for (j=nx;j<2*nx;j++) {
            filts.Px[i+j*filts.nx]=filts.Px[j+i*filts.nx]=T[i+(j-nx)*nx];
        }
    matmul("NNT",nx,nx,nx,1.0,ins->F,ins->P,0.0,T);
    for (i=0;i<nx;i++) {
        for (j=0;j<nx;j++) filts.Px[i+j*filts.nx]=T[i+j*nx];
    }
    filts.Px[2*nx+2*nx*filts.nx]+=SQR(1.0);

    trace(3,"Px=\n");
    tracemat(3,filts.Px,nx,nx,12,5);

    free(T);
}
/* estimate mono-camera motion-------------------------------------------------*/
static int estmotion(const match_set *pset,const voopt_t *opt,double *dT)
{
    return estmonort(opt,pset,dT);
}
/* update camera states--------------------------------------------------------*/
static int updatecamera(vostate_t *vo,const insopt_t *opt,const double *dT,
                        gtime_t time)
{
    double Tp[16],T[16],TT[16],rp[3],dt;

    if ((dt=timediff(time,vo->time))>10.0) {
        trace(2,"update camera pose failed\n");
        return 0;
    }
    if (fabs(dt)<1E-6) {
        trace(2,"duplicate camera update\n");
        return 0;
    }
    matcpy(rp,vo->rc,1,3);
    rt2tf(vo->Cce,vo->rc,Tp);

    matcpy(TT,dT,4,4);
    matinv(TT,4);
    matmul("NN",4,4,4,1.0,Tp,TT,0.0,T);

    tf2rt(T,vo->Cce,vo->rc);
    matmul("NN",4,4,4,1.0,vo->T,TT,0.0,T);
    matcpy(vo->T,T,4,4);

    vo->time=time;
    return 1;
}
/* update all track data------------------------------------------------------*/
static int updatetrack(const insopt_t *opt,const img_t *img)
{
    match_set *pset=&matchs.mp_dense;

#if USE_BUCKET_FEATS
    pset=&matchs.mp_bucket;
#endif
    /* update track data */
    if (!match2track(pset,matchs.pt,matchs.time,img->id,&matchs.Ip,&matchs.Ic,
                     &opt->voopt,&tracks)) {

        trace(2,"update track fail\n");
        return 0;
    }
    return 1;
}
/* ins states convert to camera states-----------------------------------------*/
static void ins2camera(const insstate_t *ins,vostate_t *vo,int flag)
{
    double T[3];
    int i;
    matmul("NT",3,3,3,1.0,ins->Cbe,ins->Cbc,0.0,vo->Cce);
    matmul("NN",3,1,3,1.0,ins->Cbe,ins->lbc,0.0,T);
    for (i=0;i<3;i++) {
        vo->rc[i]=ins->re[i]+T[i];
    }
    if (flag) {
        matcpy(vo->Cce0,vo->Cce,3,3);
        matcpy(vo->rc0 ,vo->rc ,1,3);
    }
    vo->time=ins->time;
}
/* ins states convert to camera states-----------------------------------------*/
static void camera2ins(const double *Cce,const double *rc,const double *Cbc,
                       const double *lbc,double *Cbe,double *re)
{
    double T[3];
    int i;
    matmul("NN",3,3,3,1.0,Cce,Cbc,0.0,Cbe);
    matmul("NN",3,1,3,1.0,Cbe,lbc,0.0,T);
    for (i=0;i<3&&re;i++) re[i]=rc[i]-T[i];
}
/* update track data-----------------------------------------------------------*/
static int updatealltrack(const insstate_t *ins,const insopt_t *opt,const img_t *img,
                          vostate_t *vo)
{
    /* first epoch to initial */
    if (vo->time.time==0) ins2camera(ins,vo,1);
    if (filts.flag==0) {

        initfilt(ins);
        setfiltstate(ins,vo,0);
        setfiltcov(ins);

        filts.flag=1;
        filts.n=0;
    }
    /* match feature points */
    if (img==NULL||img->feat==NULL||matchfeats(&matchs,img)<=0) {
        trace(2,"no feature points measurement data\n");
        return 0;
    }
    /* update track data */
    if (!updatetrack(opt,img)) return 0;

    /* estimate mono-camera motion */
    vo->status=estmotion(&matchs.mp_bucket,&opt->voopt,vo->dT);
    if (!vo->status) {

        trace(2,"estimate motion fail\n");
        return 0;
    }
    /* update camera states */
    updatecamera(vo,opt,vo->dT,img->time);

    /* filter workspace */
    setfiltstate(ins,vo,1);
    return 1;
}
/* propagate covariance matrix of ins-camera pose in sliding windows----------*/
static int propagate(insstate_t *ins,const insopt_t *opt)
{
    setfiltcov(ins); return 1;
}
/* check estimated states----------------------------------------------------*/
static int chkest_state(const double *x,const insopt_t *opt)
{
    static int iba=xiBa(opt),nba=xnBa(opt);
    static int ibg=xiBg(opt),nbg=xnBg(opt);
    int flag=0;

    /* check estimated states */
    if (     x[  0]==DISFLAG&&norm(x+  0,3)>15.0*D2R) flag|=1;
    if (nba&&x[iba]==DISFLAG&&norm(x+iba,3)>1E5*Mg2M) flag|=1;
    if (nbg&&x[ibg]==DISFLAG&&norm(x+ibg,3)>15.0*D2R) flag|=1;
    if (flag) {
        trace(2,"too large estimated state error\n");
        return 0;
    }
    return 1;
}
/* vo aid filter of attitude--------------------------------------------------*/
static int voflt_att(insstate_t *ins,const insopt_t *opt,const vostate_t *vo)
{
    double Ci1[9],Ci2[9];
    double dCc[9],dCi[9],phi[3],phiz[3],domgda1[9],domgda2[9],Jso3[9],r[3];
    double omg,omgz;
    double *H,*v,*R,*x;
    int i,j,ia,nx=2*ins->nx+1,nv,na,info=0,ns=ins->nx;

    trace(3,"voflt_att:\n");

    camera2ins(filts.Cce[0],filts.rc[0],ins->Cbc,ins->lbc,Ci1,NULL);
    camera2ins(filts.Cce[1],filts.rc[1],ins->Cbc,ins->lbc,Ci2,NULL);

    matmul("TN",3,3,3,0.0,filts.Cbe[0],filts.Cbe[1],0.0,dCi);
    matmul("TN",3,3,3,0.0,Ci2,Ci1,0.0,dCc);

    so3_log(dCi,phi,&omg);
    so3_log(dCc,phiz,&omgz);
    so3_jac(phi,NULL,Jso3);
    if (fabs(omg-omgz)>MAXRES_POSE) {
        trace(2,"large residual\n");
        return 0;
    }
    matmul("NN",3,3,3,1.0,Jso3,filts.Cbe[1],0.0,domgda1);
    matcpy(domgda2,Jso3,3,3);

    v=zeros(1,3); H=zeros(3,nx);
    R=zeros(3,3);
    x=zeros(1,nx);

    ia=xiA(opt); na=xnA(opt);

    for (i=0,nv=0;i<3;i++) {
        if (fabs(v[nv]=phi[i]-phiz[i])>MAXRES_POSE) {
            continue;
        }
        if (H) {
            for (j=ia;j<ia+na;j++) H[j+nv*nx]=-domgda2[i+(j-ia)*3];
            for (j=ia+ns;j<ia+ns+na;j++) {
                H[j+nv*nx]=domgda1[i+(j-ia-nx)*3];
            }
        }
        r[nv++]=SQR(VARPOSE);
    }
    if (v&&nv) {
        trace(3,"v=\n"); tracemat(3,v,3,1,12,6);
    }
    if (H&&nv) {
        trace(3,"H=\n"); tracemat(3,H,nx,nv,12,6);
    }
    if (R&&nv) {
        for (i=0;i<nv;i++) {
            R[i+i*nv]=r[i];
        }
        trace(3,"R=\n");
        tracemat(3,R,nv,nv,12,6);
    }
    if (nv<=0) {
        trace(2,"pose fusion filter fail\n");
        free(x); free(v);
        free(H); free(R);
        return 0;
    }
    if (!filter(x,filts.Px,H,v,R,nx,nv)) {

        if (!chkest_state(x,opt)) goto exit;

        for (i=0;i<ins->nx;i++) {
            for (j=0;j<ins->nx;j++) ins->P[i+j*ins->nx]=filts.Px[i+j*nx];
        }
        clp(ins,opt,x);
        ins->stat=INSS_VO;
    }
    else {
        trace(2,"filter fail\n");
    }
exit:
    free(x); free(v);
    free(H); free(R);
    return info;
}
/* vo aid filter of position--------------------------------------------------*/
static int voflt_pos(insstate_t *ins,const insopt_t *opt,const vostate_t *vo)
{
    double drc[3],dri[3],dr[3],r[3],ri1[3],ri2[3],Cbe[9];
    double W[9],Jp1[9],Jp2[9],Ja[9];
    double *H,*v,*R,*x;
    int i,j,nx=2*ins->nx+1,ip,np,ia,na,nv,info=0;
    int ns=ins->nx;

    trace(3,"voflt_pos:\n");

    camera2ins(filts.Cce[0],filts.rc[0],ins->Cbc,ins->lbc,Cbe ,ri1);
    camera2ins(filts.Cce[1],filts.rc[1],ins->Cbc,ins->lbc,NULL,ri2);

    for (i=0;i<3;i++) dr[i]=ri2[i]-ri1[i];
    matmul("TN",3,1,3,1.0,Cbe,dr,0.0,drc);

    for (i=0;i<3;i++) dr[i]=filts.re[1][i]-filts.re[0][i];
    matmul("TN",3,1,3,1.0,filts.Cbe[0],dr,0.0,dri);

    v=zeros(1,3); H=zeros(3,nx);
    R=zeros(3,3);
    x=zeros(1,nx);

    ip=xiP(opt); np=xnP(opt);
    ia=xiA(opt); na=xnA(opt);

    skewsym3(dr,W);
    matmul("TN",3,3,3,1.0,filts.Cbe[0],W,0.0,Ja);

    matt(filts.Cbe[1],3,3,Jp2);
    matt(filts.Cbe[0],3,3,Jp1);
    for (i=0,nv=0;i<3;i++) {

        if (fabs(v[nv]=dri[i]-drc[i])>MAXRES_POS) continue;
        if (H) {
            for (j=ip+ns;j<ip+ns+np;j++) H[j+nv*nx]=Jp1[i+(j-ip-ns)*3];
            for (j=ia+ns;j<ia+ns+na;j++) H[j+nv*nx]= Ja[i+(j-ia-ns)*3];
            for (j=ip;j<ip+np;j++) H[j+nv*nx]=-Jp2[i+(j-ip)*3];

            H[2*ns+nv*nx]=dri[i];
        }
        r[nv++]=SQR(VARPOS);
    }
    if (v&&nv) {
        trace(3,"v=\n"); tracemat(3,v,3,1,12,6);
    }
    if (H&&nv) {
        trace(3,"H=\n");
        tracemat(3,H,nx,nv,12,6);
    }
    if (R&&nv) {
        for (i=0;i<nv;i++) {
            R[i+i*nv]=r[i];
        }
        trace(3,"R=\n");
        tracemat(3,R,nv,nv,12,6);
    }
    if (nv<=0) {
        trace(2,"pose fusion filter fail\n");
        free(x); free(v);
        free(H); free(R);
        return 0;
    }
    if (!filter(x,filts.Px,H,v,R,nx,nv)) {

        if (!chkest_state(x,opt)) goto exit;

        for (i=0;i<ins->nx;i++) {
            for (j=0;j<ins->nx;j++) ins->P[i+j*ins->nx]=filts.Px[i+j*nx];
        }
        clp(ins,opt,x);
        filts.scale=x[2*ns];

        ins->stat=INSS_VO;
    }
    else {
        trace(2,"filter fail\n");
    }
exit:
    free(x); free(v);
    free(H); free(R);
    return info;
}
/* adjust covariance of filter------------------------------------------------*/
static void adjcovmatrix(insstate_t *ins,const insopt_t *opt,const vostate_t *vo)
{
    filts.flag=0;   
}
/* update.--------------------------------------------------------------------*/
static int updateall(insstate_t *ins,const insopt_t *opt,const img_t *img)
{
    int info=0;

    trace(3,"updateall:\n");

    if (!updatealltrack(ins,opt,img,&ins->vo)) {
        trace(2,"update track data fail\n");
        return 0;
    }
    if (filts.n>=MIN_TRACK_LEN) {

        /*vo aid filter */
        info|=voflt_att(ins,opt,&ins->vo);
        info|=voflt_pos(ins,opt,&ins->vo);

        /* reset filter */
        adjcovmatrix(ins,opt,&ins->vo);
    }
    return info;
}
/* using visual odometry to aid ins/gnss pose estimating-----------------------
 * args:    insopt_t *opt   I   ins options
 *          insstate_t *ins IO  ins states
 *          imud_t *imu     I   imu measurement data
 *          img_t *img      I   image measurement data
 *          int flag        I   update flag
 * return: status (1: ok, 0: fail)
 * ----------------------------------------------------------------------------*/
extern int voigposlc(const insopt_t *opt,insstate_t *ins,const imud_t *imu,
                     const img_t *img,int flag)
{
    trace(3,"voigposlc: time=%s\n",time_str(imu->time,4));

    switch (flag) {
        case 0: return propagate(ins,opt);
        case 1: return propagate(ins,opt)&&updateall(ins,opt,img);
        default: {
            trace(2,"not support mode\n");
        }
    }
    return 0;
}


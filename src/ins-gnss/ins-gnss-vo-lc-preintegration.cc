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
#define MIN_TRACK_LEN           3             /* min length of tracking data */
#define MAX_TRACK_LEN           20            /* max length of tracking data */
#define MAXRES_POSE             30.0*D2R      /* max residual for pose filter (rad) */
#define MAXRES_POS              30.00         /* max residual for position filter (m) */
#define VARPOSE                 1.0*D2R       /* variance of pose measurement (rad) */
#define VARPOS                  0.1           /* variance of position measurement (m) */
#define SWAP(type,x,y)          {type tmp; tmp=x; x=y; y=tmp;}

typedef struct {                              /* store ins states in precious epoch */
    gtime_t time;                             /* time of ins states */
    double Cbe[9],re[3],ve[3];                /* ins attitude/position/velocity */
    double fb[3],omgb[3];                     /* corrected specific-force (b-frame)/angular rate (b-frame) */
    double dt;
} inss_t;

typedef struct {                              /* filter workspace */
    double Cbe[2][9],re[2][3],ve[2][3];       /* ins states in precious and current epoch */
    double Cce[2][9],rc[2][3];                /* camera states in precious and current epoch */
    double *Px;                               /* error states covariance matrix */
    double scale;                             /* mono visual odometry scale */
    inss_t *insdata;                          /* ins states data */
    int nx;                                   /* number of error states */
    int n;                                    /* number of image frames */
    int ni,nimax;                             /* number of ins states in precious */
    int flag;                                 /* flag of filter */
} filt_t;

/* global variables------------------------------------------------------------*/
static track_t  tracks={0};                   /* all tracking feature points data */
static match_t  matchs={0};                   /* match feature points data */
static filt_t   filts ={0};                   /* vo aid filter workspace */

/* SO3 jacobian matrix--------------------------------------------------------*/
static void so3jac(const double *phi,double *Jri)
{
    double I[9]={1,0,0,0,1,0,0,0,1},W[9],W2[9];
    double n=norm(phi,3);
    int i;

    if (n<=1E-8) {seteye(Jri,3); return;}
    skewsym3(phi,W);
    matmul("NN",3,3,3,1.0,W,W,0.0,W2);

    for (i=0;i<9;i++) {
        Jri[i]=I[i]+0.5*W[i]+(1.0/SQR(n)+(1.0+cos(n))/(2.0*n*sin(n)))*W2[i];
    }
}
/* trace ins states ----------------------------------------------------------*/
extern void traceinss(int level, const double *Cbe,const double *re,
                      const double *ve,gtime_t time)
{
    double pos[3],Cne[9],Cnb[9],rpy[3],vel[3];
    char s[64];
    time2str(time,s,3);
    trace(level,"time  =%s\n",s);
    ecef2pos(re,pos);
    ned2xyz(pos,Cne);
    matmul3("TN",Cbe,Cne,Cnb);
    dcm2rpy(Cnb,rpy);
    trace(level,"attn =%8.5f %8.5f %8.5f\n",rpy[0]*R2D,rpy[1]*R2D,rpy[2]*R2D);

    if (ve) {
        matmul3v("T",Cne,ve,vel);
        trace(level,"veln =%8.5f %8.5f %8.5f\n",
              vel[0],vel[1],vel[2]); /* n-frame */

        matmul3v("T",Cbe,ve,vel);
        trace(level,"velb =%8.5f %8.5f %8.5f\n",
              vel[0],vel[1],vel[2]); /* b-frame */
    }
}
/* add ins states-------------------------------------------------------------*/
static int addinstate(const insstate_t *ins,filt_t *filt)
{
    inss_t inss,*p;

    /* copy ins states to buffer data */
    matcpy(inss.fb  ,ins->fb  ,3,1);
    matcpy(inss.omgb,ins->omgb,3,1);

    matcpy(inss.Cbe,ins->Cbe,3,3);
    matcpy(inss.re ,ins->re ,3,1);
    matcpy(inss.ve ,ins->ve ,3,1);

    /* time of ins states */
    inss.time=ins->time;
    inss.dt=ins->dt;

    /* expend buffer data */
    if (filt->nimax<=filt->ni) {
        if (filt->nimax<=0) filt->nimax=16*2; else filt->nimax*=2;
        if (!(p=(inss_t*)realloc(filt->insdata,sizeof(inss_t)*filt->nimax))) {

            /* add ins fail */
            filt->ni=filt->nimax=0;
            free(filt->insdata);
            filt->insdata=NULL;
            return -1;
        }
        filt->insdata=p;
    }
    /* add ins states */
    filt->insdata[filt->ni++]=inss;
    return 1;
}
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
       filts.nx=ins->nx+ins->nx;
       filts.Px=zeros(filts.nx,filts.nx);
       filts.scale=1.0;
   }
}
/* set filter ins states-------------------------------------------------------*/
static void setfiltinsstat(const insstate_t *ins,int pos)
{
    matcpy(filts.Cbe[pos],ins->Cbe,3,3);
    matcpy(filts.re [pos],ins->re ,3,1);
    matcpy(filts.ve [pos],ins->ve ,3,1);
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
/* estimate mono-camera motion-------------------------------------------------*/
static int estmotion(const match_set *pset,const voopt_t *opt,double *dT)
{
    return estmonort(opt,pset,dT,NULL);
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
    int i,j,nx=ins->nx;
    double *T,*Pp;

    if (filts.nx==0||filts.Px==NULL) return 1;

    T =zeros(nx,nx);
    Pp=zeros(nx,nx);
    for (i=0;i<nx;i++) {
        for (j=0;j<nx;j++) Pp[i+j*nx]=filts.Px[i+nx+(j+nx)*filts.nx];
    }
    /* covariance of current ins states */
    matmul33("NNT",ins->F,Pp,ins->F,nx,nx,nx,nx,T);
    for (i=0;i<nx;i++) {
        for (j=0;j<nx;j++) filts.Px[i+j*filts.nx]=T[i+j*nx];
    }
    /* covariance of precious and current ins states */
    matmul("NN",nx,nx,nx,1.0,ins->F,ins->P,0.0,T);
    for (i=0;i<nx;i++) {
        for (j=nx;j<2*nx;j++) {
            filts.Px[i+j*filts.nx]=filts.Px[j+i*filts.nx]=T[i+(j-nx)*nx];
        }
    }
    trace(3,"Px=\n");
    tracemat(3,filts.Px,filts.nx,filts.nx,15,8);

    /* add a ins states */
    addinstate(ins,&filts);

    free(T);
    return 1;
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
/* jabocians of attitude error------------------------------------------------*/
static void jacob_att(const double *z,const double *Ri,const double *Rj,
                      double *Jai,double *Jaj)
{
    double Jr[9];

    so3jac(z,Jr);
    matmul("NT",3,3,3,-1.0,Jr,Rj,0.0,Jai);
    matmul("NT",3,3,3, 1.0,Jr,Rj,0.0,Jaj);
}
/* jacobians of bg------------------------------------------------------------*/
static void jacob_bgk(const double *Rk_1,const double *Rj,const double *omgb,
                      const double dt,double *Jbgk)
{
    double dR[9],phi[3],Jr[9];

    matmul("TN",3,3,3,1.0,Rk_1,Rj,0.0,dR);
    phi[0]=omgb[0]*dt;
    phi[1]=omgb[1]*dt;
    phi[2]=omgb[2]*dt;

    so3jac(phi,Jr);
    matmul("NN",3,3,3,dt,dR,Jr,0.0,Jbgk);
}
static void jacob_bg(const double *z,const double dt,double *Jbg)
{
    inss_t *pins=filts.insdata;
    double Jbgk[9],Jbgi[9]={0},Jr[9],*Rk_1,*Rj,Rz[9];
    int i,j;

    for (Rj=pins[filts.ni-1].Cbe,i=1;i<filts.ni;i++) {
        Rk_1=pins[i].Cbe;
        jacob_bgk(Rk_1,Rj,pins[i-1].omgb,dt,Jbgk);
        for (j=0;j<9;j++) {
            Jbgi[j]+=Jbgk[j];
        }
    }
    so3jac(z,Jr);
    so3_exp(z,Rz);
    matmul33("NTN",Jr,Rz,Jbgi,3,3,3,3,Jbg);
}
/* correction of ins states---------------------------------------------------*/
static void correction(const double *dx,const insopt_t *opt,insstate_t *ins)
{
    int i,iba,ibg;

    iba=xiBa(opt); ibg=xiBg(opt);
    corratt(dx,ins->Cbe);

    /* close-loop velocity and position correction */
    for (i=0;i<3;i++) {
        ins->ve[i]-=dx[xiV(opt)+i];
        ins->re[i]-=dx[xiP(opt)+i];
    }
    /* close-loop accl and gyro bias */
    for (i=0;i<3;i++) {
        ins->ba[i]+=dx[iba+ins->nx+i];
        ins->bg[i]+=dx[ibg+ins->nx+i];
    }
}
/* vo aid filter of attitude--------------------------------------------------*/
static int voflt_att(insstate_t *ins,const insopt_t *opt,const vostate_t *vo)
{
    double Ci1[9],Ci2[9],ri1[3],ri2[3];
    double dCc[9],dCi[9],dC[9],r[3],z[3],Jai[9],Jaj[9],Jbg[9];
    double *H,*v,*R,*x;
    int i,j,nx=filts.nx,nv,ns=ins->nx;
    int na,ia1,ia2,ibg,nbg;

    trace(3,"voflt_att:\n");

    camera2ins(filts.Cce[0],filts.rc[0],ins->Cbc,ins->lbc,Ci1,ri1);
    camera2ins(filts.Cce[1],filts.rc[1],ins->Cbc,ins->lbc,Ci2,ri2);

    traceinss(3,Ci1,ri1,NULL,ins->time);

    matmul("TN",3,3,3,1.0,filts.Cbe[0],filts.Cbe[1],0.0,dCi);
    matmul("TN",3,3,3,1.0,Ci1,Ci2,0.0,dCc);

    matmul("TN",3,3,3,1.0,dCi,dCc,0.0,dC);
    so3_log(dC,z,NULL);

    H=zeros(3,nx); v=zeros(3,1);
    R=zeros(3,3); x=zeros(nx,1);

    ia1=ns+xiA (opt); ia2=xiA (opt); na=xnA(opt);
    ibg=ns+xiBg(opt); nbg=xnBg(opt);

    jacob_att(z,filts.Cbe[0],filts.Cbe[1],Jai,Jaj);
    jacob_bg (z,ins->dt,Jbg);

    for (nv=0,i=0;i<3;i++) {
        if (fabs(v[nv]=z[i])>MAXRES_POSE) continue;
#if 1
        for (j=ia1;j<ia1+na;j++) H[j+nv*nx]=Jai[i+(j-ia1)*3];
        for (j=ia2;j<ia2+na;j++) H[j+nv*nx]=Jaj[i+(j-ia2)*3];
#endif
        for (j=ibg;j<ibg+nbg;j++) {
            H[j+nv*nx]=-Jbg[i+(j-ibg)*3];
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

        trace(3,"dx=\n");
        tracemat(3,x,nx,1,12,6);

        for (i=0;i<ins->nx;i++) {
            for (j=0;j<ins->nx;j++) ins->P[i+j*ins->nx]=filts.Px[i+j*nx];
        }
        correction(x,opt,ins);

        trace(3,"Px=\n");
        tracemat(3,filts.Px,filts.nx,filts.nx,12,6);
    }
    else {
        trace(2,"filter fail\n");
    }
exit:
    free(x); free(v);
    free(H); free(R);
    return 1;
}
/* jacobians of attitude wrt bg-----------------------------------------------*/
static void jacob_dadbg(int p,int q,double *Jbg)
{
    inss_t *pins=filts.insdata;
    double Jbgk[9],*Rk_1,*Rj,dR[9],T[9];
    int i,j;

    setzero(Jbg,3,3);

    for (Rj=pins[q].Cbe,i=p+1;i<q;i++) {
        Rk_1=pins[i].Cbe;
        jacob_bgk(Rk_1,Rj,pins[i-1].omgb,pins[i-1].dt,Jbgk);

        matmul("TN",3,3,3,1.0,Rk_1,Rj,0.0,dR);
        matmul("NN",3,3,3,1.0,dR,Jbgk,0.0,T);
        for (j=0;j<9;j++) Jbg[j]+=T[j];
    }
}
/* jacobians of i-position----------------------------------------------------*/
static void jacob_posi(const double *Cbe,double *Jpi)
{
    int i,j;
    for (i=0;i<3;i++) for (j=0;j<3;j++) Jpi[i+j*3]=-Cbe[j+i*3];
}
/* jacobians of j-position ---------------------------------------------------*/
static void jacob_posj(const double *Cbe,double *Jpj)
{
    int i,j;
    for (i=0;i<3;i++) for (j=0;j<3;j++) Jpj[i+j*3]=Cbe[j+i*3];
}
/* jacobians of i-velocity----------------------------------------------------*/
static void jacob_veli(const double *Cbe,const double dt,double *Jvi)
{
    int i,j;
    for (i=0;i<3;i++) for (j=0;j<3;j++) {
            Jvi[i+j*3]=-Cbe[j+i*3]*dt;
        }
}
/* jacobians of position wrt i-attitude---------------------------------------*/
static void jacob_dpdai(const double *Cbe,const double *ri,const double *rj,
                        const double *vi,const double dt,
                        double *Jdpdai)
{
    double ge[3],dp[3],W[9];
    int i;

    pregrav(ri,ge);
    for (i=0;i<3;i++) dp[i]=rj[i]-ri[i]-vi[i]*dt-0.5*ge[i]*dt*dt;
    skewsym3(dp,W);
    matmul("TN",3,3,3,1.0,Cbe,W,0.0,Jdpdai);
}
/* jacobians of velocity wrt ba-----------------------------------------------*/
static void jacob_dvdba(const double dt,int i,int j,double *Jdvdba)
{
    double *Rk,*Ri,dR[9];
    int k,p;

    setzero(Jdvdba,3,3);

    for (Ri=filts.insdata[i].Cbe,k=i;k<=j;k++) {
        Rk=filts.insdata[k].Cbe;
        matmul("TN",3,3,3,1.0,Ri,Rk,0.0,dR);

        for (p=0;p<9;p++) {
            Jdvdba[p]+=dR[p]*dt;
        }
    }
}
/* jacobians of velocity wrt bg-----------------------------------------------*/
static void jacob_dvdbg(const double dt,int i,int j,double *Jdvdbg)
{
    double *Rk,*Ri,dR[9],W[9],Jdadbg[9];
    double T[9];
    int k,p;

    setzero(Jdvdbg,3,3);

    for (Ri=filts.insdata[i].Cbe,k=i;k<j;k++) {
        Rk=filts.insdata[k].Cbe;
        matmul("TN",3,3,3,1.0,Ri,Rk,0.0,dR);
        skewsym3(filts.insdata[k].fb,W);

        jacob_dadbg(i,k,Jdadbg);

        matmul33("NNN",dR,W,Jdadbg,3,3,3,3,T);
        for (p=0;p<9;p++) {
            Jdvdbg[p]+=T[p]*dt;
        }
    }
}
/* jacobians of position wrt ba-----------------------------------------------*/
static void jacob_dpdba(const double dt,double *Jdpdba)
{
    double Jdvdba[9],dR[9],*Rk,*Ri;
    int i,j;

    setzero(Jdpdba,3,3);

    for (Ri=filts.insdata[0].Cbe,i=0;i<filts.ni-1;i++) {
        Rk=filts.insdata[i].Cbe;
        matmul("TN",3,3,3,1.0,Ri,Rk,0.0,dR);

        jacob_dvdba(dt,0,i,Jdvdba);
        for (j=0;j<9;j++) {
            Jdpdba[j]+=Jdvdba[j]*dt-0.5*dR[j]*dt*dt;
        }
    }
}
/* jacobians of position wrt bg-----------------------------------------------*/
static void jacob_dpdbg(const double dt,double *Jdpdbg)
{
    double Jdvdbg[9],Jdadbg[9],dR[9],*Rk,*Ri;
    double W[9],T[9];
    int k,j;

    setzero(Jdpdbg,3,3);

    for (Ri=filts.insdata[0].Cbe,k=0;k<filts.ni-1;k++) {
        Rk=filts.insdata[k].Cbe;
        matmul("TN",3,3,3,1.0,Ri,Rk,0.0,dR);
        jacob_dvdbg(dt,0,k,Jdvdbg);

        jacob_dadbg(0,k,Jdadbg);

        skewsym3(filts.insdata[k].fb,W);
        matmul33("NNN",dR,W,Jdadbg,3,3,3,3,T);
        for (j=0;j<9;j++) {
            Jdpdbg[j]+=Jdvdbg[j]*dt-0.5*T[j]*dt*dt;
        }
    }
}
/* residual of increment position---------------------------------------------*/
static void pos_residual(const insstate_t *ins,const double dt,double *z)
{
    double Cbei[9],Cbej[9],ri[3],rj[3],dp[3],dpz[3],ge[3],vi[3];
    double dpc[3],dpi[3];
    int i;

    camera2ins(filts.Cce[0],filts.rc[0],ins->Cbc,ins->lbc,Cbei,ri);
    camera2ins(filts.Cce[1],filts.rc[1],ins->Cbc,ins->lbc,Cbej,rj);

    traceinss(3,Cbei,ri,NULL,ins->time);

    vi[0]=filts.ve[0][0];
    vi[1]=filts.ve[0][1];
    vi[2]=filts.ve[0][2];

    pregrav(ri,ge);
    for (i=0;i<3;i++) {
        dpz[i]=rj[i]-ri[i]-vi[i]*dt-0.5*ge[i]*dt*dt;
    }
    matcpy(ri,filts.re[0],1,3);
    matcpy(rj,filts.re[1],1,3);
    matcpy(vi,filts.ve[0],1,3);
    for (i=0;i<3;i++) {
        dp[i]=rj[i]-ri[i]-vi[i]*dt-0.5*ge[i]*dt*dt;
    }
    matmul("TN",3,1,3,1.0,Cbei,dpz,0.0,dpc);
    matmul("TN",3,1,3,1.0,Cbei,dp ,0.0,dpi);
    z[0]=dpi[0]-dpc[0];
    z[1]=dpi[1]-dpc[1];
    z[2]=dpi[2]-dpc[2];
}
/* vo aid filter of position--------------------------------------------------*/
static int voflt_pos(insstate_t *ins,const insopt_t *opt,const vostate_t *vo)
{
    double Jpi[9],Jpj[9],z[3],r[3];
    double Jdpdba[9],Jdpdbg[9];
    double *v,*H,*R,*x;
    int i,j,ipi,ipj,np,ibg,nbg,iba,nba,nv;
    int ns=ins->nx,nx=filts.nx;

    trace(3,"voflt_pos:\n");

    jacob_dpdba(ins->dt,Jdpdba); jacob_dpdbg(ins->dt,Jdpdbg);
    jacob_posi(filts.Cbe[0],Jpi);
    jacob_posj(filts.Cbe[0],Jpj);

    pos_residual(ins,ins->dt,z);

    ipi=ns+xiP (opt); ipj=xiP (opt); np=xnP(opt);
    ibg=ns+xiBg(opt); nbg=xnBg(opt);
    iba=ns+xiBa(opt); nba=xnBa(opt);

    H=zeros(3,nx); v=zeros(3, 1);
    R=zeros(3, 3); x=zeros(nx,1);

    for (nv=0,i=0;i<3;i++) {
        if (fabs(v[nv]=z[i])>MAXRES_POS) continue;

        for (j=iba;j<iba+nba;j++) H[j+nv*nx]=Jdpdba[i+(j-iba)*3];
        for (j=ibg;j<ibg+nbg;j++) H[j+nv*nx]=Jdpdbg[i+(j-ibg)*3];
#if 1
        for (j=ipi;j<ipi+np;j++) H[j+nv*nx]=Jpi[i+(j-ipi)*3];
        for (j=ipj;j<ipj+np;j++) H[j+nv*nx]=Jpj[i+(j-ipj)*3];
#endif
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
        trace(2,"position fusion filter fail\n");
        free(x); free(v);
        free(H); free(R);
        return 0;
    }
    if (!filter(x,filts.Px,H,v,R,nx,nv)) {
        if (!chkest_state(x,opt)) goto exit;

        trace(3,"dx=\n");
        tracemat(3,x,nx,1,12,6);

        for (i=0;i<ins->nx;i++) {
            for (j=0;j<ins->nx;j++) ins->P[i+j*ins->nx]=filts.Px[i+j*nx];
        }
        correction(x,opt,ins);
        
        trace(3,"Px=\n");
        tracemat(3,filts.Px,filts.nx,filts.nx,12,6);
    }
    else {
        trace(2,"filter fail\n");
    }
exit:
    free(x); free(v);
    free(H); free(R);
    return 1;
}
/* initial filter workspace---------------------------------------------------*/
static void initfiltws(insstate_t *ins)
{
    int i,j,ns=ins->nx;
    
    if (filts.nx==0) initfilt(ins);
    setfiltstate(ins,&ins->vo,0);
    for (i=ns;i<ns+ns;i++) {
        for (j=ns;j<ns+ns;j++) filts.Px[i+j*filts.nx]=ins->P[i-ns+(j-ns)*ns];
    }
    trace(3,"filt.Px=\n");
    tracemat(3,filts.Px,filts.nx,filts.nx,12,6);

    filts.flag=1; filts.n=0;
    filts.ni=0;
    addinstate(ins,&filts);
}
/* update.--------------------------------------------------------------------*/
static int updateall(insstate_t *ins,const insopt_t *opt,const img_t *img)
{
    trace(3,"updateall:\n");
    int info=0;

    /* first epoch to initial */
    if (ins->vo.time.time==0) ins2camera(ins,&ins->vo,1);
    if (filts.flag==0) {

        initfiltws(ins);
    }
    if (!updatealltrack(ins,opt,img,&ins->vo)) return 0;
    if (filts.n>=MIN_TRACK_LEN) {

        /* vo aid filter */
        info|=voflt_att(ins,opt,&ins->vo);
        info|=voflt_pos(ins,opt,&ins->vo);
        if (info) {
            ins->stat=INSS_VO;
        }
        /* reset camera states */
        ins2camera(ins,&ins->vo,1);

        /* reset filter */
        filts.flag=0;
    }
    /* re-initial filter */
    if (filts.flag==0) initfiltws(ins);
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
extern int voigposlcpr(const insopt_t *opt,insstate_t *ins,const imud_t *imu,
                       const img_t *img,int flag)
{
    trace(3,"voigposlcpre: time=%s\n",time_str(imu->time,4));

    switch (flag) {
        case 0: return propagate(ins,opt);
        case 1: return propagate(ins,opt)&&updateall(ins,opt,img);
        default: {
            trace(2,"not support mode\n");
        }
    }
    return 0;
}


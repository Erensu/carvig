/*------------------------------------------------------------------------------
 * ins-gnss-vo.cc : ins-gnss-vo coupled common functions
 *
 * reference :
 *    [1] P.D.Groves, Principles of GNSS, Intertial, and Multisensor Integrated
 *        Navigation System, Artech House, 2008
 *    [2] Tedaldi D, Pretto A, Menegatti E. A robust and easy to implement method
 *        for IMU calibration without external equipments,2014.
 *    [3] Skog I, Handel P. Effects of time synchronization errors in GNSS-aided
 *        INS 2008.
 *    [4] Shin E H. Accuracy Improvement of Low Cost INS/GPS for Land Applications
 *    [5] Li M, Mourikis A I. High-precision, consistent EKF-based visual–inertial
 *        odometry[J].International Journal of Robotics Research,2013,32(6):690-711.
 *    [6] Monocular Visual Inertial Odometry on a Mobile Device.
 *    [7] Mourikis A / , Roumeliotis S / . A Multi-State Constraint Kalman Filter
 *        for Vision-aided Inertial Navigation[C]// IEEE International Conference
 *        on Robotics & Automation. IEEE, 2007.
 *    [8] Forster C , Carlone L , Dellaert F , et al. On-Manifold Preintegration
 *        for Real-Time Visual-Inertial Odometry[J]. IEEE Transactions on Robotics,
 *        2015, 33(1):1-21.
 *    [9] Observability-constrained vision-aided inertial navigation.
 *    [10] Li M ,Mourikis A I. Improving the accuracy of EKF-based visual-inertial
 *         odometry[C]// IEEE International Conference on Robotics & Automation.
 *         IEEE, 2012.
 *    [11] Li M ,Mourikis A I. High-precision, consistent EKF-based visual-inertial
 *         odometry[M].Sage Publications, Inc. 2013.
 *
 * version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
 * history : 2018/12/02 1.0 new
 *-----------------------------------------------------------------------------*/
#include "carvig.h"

/* constants-------------------------------------------------------------------*/
#define MIN_TRACK_LEN               2         /* min length of tracking data */
#define MAX_TRACK_LEN               20        /* max length of tracking data */
#define VAR_FEAT                    SQR(3.0)  /* variance of feature point in image coordinates */
#define SWAP(type,x,y)              {type tmp; tmp=x; x=y; y=tmp;}
#define MAXITER                     10        /* max iteration for Gauss Newton optimization */
#define MAX_FEAT_RES                5.0       /* max residual for a feature point measurement */
#define MAX_GNCOST_NORM             0.01      /* set to inf to allow any triangulation, no matter how bad */
#define OUT_DETECT                  1         /* Mahalanobis gating test for the residual */
#define USE_BUCKET_FEATS            1         /* use bucket matched feature points data */

/* type definitions ----------------------------------------------------------*/
typedef struct hashtable {                    /* hash table type */
    int id;                                   /* id of track feature point */
    const trackd_t *ptrk;                     /* pointer of this track feature point */
    UT_hash_handle hh;                        /* makes this structure hashable */
} hashtable_t;

typedef struct cams {                         /* camera pose states type */
    gtime_t time;                             /* timestamp of ins states */
    long int find;                            /* frame index corresponding to camera */
    unsigned char flag;                       /* flag of camera pose (0:initial, 1: updated, 2:bad) */
    double re[3],ve[3],Cce[9];                /* camera position/velocity/attitude states in ecef */
    hashtable *trackfeat;                     /* tracking feature points */
} cams_t;

typedef struct vofilt {                       /* vo filter workspace type */
    int n,nmax;                               /* number and max number of history ins states */
    int nx;                                   /* number error states of current sliding windows */
    double *Px;                               /* current ins states covariance included camera pose in sliding window */
    cams_t *data;                             /* all tracked camera states */
} vofilt_t;

/* global variables------------------------------------------------------------*/
static vofilt_t vofilt={0};                   /* vo aid filter workspace */
static track_t  tracks={0};                   /* all tracking feature points data */
static match_t  matchs={0};                   /* match feature points data */
static const insstate_t *inss=NULL;           /* pointer to ins states */
static const insopt_t   *opts=NULL;           /* pointer to ins options */

static int (*featpos_funcptr)(const cams_t *cam1,const cams_t *cam2,const insstate_t *ins,
                              const double *uv1,
                              const double *uv2,
                              const double *K,
                              double *pf);    /* function pointer of feature position compute function */

/* add a track feature to hash table------------------------------------------*/
static void hash_add_feature(struct hashtable **ht,const int id,const trackd_t *ptk)
{
    struct hashtable *s=NULL;
    HASH_FIND_INT(*ht,&id,s);  /* id already in the hash? */
    if (s==NULL) {

        s=(struct hashtable *)malloc(sizeof *s);
        s->id=id;
        s->ptrk=ptk;
        HASH_ADD_INT(*ht,id,s);  /* id: name of key field */
    }
}
/* find a track feature from hash table---------------------------------------*/
static hashtable* hash_find_feature(const struct hashtable *ht,const int id)
{
    struct hashtable *s=NULL; 

    HASH_FIND_INT(ht,&id,s);  /* s: output pointer */
    return s;
}
/* remove element from hash table---------------------------------------------*/
static void hash_delete_feature(struct hashtable **ht,const int id)
{
    struct hashtable *s=NULL;
    if (!(s=hash_find_feature(*ht,id))) {
        trace(2,"no feature to delete\n");
        return;
    }
    HASH_DEL(*ht,s);  /* user: pointer to deletee */
    free(s);
}
/* delete hash table----------------------------------------------------------*/
static void hash_destroy(struct hashtable **ht)
{
    struct hashtable *current,*tmp;

    HASH_ITER(hh,*ht,current,tmp) {
        HASH_DEL(*ht,current);  /* delete; users advances to next */
        free(current);          /* optional- if you want to free  */
    }
}
/* create a hash table--------------------------------------------------------*/
static struct hashtable *hashtable()
{
    return NULL;
}
/* counts of elements in hash table-------------------------------------------*/
static int hash_counts(const struct hashtable *ht)
{
    return HASH_COUNT(ht);
}
/* initial vo-aid-------------------------------------------------------------*/
extern void initvoaid(insopt_t *opt)
{
    trace(3,"initvoaid:\n");

    opt->voopt.match.f =opt->voopt.calib.f;
    opt->voopt.match.fu=opt->voopt.calib.fu;
    opt->voopt.match.fv=opt->voopt.calib.fv;

    opt->voopt.match.cu=opt->voopt.calib.cu;
    opt->voopt.match.cv=opt->voopt.calib.cv;

    init_match(&matchs,&opt->voopt.match);
    inittrackimgbuf(&opt->voopt);

    vofilt.nx=xnX(opt);
    vofilt.Px=zeros(vofilt.nx,vofilt.nx); vofilt.n=vofilt.nmax=0;
    return;
}
/* free vo-aid----------------------------------------------------------------*/
extern void freevoaid()
{
    trace(3,"freevoaid:\n");
    int i;

    freetrackset(&tracks);
    free_match(&matchs);
    if (vofilt.Px) free(vofilt.Px);

    for (i=0;i<vofilt.n;i++) {
        hash_destroy(&vofilt.data[i].trackfeat);
    }
    if (vofilt.data) {
        free(vofilt.data);
    }
    freetrackimgbuf(); return;
}
/* resize matrix--------------------------------------------------------------
 * args:    double *A  IO  resized matrix
 *          int n,m    I   size of matrix before resize
 *          int p,q    I   size of matrix after resize
 * return : 1 (ok) or 0 (fail)
 * ---------------------------------------------------------------------------*/
extern int resize(double **A,int m,int n,int p,int q)
{
    trace(3,"resize:\n");
    double *Ap=zeros(p,q);

    matcpy(Ap,*A,m,n);
    free(*A); *A=Ap;
    return 1;
}
/* get camera pose using current ins stats------------------------------------*/
static void campose(const insstate_t *ins,const insopt_t *opt,cams_t *cams,
                    double *Jp,double *Ja,double *Jv)
{
    double T[3],we[3]={0,0,-OMGE},web[3];
    double W[9],W1[9],W2[9];
    int i,j;
    static double I[9]={1,0,0,0,1,0,0,0,1};

    static int ip=xiP(opt),ia=xiA(opt),ilc=xiCl (opt);
    static int np=xnP(opt),na=xnA(opt),nlc=xnCla(opt);
    static int iv=xiV(opt),nv=xnV(opt);
    static int ibg=xiBg(opt),nbg=xnBg(opt);
    
    /* compute camera pose and velocity */
    if (cams) {
        matmul("NT",3,3,3,1.0,ins->Cbe,ins->Cbc,0.0,cams->Cce);
        matmul("NN",3,1,3,1.0,ins->Cbe,ins->lbc,0.0,T);
        for (i=0;i<3;i++) {
            cams->re[i]=ins->re[i]+T[i];
        }
        matmul("NN",3,1,3,1.0,ins->Cbe,we,0.0,T);
        for (i=0;i<3;i++) {
            T[i]=T[i]+ins->omgb[i];
        }
        skewsym3(T,W);
        matmul("NN",3,1,3,1.0,W,ins->lbc,0.0,T);
        matmul("NN",3,1,3,1.0,ins->Cbe,T,0.0,cams->ve);
        for (i=0;i<3;i++) {
            cams->ve[i]=cams->ve[i]+ins->ve[i];
        }
        cams->time=ins->time;
        cams->flag=0; /* initial flag */
    }
    /* jacobian for camera position */
    if (Jp) {

        matmul("NN",3,3,3,1.0,ins->Cbe,ins->lbc,0.0,T);
        skewsym3(T,W);

        /* jacobians wrt. attitude */
        for (i=0;i<3;i++) {
            for (j=ia;j<na+ia;j++) Jp[i+j*3]=W[i+(j-ia)*3];
        }
        /* jacobians wrt. ins position */
        for (i=0;i<3;i++) {
            for (j=ip;j<np+ip;j++) Jp[i+j*3]=I[i+(j-ip)*3];
        }
        /* jacobians wrt. camera lever arm  */
        if (opt->estcaml) {
            for (i=0;i<3;i++) for (j=ilc;j<ilc+nlc;j++) Jp[i+j*3]=ins->Cbe[i+(j-ilc)*3];
        }
        trace(3,"Jp=\n");
        tracemat(3,Jp,3,ins->nx,12,6);
    }
    /* jacobian for camera attitude */
    if (Ja) {
        for (i=0;i<3;i++) for (j=ia;j<ia+na;j++) {
                Ja[i+j*3]=I[i+(j-ia)*3];
            }
        trace(3,"Ja=\n");
        tracemat(3,Ja,3,ins->nx,12,6);
    }
    /* jacobian for camera velocity */
    if (Jv) {
        matmul("TN",3,1,3,1.0,ins->Cbe,we,0.0,T);
        for (i=0;i<3;i++) {
            web[i]=ins->omgb[i]+T[i];
        }
        skewsym3(web,W);
        matmul33("NNN",ins->Cbe,W,ins->lbc,3,3,3,1,T);
        skewsym3(T,W);

        matmul("NN",3,1,3,1.0,ins->Cbe,ins->lbc,0.0,T);
        skewsym3(T,W1);
        matmul("NN",3,3,3,1.0,W1,Omge,0.0,W2);

        /* wrt. attitude error state */
        for (i=0;i<3;i++) {
            for (j=ia;j<ia+na;j++) Jv[i+j*3]=W[i+(j-ia)*3]+W2[i+(j-ia)*3];
        }
        /* wrt. velocity */
        for (i=0;i<3;i++) {
            for (j=iv;j<iv+nv;j++) Jv[i+j*3]=I[i+(j-iv)*3];
        }
        /* wrt. bg */
        skewsym3(ins->lbc,W);
        matmul("NN",3,3,3,-1.0,ins->Cbe,W,0.0,W1);
        for (i=0;i<3;i++) {
            for (j=ibg;j<ibg+nbg;j++) Jv[i+j*3]=W1[i+(j-ibg)*3];
        }
        /* wrt. camera lever arm */
        if (opt->estcaml) {
            skewsym3(web,W);
            matmul("NN",3,3,3,1.0,ins->Cbe,W,0.0,W1);
            for (i=0;i<3;i++) {
                for (j=ilc;j<ilc+nlc;j++) Jv[i+j*3]=W1[i+(j-ilc)*3];
            }
        }
        trace(3,"Jv=\n");
        tracemat(3,Jv,3,ins->nx,12,6);
    }
}
/* add a ins state to filter--------------------------------------------------*/
static int addinstat(vofilt_t *filt,const insstate_t *in,const insopt_t *opt,
                     int cur_fid)
{
    cams_t cams,*p;
    trace(3,"addinstat: time=%s\n",time_str(in->time,4));

    /* compute camera pose */
    campose(in,opt,&cams,NULL,NULL,NULL);

    /* create hash table */
    cams.trackfeat=hashtable();
    cams.find=cur_fid;

    /* add ins states */
    if (filt->nmax<=filt->n) {
        if (filt->nmax<=0) filt->nmax=16*2; else filt->nmax*=2;
        if (!(p=(cams_t*)realloc(filt->data,sizeof(cams_t)*filt->nmax))) {

            /* add ins fail */
            filt->n=filt->nmax=0; free(filt->data);
            filt->data=NULL;
            return -1;
        }
        filt->data=p;
        int nxc=filt->nmax*9+in->nx;
        if (!resize(&filt->Px,filt->nx,filt->nx,nxc,nxc)) {
            trace(2,"Py resize matrix fail\n");

            free(filt->Px); filt->Px=NULL;
            return -1;
        }
    }
    filt->data[filt->n++]=cams;
    return 1;
}
/* propagate covariance matrix of ins-camera pose in sliding windows----------*/
static int propagate(vofilt_t *filt,const insstate_t *ins)
{
    int i,j,nx=filt->nx;
    double *T,*Pp;

    trace(3,"propagate: nx=%d\n",filt->nx);

    if (filt->n==0||filt->nx==ins->nx) {
        for (i=0;i<nx;i++) for (j=0;j<nx;j++) filt->Px[i+j*nx]=ins->P[i+j*nx];
        return 1;
    }
    Pp=mat(nx,nx); T=mat(nx,nx);
    for (i=0;i<nx-ins->nx;i++) {
        for (j=0;j<ins->nx;j++) T[i+j*(nx-ins->nx)]=filt->Px[ins->nx+i+j*nx];
    }
    trace(3,"P(0)=\n");
    tracemat(3,filt->Px,filt->nx,filt->nx,12,8);

    matmul("NT",nx-ins->nx,ins->nx,ins->nx,1.0,T,ins->F,0.0,Pp);
    for (i=0;i<nx-ins->nx;i++) {
        for (j=0;j<ins->nx;j++) {
            filt->Px[(ins->nx+i)*nx+j]=filt->Px[ins->nx+i+j*nx]=Pp[i+j*(nx-ins->nx)];
        }
    }
    for (i=0;i<ins->nx;i++) {
        for (j=0;j<ins->nx;j++) {
            filt->Px[i+j*nx]=ins->P[i+j*ins->nx];
        }
    }
    trace(3,"P(1)=\n");
    tracemat(3,filt->Px,filt->nx,filt->nx,12,8);

    free(T); free(Pp);
    return 1;
}
/* add a ins state to vo-filter workspace-------------------------------------*/
static void augstates(vofilt_t *filt,const insstate_t *ins,const insopt_t *opt,
                      int cur_fid)
{
    double *J,*Jp,*Jv,*Ja,*P,*I;
    int i,j,nx=ins->nx;
    int nxp=filt->nx,nxc=filt->nx+9;

    trace(3,"augstates:\n");

    /* first camera state */
    if (filt->n==0||filt->nx==ins->nx) {
        for (i=0;i<nx;i++) for (j=0;j<nx;j++) filt->Px[i+j*nx]=ins->P[i+j*nx];
    }
    /* add a camera pose to sliding window */
    if (addinstat(filt,ins,opt,cur_fid)<0) {
        trace(2,"add camera pose fail\n");
        return;
    }
    Jp=zeros(3,nx); Jv=zeros(3,nx); Ja=zeros(3,nx); P=zeros(nxc,nxc);
    J=zeros(nxc,nxp); I=eye(nxp);

    /* jacobians */
    campose(ins,opt,NULL,Jp,Ja,Jv);

    asi_blk_mat(J,nxc,nxp,Ja,3,nx,nxp+0,0);
    asi_blk_mat(J,nxc,nxp,Jv,3,nx,nxp+3,0);
    asi_blk_mat(J,nxc,nxp,Jp,3,nx,nxp+6,0);
    asi_blk_mat(J,nxc,nxp,I,nxp,nxp,0,0);

    matmul33("NNT",J,filt->Px,J,nxc,nxp,nxp,nxc,P);

    trace(3,"J=\n");
    tracemat(3,J,nxc,nxp,12,6);

    trace(3,"P(0)=\n");
    tracemat(3,filt->Px,nxp,nxp,12,6);

    trace(3,"P(1)=\n");
    tracemat(3,P,nxc,nxc,12,6);
    
    /* augment covariance matrix */
    for (i=0;i<nxc;i++) {
        for (j=0;j<nxc;j++) filt->Px[i+j*nxc]=P[i+j*nxc];
    }
    /* update nx of `Px' */
    filt->nx=nxc;

    free(J); free(Jp); free(Ja); free(Jv);
    free(P); free(I);
}
/* handle of remove camera pose-----------------------------------------------*/
static int handlerm(vofilt_t *filt,const insopt_t *opt,int idx)
{
    int i,j,*index,k=0;

    index=imat(filt->nx,1);

    for (i=idx;i<filt->n-1;i++) filt->data[i]=filt->data[i+1];
    for (i=0;i<filt->n;i++) {
        if (idx==i) continue; index[k++]=i;
    }
    for (i=0;i<k;i++) {
        for (j=0;j<k;j++) {
            filt->Px[i+j*k]=filt->Px[index[i]+index[j]*filt->n];
        }
    }
    free(index);
    filt->n--; return k;
}
/* remove camera pose from sliding window-------------------------------------*/
static int rmcamera(vofilt_t *filt,const insopt_t *opt)
{
    int i;
    trace(3,"rmcamera:\n");

    for (i=0;i<filt->n;i++) {
        if (hash_counts(filt->data[i].trackfeat)||filt->data[i].flag!=2) continue;
        hash_destroy(&filt->data[i].trackfeat);

        /* remove camera */
        if (handlerm(filt,opt,i)!=filt->n) continue;
    }
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
/* add current tracking feature points to current camera----------------------*/
static int addcurrenttrack(cams_t *cam)
{
    int i,j;

    trace(3,"add new tracks: %d\n",tracks.nnew);
    trace(3,"update tracks: %d\n",tracks.nupd);

    /* add new track to camera tracking-feature list */
    for (i=0;i<tracks.nnew;i++) {
        j=tracks.newtrack[i];
        hash_add_feature(&cam->trackfeat,tracks.data[j].uid,&tracks.data[j]);
    }
    /* update old track */
    for (i=0;i<tracks.nupd;i++) {
        j=tracks.updtrack[i];
        hash_add_feature(&cam->trackfeat,tracks.data[j].uid,
                         &tracks.data[j]);
    }
    return tracks.nnew||tracks.nupd;
}
/* get camera pointer through given `uid'-------------------------------------*/
static cams_t* getcampointor(int uid)
{
    cams_t *pcam=NULL;
    int i;
    for (i=0;i<vofilt.n;i++) {
        if (vofilt.data[i].find==uid) {pcam=&vofilt.data[i]; break;}
    }
    return pcam;
}
/* update current camera-frame track data-------------------------------------*/
static int updatecamtrack(const insopt_t *opt,gtime_t time,int pre_id,int cur_ind)
{
    cams_t *cam;
    trace(3,"updatecamtrack: \n");

    if ((cam=getcampointor(pre_id))&&hash_counts(cam->trackfeat)==0) {

        /* if no matched feature of precious camera
         * then we fill in with current tracks
         * */
        addcurrenttrack(cam);
    }
    /* add current track features */
    if ((cam=getcampointor(cur_ind))&&addcurrenttrack(cam)) return 1;
    return 0;
}
/* distort feature point------------------------------------------------------*/
static void distortfeat(const insstate_t *ins,const double *pfn,double *pfnd)
{
    cam_t camp={0};
    camp.k1=ins->k1; camp.k2=ins->k2;
    camp.p1=ins->p1;
    camp.p2=ins->p2;
    distortradtan(&camp,pfn,pfnd,NULL);
}
/* calculates measurement variance for a feature point------------------------*/
static void featureR(const trackd *ftrack,const feature *feat,const cams_t *cam,
                     const insstate_t *ins,double *R)
{
    trace(3,"featureR:\n");
#if 0
    /* other method to determinate feature variance */
#else
    /* default feature variance */
    R[0]=VAR_FEAT/SQR(ins->fx);
    R[3]=VAR_FEAT/SQR(ins->fy);
#endif
}
/* calculates measurement residual for a feature point------------------------*/
static void featurev(const trackd *ftrack,const feature *feat,const cams_t *cam,
                     const insstate_t *ins,double *v)
{
    double pf[3],pfn[3],pfnd[3];
    double dp[3],uv[2];
    int i;

    trace(3,"featurev:\n");

    /* feature position in c-frame */
    for (i=0;i<3;i++) dp[i]=ftrack->xyz[i]-cam->re[i];
    matmul3v("T",cam->Cce,dp,pf);

    /* distorts feature point */
    for (i=0;i<3;i++) pfn[i]=pf[i]/pf[2];
    distortfeat(ins,pfn,pfnd);

    /* project into image */
    uv[0]=ins->fx*pfnd[0]+ins->ox;
    uv[1]=ins->fy*pfnd[1]+ins->oy;

    if (v) {
        v[0]=uv[0]-feat->u; v[1]=uv[1]-feat->v;
        trace(3,"v=\n");
        tracemat(3,v,1,2,12,6);
    }
}
/* calculates measurement jacobians matrix for a feature point----------------*/
static void featureH(const trackd_t *trackf,const cams_t *cam,
                     const insstate_t *ins,double *Hf,double *Hxc,
                     double *Hfo,
                     double *Hkp)
{
    double pf[3],dp[3],pfn[3];
    double J1[9],J2[6]={0},J3[4]={0},J4[9],J5[9],J6[4],W[9],pfnd[3];
    double Jac[9],Jpc[9];
    double T[9],r2;
    int i,j;

    trace(3,"featureH:\n");

    /* feature position in c-frame */
    for (i=0;i<3;i++) dp[i]=trackf->xyz[i]-cam->re[i];
    matmul3v("T",cam->Cce,dp,pf);

    /* distorts */
    for (i=0;i<3;i++) pfn[i]=pf[i]/pf[2];
    distortfeat(ins,pfn,pfnd);

    /* Hf: 2*3 for a feature point */
    matcpy(J1,cam->Cce,3,3);

    J2[0]=1.0/SQR(pf[2]); J2[3]=1.0/SQR(pf[2]);
    J2[4]=-pf[0];
    J2[5]=-pf[1];

    J3[0]=ins->fx; J3[3]=ins->fy;
    if (Hf) {
        matmul33("NNN",J3,J6,J2,2,2,2,3,W);
        matmul("NN",2,3,3,1.0,W,J1,0.0,Hf);

        trace(3,"Hf=\n"); tracemat(3,Hf,2,3,12,6);
    }
    /* Hxc: 2*9 for a camera pose */
    skewsym3(dp,W);
    matmul("TN",3,3,3,1.0,cam->Cce,W,0.0,J4);
    matmul33("NNN",J3,J6,J2,2,3,3,3,T);
    matmul("NN",2,3,3,1.0,T,J4,0.0,Jac);

    matt(cam->Cce,3,3,J5); for (i=0;i<9;i++) J5[i]=-J5[i];
    matmul33("NNN",J3,J2,J5,2,3,3,3,Jpc);
    if (Hxc) {
        for (i=0;i<2;i++) {
            for (j=0;j<3;j++) Hxc[i+j*9]=Jac[i+(j-0)*3];
            for (j=6;j<9;j++) Hxc[i+j*9]=Jpc[i+(j-6)*3];
        }
        trace(3,"Hxc=\n"); tracemat(3,Hxc,2,9,12,6);
    }
    /* Hfo: 2*4 for a camera */
    if (Hfo) {
        Hfo[0]=pfnd[0]; Hfo[4]=1.0;
        Hfo[3]=pfnd[1]; Hfo[7]=1.0;

        trace(3,"Hfo=\n"); tracemat(3,Hfo,2,4,12,6);
    }
    /* Hkp: 2*4 for a camera */
    if (Hkp) {
        r2=SQR(pfn[0])+SQR(pfn[1]);

        Hkp[0]=pfn[0]*r2; Hkp[2]=pfn[0]*SQR(r2);
        Hkp[4]=2.0*pfn[0]*pfn[1];
        Hkp[6]=r2+2.0*SQR(pfn[0]);

        Hkp[1]=pfn[1]*r2; Hkp[3]=pfn[1]*SQR(r2);
        Hkp[5]=r2+2.0*SQR(pfn[1]);
        Hkp[7]=2.0*pfn[0]*pfn[1];

        trace(3,"Hkp=\n"); tracemat(3,Hkp,2,4,12,6);
    }
}
/* find all camera poses where feature points are observed--------------------*/
static int findallcamera(const trackd_t *feat,const insopt_t *opt,gtime_t time,
                         int *index)
{
    int i,k,j;
    for (i=0,k=0;i<vofilt.n;i++) {
        if (!hash_find_feature(vofilt.data[i].trackfeat,feat->uid)) continue;
        index[k++]=i;
    }
    for (i=0;i<feat->n;i++) for (j=0;j<k;j++) {
        if (fabs(timediff(feat->data[i].time,vofilt.data[index[j]].time))<1E-6) {
            break;
        }
        SWAP(int,index[i],index[j]);
    }
    return k;
}
/* get camera calibration matrix----------------------------------------------*/
static void getK(const insstate_t *ins,double *K)
{
    /* camera calibration matrix */
    K[0]=ins->fx; K[4]=ins->fy;
    K[6]=ins->ox; K[7]=ins->oy;
    K[8]=1.0;
}
/* rotation parameters and translation parameters convert to transformation---
 * arg   : double *R  I  rotation matrix
 *         double *t  I  translation matrix
 *         double *T  O  transformation matrix
 * return: none
 * --------------------------------------------------------------------------*/
extern void rt2tf(const double *R,const double *t,double *T)
{
    setzero(T,4,4);
    T[0]=R[0]; T[4]=R[3]; T[ 8]=R[6]; T[12]=t[0];
    T[1]=R[1]; T[5]=R[4]; T[ 9]=R[7]; T[13]=t[1];
    T[2]=R[2]; T[6]=R[5]; T[10]=R[8]; T[14]=t[2]; T[15]=1.0;
}
/* transform matrix convert to rotation and translation parameters-----------*/
extern void tf2rt(const double *T,double *R,double *t)
{
    if (R) {
        seteye(R,3);
        R[0]=T[0]; R[3]=T[4]; R[6]=T[8 ];
        R[1]=T[1]; R[4]=T[5]; R[7]=T[9 ];
        R[2]=T[2]; R[5]=T[6]; R[8]=T[10];
    }
    if (t) {
        setzero(t,1,3);
        t[0]=T[12]; t[1]=T[13]; t[2]=T[14];
    }
}
/* compute relative transformation between two frame -------------------------*/
static void reltf(gtime_t ts,gtime_t te,const cams_t *cam1,const cams_t *cam2,
                  double *C,double *t)
{

}
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
/* bearing vector of feature point--------------------------------------------*/
static int bearingvector(const cams_t *cam,const insstate_t *ins,const double *uv,
                         const double *K,double *bv)
{
    double img_p[2],uimg_p[2],n; int i;
    cam_t camp;

    trace(2,"bearingvector:\n");

    /* camera calibration parameters */
    camp.k1=ins->k1; camp.k2=ins->k2;
    camp.p1=ins->p1; camp.p2=ins->p2;

    /* unscale and center */
    img_p[0]=(uv[0]-K[6])/K[0];
    img_p[1]=(uv[1]-K[7])/K[4];

    /* undistort */
    if (!undistortradtan(&camp,img_p,uimg_p,NULL)) {

        trace(2,"undistort fail\n");
        return 0;
    }
    /* project 1 into z direction and normalization*/
    bv[0]=uimg_p[0];
    bv[1]=uimg_p[1];
    bv[2]=1.0;
    return 1;
}
/* initial feature point position using bearing vector------------------------*/
static int initfeatpos_bearvetcor(const cams_t *cam1,const cams_t *cam2,const insstate_t *ins,
                                  const double *uv1,const double *uv2,
                                  const double *K,double *pf)
{
    double T1[16],T2[16],dT[16],R[9],t[3],A[16];
    double bv1[3],bv2[3],P1[12],P2[12];
    int i;

    trace(3,"initfeatpos_bearvetcor:\n");

    rt2tf(cam1->Cce,cam1->re,T1); rt2tf(cam2->Cce,cam2->re,T2);
    if (!matinv(T1,4)) {
        matmul("NN",4,4,4,1.0,T1,T2,0.0,dT);
    }
    else {
        seteye(dT,4);
    }
    matinv(dT,4); tf2rt(dT,R,t);

    /* bearing vector: cam1->uv1, cam2->uv2 */
    bearingvector(cam1,ins,uv1,K,bv1);
    bearingvector(cam2,ins,uv2,K,bv2);

    /* camera project matrix */
    prjmatrix(R,t,K,P2); seteye(R,3); setzero(t,1,3);
    prjmatrix(R,t,K,P1);

    for (i=0;i<4;i++) {
        A[0+4*i]=P1[2+3*i]*bv1[0]-P1[0+3*i]*bv1[2];
        A[1+4*i]=P1[2+3*i]*bv1[1]-P1[1+3*i]*bv1[2];
        A[2+4*i]=P2[2+3*i]*bv2[0]-P2[0+3*i]*bv2[2];
        A[3+4*i]=P2[2+3*i]*bv2[1]-P2[1+3*i]*bv2[2];
    }
    if (!svd(A,4,4,NULL,NULL,dT)) {
        trace(2,"Jacobians svd fail\n");
        return 0;
    }
    /* return false if this point is at infinity */
    if (fabs(dT[3+3*4])<1E-10) {
        trace(2,"feature point is at infinity\n");
        return 0;
    }
    for (i=0;i<3;i++) pf[i]=dT[i+3*4]/dT[3+3*4];
    trace(3,"feature position: %8.4lf  %8.4lf  %8.4lf\n",pf[0],pf[1],pf[2]);
    return 1;
}
/* initial feature point position using fast non-linear approximation---------*/
static int initfeatpos_nolinear(const cams_t *cam1,const cams_t *cam2,const insstate_t *ins,
                                const double *uv1,const double *uv2,
                                const double *K,double *pf)
{
    double A[4],b[2],R[9],t[3],b1[3],b2[3],bt[3],lambda[2];
    double T1[16],T2[16],dT[16],xm[3],xn[3];
    int i;

    trace(3,"initfeatpos_nolinear:\n");

    /* transformation */
    rt2tf(cam1->Cce,cam1->re,T1); rt2tf(cam2->Cce,cam2->re,T2);
    if (!matinv(T1,4)) {
        matmul("NN",4,4,4,1.0,T1,T2,0.0,dT);
    }
    else seteye(dT,4); tf2rt(dT,R,t);

    /* bearing vector: cam1->uv1, cam2->uv2 */
    bearingvector(cam1,ins,uv1,K,b1);
    bearingvector(cam2,ins,uv2,K,b2);

    matmul3("N",R,b2,bt);

    b[0]=dot(t,b1,3);
    b[1]=dot(t,bt,3);
    
    A[0]= dot(b1,b1,3); A[1]=dot(b1,bt,3); A[2]=-A[1];
    A[3]=-dot(bt,bt,3);

    if (matinv(A,2)) return 0;
    matmul("NN",2,1,2,1.0,A,b,0.0,lambda);

    for (i=0;i<3;i++) {
        xm[i]=lambda[0]*b1[i]; xn[i]=t[i]+lambda[1]*bt[i];
        pf[i]=(xm[i]+xn[i])/2.0;
    }
    trace(3,"feature position: %8.4lf  %8.4lf  %8.4lf\n",pf[0],pf[1],pf[2]);
    return 1;
}
/* initial feature point position using chieral method------------------------*/
static int initfeatpos_chieral(const cams_t *cam1,const cams_t *cam2,const insstate_t *ins,
                               const double *uv1,const double *uv2,
                               const double *K,double *pf)
{
    double C[9],t[3],P1[12],P2[12],T1[16],T2[16],dT[16];
    double J[16],V[16],x1[3],x2[3],X[4];
    int j;

    trace(3,"initfeatpos_chieral:\n");

    /* relative camera transformation */
    rt2tf(cam1->Cce,cam1->re,T1); rt2tf(cam2->Cce,cam2->re,T2);
    if (!matinv(T1,4)) {
        matmul("NN",4,4,4,1.0,T1,T2,0.0,dT);
    }
    else {
        seteye(dT,4);
    }
    tracemat(3,dT,4,4,12,6);
    matinv(dT,4); tf2rt(dT,C,t);

    /* camera project matrix */
    prjmatrix(C,t,K,P2); seteye(C,3); setzero(t,1,3);
    prjmatrix(C,t,K,P1);
    for (j=0;j<4;j++) {
        J[0+4*j]=P1[2+3*j]*uv1[0]-P1[0+3*j];
        J[1+4*j]=P1[2+3*j]*uv1[1]-P1[1+3*j];
        J[2+4*j]=P2[2+3*j]*uv2[0]-P2[0+3*j];
        J[3+4*j]=P2[2+3*j]*uv2[1]-P2[1+3*j];
    }
    if (!svd(J,4,4,NULL,NULL,V)) return 0;

    /* return false if this point is at infinity */
    if (fabs(V[3+3*4])<1E-10) {
        trace(2,"feature point is at infinity\n");
        return 0;
    }
    for (j=0;j<4;j++) X[j]=V[j+3*4];
    matmul("NN",3,1,4,1.0,P1,X,0.0,x1);
    matmul("NN",3,1,4,1.0,P2,X,0.0,x2);

    if (x1[2]*X[3]<0.0||x2[2]*X[3]<0.0) return 0;
    for (j=0;j<3;j++) {
        pf[j]=V[j+3*4]/V[3+3*4];
    }
    trace(3,"feature position: %8.4lf  %8.4lf  %8.4lf\n",pf[0],pf[1],pf[2]);
    return 1;
}
/* initial feature point position --------------------------------------------*/
static int initfeatpos(const cams_t *cam1,const cams_t *cam2,const insstate_t *ins,
                       const double *uv1,
                       const double *uv2,const double *K,double *pf,
                       int method)
{
    trace(2,"initfeatpos: method=%d\n",method);

    switch (method) {
        case 0: featpos_funcptr=initfeatpos_chieral;    break;
        case 1: featpos_funcptr=initfeatpos_nolinear;   break;
        case 2: featpos_funcptr=initfeatpos_bearvetcor; break;
        default:
            featpos_funcptr=initfeatpos_chieral; break;
    }
    return featpos_funcptr(cam1,cam2,ins,uv1,uv2,K,pf);
}
/* estimate feature point position in ecef------------------------------------*/
static int estfeatpos(trackd_t *feat,const insopt_t *opt,const insstate_t *ins,
                      gtime_t time,int *index,double *pf)
{
    static double uv1[2],uv2[2];
    double K[9]={0};
    cams_t *cam1,*cam2;
    int k=0;

    trace(3,"estfeatpos:\n");

    if (!(k=findallcamera(feat,opt,time,index))||(k!=feat->n)) {
        trace(2,"no found camera observed this feature\n");
        return 0;
    }
    if (k<MIN_TRACK_LEN) return 0;

    /* camera pointer */
    cam1=&vofilt.data[index[  0]];
    cam2=&vofilt.data[index[k-1]];
    
    /* feature point measurement data */
    uv1[0]=feat->data[        0].u; uv1[1]=feat->data[0].v;
    uv2[0]=feat->data[feat->n-1].u;
    uv2[1]=feat->data[feat->n-1].v;
    getK(ins,K);

    /* draw track feature for debug */
    drawtrackd(feat,&opt->voopt);

    /* initial feature point position */
    if (!initfeatpos(cam1,cam2,ins,uv1,uv2,K,pf,0)) {
        trace(2,"initial feature position fail\n");
        return 0;
    }
    feat->flag=TRACK_INIT; /* initial flag */
    return 1;
}
/* get feature point measurement data-----------------------------------------*/
static int getfeatmeas(const feature *feat,const insstate_t *ins,double *pfn)
{
    double uv[3],K[9]={0},pf[3];
    cam_t cam;

    /* camera calibration parameters */
    cam.k1=ins->k1; cam.k2=ins->k2;
    cam.p1=ins->p1; cam.p2=ins->p2;
    getK(ins,K);

    uv[0]=feat->u; uv[1]=feat->v; pfn[2]=uv[2]=1.0;

    if (matinv(K,3)) return 0;
    matmul3v("N",K,uv,pf);

    /* undistort feature point */
    undistortradtan(&cam,pf,pfn,NULL);
    return 1;
}
/* jacobians of measurement wrt. (alpha,beta,rho)-----------------------------*/
static void jacobian(const double *C,const double *t,const double *h,double *J)
{
    J[0]=-C[0]/h[2]+(h[0]/SQR(h[2]))*C[2]; /* for alpha. */
    J[1]=-C[1]/h[2]+(h[1]/SQR(h[2]))*C[2];

    J[2]=-C[3]/h[2]+(h[0]/SQR(h[2]))*C[5]; /* for beta. */
    J[3]=-C[4]/h[2]+(h[1]/SQR(h[2]))*C[5];

    J[4]=-t[0]/h[2]+(h[0]/SQR(h[2]))*t[2]; /* for rho. */
    J[5]=-t[1]/h[2]+(h[1]/SQR(h[2]))*t[2];
}
/* iteration for Gauss Newton optimization------------------------------------*/
static int itergnop(const int* index,int k,trackd_t *feat,const insopt_t *opt,
                    const insstate_t *ins,const double *x,double *Jprev,
                    double *dx,double *Jderiv)
{
    double *E,*W,*v,C[9],t[3],zhat[3],T[3],h[3],J[6],*w;
    double Jnew,*EWE,*Wi,*Ev;
    cams_t *cam1,*cam2;
    int i,j,l,n;

    trace(3,"itergnop:\n");

    W =zeros(2*k,2*k); E=zeros(2*k,3); v=zeros(2*k,1);
    Wi=zeros(2*k,2*k);
    w =zeros(2*k,  1);

    for (n=0,j=0,cam1=&vofilt.data[index[0]];j<k;j++) {

        /* relative camera transformation */
        cam2=&vofilt.data[index[j]];
        reltf(cam1->time,cam2->time,cam1,cam2,C,t);

        /* for the weight matrix */
        w[2*n+0]=VAR_FEAT/SQR(ins->fx);
        w[2*n+1]=VAR_FEAT/SQR(ins->fy);

        /* feature observation */
        if (!getfeatmeas(&feat->data[j],ins,zhat)) continue;
        T[0]=x[0]; T[1]=x[1]; T[2]=1.0;

        matmul3v("N",C,T,h);
        for (i=0;i<3;i++) h[i]=T[i]+x[2]*t[i];

        /* form the error vector */
        v[2*n+0]=zhat[0]-h[0]/h[2]; v[2*n+1]=zhat[1]-h[1]/h[2];

        /* form the jacobians matrix */
        jacobian(C,t,h,J);

        for (i=0;i<2;i++) {
            for (l=0;l<3;l++) E[2*n+i+l*3]=J[i+l*3];
        }
        n++; /* number of valid observation */
    }
    /* check is valid */
    if (!n) {
        free(v); free(W); free(E); free(Wi);
        return 0;
    }
    for (j=0;j<n;j++) {
        W[2*j+2*j*(2*n)]=w[2*j]; W[2*j+1+(2*j+1)*(2*n)]=w[2*j+1];
    }
    trace(3,"W=\n"); tracemat(3,W,2*n,2*n,12,6);
    trace(3,"v=\n"); tracemat(3,v,2*n,1,12,6);
    trace(3,"E=\n"); tracemat(3,E,2*n,3,12,6);

    /* calculate the cost function */
    matcpy(Wi,W,2*n,2*n);
    if (matinv(Wi,2*n)) {
        free(v); free(W); free(E); free(Wi); return 0;
    }
    matmul33("TNN",v,Wi,v,1,2*n,2*n,1,&Jnew);

    /* solve */
    EWE=mat(3,3); Ev=mat(3,1);
    matmul33("TNN",E,Wi,E,3,2*n,2*n,3,EWE);
    matmul33("TNN",E,Wi,v,3,2*n,2*n,1,Ev);
    if (matinv(EWE,3)) {
        free(v); free(W); free(E); free(Wi);
        free(EWE); free(Ev);
        return 0;
    }
    matmul3("N",EWE,Ev,dx);
    for (i=0;i<3;i++) dx[i]=-dx[i];

    /* calculate the cost */
    *Jderiv=fabs((0.5*Jnew-*Jprev)/(0.5*Jnew)); *Jprev=0.5*Jnew;

    free( v); free(  W); free( E);
    free(Wi); free(EWE); free(Ev);
    return 1;
}
/* calculate the position estimate of the feature using Gauss Newton optimization
 * ---------------------------------------------------------------------------*/
static int calcgnposest(trackd_t *feat,const insopt_t *opt,const insstate_t *ins,
                        gtime_t time,int *index)
{
    int i,k=0;
    double Jprev=1E9,Jderiv,x[3],dx[3],pf[3],Jn;

    trace(3,"calcgnposest:\n");

    /* feature position in first camera frame */
    if (!(k=estfeatpos(feat,opt,ins,time,index,pf))) return 0;

    /* re-parameter for feature position */
    for (i=0;i<3;i++) x[i]=pf[i]/pf[2];

    for (i=0;i<MAXITER;i++) {

        /* iteration for optimization */
        if (!itergnop(index,k,feat,opt,ins,x,&Jprev,dx,&Jderiv)) continue;

        if (Jderiv<0.01) break;
        else {
            for (i=0;i<3;i++) x[i]+=dx[i];
        }
    }
    if (i==MAXITER) {
        trace(2,"optimization fail\n");
        return 0;
    }
    dx[0]=x[0]; dx[1]=x[1];
    dx[2]= 1.0;
    matmul("NN",3,1,3,1.0/x[2],vofilt.data[index[0]].Cce,
           dx,0.0,pf);

    /* update feature */
    for (i=0;i<3;i++) {
        feat->xyz[i]=pf[i]+vofilt.data[index[0]].re[i];
    }
    Jn=Jprev/SQR(k);
    return k&&Jn<MAX_GNCOST_NORM;
}
/* measurement update for a feature point-------------------------------------*/
static int featmeas(trackd_t *trk,const insstate_t *ins,const insopt_t *opt,
                    gtime_t time,double *H,double *v,double *R,double *Hf,int *index,int *ki)
{
    int i,j,k,nv,nx=vofilt.nx,ixs;
    double *Hc,*Hkp,*Hfo,vc[2],Rc[4],*r,*Hfb;
    cams_t *pcam=NULL;

    static int ifo=xiCfo(opt),nfo=xnCfo(opt);
    static int ikp=xiCkp(opt),nkp=xnCkp(opt);
    
    /*  compute feature position*/
    if (!(*ki=k=calcgnposest(trk,opt,ins,time,index))) {
        return 0;
    }
    Hc =zeros(2,9); Hkp=zeros(2,4);
    Hfo=zeros(2,4); Hfb=zeros(2,3);
    r=zeros(2*k,1);

    /* for every camera pose */
    for (i=0,nv=0,*ki=0,pcam=vofilt.data;i<k&&pcam;i++) {

        /* residual vector */
        featurev(trk,&trk->data[i],&pcam[index[i]],ins,
                 vc);

        /* jacobians matrix */
        featureH(trk,&pcam[index[i]],ins,
                 Hfb,Hc,Hfo,Hkp);

        /* measurement variance */
        featureR(trk,&trk->data[i],&pcam[index[i]],ins,
                 Rc);

        /* index of camera pose */
        ixs=ins->nx+index[i]*9;
        if (H) {
            /* for fx,fy,ox,oy */
            for (j=0;j<4&&nfo;j++) {
                H[(2*nv+0)+(ifo+j)*nx]=Hfo[4*j+0];
                H[(2*nv+1)+(ifo+j)*nx]=Hfo[4*j+1];
            }
            /* for k1,k2,p1,p2 */
            for (j=0;j<4&&nkp;j++) {
                H[(2*nv+0)+(ikp+j)*nx]=Hkp[4*j+0];
                H[(2*nv+1)+(ikp+j)*nx]=Hkp[4*j+1];
            }
            /* for camera pose */
            for (j=0;j<9;j++) {
                H[(2*nv+0)+(ixs+j)*nx]=Hc[9*j+0]; /* col-major=>row-major*/
                H[(2*nv+1)+(ixs+j)*nx]=Hc[9*j+1];
            }
        }
        if (v) {
            v[2*nv+0]=vc[0]; v[2*nv+1]=vc[1];
        }
        r[2*nv]=Rc[0]; r[2*nv+1]=Rc[3];
        if (Hf) {
            for (j=0;j<3;j++) {
                Hf[2*nv+0+j*2]=Hfb[0+j*2]; /* for feature position */
                Hf[2*nv+1+j*2]=Hfb[1+j*2];
            }
        }
        if (fabs(v[2*nv+0])>MAX_FEAT_RES||fabs(v[2*nv+1])>MAX_FEAT_RES) continue;
        index[(*ki)++]=index[i];
        nv++;
    }
    if (nv<=0) {
        trace(2,"no feature point measurement\n");
        free(Hc ); free(Hfo);
        free(Hkp); free(r);
        return 0;
    }
    for (i=0;i<nv&&R;i++) {
        R[2*i+0+(2*i+0)*2*nv]=r[2*i+0];
        R[2*i+1+(2*i+1)*2*nv]=r[2*i+1];
    }
    trace(3,"R=\n"); tracemat(3,R,2*nv,2*nv,12,6);
    trace(3,"v=\n"); tracemat(3,v,2*nv,1,12,6);

    trace(3,"H=\n");
    tracemat(3,H,nx,2*nv,12,6);

    free(Hc ); free(Hfo);
    free(Hkp); free(r);
    return 2*nv;
}
/* use left nullspace of H to `H/v/R'-----------------------------------------*/
static int nulltrick(const double *H,const double *v,const double *R,
                     int nx,int nv,
                     const double *Hf, int m,int n,
                     double *Ho,double *Ro,double *vo,int *p,int *q)
{
    double *N,*Hft;
    int i,j;

    N=mat(m,n); Hft=mat(m,n); matt(Hf,m,n,Hft);

    if (!null(Hft,n,m,N,&i,&j)) {
        return 0;
    }
    matmul("TN",j,nx,i,1.0,N,H,0.0,Ho);
    matmul("TN",j,nv,i,1.0,N,v,0.0,vo);

    matmul33("TNN",N,R,N,j,nv,nv,j,Ro);

    *p=i; *q=j; free(N); free(Hft);
    return 1;
}
/* update camera pose and ins states------------------------------------------*/
static void updatestat(double *dx,insstate_t *ins,const insopt_t *opt,
                       const int* index,int k)
{
    static int ifo=xiCfo(opt),nfo=xnCfo(opt);
    static int ikp=xiCkp(opt),nkp=xnCkp(opt);
    static int ncl=xnCl(opt);
    int i,j,ns=ins->nx;

    /* update ins states */
    clp(ins,opt,dx);

    /* update camera states */
    for (i=0;i<k;i++) {
        corratt(dx+ns+9*index[i],vofilt.data[index[i]].Cce);

        for (j=0;j<3;j++) {
            vofilt.data[index[i]].re[j]-=dx[ns+9*index[i]+6];
            vofilt.data[index[i]].ve[j]-=dx[ns+9*index[i]+3];
        }
        vofilt.data[index[i]].flag=1;
    }
    /* update camera calibration parameters */
    if (nfo) {
        ins->fx=ins->fx-dx[ifo+0]; ins->fy=ins->fy-dx[ifo+1];
        ins->ox=ins->ox-dx[ifo+2];
        ins->oy=ins->oy-dx[ifo+3];
    }
    if (nkp) {
        ins->k1=ins->k1-dx[ikp+0];
        ins->k2=ins->k2-dx[ikp+1];
        ins->p1=ins->p1-dx[ikp+2];
        ins->p2=ins->p2-dx[ikp+3];
    }
    for (i=0;i<ncl;i++) dx[i]=0.0;
}
/* outlier detection for feature point measurement----------------------------*/
static int outdetect(const double *P,const double *Ho,const double *Ro,
                     const double *vo,int nv,int nx)
{
    double *HPHT,*HRP,r;
    int i;

    HPHT=mat(nv,nv); HRP=mat(nv,nv);
    matmul33("NNT",Ho,P,Ho,nv,nx,nx,nv,HPHT);
    for (i=0;i<nv*nv;i++) {
        HRP[i]=Ro[i]+HPHT[i];
    }
    if (matinv(HRP,nv)) {
        free(HPHT); free(HRP); return 0;
    }
    matmul33("TNN",vo,HRP,vo,1,nv,nv,1,&r);
    free(HPHT); free(HRP);
    return r>chisqr[nv];
}
/* feature point measurement update-------------------------------------------*/
static int updatefeatmeas(const insopt_t *opt,insstate_t *ins,gtime_t time)
{
    double *H,*R,*v,*Ho,*Ro,*vo,*Hf,*x,*Hot;
    int i,nv,nx=vofilt.nx,n=vofilt.n,index[MAX_TRACK_LEN],ki=0,j,flag=1;
    int p,q;

    trace(3,"updatefeatmeas: time=%s\n",time_str(time,4));

    H =zeros(n*2,nx); R=zeros(n*2,n*2);
    v =zeros(n*2, 1);
    Hf=zeros(n*2, 3);

    Ho=mat(n*2,nx); Ro=mat(n*2,n*2);
    vo=mat(n*2, 1);
    x=zeros(nx,1); Hot=mat(2*n,nx);

    for (nv=0,i=0;i<tracks.n;i++) {
        if (tracks.data[i].flag!=2) continue;
        if (!(nv=featmeas(&tracks.data[i],ins,opt,time,H,v,R,Hf,index,&ki))) flag=0;

        /* do null-space trick */
        if (flag&&!nulltrick(H,v,R,nx,nv,Hf,nv,3,Ho,Ro,vo,&p,&q)) flag=0;

#if OUT_DETECT
        /* outlier detection */
        if (flag&&outdetect(vofilt.Px,Ho,Ro,vo,nv,nx)) flag=0;
#endif
        /* remove feature */
        for (j=0;j<vofilt.n;j++) {
            if (!hash_find_feature(vofilt.data[j].trackfeat,tracks.data[i].uid)) continue;
            hash_delete_feature(&vofilt.data[j].trackfeat,tracks.data[i].uid);
        }
        /* remove match index in hash table */
        rmmatchindex(tracks.data[i].last_idx);

        tracks.data[i].last_idx=-1;
        if (flag==0) continue;

        /* ekf filter */
        matt(Ho,q,nx,Hot);
        if (filter(x,vofilt.Px,Hot,vo,Ro,nx,q)) {

            trace(2,"filter error\n");
            flag++; continue;
        }
        /* update camera pose and ins states */
        updatestat(x,ins,opt,index,ki);

        tracks.data[i].flag=TRACK_FILTER;
    }
    free(H ); free(v ); free(R ); free(x);
    free(Hf); free(Ho); free(vo);
    free(Ro); free(Hot);
    return nv;
}
/* update.--------------------------------------------------------------------*/
static int updateall(insstate_t *ins,const insopt_t *opt,const img_t *img)
{
    trace(3,"updateall:\n");

    /* augment states with a new camera pose */
    augstates(&vofilt,ins,opt,img->id);

    /* match feature points */
    if (img==NULL||matchfeats(&matchs,img)<=0) {

        trace(2,"match feature points fail\n");
        return 0;
    }
    /* update all track data */
    if (!updatetrack(opt,img)) return 0;
    if (!updatecamtrack(opt,img->time,matchs.Ip.id,img->id)) return 0;

    /* update camera/ins states */
    updatefeatmeas(opt,ins,img->time);

    /* remove redundant camera */
    rmcamera(&vofilt,opt);
    return 1;
}
/* using visual odometry to aid ins/gnss pose estimating-----------------------
 * args:    insopt_t *opt   I   ins options
 *          insstate_t *ins IO  ins states
 *          imud_t *imu     I   imu measurement data
 *          img_t *img      I   image measurement data
 *          int flag        I   flag-0: propagate ins-camera pose covariance
 *                              flag-1: feature points measurement update
 * return: status (1: ok, 0: fail)
 * ----------------------------------------------------------------------------*/
extern int voigpos(const insopt_t *opt,insstate_t *ins,const imud_t *imu,
                   const img_t *img,int flag)
{
    trace(3,"voigpos: time=%s\n",time_str(imu->time,4));

    switch (flag) {
        case 0: return propagate(&vofilt,ins); 
        case 1: return propagate(&vofilt,ins)&&updateall(ins,opt,img);
    }
    if (inss==NULL) inss=ins;
    if (opts==NULL) opts=opt;
    return 0;
}
/* predict image point using homography---------------------------------------
 * args:    insopt_t *opt   I   ins options
 *          double *R,*t    I   rotation and translation between frames
 *          double *K       I   camera calibration parameters
 *          double *uv      I   image coordinate of precious frame
 *          double *pre_uv  O   predict image coordinate in current frame
 * return: status (1: ok, 0: fail)
 * ---------------------------------------------------------------------------*/
extern int predictfeat(const double *R,const double *t,const double *K,
                       const double *uv,double *uvp)
{
    double H[9],Ki[9],p[3],pp[3],T[16];
    double C[9],r[9];

    trace(3,"predictfeatH:\n");

    /* pose of current relative to precious */
    rt2tf(R,t,T); if (matinv(T,4)) return 0;
    tf2rt(T,C,r);

    /* homography by rotation */
    matcpy(Ki,K,3,3); if (matinv(Ki,3)) return 0;
    matmul33("NNN",K,C,Ki,3,3,3,3,H);

    trace(3,"H=\n"); tracemat(3,H,3,3,12,6);

    /* predict image coordinate */
    p[0]=uv[0];
    p[1]=uv[1];
    p[2]=1.0;
    matmul("NN",3,1,3,1.0,H,p,0.0,pp);

    uvp[0]=pp[0]/pp[2]; uvp[1]=pp[1]/pp[2];
    return 1;
}
/* compute coarse feature point position--------------------------------------
 * args:    trackd_t *feat  IO  tracking feature point
 *          gtime_t　time   I   current time
 *          double *pf      O   feature point position relative to first frame
 * return: status (1: ok, 0: fail)
 * ---------------------------------------------------------------------------*/
extern int getfeaturepos(trackd_t *feat,gtime_t time,double *pf)
{
    int index[MAX_TRACK_LEN];
    return estfeatpos(feat,opts,inss,time,index,pf);
}
/* get camera pose which observe given feature point--------------------------
 * args:    gtime_t  time  I  current time
 *          double *R      O  rotation matrix
 *          double *t      O  translation matrix
 * return: status (1: ok, 0: fail)
 * ---------------------------------------------------------------------------*/
extern int getcamerapose(gtime_t time,double *R,double *t)
{
    int i;
    for (i=0;i<vofilt.n;i++) {
        if (fabs(timediff(time,vofilt.data[i].time))>1E-5) continue;
        if (R) matcpy(R,vofilt.data[i].Cce,3,3);
        if (t) matcpy(t,vofilt.data[i].re ,3,1);
        return 1;
    }
    return 0;
}





#include "carvig.h"

typedef struct pose {
    double sow;
    double rr[3],Cbe[9],vn[3];
} pose_t;

typedef struct poseall {
    pose_t *data;
    int n,nmax;
} poseall_t;

#define add_data_func_delc(data_type,add_data_type)            \
    static int add_##data_type##_##add_data_type(              \
                             data_type *ivec,                  \
                             const add_data_type* n)           \
    {                                                          \
        add_data_type *obs_data;                               \
        if (ivec->nmax<=ivec->n) {                             \
            if (ivec->nmax<=0) ivec->nmax=10*2;                \
            else ivec->nmax*=2;                                \
            if (!(obs_data=(add_data_type*)realloc(ivec->data, \
                sizeof(add_data_type)*ivec->nmax))) {          \
                free(ivec->data);                              \
                ivec->data=NULL; ivec->n=ivec->nmax=0;         \
            }                                                  \
            ivec->data=obs_data;                               \
        }                                                      \
        ivec->data[ivec->n++]=*n;                              \
        return 1;                                              \
    }                                                          \

add_data_func_delc(poseall,pose)

static void getatt1(const double *re,const double *Cbe,double *rpy)
{
    double llh[3],C[9],Cnb[9];

    ecef2pos(re,llh);
    ned2xyz(llh,C);
    matmul("TN",3,3,3,1.0,Cbe,C,0.0,Cnb);
    dcm2rpy(Cnb,rpy);
}

int main()
{
    const char *file0="/media/sujinglan/Windows/Users/sujinglan/Desktop/gps-ins-cam/20190308-fb/20190308/pro/Loose/pro.txt";
    const char *file1="/media/sujinglan/Windows/Users/sujinglan/Desktop/gps-ins-cam/out_vf";
    const char *ofile="./out";
    const char *ofile_vo="./out_vo";
    char buff[1024];
    FILE *fp0=fopen(file0,"r");
    FILE *fp1=fopen(file1,"r");
    FILE *fpo=fopen(ofile,"w");
    FILE *fpo_vo=fopen(ofile_vo,"w");

    int count=0;

    poseall_t pose={0};

    while (fgets(buff,sizeof(buff),fp0)) {

        pose_t p={0};
        double rpy[3],Cnb[9],Cne[9],pos[3];

        if (sscanf(buff,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&p.sow,p.rr,p.rr+1,p.rr+2,
                   rpy,rpy+1,rpy+2,p.vn,p.vn+1,p.vn+2)<10) {
            continue;
        }
        rpy[0]*=D2R;
        rpy[1]*=D2R;
        rpy[2]*=D2R;
        rpy2dcm(rpy,Cnb); ecef2pos(p.rr,pos);
        ned2xyz(pos,Cne);
        matmul("NT",3,3,3,1.0,Cne,Cnb,0.0,p.Cbe);

        add_poseall_pose(&pose,&p);

        count++;
        if (count>=90000) break;
    }
    double tmp[2+6+6];
    pose_t *pp,*pc;

    while (fgets(buff,sizeof(buff),fp1)) {
        if (sscanf(buff,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf\n",
                   &tmp[0],&tmp[1],
                   &tmp[2],&tmp[3],&tmp[4],&tmp[5],&tmp[6],&tmp[7],
                   &tmp[8],&tmp[9],&tmp[10],&tmp[11],&tmp[12],&tmp[13])<14) {
            continue;
        }
        pp=NULL;
        pc=NULL;
        for (int i=0;i<pose.n;i++) {
            if (fabs(pose.data[i].sow-tmp[0])<DTTOL) {
                pp=&pose.data[i];
                break;
            }
        }
        for (int i=0;i<pose.n;i++) {
            if (fabs(pose.data[i].sow-tmp[1])<DTTOL) {
                pc=&pose.data[i];
                break;
            }
        }
        if (pp==NULL||pc==NULL) {
            continue;
        }

        double dC[9],phi[3],dr[3],drr[3];

        matmul("TN",3,3,3,1.0,pp->Cbe,pc->Cbe,0.0,dC);
        so3_log(dC,phi,NULL);

        drr[0]=pc->rr[0]-pp->rr[0];
        drr[1]=pc->rr[1]-pp->rr[1];
        drr[2]=pc->rr[2]-pp->rr[2];
        matmul("TN",3,1,3,1.0,pp->Cbe,drr,0.0,dr);

        fprintf(fpo,"%6.4lf  %6.4lf  %10.6lf  %10.6lf  %10.6lf  %10.6lf  %10.6lf  %10.6lf  %10.6lf  %10.6lf  %10.6lf\n",
                tmp[0],tmp[1],
                phi[0],phi[1],phi[2],
                dr[0],dr[1],dr[2],
                pp->vn[0],pp->vn[1],pp->vn[2]);
        fprintf(fpo_vo,"%s",buff);
        fflush(fpo_vo);
        fflush(fpo);
    }
    return 0;
}
/*------------------------------------------------------------------------------
 * ins-track.cc : feature points tracking functions
 *
 * reference :
 *    [1] P.D.Groves, Principles of GNSS, Intertial, and Multisensor Integrated
 *        Navigation System, Artech House, 2008
 *    [2] Bruce D. Lucas and Takeo Kanade. An Iterative Image Registration Technique
 *        with an Application to Stereo Vision. International Joint Conference on
 *        Artificial Intelligence, pages 674-679, 1981.
 *    [3] Carlo Tomasi and Takeo Kanade. Detection and Tracking of Point Features.
 *        Carnegie Mellon University Technical Report CMU-CS-91-132, April 1991.
 *    [4] Jianbo Shi and Carlo Tomasi. Good Features to Track. IEEE Conference on
 *        Computer Vision and Pattern Recognition, pages 593-600, 1994.
 *    [5] Stan Birchfield. Derivation of Kanade-Lucas-Tomasi Tracking Equation.
 *        Unpublished, January 1997.
 *
 * version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
 * history : 2018/10/25 1.0 new
 *-----------------------------------------------------------------------------*/
#include <carvig.h>

/* constants------------------------------------------------------------------*/
#define REFINE          1         /* tracking refine feature points */
#define OUTPUT_PPM      1         /* output tracking data to ppm file */

/* global variables-----------------------------------------------------------*/
static long int id_seed=1;        /* generate a new feature id */

/* initial track---------------------------------------------------------------
 * args:    trackd_t *data  IO  track set data
 *          voopt_t *opt    I   track options
 * return: status (1: ok,0: fail)
 * ----------------------------------------------------------------------------*/
extern int inittrack(trackd_t *data,const voopt_t *opt)
{
    gtime_t t0={0};
    trace(3,"inittrack:\n");

    data->ts=data->te=t0;
    data->n=data->nmax=0; data->uid=id_seed++; sprintf(data->name,"%ld",data->uid);
    data->data=NULL;
    data->I   =NULL;
    return 1;
}
/* find a track in tracking set data------------------------------------------*/
static int findtrack(const track_t *track,const match_point *mp)
{
    int i;
    for (i=0;i<track->n;i++) {
        if (mp->ip==track->data[i].last_idx) return i;
    }
    return track->n;
}
/* copy image data------------------------------------------------------------*/
static int copyimg(img_t *out,const img_t *in)
{
    out->data=(unsigned char*)malloc(sizeof(unsigned char)*in->w*in->h);
    memcpy(out->data,in->data,sizeof(unsigned char)*in->w*in->h);
    out->h=in->h;
    out->w=in->w;
    out->time=in->time;
}
/* add new feature and image data to track------------------------------------*/
static int addnewfeatimg(trackd_t *track,const feature *feat,const img_t *img)
{
    img_t *img_data; int i;
    feature *obs_data;

    if (track->nmax<=track->n) {
        if (track->nmax<=0) track->nmax=1024; else track->nmax*=2;
        if (!(obs_data=(feature *)realloc(track->data,sizeof(feature)*track->nmax))) {
            free(track->data); track->data=NULL;

            track->n=track->nmax=0;
            return -1;
        }
        track->data=obs_data;
#if TRACE_TRACK
        if (!(img_data=(img_t*)realloc(track->I,sizeof(img_t)*track->nmax))) {
            for (i=0;i<track->n;i++) free(track->I[i].data); track->I[i].data=NULL;

            free(track->I); track->I=NULL;
            track->n=track->nmax=0;
            return -1;
        }
        track->I=img_data;
#endif
    }
#if TRACE_TRACK
    copyimg(&track->I[track->n],img);
#endif
    track->data[track->n++]=*feat;
    return 1;
}
/* add new track to list------------------------------------------------------*/
static int addnewtrack(track_t *track,const trackd_t *data)
{
    trackd_t *obs_data;

    if (track->nmax<=track->n) {
        if (track->nmax<=0) track->nmax=16*2; else track->nmax*=2;
        if (!(obs_data=(trackd_t *)realloc(track->data,sizeof(trackd_t)*track->nmax))) {

            free(track->data);
            track->data=NULL;
            track->n=track->nmax=0; return -1;
        }
        track->data=obs_data;
    }
    track->data[track->n++]=*data;
    return 1;
}
/* create a new track---------------------------------------------------------*/
static int newtrack(const match_point_t *mp,gtime_t tp,gtime_t tc,int curr_frame,
                    const voopt_t *opt,track_t *track,
                    const img_t *pimg,const img_t *cimg)
{
    /* new track. */
    feature fp,fc; trackd_t ntrack;

    fp.time=tp; fp.u=mp->up; fp.v=mp->vp;
    fp.valid=1; fp.status=FEAT_CREATE;

    fc.time=tc; fc.u=mp->uc; fc.v=mp->vc;
    fc.valid=1; fc.status=FEAT_CREATE;

    /* add feature. */
    inittrack(&ntrack,opt);
    addnewfeatimg(&ntrack,&fp,pimg); addnewfeatimg(&ntrack,&fc,cimg);

    /* update frame id */
    ntrack.first_frame=curr_frame-1;
    ntrack.last_frame =curr_frame;
    ntrack.last_idx   =mp->ic;

    /* timestamp */
    ntrack.ts=tp; ntrack.te=tc;

    ntrack.flag=1;

    /* add new track */
    if (addnewtrack(track,&ntrack)<=0) {
        trace(2,"add new track fail\n");
        return 0;
    }
    return track->n-1;
}
/* matched feature points set data convert to track set-----------------------
 * args:  match_set *set   I  matched feature points set data
 *        gtime_t tp,tc    I  precious/current time
 *        int curr_frame   I  uid of current frame
 *        img_t *pimg,cimg I  precious/current image data (NULL: no input)
 *        voopt_t *opt     I  options
 *        track_t *track   O  tracking data
 * return: status (1: ok,0: fail)
 * ---------------------------------------------------------------------------*/
extern int match2track(const match_set *mset,gtime_t tp,gtime_t tc,int curr_frame,
                       const img_t *pimg,const img_t *cimg,
                       const voopt_t *opt,track_t *track)
{
    register int i,idx=0;
    feature feat;

    track->nnew=0; track->nupd=0;

    trace(3,"match2track:\n");

    /* initial flag of all track */
    for (i=0;i<track->n;i++) track->data[i].flag=2;

    /* initial track */
    if (track->n==0) {
        for (i=0;i<mset->n;i++) {

            /* create a new track */
            newtrack(&mset->data[i],tp,tc,curr_frame,opt,track,pimg,cimg);

            /* update index */
            track->newtrack[track->nnew]=i;
            track->nnew++;
        }
        return 1;
    }
    for (i=0;i<mset->n;i++) {
        if (findtrack(track,&mset->data[i])<0) {

            /* create a new track */
            if ((idx=newtrack(&mset->data[i],tp,tc,curr_frame,opt,track,pimg,cimg)<0)) return 0;

            /* update index */
            track->newtrack[track->nnew]=idx;
            track->nnew++;
            continue;
        }
#if REFINE
        track->data[idx].data[track->data[idx].n-1].u=mset->data[i].up;
        track->data[idx].data[track->data[idx].n-1].v=mset->data[i].vp;
#endif
        /* update track */
        feat.u=mset->data[i].uc;
        feat.v=mset->data[i].vc;

        feat.time=tc;
        feat.valid=1; feat.status=FEAT_CREATE;
        addnewfeatimg(&track->data[idx],&feat,cimg);

        /* track flag. */
        track->data[idx].last_idx  =mset->data[i].ic;
        track->data[idx].last_frame=curr_frame;

        /* timestamp */
        track->data[idx].te=tc;

        /* update index */
        track->updtrack[track->nupd]=idx;
        track->nupd++;

        track->data[idx].flag=0;
    }
    return track->n>0;
}
/* free tracking--------------------------------------------------------------*/
extern void freetrack(trackd_t *track)
{
    trace(3,"freetrack:\n");

    int i;
    if (track->data) {
        free(track->data); track->data=NULL;
    }
    if (track->I) {
        for (i=0;i<track->n;i++) free(track->I[i].data);
    }
    track->n=track->nmax=0;
    track->first_frame=0;
    track->last_frame =0; track->last_idx=0; track->uid=0;
}
/* free tracking set data ----------------------------------------------------*/
extern void freetrackset(track_t *track)
{
    int i;
    trace(3,"freetrackset:\n");

    if (track->data) {
        for (i=0;i<track->n;i++) freetrack(&track->data[i]);
    }
    track->n=track->nmax=0;
}
/* write feature points to PPM file-------------------------------------------
 * args:    feature *featurelist    I  feature points list
 *          int nfeature            I  number of feature points
 *          unsigned char *greyimg  I  image data
 *          int ncols,nrows         I  size of image
 *          char *filename          I  ppm file path
 * return: none
 * ---------------------------------------------------------------------------*/
extern void feat2ppm(feature *featurelist, int nfeature,unsigned char *greyimg,
                     int ncols,int nrows,char *filename)
{
    trace(3,"feat2ppm:\n");

    int nbytes=ncols*nrows*sizeof(char);
    unsigned char *redimg,*grnimg,*bluimg;
    int offset;
    int x,y,xx,yy;
    int i;

    /* Allocate memory for component images */
    redimg=(unsigned char*)malloc(nbytes);
    grnimg=(unsigned char*)malloc(nbytes);
    bluimg=(unsigned char*)malloc(nbytes);
    if (redimg==NULL||grnimg==NULL||bluimg==NULL) return;

    memcpy(redimg,greyimg,nbytes);
    memcpy(grnimg,greyimg,nbytes);
    memcpy(bluimg,greyimg,nbytes);

    /* Overlay features in red */
    for (i=0;i<nfeature;i++) {
        x=(int)(featurelist[i].u+0.5);
        y=(int)(featurelist[i].v+0.5);
        for (yy=y-1;yy<=y+1;yy++) for (xx=x-1;xx<=x+1;xx++) {
            if (xx>=0&&yy>=0&&xx<ncols&&yy<nrows)  {
                offset=yy*ncols+xx;
                *(redimg+offset)=255;
                *(grnimg+offset)=0;
                *(bluimg+offset)=0;
            }
        }
    }
    /* write to PPM file */
    ppmWriteFileRGB(filename,redimg,grnimg,bluimg,ncols,nrows);

    /* free memory */
    free(redimg);
    free(grnimg);
    free(bluimg);
}
/* trace tracking data--------------------------------------------------------*/
extern void tracetrack(const track_t *track)
{
    char dir[126],path[126];
    trackd_t *pta;
    int i,j;

    trace(3,"tracetrack: n=%d\n",track->n);

    for (i=0;i<track->n;i++) {
        pta=&track->data[i]; trace(3,"\ntrack-uid=%4d:  \n");
        for (j=0;j<pta->n;j++) {

            /* [time  frame_id  (u,v)] */
            trace(3,"    %s  %4d  (%8.3lf  %8.3lf)  \n",
                  time_str(pta->data[j].time,4),
                  pta->first_frame+j,pta->data[j].u,pta->data[j].v);
        }
        /* output to PPM file */
#if OUTPUT_PPM&TRACE_TRACK
        sprintf(dir,"/media/sujinglan/Files/carvig-debug/test_track/track_%ld",pta->uid);

        /* create new dir. */
        if (access(dir,F_OK)==-1) {
            if (mkdir(dir,0777)!=0) continue;
        }
        for (j=0;j<pta->n;j++) {

            sprintf(path,"%s/%d.ppm",dir,j);
            feat2ppm(&pta->data[j],1,pta->I[j].data,pta->I[j].w,
                     pta->I[j].h,path);
        }
#endif
    }
}
/* get track given track-uid--------------------------------------------------
 * args:    track_t *track  I  tracking data
 *          int uid         I  track uid need to find
 * return: pointer of track (NULL: no found)
 * ---------------------------------------------------------------------------*/
extern trackd_t *gettrack(const track_t *track,int uid)
{
    trace(3,"gettrack:\n");
    int i;
    for (i=0;i<track->n;i++) {
        if (track->data[i].uid==uid) return &track->data[i];
    }
    return NULL;
}
/* check feature point whether is out of view or track lost-------------------
 * args:    track_t *track  I  tracking data
 *          voopt_t *opt    I  visual odomtery options
 *          int id          I  feature point id
 *          gtime_t time    I  tracking timestamp
 * return:  status (1: out of view,0: in view,-1: no found)
 * ---------------------------------------------------------------------------*/
extern int outofview(const track_t *track,const voopt_t *opt,int id,gtime_t time)
{
    trackd_t *pt=NULL;
    int i;

    for (i=0;i<track->n;i++) {
        if (track->data[i].uid==id) {pt=&track->data[i]; break;}
    }
    if (!pt) return -1;
    for (i=0;i<pt->n;i++) {
        if (fabs(timediff(time,pt->data[i].time))<1E-6) return 1;
    }
    return 0;
}




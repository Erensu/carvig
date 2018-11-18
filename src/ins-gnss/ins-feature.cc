/*------------------------------------------------------------------------------
 * ins-feature.cc : klt feature point track functions
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
 * history : 2017/10/17 1.0 new
 *-----------------------------------------------------------------------------*/
#include <carvig.h>
#include <vision.h>

/* type definitions ----------------------------------------------------------*/
typedef struct  {
    int width;
    float data[MAX_KERNEL_WIDTH];
} ConvolutionKernel;
typedef enum {SELECTING_ALL, REPLACING_SOME} selectionMode;

/* global variable -----------------------------------------------------------*/
static const int mindist                  =10;
static const int window_size              =7;
static const int min_eigenvalue           =1;
static const float min_determinant        =0.01;
static const float min_displacement       =0.1;
static const int max_iterations           =10;
static const float max_residue            =10.0;
static const float grad_sigma             =1.0;
static const float smooth_sigma_fact      =0.1;
static const float pyramid_sigma_fact     =0.9;
static const float step_factor            =1.0;
static const KLT_BOOL sequentialMode      =false;
static const KLT_BOOL lighting_insensitive=false;

/* for affine mapping*/
static const int affineConsistencyCheck    =-1;
static const int affine_window_size        =15;
static const int affine_max_iterations     =10;
static const float affine_max_residue      =10.0;
static const float affine_min_displacement =0.02;
static const float affine_max_displacement_differ=1.5;

static const KLT_BOOL smoothBeforeSelecting=true;
static const KLT_BOOL writeInternalImages  =false;
static const int search_range              =15;
static const int nSkippedPixels            =0;
static float sigma_last                    =-10.0;

/* Kernels */
static ConvolutionKernel gauss_kernel;
static ConvolutionKernel gaussderiv_kernel;

/*----------------------------------------------------------------------------
 * _quicksort
 * Replacement for qsort().  Computing time is decreased by taking
 * advantage of specific knowledge of our array (that there are
 * three ints associated with each point).
 *
 * This routine generously provided by
 *      Manolis Lourakis <lourakis@csi.forth.gr>
 *
 * NOTE: The results of this function may be slightly different from
 * those of qsort().  This is due to the fact that different sort
 * algorithms have different behaviours when sorting numbers with the
 * same value: Some leave them in the same relative positions in the
 * array, while others change their relative positions. For example,
 * if you have the array [c d b1 a b2] with b1=b2, it may be sorted as
 * [a b1 b2 c d] or [a b2 b1 c d].
 *----------------------------------------------------------------------------*/
#define SWAP3(list, i, j)               \
{                                       \
     register int *pi, *pj, tmp;        \
     pi=list+3*(i); pj=list+3*(j);      \
                                        \
     tmp=*pi;                           \
     *pi++=*pj;                         \
     *pj++=tmp;                         \
                                        \
     tmp=*pi;                           \
     *pi++=*pj;                         \
     *pj++=tmp;                         \
                                        \
     tmp=*pi;                           \
     *pi=*pj;                           \
     *pj=tmp;                           \
}
static void quicksort(int *pointlist, int n)
{
    unsigned int i,j,ln,rn;

    while (n>1) {
        SWAP3(pointlist,0,n/2);
        for (i=0,j = n;;) {
            do --j;
            while (pointlist[3*j+2]<pointlist[2]);
            do
                ++i;
            while (i<j&&pointlist[3*i+2]>pointlist[2]);
            if (i>=j) break;
            SWAP3(pointlist,i,j);
        }
        SWAP3(pointlist,j,0);
        ln=j;
        rn=n-++j;
        if (ln<rn) {
            quicksort(pointlist,ln);
            pointlist+=3*j;
            n=rn;
        }
        else {
            quicksort(pointlist+3*j,rn);
            n=ln;
        }
    }
}
#undef SWAP3
/*----------------------------------------------------------------------------
 * _comparePoints
 *
 * Used by qsort (in _KLTSelectGoodFeatures) to determine
 * which feature is better.
 * By switching the '>' with the '<', qsort is fooled into sorting
 * in descending order.
 *----------------------------------------------------------------------------*/
#ifdef KLT_USE_QSORT
static int _comparePoints(const void *a, const void *b)
{
    int v1=*(((int *)a)+2);
    int v2=*(((int *)b)+2);

    if      (v1>v2)  return -1;
    else if (v1<v2)  return 1;
    else return 0;
}
#endif
/*----------------------------------------------------------------------------*/
static void fillFeaturemap(int x, int y,unsigned char *featuremap,
                           int mindist,int ncols,int nrows)
{
    register int ix, iy;

    for (iy=y-mindist;iy<=y+mindist;iy++)
        for (ix=x-mindist;ix<=x+mindist;ix++) {
            if (ix>=0&&ix<ncols&&iy>=0&&iy<nrows) featuremap[iy*ncols+ix]=1;
        }
}
/*----------------------------------------------------------------------------
 * enforceMinimumDistance
 *
 * Removes features that are within close proximity to better features.
 *
 * INPUTS
 * featurelist:  A list of features.  The nFeatures property
 *               is used.
 *
 * OUTPUTS
 * featurelist:  Is overwritten.  Nearby "redundant" features are removed.
 *               Writes -1's into the remaining elements.
 *
 * RETURNS
 * The number of remaining features.
 *----------------------------------------------------------------------------*/
static void enforceMinimumDistance(
        int *pointlist,              /* featurepoints */
        int npoints,                 /* number of featurepoints */
        pfeaturelist_t featurelist,  /* features */
        int ncols, int nrows,        /* size of images */
        int mindist,                 /* min. dist b/w features */
        int min_eigenvalue,          /* min. eigenvalue */
        KLT_BOOL overwriteAllFeatures)
{
    register int indx;         /* Index into features */
    register int x,y,val;      /* Location and trackability of pixel under consideration */
    unsigned char *featuremap; /* Boolean array recording proximity of features */
    int *ptr;

    /* Cannot add features with an eigenvalue less than one */
    if (min_eigenvalue<1) min_eigenvalue=1;

    /* Allocate memory for feature map and clear it */
    featuremap=(unsigned char *)malloc(ncols*nrows*sizeof(unsigned char));
    memset(featuremap,0,ncols*nrows);

    /* Necessary because code below works with (mindist-1) */
    mindist--;

    /* If we are keeping all old good features, then add them to the featuremap */
    if (!overwriteAllFeatures)
        for (indx=0;indx<featurelist->n;indx++)
            if (featurelist->feature[indx]->valid>=0) {
                x=(int)featurelist->feature[indx]->x;
                y=(int)featurelist->feature[indx]->y;
                fillFeaturemap(x,y,featuremap,mindist,ncols,nrows);
            }

    /* For each feature point, in descending order of importance, do ... */
    ptr =pointlist;
    indx=0;
    while (1) {

        /* If we can't add all the points, then fill in the rest
           of the featurelist with -1's */
        if (ptr>=pointlist+3*npoints) {
            while (indx<featurelist->n) {
                if (overwriteAllFeatures||
                    featurelist->feature[indx]->valid<0) {
                    featurelist->feature[indx]->x=-1;
                    featurelist->feature[indx]->y=-1;
                    featurelist->feature[indx]->valid  =KLT_NOT_FOUND;
                    featurelist->feature[indx]->aff_img=NULL;
                    featurelist->feature[indx]->aff_img_gradx=NULL;
                    featurelist->feature[indx]->aff_img_grady=NULL;
                    featurelist->feature[indx]->aff_x=-1.0f;
                    featurelist->feature[indx]->aff_y=-1.0f;
                    featurelist->feature[indx]->aff_Axx=1.0;
                    featurelist->feature[indx]->aff_Ayx=0.0;
                    featurelist->feature[indx]->aff_Axy=0.0;
                    featurelist->feature[indx]->aff_Ayy=1.0;
                }
                indx++;
            }
            break;
        }
        x  =*ptr++;
        y  =*ptr++;
        val=*ptr++;

        /* Ensure that feature is in-bounds */
        assert(x>=0);
        assert(x<ncols);
        assert(y>=0);
        assert(y<nrows);

        while (!overwriteAllFeatures &&indx<featurelist->n
               &&featurelist->feature[indx]->valid>=0)
            indx++;

        if (indx>=featurelist->n) break;

        /* If no neighbor has been selected, and if the minimum
           eigenvalue is large enough, then add feature to the current list */
        if (!featuremap[y*ncols+x]&&val>=min_eigenvalue)  {
            featurelist->feature[indx]->x    =(float)x;
            featurelist->feature[indx]->y    =(float)y;
            featurelist->feature[indx]->valid=(int)val;
            featurelist->feature[indx]->aff_img=NULL;
            featurelist->feature[indx]->aff_img_gradx=NULL;
            featurelist->feature[indx]->aff_img_grady=NULL;
            featurelist->feature[indx]->aff_x=-1.0f;
            featurelist->feature[indx]->aff_y=-1.0f;
            featurelist->feature[indx]->aff_Axx=1.0;
            featurelist->feature[indx]->aff_Ayx=0.0;
            featurelist->feature[indx]->aff_Axy=0.0;
            featurelist->feature[indx]->aff_Ayy=1.0;
            indx++;

            /* Fill in surrounding region of feature map, but
               make sure that pixels are in-bounds
               */
            fillFeaturemap(x,y,featuremap,mindist,ncols,nrows);
        }
    }
    /* Free feature map  */
    free(featuremap);
}
/* createArray2D--------------------------------------------------------------*/
static void** createArray2D(int ncols, int nrows, int nbytes)
{
    char **tt;
    register int i;

    tt=(char **)malloc(nrows*sizeof(void*)+ncols*nrows*nbytes);
    if (tt==NULL) {
        trace(2,"createArray2D) Out of memory");
    }
    for (i=0;i<nrows;i++) {
        tt[i]=((char *)tt)+(nrows*sizeof(void *)+i*ncols*nbytes);
    }
    return (void**)tt;
}
/* computeKernels-------------------------------------------------------------*/
static void computeKernels(float sigma,ConvolutionKernel *gauss,
                           ConvolutionKernel *gaussderiv)
{
    const float factor=0.01;   /* for truncating tail */
    register int i;

    trace(3,"_computeKernels:\n");

    assert(MAX_KERNEL_WIDTH%2==1);
    assert(sigma>=0.0);

    /* Compute kernels, and automatically determine widths */
    {
        const int hw=MAX_KERNEL_WIDTH/2;
        float max_gauss=1.0,max_gaussderiv=(float)(sigma*exp(-0.5));

        /* Compute gauss and deriv */
        for (i=-hw;i<=hw;i++) {
            gauss->data[i+hw]     =(float)exp(-i*i/(2*sigma*sigma));
            gaussderiv->data[i+hw]=-i*gauss->data[i+hw];
        }
        /* Compute widths */
        gauss->width=MAX_KERNEL_WIDTH;
        for (i=-hw;fabs(gauss->data[i+hw]/max_gauss)<factor;
             i++,gauss->width-=2);

        gaussderiv->width=MAX_KERNEL_WIDTH;
        for (i=-hw;fabs(gaussderiv->data[i+hw]/max_gaussderiv)<factor;
             i++,gaussderiv->width-=2);

        if (gauss->width==MAX_KERNEL_WIDTH||gaussderiv->width==MAX_KERNEL_WIDTH) {

            trace(2,"(_computeKernels) MAX_KERNEL_WIDTH %d is too small for "
                    "a sigma of %f",MAX_KERNEL_WIDTH,sigma);
        }
    }
    /* Shift if width less than MAX_KERNEL_WIDTH */
    for (i=0;i<gauss->width;i++) {
        gauss->data[i]=gauss->data[i+(MAX_KERNEL_WIDTH-gauss->width)/2];
    }
    for (i=0;i<gaussderiv->width;i++) {
        gaussderiv->data[i]=gaussderiv->data[i+(MAX_KERNEL_WIDTH-gaussderiv->width)/2];
    }
    /* Normalize gauss and deriv */
    {
        const int hw=gaussderiv->width/2;
        float den;

        den=0.0;
        for (i=0;i<gauss->width;i++) den+=gauss->data[i];
        for (i=0;i<gauss->width;i++) gauss->data[i]/=den;
        den=0.0;
        for (i=-hw;i<=hw;i++) den-=i*gaussderiv->data[i+hw];
        for (i=-hw;i<=hw;i++) gaussderiv->data[i+hw]/=den;
    }
    sigma_last=sigma;
    return;
}
/*----------------------------------------------------------------------------*/
static void KLTGetKernelWidths(float sigma,int *gauss_width,
                               int *gaussderiv_width)
{
    trace(3,"_KLTGetKernelWidths:\n");

    computeKernels(sigma,&gauss_kernel,&gaussderiv_kernel);
    *gauss_width=gauss_kernel.width;
    *gaussderiv_width=gaussderiv_kernel.width;
}
/*----------------------------------------------------------------------------*/
static float KLTComputeSmoothSigma(tracking_context_t *tc)
{
    return tc->smooth_sigma_fact*MAX(tc->window_width,tc->window_height);
}
/*----------------------------------------------------------------------------*/
static float pyramidSigma(tracking_context_t* tc)
{
    return tc->pyramid_sigma_fact*tc->subsampling;
}
/*----------------------------------------------------------------------------*/
static void klt_UpdateTCBorder(tracking_context_t* tc)
{
    register float val;
    register int pyramid_gauss_hw;
    register int smooth_gauss_hw;
    register int gauss_width, gaussderiv_width;
    register int num_levels = tc->nPyramidLevels;
    register int n_invalid_pixels;
    register int window_hw;
    register int ss = tc->subsampling;
    register int ss_power;
    register int border;
    register int i;

    trace(3,"klt_UpdateTCBorder:\n");

    /* Check window size (and correct if necessary) */
    if (tc->window_width%2!=1) {
        tc->window_width=tc->window_width+1;
        trace(3,"(KLTUpdateTCBorder) Window width must be odd.  "
                "Changing to %d.\n",tc->window_width);
    }
    if (tc->window_height%2!=1) {
        tc->window_height=tc->window_height+1;
        trace(3,"(KLTUpdateTCBorder) Window height must be odd.  "
                "Changing to %d.\n",tc->window_height);
    }
    if (tc->window_width<3) {
        tc->window_width=3;
        trace(3,"(KLTUpdateTCBorder) Window width must be at least three.\n"
                "Changing to %d.\n",tc->window_width);
    }
    if (tc->window_height<3) {
        tc->window_height=3;
        trace(3,"(KLTUpdateTCBorder) Window height must be at least three.\n"
                "Changing to %d.\n",tc->window_height);
    }
    window_hw=MAX(tc->window_width,tc->window_height)/2;

    /* Find widths of convolution windows */
    KLTGetKernelWidths(KLTComputeSmoothSigma(tc),&gauss_width,&gaussderiv_width);
    smooth_gauss_hw=gauss_width/2;

    KLTGetKernelWidths(pyramidSigma(tc),&gauss_width,&gaussderiv_width);
    pyramid_gauss_hw=gauss_width/2;

    /* Compute the # of invalid pixels at each level of the pyramid.
       n_invalid_pixels is computed with respect to the ith level
       of the pyramid.  So, e.g., if n_invalid_pixels = 5 after
       the first iteration, then there are 5 invalid pixels in
       level 1, which translated means 5*subsampling invalid pixels
       in the original level 0.
       */
    n_invalid_pixels=smooth_gauss_hw;
    for (i=1;i<num_levels;i++)  {
        val=((float) n_invalid_pixels+pyramid_gauss_hw)/ss;
        n_invalid_pixels=(int)(val+0.99);  /* Round up */
    }
    /* ss_power = ss^(num_levels-1) */
    ss_power=1;
    for (i=1;i<num_levels;i++) ss_power*=ss;

    /* Compute border by translating invalid pixels back into */
    /* original image */
    border=(n_invalid_pixels+window_hw)*ss_power;

    tc->borderx=border;
    tc->bordery=border;
    return;
}
/* KLT change TCP-Pyram ID----------------------------------------------------
 * args:    tracking_context_t *tc  IO  klt tracking context
 *          int search_range        I   klt search range
 * return: none
 * ---------------------------------------------------------------------------*/
static void klt_changeTCPyramid(tracking_context_t* tc,int search_range)
{
    register float window_halfwidth,subsampling;

    trace(3,"klt_changeTCPyramid:\n");

    /* check window size (and correct if necessary) */
    if (tc->window_width%2!=1) {
        tc->window_width=tc->window_width+1;
        trace(2,"(KLTChangeTCPyramid) Window width must be odd.Changing to %d.\n",
              tc->window_width);
    }
    if (tc->window_height%2!=1) {
        tc->window_height=tc->window_height+1;
        trace(2,"(KLTChangeTCPyramid) Window height must be odd.Changing to %d.\n",
              tc->window_height);
    }
    if (tc->window_width<3) {
        tc->window_width=3;
        trace(2,"(KLTChangeTCPyramid) Window width must be at least three."
                      " Changing to %d.\n",
              tc->window_width);
    }
    if (tc->window_height<3) {
        tc->window_height=3;
        trace(2,"(KLTChangeTCPyramid) Window height must be at least three.\n"
                "Changing to %d.\n",tc->window_height);
    }
    window_halfwidth=MIN(tc->window_width,tc->window_height)/2.0f;

    subsampling=((float)search_range)/window_halfwidth;

    if (subsampling<1.0) {		/* 1.0 = 0+1 */
        tc->nPyramidLevels=1;
    }
    else if (subsampling<=3.0) {	/* 3.0 = 2+1 */
        tc->nPyramidLevels=2;
        tc->subsampling=2;
    }
    else if (subsampling<= 5.0) {	/* 5.0 = 4+1 */
        tc->nPyramidLevels=2;
        tc->subsampling=4;
    }
    else if (subsampling<=9.0) {	/* 9.0 = 8+1 */
        tc->nPyramidLevels=2;
        tc->subsampling = 8;
    }
    else {
        /* The following lines are derived from the formula:
           search_range =
           window_halfwidth * \sum_{i=0}^{nPyramidLevels-1} 8^i,
           which is the same as:
           search_range =
           window_halfwidth * (8^nPyramidLevels - 1)/(8 - 1).
           Then, the value is rounded up to the nearest integer.
           */
        float val=(float)(log(7.0*subsampling+1.0)/log(8.0));
        tc->nPyramidLevels=(int)(val+0.99);
        tc->subsampling=8;
    }
    return;
}
/* create klt tracking context------------------------------------------------*/
extern tracking_context_t* klt_create_tracking_context()
{
    tracking_context_t *tc;

    trace(3,"klt_create_tracking_context:\n");

    /* allocate memory */
    tc=(tracking_context_t*)malloc(sizeof(tracking_context_t));

    /* set values to default values */
    tc->mindist              =mindist;
    tc->window_width         =window_size;
    tc->window_height        =window_size;
    tc->sequentialMode       =sequentialMode;
    tc->smoothBeforeSelecting=smoothBeforeSelecting;
    tc->writeInternalImages  =writeInternalImages;
    tc->lighting_insensitive =lighting_insensitive;
    tc->min_eigenvalue       =min_eigenvalue;
    tc->min_determinant      =min_determinant;
    tc->max_iterations       =max_iterations;
    tc->min_displacement     =min_displacement;
    tc->max_residue          =max_residue;
    tc->grad_sigma           =grad_sigma;
    tc->smooth_sigma_fact    =smooth_sigma_fact;
    tc->pyramid_sigma_fact   =pyramid_sigma_fact;
    tc->step_factor          =step_factor;
    tc->nSkippedPixels       =nSkippedPixels;
    tc->pyramid_last         =NULL;
    tc->pyramid_last_gradx   =NULL;
    tc->pyramid_last_grady   =NULL;

    /* for affine mapping */
    tc->affineConsistencyCheck        =affineConsistencyCheck;
    tc->affine_window_width           =affine_window_size;
    tc->affine_window_height          =affine_window_size;
    tc->affine_max_iterations         =affine_max_iterations;
    tc->affine_max_residue            =affine_max_residue;
    tc->affine_min_displacement       =affine_min_displacement;
    tc->affine_max_displacement_differ=affine_max_displacement_differ;

    /* Change nPyramidLevels and subsampling */
    klt_changeTCPyramid(tc,search_range);

    /* Update border, which is dependent upon  */
    /* smooth_sigma_fact, pyramid_sigma_fact, window_size, and subsampling */
    klt_UpdateTCBorder(tc);
    return tc;
}
/* create feature points list-------------------------------------------------
 * args:    int nFeatures  I  number of feature points
 * return: pointer of feature points list
 * ---------------------------------------------------------------------------*/
extern featurelist_t* klt_create_featurelist(int nFeatures)
{
    featurelist_t* fl;
    pfeature_t first;
    int nbytes=sizeof(featurelist_t)+nFeatures*sizeof(pfeature_t)+
            nFeatures*sizeof(feature_t);
    register int i;

    trace(3,"klt_create_featurelist:\n");

    /* allocate memory for feature list */
    fl=(featurelist_t*)malloc(nbytes);

    /* set parameters */
    fl->n=nFeatures;

    /* set pointers */
    fl->feature=(pfeature_t*)(fl+1);
    first=(pfeature_t)(fl->feature+nFeatures);
    for (i=0;i<nFeatures;i++) {
        fl->feature[i]=first+i;
        fl->feature[i]->aff_img=NULL;  /* initialization fixed by Sinisa Segvic */
        fl->feature[i]->aff_img_gradx=NULL;
        fl->feature[i]->aff_img_grady=NULL;
    }
    /* return feature list */
    return fl;
}
/* create feature frame-------------------------------------------------------
 * args:    int nFrames  I  number of feature frames
 * return: pointer of feature frames
 * ---------------------------------------------------------------------------*/
extern feature_history_t* klt_create_feature_history(int nFrames)
{
    pfeature_history_t fh;
    pfeature_t first;
    int nbytes=sizeof(feature_history_t)+nFrames*sizeof(pfeature_t)+
            nFrames*sizeof(feature_t);
    register int i;

    trace(3,"klt_create_feature_history:\n");

    /* Allocate memory for feature history */
    fh=(pfeature_history_t)malloc(nbytes);

    /* Set parameters */
    fh->nFrames=nFrames;

    /* Set pointers */
    fh->feature=(pfeature_t *)(fh+1);
    first=(pfeature_t)(fh->feature+nFrames);
    for (i=0;i<nFrames;i++) {
        fh->feature[i]=first+i;
    }
    /* return feature history */
    return fh;
}
/* create feature point table-------------------------------------------------
 * args:     int nFrames    I  number of frames
 *           int nFeatures  I  number of feature points
 * return: pointer of feature points table
 * ---------------------------------------------------------------------------*/
extern feature_table_t* klt_create_feature_table(int nFrames,int nFeatures)
{
    pfeature_table_t ft;
    pfeature_t first;
    int nbytes=sizeof(feature_table_t);
    register int i,j;

    trace(3,"klt_create_feature_table:\n");

    /* Allocate memory for feature history */
    ft=(pfeature_table_t)malloc(nbytes);

    /* Set parameters */
    ft->nFrames  =nFrames;
    ft->nFeatures=nFeatures;

    /* Set pointers */
    ft->feature = (pfeature_t **)createArray2D(nFrames,nFeatures,sizeof(pfeature_t));
    first=(pfeature_t)malloc(nFrames*nFeatures*sizeof(feature_t));
    for (j=0;j<nFeatures;j++) {
        for (i=0;i<nFrames;i++) ft->feature[j][i]=first+j*nFrames+i;
    }
    /* return feature table */
    return ft;
}
/* given a pointer to image data (probably unsigned chars), copy data to a
 * float image.
 * args:    unsigned char* img    I  input image data
 *          int nclos,nrows       I  size of image
 *          floatimg_t *floatimg  O  output float image data
 * return: none
 * ---------------------------------------------------------------------------*/
extern void tofloatImage(unsigned char *img,int ncols, int nrows,
                         floatimg_t* floatimg)
{
    trace(3,"tofloatImage:\n");

    unsigned char *ptrend=img+ncols*nrows;
    float *ptrout=floatimg->data;

    /* Output image must be large enough to hold result */
    assert(floatimg->ncols>=ncols);
    assert(floatimg->nrows>=nrows);

    floatimg->ncols=ncols;
    floatimg->nrows=nrows;

    while (img<ptrend) *ptrout++=(float)*img++;
}
/* create a float image-------------------------------------------------------
 * args:    int ncol,nrow  I  size of image data
 * return: pointer of image data
 * ---------------------------------------------------------------------------*/
extern floatimg_t* create_floatimg(int ncols,int nrows)
{
    floatimg_t* floatimg;
    int nbytes=sizeof(floatimg_t)+ncols*nrows*sizeof(float);

    floatimg=(floatimg_t*)malloc(nbytes);
    if (floatimg==NULL) {
        trace(2,"(create_floatimg) Out of memory");
    }
    floatimg->ncols=ncols;
    floatimg->nrows=nrows;
    floatimg->data =(float*)(floatimg+1);

    return floatimg;
}
/* free float image data------------------------------------------------------
 * args:    floatimg_t *floatimg  IO  input/output float image data
 * return: none
 * ---------------------------------------------------------------------------*/
extern void free_floatimg(floatimg_t* floatimg)
{
    trace(3,"free_floatimg:\n");
    free(floatimg); return;
}
/*----------------------------------------------------------------------------*/
static void convolveImageHoriz(floatimg_t* imgin,ConvolutionKernel kernel,
                               floatimg_t* imgout)
{
    float *ptrrow=imgin->data;           /* Points to row's first pixel */
    register float *ptrout=imgout->data, /* Points to next output pixel */
            *ppp;

    register float sum;
    register int radius=kernel.width/2;
    register int ncols =imgin->ncols,nrows=imgin->nrows;
    register int i,j,k;

    /* Kernel width must be odd */
    assert(kernel.width%2==1);

    /* Must read from and write to different images */
    assert(imgin!=imgout);

    /* Output image must be large enough to hold result */
    assert(imgout->ncols>=imgin->ncols);
    assert(imgout->nrows>=imgin->nrows);

    /* For each row, do ... */
    for (j=0;j<nrows;j++)  {

        /* Zero leftmost columns */
        for (i=0;i<radius;i++) *ptrout++=0.0;

        /* Convolve middle columns with kernel */
        for (;i<ncols-radius;i++)  {
            ppp=ptrrow+i-radius;
            sum=0.0;
            for (k=kernel.width-1;k>=0;k--) sum+=*ppp++*kernel.data[k];
            *ptrout++ = sum;
        }
        /* Zero rightmost columns */
        for (;i<ncols;i++) *ptrout++=0.0;
        ptrrow+=ncols;
    }
    return;
}
/*----------------------------------------------------------------------------*/
static void convolveImageVert(floatimg_t* imgin,ConvolutionKernel kernel,
                              floatimg_t* imgout)
{
    float *ptrcol=imgin->data;            /* Points to row's first pixel */
    register float *ptrout=imgout->data,  /* Points to next output pixel */
            *ppp;
    register float sum;
    register int radius=kernel.width/2;
    register int ncols =imgin->ncols,nrows=imgin->nrows;
    register int i,j,k;

    /* Kernel width must be odd */
    assert(kernel.width%2==1);

    /* Must read from and write to different images */
    assert(imgin!=imgout);

    /* Output image must be large enough to hold result */
    assert(imgout->ncols>=imgin->ncols);
    assert(imgout->nrows>=imgin->nrows);

    /* For each column, do ... */
    for (i=0;i<ncols;i++)  {

        /* Zero topmost rows */
        for (j=0;j<radius;j++) {
            *ptrout=0.0; ptrout+=ncols;
        }
        /* Convolve middle rows with kernel */
        for (;j<nrows-radius;j++) {
            ppp=ptrcol+ncols*(j-radius);
            sum=0.0;
            for (k=kernel.width-1;k>=0;k--) {
                sum+=*ppp*kernel.data[k];
                ppp+=ncols;
            }
            *ptrout=sum; ptrout+=ncols;
        }
        /* Zero bottommost rows */
        for (;j<nrows;j++) {
            *ptrout=0.0; ptrout+=ncols;
        }
        ptrcol++;
        ptrout-=nrows*ncols-1;
    }
    return;
}
/*----------------------------------------------------------------------------*/
static void convolveSeparate(floatimg_t* imgin,
                             ConvolutionKernel horiz_kernel,
                             ConvolutionKernel vert_kernel,
                             floatimg_t* imgout)
{
    /* Create temporary image */
    floatimg_t* tmpimg;
    tmpimg=create_floatimg(imgin->ncols, imgin->nrows);

    /* do convolution */
    convolveImageHoriz(imgin,horiz_kernel,tmpimg);
    convolveImageVert (tmpimg,vert_kernel,imgout);

    /* free memory */
    free_floatimg(tmpimg);
    return;
}
/* compute  of image----------------------------------------------------------
 * args:    floatimg_t* img    I  input float image data
 *          float sigma        I  sigma for process image data
 *          floatimg_t* gradx  O  x-gradients of input image data
 *          floatimg_t* grady  O  y-gradients of input image data
 * return: none
 * ---------------------------------------------------------------------------*/
extern void compute_gradients(floatimg_t *img, float sigma, floatimg_t *gradx,
                              floatimg_t *grady)
{
    /* Output images must be large enough to hold result */
    assert(gradx->ncols>=img->ncols);
    assert(gradx->nrows>=img->nrows);
    assert(grady->ncols>=img->ncols);
    assert(grady->nrows>=img->nrows);

    /* Compute kernels, if necessary */
    if (fabs(sigma-sigma_last)>0.05) {
        computeKernels(sigma,&gauss_kernel,&gaussderiv_kernel);
    }
    convolveSeparate(img,gaussderiv_kernel,gauss_kernel,gradx);
    convolveSeparate(img,gauss_kernel,gaussderiv_kernel,grady);
    return;
}
/* free pyramid---------------------------------------------------------------
 * args:    pyramid_t*  IO  pyramid
 * return: none
 * ---------------------------------------------------------------------------*/
extern void free_pyramid(pyramid_t *pyramid)
{
    register int i;

    /* free images */
    for (i=0;i<pyramid->nLevels;i++) {
        free_floatimg(pyramid->img[i]);
    }
    /* free structure */
    free(pyramid); return;
}
/* free KLT tracking context--------------------------------------------------
 * args:    tracking_context_t* context  IO  input track context
 * return: none
 * ---------------------------------------------------------------------------*/
extern void klt_free_track_context(tracking_context_t *tc)
{
    trace(3,"klt_free_track_context:\n");

    if (tc->pyramid_last)  {
        free_pyramid((ppyramid_t)tc->pyramid_last);
    }
    if (tc->pyramid_last_gradx) {
        free_pyramid((ppyramid_t)tc->pyramid_last_gradx);
    }
    if (tc->pyramid_last_grady) {
        free_pyramid((ppyramid_t)tc->pyramid_last_grady);
    }
    free(tc);
    return;
}
/* free feature points list---------------------------------------------------
 * args:    featurelist_t *flist  IO  input feature list
 * return: none
 * ---------------------------------------------------------------------------*/
extern void klt_free_featurelist(featurelist_t *flist)
{
    /* for affine mapping */
    int indx;

    trace(3,"klt_free_featurelist:\n");

    for (indx=0;indx<flist->n;indx++) {
        /* free image and gradient  */
        free_floatimg(flist->feature[indx]->aff_img);
        free_floatimg(flist->feature[indx]->aff_img_gradx);
        free_floatimg(flist->feature[indx]->aff_img_grady);
        flist->feature[indx]->aff_img = NULL;
        flist->feature[indx]->aff_img_gradx = NULL;
        flist->feature[indx]->aff_img_grady = NULL;
    }
    free(flist);
    return;
}
/* free feature points frame---------------------------------------------------
 * args:    feature_history_t *fhist  IO  input feature frame
 * return: none
 * ---------------------------------------------------------------------------*/
extern void klt_free_featurehist(feature_history_t *fhist)
{
    trace(3,"klt_free_featurehist:\n");
    if (fhist) free(fhist);
}
/* free feature points table---------------------------------------------------
 * args:    feature_table_t *table  IO  input feature points table
 * return: none
 * ---------------------------------------------------------------------------*/
extern void klt_free_featuretable(feature_table_t *table)
{
    trace(3,"klt_free_featuretable:\n");
    if (table->feature) {
        free(table->feature[0][0]);
        free(table->feature);
        free(table);
    }
    return;
}
/*----------------------------------------------------------------------------
 * Given the three distinct elements of the symmetric 2x2 matrix
 *                     [gxx gxy]
 *                     [gxy gyy],
 * Returns the minimum eigenvalue of the matrix.
 *----------------------------------------------------------------------------*/
static float minEigenvalue(float gxx,float gxy,float gyy)
{
    return (float)((gxx+gyy-sqrt((gxx-gyy)*(gxx-gyy)+4*gxy*gxy))/2.0);
}
/*----------------------------------------------------------------------------*/
static void sortPointList(int *pointlist,int npoints)
{
#ifdef KLT_USE_QSORT
    qsort(pointlist,npoints,3*sizeof(int),_comparePoints);
#else
    quicksort(pointlist,npoints);
#endif
}
/*----------------------------------------------------------------------------*/
static void KLTSelectGoodFeatures(tracking_context_t* tc, unsigned char *img,
                                  int ncols,int nrows,
                                  pfeaturelist_t featurelist,
                                  selectionMode mode)
{
    pfloatimg_t floatimg,gradx,grady;
    int window_hw,window_hh;
    int *pointlist;
    int npoints=0;
    KLT_BOOL overwriteAllFeatures=(mode==SELECTING_ALL)?
                                   true:false;
    KLT_BOOL floatimages_created=false;

    /* Check window size (and correct if necessary) */
    if (tc->window_width%2!=1) {
        tc->window_width=tc->window_width+1;
        trace(2,"Tracking context's window width must be odd.  "
                           "Changing to %d.\n",tc->window_width);
    }
    if (tc->window_height%2!=1) {
        tc->window_height=tc->window_height+1;
        trace(2,"Tracking context's window height must be odd.  "
                           "Changing to %d.\n",tc->window_height);
    }
    if (tc->window_width<3) {
        tc->window_width=3;
        trace(2,"Tracking context's window width must be at least three.\n"
                           "Changing to %d.\n",
              tc->window_width);
    }
    if (tc->window_height<3) {
        tc->window_height=3;
        trace(2,"Tracking context's window height must be at least three.\n"
                           "Changing to %d.\n",
                   tc->window_height);
    }
    window_hw=tc->window_width/2;
    window_hh=tc->window_height/2;

    /* Create pointlist, which is a simplified version of a featurelist, */
    /* for speed.  Contains only integer locations and values. */
    pointlist=(int *)malloc(ncols*nrows*3*sizeof(int));

    /* Create temporary images, etc. */
    if (mode==REPLACING_SOME&&
        tc->sequentialMode && tc->pyramid_last!=NULL)  {
        floatimg=((ppyramid_t) tc->pyramid_last)->img[0];
        gradx=((ppyramid_t)tc->pyramid_last_gradx)->img[0];
        grady=((ppyramid_t)tc->pyramid_last_grady)->img[0];
        assert(gradx!=NULL);
        assert(grady!=NULL);
    }
    else {
        floatimages_created=true;
        floatimg=create_floatimg(ncols,nrows);
        gradx   =create_floatimg(ncols,nrows);
        grady   =create_floatimg(ncols,nrows);
        if (tc->smoothBeforeSelecting) {
            pfloatimg_t tmpimg;
            tmpimg=create_floatimg(ncols,nrows);
            tofloatImage(img,ncols,nrows,tmpimg);
            compute_smoothedimg(tmpimg,ComputeSmoothSigma(tc),floatimg);
            free_floatimg(tmpimg);
        }
        else tofloatImage(img,ncols,nrows,floatimg);

        /* Compute gradient of image in x and y direction */
        compute_gradients(floatimg,tc->grad_sigma,gradx,grady);
    }
    /* write internal images */
    if (tc->writeInternalImages)  {
        writeFloatImageToPGM(floatimg,(char*)"kltimg_sgfrlf.pgm");
        writeFloatImageToPGM(gradx,(char*)"kltimg_sgfrlf_gx.pgm");
        writeFloatImageToPGM(grady,(char*)"kltimg_sgfrlf_gy.pgm");
    }
    /* compute trackability of each image pixel as the minimum
       of the two eigenvalues of the Z matrix */
    {
        register float gx, gy;
        register float gxx, gxy, gyy;
        register int xx, yy;
        register int *ptr;
        float val;
        unsigned int limit=1;
        int borderx=tc->borderx;	/* Must not touch cols */
        int bordery=tc->bordery;	/* lost by convolution */
        int x, y;
        int i;

        if (borderx<window_hw) borderx=window_hw;
        if (bordery<window_hh) bordery=window_hh;

        /* Find largest value of an int */
        for (i=0;i<sizeof(int);i++) limit*=256;
        limit=limit/2-1;

        /* For most of the pixels in the image, do ... */
        ptr=pointlist;
        for (y=bordery;y<nrows-bordery;y+=tc->nSkippedPixels+1)
            for (x=borderx;x<ncols-borderx;x+=tc->nSkippedPixels+1)  {

                /* Sum the gradients in the surrounding window */
                gxx=0; gxy=0; gyy=0;
                for (yy=y-window_hh;yy<=y+window_hh;yy++)
                    for (xx=x-window_hw;xx<=x+window_hw;xx++)  {
                        gx=*(gradx->data+ncols*yy+xx);
                        gy=*(grady->data+ncols*yy+xx);
                        gxx+=gx*gx;
                        gxy+=gx*gy;
                        gyy+=gy*gy;
                    }
                /* Store the trackability of the pixel as the minimum
                   of the two eigenvalues */
                *ptr++=x;
                *ptr++=y;
                val=minEigenvalue(gxx,gxy,gyy);
                if (val>limit) {
                    trace(2,"minimum eigenvalue %f is "
                            "greater than the capacity of an int; setting "
                            "to maximum value",val);
                    val=(float)limit;
                }
                *ptr++=(int)val;
                npoints++;
            }
    }
    /* Sort the features  */
    sortPointList(pointlist,npoints);

    /* Check tc->mindist */
    if (tc->mindist<0)  {
        trace(2,"Tracking context field tc->mindist "
                "is negative (%d); setting to zero",
              tc->mindist);
        tc->mindist=0;
    }
    /* Enforce minimum distance between features */
    enforceMinimumDistance(
            pointlist,
            npoints,
            featurelist,
            ncols,nrows,
            tc->mindist,
            tc->min_eigenvalue,
            overwriteAllFeatures);

    /* Free memory */
    free(pointlist);
    if (floatimages_created)  {
        free_floatimg(floatimg);
        free_floatimg(gradx);
        free_floatimg(grady);
    }
    return;
}
/* select good feature points from tracking context---------------------------
 * args:    tracking_context_t* tc  IO  tracking context
 *          unsigned char* img      I   input image data
 *          int ncols,nrows         I   size of image
 *          featurelist_t* fl       O   feature points list
 * return: none
 * ---------------------------------------------------------------------------*/
extern void klt_select_goodfeatures(tracking_context_t* tc,unsigned char *img,
                                    int ncols,int nrows,
                                    featurelist_t *fl)
{
    trace(3,"klt_select_goodfeatures:\n");
    KLTSelectGoodFeatures(tc,img,ncols,nrows,fl,SELECTING_ALL);
}
/* image smooth---------------------------------------------------------------
 * args:    pfloatimg_t *img    I  input image data
 *          float sigma         I  sigma for smoothing image
 *          pfloatimg_t smooth  O  smoothed image
 * return: none
 * ---------------------------------------------------------------------------*/
extern void compute_smoothedimg(pfloatimg_t img,float sigma,
                                pfloatimg_t smooth)
{
    trace(3,"compute_smoothedimg:\n");

    /* output image must be large enough to hold result */
    assert(smooth->ncols>=img->ncols);
    assert(smooth->nrows>=img->nrows);

    /* compute kernel, if necessary; gauss_deriv is not used */
    if (fabs(sigma-sigma_last)>0.05) {
        computeKernels(sigma,&gauss_kernel,&gaussderiv_kernel);
    }
    convolveSeparate(img,gauss_kernel,
                     gauss_kernel,smooth);
    return;
}
/* get smooth sigma value-----------------------------------------------------
 * args:    tracking_context_t* tc  IO  tracking context
 * return: sigma value
 * ---------------------------------------------------------------------------*/
extern float ComputeSmoothSigma(tracking_context_t* tc)
{
    return (tc->smooth_sigma_fact*MAX(tc->window_width,tc->window_height));
}
/*----------------------------------------------------------------------------*/
extern int klt_count_remaining_features(pfeaturelist_t fl)
{
    int count=0;
    int i;

    for (i=0;i<fl->n;i++) {
        if (fl->feature[i]->valid>=0) count++;
    }
    return count;
}
/*----------------------------------------------------------------------------
 * KLTReplaceLostFeatures
 *
 * Main routine, visible to the outside.  Replaces the lost features
 * in an image.
 *
 * INPUTS
 * tc:	Contains parameters used in computation (size of image,
 *        size of window, min distance b/w features, sigma to compute
 *        image gradients, # of features desired).
 * img:	Pointer to the data of an image (probably unsigned chars).
 *
 * OUTPUTS
 * features:	List of features.  The member nFeatures is computed.
 *----------------------------------------------------------------------------*/
extern void klt_replace_lost_features(tracking_context_t *tc,unsigned char *img,
                                     int ncols,int nrows,pfeaturelist_t fl)
{
    int nLostFeatures=fl->n-klt_count_remaining_features(fl);

    trace(3,"klt_replace_lost_features:\n");

    /* If there are any lost features, replace them */
    if (nLostFeatures>0) {
        KLTSelectGoodFeatures(tc,img,ncols,nrows,fl,REPLACING_SOME);
    }
}
/* create a pyramid-----------------------------------------------------------
 * args:    int ncols,nrows  I  size of image
 *          int sunsampling  I  step of sub-sample
 *          int nlevels      I  number of pyramid levels
 * return: pointer to a pyramid
 * ---------------------------------------------------------------------------*/
extern ppyramid_t klt_create_pyramid(int ncols,int nrows,int subsampling,
                                     int nlevels)
{
    ppyramid_t pyramid;
    int nbytes=sizeof(pyramid_t) +
               nlevels*sizeof(pfloatimg_t *)+
               nlevels*sizeof(int)+
               nlevels*sizeof(int);
    int i;

    if (subsampling!=2&&subsampling!=4&&
        subsampling!=8&&subsampling!=16&&subsampling!=32) {
        trace(2,"Pyramid's subsampling must be either 2, 4, 8, 16, or 32\n");
    }
    /* Allocate memory for structure and set parameters */
    pyramid=(ppyramid_t)malloc(nbytes);
    if (pyramid==NULL) {
        trace(2,"(klt_create_pyramid) Out of memory");
    }
    /* set parameters */
    pyramid->subsampling=subsampling;
    pyramid->nLevels=nlevels;
    pyramid->img  =(pfloatimg_t*)(pyramid+1);
    pyramid->ncols=(int*)(pyramid->img+nlevels);
    pyramid->nrows=(int*)(pyramid->ncols+nlevels);

    /* Allocate memory for each level of pyramid and assign pointers */
    for (i=0;i<nlevels;i++)  {
        pyramid->img[i]  =create_floatimg(ncols,nrows);
        pyramid->ncols[i]=ncols;
        pyramid->nrows[i]=nrows;
        ncols/=subsampling;
        nrows/=subsampling;
    }
    return pyramid;
}
/* compute pyramid of image---------------------------------------------------
 * args:    pfloatimg_t *img     IO  image float data
 *          ppyramid_t *pyramid  IO  pyramid of image data
 *          float sigma_fact     I   sigma fact for processing image data
 * return: none
 * ---------------------------------------------------------------------------*/
extern void klt_compute_pyramid(pfloatimg_t img,ppyramid_t pyramid,
                                float sigma_fact)
{
    pfloatimg_t currimg,tmpimg;
    int ncols=img->ncols,nrows=img->nrows;
    int subsampling=pyramid->subsampling;
    int subhalf=subsampling/2;
    float sigma=subsampling*sigma_fact;  /* empirically determined */
    int oldncols;
    int i,x,y;

    if (subsampling!=2&&subsampling!=4&&subsampling!=8&&subsampling!=16&&subsampling!=32){
        trace(2,"(_KLTComputePyramid)  Pyramid's subsampling must "
                "be either 2, 4, 8, 16, or 32");
        return;
    }
    assert(pyramid->ncols[0]==img->ncols);
    assert(pyramid->nrows[0]==img->nrows);

    /* Copy original image to level 0 of pyramid */
    memcpy(pyramid->img[0]->data,img->data,ncols*nrows*sizeof(float));

    currimg=img;
    for (i=1;i<pyramid->nLevels;i++) {
        tmpimg=create_floatimg(ncols,nrows);
        compute_smoothedimg(currimg,sigma,tmpimg);

        /* Subsample */
        oldncols=ncols;
        ncols/=subsampling; nrows/=subsampling;
        for (y=0;y<nrows;y++)
            for (x=0;x<ncols;x++) {
                pyramid->img[i]->data[y*ncols+x]=tmpimg->data[(subsampling*y+subhalf)*oldncols+(subsampling*x+subhalf)];
            }
        /* Reassign current image */
        currimg=pyramid->img[i];
        free_floatimg(tmpimg);
    }
    return;
}





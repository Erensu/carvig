/*------------------------------------------------------------------------------
* vision.h : computer vision header file
*
* version : $Revision:$ $Date:$
* history : 2018/11/17 1.0  new (first version)
*-----------------------------------------------------------------------------*/
#ifndef VISION_H
#define VISION_H
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifdef __cplusplus
extern "C"{
#endif

/* constants -----------------------------------------------------------------*/
#define KLT_BOOL int

#define KLT_TRACKED           0
#define KLT_NOT_FOUND        -1
#define KLT_SMALL_DET        -2
#define KLT_MAX_ITERATIONS   -3
#define KLT_OOB              -4
#define KLT_LARGE_RESIDUE    -5
#define MAX_KERNEL_WIDTH 	 71

#define SWAP_ME(X,Y) {temp=(X);(X)=(Y);(Y)=temp;}

/* type definitions ----------------------------------------------------------*/
typedef struct  {                /* float image data type */
    int ncols;                   /* number of image data rows */
    int nrows;                   /* number of image data cols */
    float *data;                 /* float image data */
} floatimg_t,*pfloatimg_t;

typedef struct {                 /* feature point data type */
    float x,y;                   /* image coordinate of feature point */
    int valid;                   /* valid flag for feature point */

    /* for affine mapping */
    pfloatimg_t aff_img;
    pfloatimg_t aff_img_gradx;
    pfloatimg_t aff_img_grady;
    float aff_x,aff_y;
    float aff_Axx,aff_Ayx;
    float aff_Axy,aff_Ayy;
} feature_t,*pfeature_t;

typedef struct  {                /* KLT feature track context type */
    int mindist;		    	 /* min distance b/w features */
    int window_width,window_height;
    KLT_BOOL sequentialMode;	 /* whether to save most recent image to save time */

    /* can set to TRUE manually, but don't set to */
    /* FALSE manually */
    KLT_BOOL smoothBeforeSelecting;	/* whether to smooth image before */
    KLT_BOOL writeInternalImages;	/* selecting features: whether to write internal images */
    KLT_BOOL lighting_insensitive;  /* tracking features:whether to normalize for gain and bias (not in original algorithm) */

    int min_eigenvalue;		        /* smallest eigenvalue allowed for selecting */
    float min_determinant;	        /* th for determining lost */
    float min_displacement;	        /* th for stopping tracking when pixel changes little */
    int max_iterations;		        /* th for stopping tracking when too many iterations */
    float max_residue;		        /* th for stopping tracking when residue is large */
    float grad_sigma;
    float smooth_sigma_fact;
    float pyramid_sigma_fact;
    float step_factor;              /* size of Newton steps; 2.0 comes from equations, 1.0 seems to avoid overshooting */
    int nSkippedPixels;		        /* # of pixels skipped when finding features */
    int borderx;			        /* border in which features will not be found */
    int bordery;
    int nPyramidLevels;	          	/* computed from search_ranges */
    int subsampling;		        /* 		" */

    /* for affine mapping */
    int affine_window_width, affine_window_height;
    int affineConsistencyCheck; /* whether to evaluates the consistency of features with affine mapping
                                   -1 = don't evaluates the consistency
                                    0 = evaluates the consistency of features with translation mapping
                                    1 = evaluates the consistency of features with similarity mapping
                                    2 = evaluates the consistency of features with affine mapping
                                */
    int affine_max_iterations;
    float affine_max_residue;
    float affine_min_displacement;
    float affine_max_displacement_differ; /* th for the difference between the displacement calculated
                                             by the affine tracker and the frame to frame tracker in pel
                                             */
    /* User must not touch these */
    void *pyramid_last;
    void *pyramid_last_gradx;
    void *pyramid_last_grady;
} tracking_context_t,*ptracking_context_t;

typedef struct  {                    /* feature list data type */
    int n;                           /* number of feature points */
    pfeature_t *feature;             /* feature list data */
} featurelist_t,*pfeaturelist_t;

typedef struct  {
    int nFrames;
    pfeature_t *feature;
} feature_history_t,*pfeature_history_t;

typedef struct  {
    int nFrames;
    int nFeatures;
    pfeature_t **feature;
} feature_table_t,*pfeature_table_t;

typedef struct  {
    int subsampling;
    int nLevels;
    pfloatimg_t *img;
    int *ncols,*nrows;
} pyramid_t,*ppyramid_t;

/* functions -----------------------------------------------------------------*/
extern tracking_context_t* klt_create_tracking_context();
extern featurelist_t* klt_create_featurelist(int nFeatures);
extern feature_history_t* klt_create_feature_history(int nFrames);
extern feature_table_t* klt_create_feature_table(int nFrames,int nFeatures);

extern floatimg_t* create_floatimg(int ncol,int nrow);

extern void klt_free_track_context(tracking_context_t *context);
extern void klt_free_featurelist(featurelist_t *flist);
extern void klt_free_featurehist(feature_history_t *fhist);
extern void klt_free_featuretable(feature_table_t *table);

extern void free_floatimg(floatimg_t* floatimg);
extern void free_pyramid(pyramid_t *pyramid);

extern void klt_select_goodfeatures(tracking_context_t* tc,unsigned char *img,
                                    int ncols,int nrows,
                                    featurelist_t *fl);
extern int klt_count_remaining_features(pfeaturelist_t fl);
extern void klt_replace_lost_features(tracking_context_t *tc,unsigned char *img,
                                      int ncols,int nrows,pfeaturelist_t fl);

extern void tofloatImage(unsigned char *img,int ncols, int nrows,
                         floatimg_t* floatimg);
extern void compute_gradients(floatimg_t *img, float sigma,
                              floatimg_t *gradx, floatimg_t *grady);

extern void compute_smoothedimg(pfloatimg_t img,float sigma,
                                pfloatimg_t smooth);
extern float ComputeSmoothSigma(tracking_context_t* tc);

extern void writeFloatImageToPGM(pfloatimg_t img,char *filename);

extern unsigned char* pgmReadFile(char *fname,unsigned char *img,
                                  int *ncols,int *nrows);
extern void pgmWriteFile(char *fname,unsigned char *img,
                         int ncols,int nrows);

extern void ppmWriteFileRGB(char *fname,unsigned char *redimg,
                            unsigned char *greenimg,unsigned char *blueimg,
                            int ncols,int nrows);

extern unsigned char* pgmRead(FILE *fp,unsigned char *img,
                              int *ncols, int *nrows);
extern void pgmWrite(FILE *fp,unsigned char *img, int ncols,int nrows);
extern void ppmWrite(FILE *fp,unsigned char *redimg,
                     unsigned char *greenimg,unsigned char *blueimg,
                     int ncols,int nrows);

extern ppyramid_t klt_create_pyramid(int ncols,int nrows,int subsampling,
                                     int nlevels);
extern void klt_compute_pyramid(pfloatimg_t img,ppyramid_t pyramid,
                                float sigma_fact);
extern void klt_track_features(tracking_context_t* tc,
                               unsigned char *img1,unsigned char *img2,
                               int ncols,int nrows,
                               pfeaturelist_t featurelist);

extern void writeFeatureListToPPM(pfeaturelist_t featurelist,
                                  unsigned char *greyimg,int ncols,int nrows,
                                  char *filename);

#ifdef __cplusplus
}
#endif
#endif
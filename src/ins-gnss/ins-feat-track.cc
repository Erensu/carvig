/*------------------------------------------------------------------------------
 * ins-feat-track.cc : klt feature point track functions
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

/*-----------------------------------------------------------------------------*/
static KLT_BOOL outOfBounds(float x,float y,int ncols,int nrows,int borderx,
                            int bordery)
{
    return (x<borderx||x>ncols-1-borderx||y<bordery||y>nrows-1-bordery);
}
/*----------------------------------------------------------------------------*/
static void am_getSubFloatImage(pfloatimg_t img,float x, float y,
                                pfloatimg_t window)
{
    register int hw=window->ncols/2, hh = window->nrows/2;
    int x0=(int)x;
    int y0=(int)y;
    float *windata=window->data;
    int offset;
    register int i,j;

    assert(x0-hw>=0);
    assert(y0-hh>=0);
    assert(x0+hw<=img->ncols);
    assert(y0+hh<=img->nrows);

    /* copy values */
    for (j=-hh;j<=hh;j++)
        for (i=-hw;i<=hw;i++)  {
            offset=(j+y0)*img->ncols+(i+x0);
            *windata++=*(img->data+offset);
        }
}
/*----------------------------------------------------------------------------*/
static float* allocateFloatWindow(int width,int height)
{
    float* fw;

    fw=(float*)malloc(width*height*sizeof(float));
    if (fw==NULL) {
        trace(2,"(_allocateFloatWindow) Out of memory.\n");
    }
    return fw;
}
/*----------------------------------------------------------------------------*/
static float **am_matrix(long nr, long nc)
{
    float **m;
    int a;
    m=(float**)malloc((size_t)(nr*sizeof(float*)));
    m[0]=(float*)malloc((size_t)((nr*nc)*sizeof(float)));
    for (a=1;a<nr;a++) m[a]=m[a-1]+nc;
    return m;
}
/*----------------------------------------------------------------------------*/
static void am_free_matrix(float **m)
{
    free(m[0]); free(m);
}
/*----------------------------------------------------------------------------
 * Solves the 2x2 matrix equation
 *         [gxx gxy] [dx] = [ex]
 *         [gxy gyy] [dy] = [ey]
 * for dx and dy.
 *
 * Returns KLT_TRACKED on success and KLT_SMALL_DET on failure
 *----------------------------------------------------------------------------*/
static int solveEquation(float gxx, float gxy, float gyy,float ex, float ey,
                         float small,float *dx,float *dy)
{
    float det=gxx*gyy-gxy*gxy;
    if (det<small) return KLT_SMALL_DET;

    *dx=(gyy*ex-gxy*ey)/det;
    *dy=(gxx*ey-gxy*ex)/det;
    return KLT_TRACKED;
}
/*----------------------------------------------------------------------------*/
static float sumAbsFloatWindow(float* fw,int width,int height)
{
    float sum=0.0;
    int w;

    for (;height>0;height--) {
        for (w=0;w<width;w++) sum+=(float)fabs(*fw++);
    }
    return sum;
}
/*----------------------------------------------------------------------------
 * Given a point (x,y) in an image, computes the bilinear interpolated
 * gray-level value of the point in the image.
 *----------------------------------------------------------------------------*/
static float interpolate(float x,float y,pfloatimg_t img)
{
    int xt=(int)x;  /* coordinates of top-left corner */
    int yt=(int)y;
    float ax=x-xt;
    float ay=y-yt;
    float *ptr=img->data+(img->ncols*yt)+xt;

    assert (xt>=0&&yt>=0&&xt<=img->ncols-2&&yt<=img->nrows-2);
    return ((1-ax)*(1-ay)**ptr+ax*(1-ay)**(ptr+1)+
            (1-ax)*ay**(ptr+(img->ncols))+ax*ay**(ptr+(img->ncols)+1));
}
/*----------------------------------------------------------------------------
 * aligns the gradients with the affine transformed window
 *----------------------------------------------------------------------------*/
static void am_getGradientWinAffine(
        pfloatimg_t in_gradx,
        pfloatimg_t in_grady,
        float x, float y,      /* center of window*/
        float Axx, float Ayx , float Axy, float Ayy,  /* affine mapping */
        int width, int height, /* size of window */
        float* out_gradx,      /* output */
        float* out_grady)      /* output */
{
    register int hw=width/2,hh=height/2;
    register int i,j;
    float mi,mj;

    /* Compute values */
    for (j=-hh;j<=hh;j++)
        for (i=-hw;i<=hw;i++)  {
            mi=Axx*i+Axy*j;
            mj=Ayx*i+Ayy*j;
            *out_gradx++=interpolate(x+mi,y+mj,in_gradx);
            *out_grady++=interpolate(x+mi,y+mj,in_grady);
        }
}
/*----------------------------------------------------------------------------*/
static int am_gauss_jordan_elimination(float **a, int n, float **b, int m)
{
    /* re-implemented from Numerical Recipes in C */
    int *indxc,*indxr,*ipiv;
    int i,j,k,l,ll;
    float big,dum,pivinv,temp;
    int col = 0;
    int row = 0;

    indxc=(int *)malloc((size_t) (n*sizeof(int)));
    indxr=(int *)malloc((size_t) (n*sizeof(int)));
    ipiv=(int *)malloc((size_t) (n*sizeof(int)));
    for (j=0;j<n;j++) ipiv[j]=0;
    for (i=0;i<n;i++) {
        big=0.0;
        for (j=0;j<n;j++)
            if (ipiv[j] != 1)
                for (k=0;k<n;k++) {
                    if (ipiv[k] == 0) {
                        if (fabs(a[j][k]) >= big) {
                            big= (float) fabs(a[j][k]);
                            row=j;
                            col=k;
                        }
                    } else if (ipiv[k] > 1) return KLT_SMALL_DET;
                }
        ++(ipiv[col]);
        if (row != col) {
            for (l=0;l<n;l++) SWAP_ME(a[row][l],a[col][l])
            for (l=0;l<m;l++) SWAP_ME(b[row][l],b[col][l])
        }
        indxr[i]=row;
        indxc[i]=col;
        if (a[col][col] == 0.0) return KLT_SMALL_DET;
        pivinv=1.0f/a[col][col];
        a[col][col]=1.0;
        for (l=0;l<n;l++) a[col][l] *= pivinv;
        for (l=0;l<m;l++) b[col][l] *= pivinv;
        for (ll=0;ll<n;ll++)
            if (ll != col) {
                dum=a[ll][col];
                a[ll][col]=0.0;
                for (l=0;l<n;l++) a[ll][l] -= a[col][l]*dum;
                for (l=0;l<m;l++) b[ll][l] -= b[col][l]*dum;
            }
    }
    for (l=n-1;l>=0;l--) {
        if (indxr[l] != indxc[l])
            for (k=0;k<n;k++)
                SWAP_ME(a[k][indxr[l]],a[k][indxc[l]]);
    }
    free(ipiv);
    free(indxr);
    free(indxc);

    return KLT_TRACKED;
}
/*----------------------------------------------------------------------------
 * Given two images and the window center in both images,
 * aligns the images with the window and computes the difference
 * between the two overlaid images using the affine mapping.
 *       A =  [ Axx Axy]
 *            [ Ayx Ayy]
 *----------------------------------------------------------------------------*/
static void am_computeIntensityDifferenceAffine(
        pfloatimg_t img1,        /* images */
        pfloatimg_t img2,
        float x1, float y1,      /* center of window in 1st img */
        float x2, float y2,      /* center of window in 2nd img */
        float Axx, float Ayx , float Axy, float Ayy,    /* affine mapping */
        int width, int height,   /* size of window */
        float* imgdiff)          /* output */
{
    register int hw=width/2,hh=height/2;
    float g1,g2;
    register int i,j;
    float mi,mj;

    /* Compute values */
    for (j=-hh;j<=hh;j++)
        for (i=-hw;i<=hw;i++) {
            g1=interpolate(x1+i,y1+j,img1);
            mi=Axx*i+Axy*j;
            mj=Ayx*i+Ayy*j;
            g2=interpolate(x2+mi,y2+mj,img2);
            *imgdiff++=g1-g2;
        }
}
/*----------------------------------------------------------------------------*/
static void am_compute6by6GradientMatrix(float* gradx,float* grady,
                                         int width,  int height,float **T)
{
    register int hw = width/2, hh = height/2;
    register int i, j;
    float gx, gy, gxx, gxy, gyy,  x, y, xx, xy, yy;

    /* Set values to zero */
    for (j = 0 ; j < 6 ; j++)  {
        for (i = j ; i < 6 ; i++)  {
            T[j][i] = 0.0;
        }
    }
    for (j = -hh ; j <= hh ; j++) {
        for (i = -hw ; i <= hw ; i++)  {
            gx = *gradx++;
            gy = *grady++;
            gxx = gx * gx;
            gxy = gx * gy;
            gyy = gy * gy;
            x = (float) i;
            y = (float) j;
            xx = x * x;
            xy = x * y;
            yy = y * y;

            T[0][0] += xx * gxx;
            T[0][1] += xx * gxy;
            T[0][2] += xy * gxx;
            T[0][3] += xy * gxy;
            T[0][4] += x  * gxx;
            T[0][5] += x  * gxy;

            T[1][1] += xx * gyy;
            T[1][2] += xy * gxy;
            T[1][3] += xy * gyy;
            T[1][4] += x  * gxy;
            T[1][5] += x  * gyy;

            T[2][2] += yy * gxx;
            T[2][3] += yy * gxy;
            T[2][4] += y  * gxx;
            T[2][5] += y  * gxy;

            T[3][3] += yy * gyy;
            T[3][4] += y  * gxy;
            T[3][5] += y  * gyy;

            T[4][4] += gxx;
            T[4][5] += gxy;

            T[5][5] += gyy;
        }
    }
    for (j = 0 ; j < 5 ; j++)  {
        for (i = j+1 ; i < 6 ; i++)  {
            T[i][j] = T[j][i];
        }
    }
}
/*----------------------------------------------------------------------------*/
static void am_compute4by1ErrorVector(float* imgdiff,float* gradx,float* grady,
                                      int width,int height,float **e)
{
    register int hw = width/2, hh = height/2;
    register int i, j;
    register float diff,  diffgradx,  diffgrady;

    /* Set values to zero */
    for(i = 0; i < 4; i++) e[i][0] = 0.0;

    /* Compute values */
    for (j = -hh ; j <= hh ; j++) {
        for (i = -hw ; i <= hw ; i++)  {
            diff = *imgdiff++;
            diffgradx = diff * (*gradx++);
            diffgrady = diff * (*grady++);
            e[0][0] += diffgradx * i + diffgrady * j;
            e[1][0] += diffgrady * i - diffgradx * j;
            e[2][0] += diffgradx;
            e[3][0] += diffgrady;
        }
    }
    for(i = 0; i < 4; i++) e[i][0] *= 0.5;
}
/*----------------------------------------------------------------------------*/
static void am_compute6by1ErrorVector(float* imgdiff,float* gradx,float* grady,
                                      int width,int height,float **e)
{
    register int hw = width/2, hh = height/2;
    register int i, j;
    register float diff,  diffgradx,  diffgrady;

    /* Set values to zero */
    for(i = 0; i < 6; i++) e[i][0] = 0.0;

    /* Compute values */
    for (j = -hh ; j <= hh ; j++) {
        for (i = -hw ; i <= hw ; i++)  {
            diff = *imgdiff++;
            diffgradx = diff * (*gradx++);
            diffgrady = diff * (*grady++);
            e[0][0] += diffgradx * i;
            e[1][0] += diffgrady * i;
            e[2][0] += diffgradx * j;
            e[3][0] += diffgrady * j;
            e[4][0] += diffgradx;
            e[5][0] += diffgrady;
        }
    }
    for(i = 0; i < 6; i++) e[i][0] *= 0.5;
}
/*----------------------------------------------------------------------------
 * Given two images and the window center in both images,
 * aligns the images wrt the window and computes the difference
 * between the two overlaid images; normalizes for overall gain and bias.
 *----------------------------------------------------------------------------*/
static void computeIntensityDifferenceLightingInsensitive(
        pfloatimg_t img1,       /* images */
        pfloatimg_t img2,
        float x1, float y1,     /* center of window in 1st img */
        float x2, float y2,     /* center of window in 2nd img */
        int width, int height,  /* size of window */
        float* imgdiff)         /* output */
{
    register int hw = width/2, hh = height/2;
    float g1, g2, sum1_squared = 0, sum2_squared = 0;
    register int i, j;

    float sum1 = 0, sum2 = 0;
    float mean1, mean2,alpha,belta;
    /* Compute values */
    for (j = -hh ; j <= hh ; j++)
        for (i = -hw ; i <= hw ; i++)  {
            g1 = interpolate(x1+i, y1+j, img1);
            g2 = interpolate(x2+i, y2+j, img2);
            sum1 += g1;    sum2 += g2;
            sum1_squared += g1*g1;
            sum2_squared += g2*g2;
        }
    mean1=sum1_squared/(width*height);
    mean2=sum2_squared/(width*height);
    alpha = (float) sqrt(mean1/mean2);
    mean1=sum1/(width*height);
    mean2=sum2/(width*height);
    belta = mean1-alpha*mean2;

    for (j = -hh ; j <= hh ; j++)
        for (i = -hw ; i <= hw ; i++)  {
            g1 = interpolate(x1+i, y1+j, img1);
            g2 = interpolate(x2+i, y2+j, img2);
            *imgdiff++ = g1- g2*alpha-belta;
        }
}
/*----------------------------------------------------------------------------
 * Given two gradients and the window center in both images,
 * aligns the gradients wrt the window and computes the sum of the two
 * overlaid gradients; normalizes for overall gain and bias.
 *----------------------------------------------------------------------------*/
static void computeGradientSumLightingInsensitive(
        pfloatimg_t gradx1,      /* gradient images */
        pfloatimg_t grady1,
        pfloatimg_t gradx2,
        pfloatimg_t grady2,
        pfloatimg_t img1,        /* images */
        pfloatimg_t img2,

        float x1, float y1,      /* center of window in 1st img */
        float x2, float y2,      /* center of window in 2nd img */
        int width, int height,   /* size of window */
        float* gradx,            /* output */
        float* grady)     
{
    register int hw = width/2, hh = height/2;
    float g1, g2, sum1_squared = 0, sum2_squared = 0;
    register int i, j;

    float sum1 = 0, sum2 = 0;
    float mean1, mean2, alpha;
    for (j = -hh ; j <= hh ; j++)
        for (i = -hw ; i <= hw ; i++)  {
            g1 = interpolate(x1+i, y1+j, img1);
            g2 = interpolate(x2+i, y2+j, img2);
            sum1_squared += g1;    sum2_squared += g2;
        }
    mean1 = sum1_squared/(width*height);
    mean2 = sum2_squared/(width*height);
    alpha = (float) sqrt(mean1/mean2);

    /* Compute values */
    for (j = -hh ; j <= hh ; j++)
        for (i = -hw ; i <= hw ; i++)  {
            g1 = interpolate(x1+i, y1+j, gradx1);
            g2 = interpolate(x2+i, y2+j, gradx2);
            *gradx++ = g1 + g2*alpha;
            g1 = interpolate(x1+i, y1+j, grady1);
            g2 = interpolate(x2+i, y2+j, grady2);
            *grady++ = g1+ g2*alpha;
        }
}
/*----------------------------------------------------------------------------
 * Given two images and the window center in both images,
 * aligns the images wrt the window and computes the difference
 * between the two overlaid images.
 *----------------------------------------------------------------------------*/
static void computeIntensityDifference(
        pfloatimg_t img1,       /* images */
        pfloatimg_t img2,
        float x1, float y1,     /* center of window in 1st img */
        float x2, float y2,     /* center of window in 2nd img */
        int width, int height,  /* size of window */
        float* imgdiff)         /* output */
{
    register int hw = width/2, hh = height/2;
    float g1, g2;
    register int i, j;

    /* Compute values */
    for (j = -hh ; j <= hh ; j++)
        for (i = -hw ; i <= hw ; i++)  {
            g1 = interpolate(x1+i, y1+j, img1);
            g2 = interpolate(x2+i, y2+j, img2);
            *imgdiff++ = g1 - g2;
        }
}
/*----------------------------------------------------------------------------
 * Given two gradients and the window center in both images,
 * aligns the gradients wrt the window and computes the sum of the two
 * overlaid gradients.
 *----------------------------------------------------------------------------*/
static void computeGradientSum(
        pfloatimg_t gradx1,      /* gradient images */
        pfloatimg_t grady1,
        pfloatimg_t gradx2,
        pfloatimg_t grady2,
        float x1, float y1,      /* center of window in 1st img */
        float x2, float y2,      /* center of window in 2nd img */
        int width, int height,   /* size of window */
        float* gradx,      /* output */
        float* grady)      /*   " */
{
    register int hw = width/2, hh = height/2;
    float g1, g2;
    register int i, j;

    /* Compute values */
    for (j = -hh ; j <= hh ; j++)
        for (i = -hw ; i <= hw ; i++)  {
            g1 = interpolate(x1+i, y1+j, gradx1);
            g2 = interpolate(x2+i, y2+j, gradx2);
            *gradx++ = g1 + g2;
            g1 = interpolate(x1+i, y1+j, grady1);
            g2 = interpolate(x2+i, y2+j, grady2);
            *grady++ = g1 + g2;
        }
}
/*----------------------------------------------------------------------------*/
static void _am_compute4by4GradientMatrix(float* gradx,float* grady,int width,
                                          int height,float **T)
{
    register int hw = width/2, hh = height/2;
    register int i, j;
    float gx, gy, x, y;

    /* Set values to zero */
    for (j = 0 ; j < 4 ; j++)  {
        for (i = 0 ; i < 4 ; i++)  {
            T[j][i] = 0.0;
        }
    }
    for (j = -hh ; j <= hh ; j++) {
        for (i = -hw ; i <= hw ; i++)  {
            gx = *gradx++;
            gy = *grady++;
            x = (float) i;
            y = (float) j;
            T[0][0] += (x*gx+y*gy) * (x*gx+y*gy);
            T[0][1] += (x*gx+y*gy)*(x*gy-y*gx);
            T[0][2] += (x*gx+y*gy)*gx;
            T[0][3] += (x*gx+y*gy)*gy;

            T[1][1] += (x*gy-y*gx) * (x*gy-y*gx);
            T[1][2] += (x*gy-y*gx)*gx;
            T[1][3] += (x*gy-y*gx)*gy;

            T[2][2] += gx*gx;
            T[2][3] += gx*gy;

            T[3][3] += gy*gy;
        }
    }
    for (j = 0 ; j < 3 ; j++)  {
        for (i = j+1 ; i < 4 ; i++)  {
            T[i][j] = T[j][i];
        }
    }
}
/*----------------------------------------------------------------------------*/
static void compute2by1ErrorVector(
        float* imgdiff,
        float* gradx,
        float* grady,
        int width,         /* size of window */
        int height,
        float step_factor, /* 2.0 comes from equations, 1.0 seems to avoid overshooting */
        float *ex,         /* return values */
        float *ey)
{
    register float diff;
    register int i;

    /* Compute values */
    *ex = 0;  *ey = 0;
    for (i = 0 ; i < width * height ; i++)  {
        diff = *imgdiff++;
        *ex += diff * (*gradx++);
        *ey += diff * (*grady++);
    }
    *ex *= step_factor;
    *ey *= step_factor;
}
/*----------------------------------------------------------------------------*/
static void compute2by2GradientMatrix(
        float* gradx,
        float* grady,
        int width,   /* size of window */
        int height,
        float *gxx,  /* return values */
        float *gxy,
        float *gyy)

{
    register float gx,gy;
    register int i;

    /* Compute values */
    *gxx=0.0;*gxy=0.0;*gyy=0.0;
    for (i=0;i<width*height;i++) {
        gx=*gradx++;
        gy=*grady++;
        *gxx+=gx*gx;
        *gxy+=gx*gy;
        *gyy+=gy*gy;
    }
}
static void am_compute4by4GradientMatrix(
        float* gradx,
        float* grady,
        int width,   /* size of window */
        int height,
        float **T)   /* return values */
{
    register int hw = width/2, hh = height/2;
    register int i, j;
    float gx, gy, x, y;

    /* Set values to zero */
    for (j = 0 ; j < 4 ; j++)  {
        for (i = 0 ; i < 4 ; i++)  {
            T[j][i] = 0.0;
        }
    }
    for (j = -hh ; j <= hh ; j++) {
        for (i = -hw ; i <= hw ; i++)  {
            gx = *gradx++;
            gy = *grady++;
            x = (float) i;
            y = (float) j;
            T[0][0] += (x*gx+y*gy) * (x*gx+y*gy);
            T[0][1] += (x*gx+y*gy)*(x*gy-y*gx);
            T[0][2] += (x*gx+y*gy)*gx;
            T[0][3] += (x*gx+y*gy)*gy;

            T[1][1] += (x*gy-y*gx) * (x*gy-y*gx);
            T[1][2] += (x*gy-y*gx)*gx;
            T[1][3] += (x*gy-y*gx)*gy;

            T[2][2] += gx*gx;
            T[2][3] += gx*gy;

            T[3][3] += gy*gy;
        }
    }
    for (j = 0 ; j < 3 ; j++)  {
        for (i = j+1 ; i < 4 ; i++)  {
            T[i][j] = T[j][i];
        }
    }
}
/*----------------------------------------------------------------------------
 * Tracks a feature point from the image of first occurrence to the actual image.
 *
 * RETURNS
 * KLT_SMALL_DET or KLT_LARGE_RESIDUE or KLT_OOB if feature is lost,
 * KLT_TRACKED otherwise.
 *----------------------------------------------------------------------------*/
static int am_trackFeatureAffine(
               float x1,  /* location of window in first image */
               float y1,
               float *x2, /* starting location of search in second image */
               float *y2,
               pfloatimg_t img1,
               pfloatimg_t gradx1,
               pfloatimg_t grady1,
               pfloatimg_t img2,
               pfloatimg_t gradx2,
               pfloatimg_t grady2,
               int width,           /* size of window */
               int height,
               float step_factor,   /* 2.0 comes from equations, 1.0 seems to avoid overshooting */
               int max_iterations,
               float small,         /* determinant threshold for declaring KLT_SMALL_DET */
               float th,            /* displacement threshold for stopping  */
               float th_aff,
               float max_residue,   /* residue threshold for declaring KLT_LARGE_RESIDUE */
               int lighting_insensitive,  /* whether to normalize for gain and bias */
               int affine_map,      /* whether to evaluates the consistency of features with affine mapping */
               float mdd,           /* difference between the displacements */
               float *Axx, float *Ayx,
               float *Axy, float *Ayy) /* used affine mapping */
{
    float *imgdiff,*gradx,*grady;
    float gxx,gxy,gyy,ex,ey,dx,dy;
    int iteration=0;
    int status=0;
    int hw=width/2;
    int hh=height/2;
    int nc1=img1->ncols;
    int nr1=img1->nrows;
    int nc2=img2->ncols;
    int nr2=img2->nrows;
    float **a;
    float **T;
    float one_plus_eps=1.001f;   /* To prevent rounding errors */
    float old_x2=*x2;
    float old_y2=*y2;
    KLT_BOOL convergence=false;

    trace(3,"am_trackFeatureAffine:\n");

    /* Allocate memory for windows */
    imgdiff=allocateFloatWindow(width,height);
    gradx  =allocateFloatWindow(width,height);
    grady  =allocateFloatWindow(width,height);
    T=am_matrix(6,6);
    a=am_matrix(6,1);

    /* Iteratively update the window position */
    do  {
        if(!affine_map) {
            /* pure translation tracker */

            /* If out of bounds, exit loop */
            if (x1-hw <0.0f||nc1-( x1+hw)<one_plus_eps||
                *x2-hw<0.0f||nc2-(*x2+hw)<one_plus_eps||
                y1-hh <0.0f||nr1-( y1+hh)<one_plus_eps||
                *y2-hh<0.0f||nr2-(*y2+hh)<one_plus_eps) {
                status=KLT_OOB;
                break;
            }
            /* Compute gradient and difference windows */
            if (lighting_insensitive) {
                computeIntensityDifferenceLightingInsensitive(img1, img2,x1,y1,*x2,*y2,
                                                              width,height,imgdiff);
                computeGradientSumLightingInsensitive(gradx1,grady1,gradx2,grady2,
                                                      img1, img2,x1,y1,*x2,*y2,
                                                      width,height,gradx,grady);
            } else {
                computeIntensityDifference(img1,img2,x1,y1,*x2,*y2,
                                           width,height,imgdiff);
                computeGradientSum(gradx1,grady1,gradx2,grady2,
                                   x1,y1,*x2,*y2,width,height,gradx,grady);
            }
            /* Use these windows to construct matrices */
            compute2by2GradientMatrix(gradx,grady,width,height,
                                      &gxx,&gxy,&gyy);
            compute2by1ErrorVector(imgdiff,gradx,grady,width,height,step_factor,
                                   &ex,&ey);

            /* Using matrices, solve equation for new displacement */
            status=solveEquation(gxx,gxy,gyy,ex,ey,small,&dx,&dy);

            convergence=(fabs(dx)<th&&fabs(dy)<th);

            *x2+=dx;
            *y2+=dy;

        }
        else {
            /* affine tracker */
            float ul_x=*Axx*(-hw)+*Axy*  hh+*x2;  /* upper left corner */
            float ul_y=*Ayx*(-hw)+*Ayy*  hh+*y2;
            float ll_x=*Axx*(-hw)+*Axy*(-hh)+*x2;  /* lower left corner */
            float ll_y=*Ayx*(-hw)+*Ayy*(-hh)+*y2;
            float ur_x=*Axx*hw+*Axy*  hh +*x2;  /* upper right corner */
            float ur_y=*Ayx*hw+*Ayy*  hh +*y2;
            float lr_x=*Axx*hw+*Axy*(-hh)+*x2;  /* lower right corner */
            float lr_y=*Ayx*hw+*Ayy*(-hh)+*y2;

            /* If out of bounds, exit loop */
            if ( x1-hw<0.0f||nc1-(x1+hw)<one_plus_eps||
                 y1-hh<0.0f||nr1-(y1+hh)<one_plus_eps||
                 ul_x <0.0f||nc2-(ul_x )<one_plus_eps||
                 ll_x <0.0f||nc2-(ll_x )<one_plus_eps||
                 ur_x <0.0f||nc2-(ur_x )<one_plus_eps||
                 lr_x <0.0f||nc2-(lr_x )<one_plus_eps||
                 ul_y <0.0f||nr2-(ul_y )<one_plus_eps||
                 ll_y <0.0f||nr2-(ll_y )<one_plus_eps||
                 ur_y <0.0f||nr2-(ur_y )<one_plus_eps||
                 lr_y <0.0f||nr2-(lr_y )<one_plus_eps) {
                status=KLT_OOB;
                break;
            }
            am_computeIntensityDifferenceAffine(img1,img2,x1,y1,*x2,*y2,*Axx,*Ayx,*Axy,*Ayy,
                                                width,height,imgdiff);

            am_getGradientWinAffine(gradx2,grady2,*x2,*y2,*Axx,*Ayx,*Axy,*Ayy,
                                    width, height,gradx,grady);

            switch (affine_map) {
                case 1:
                    am_compute4by1ErrorVector(imgdiff, gradx,grady,width, height,a);
                    am_compute4by4GradientMatrix(gradx,grady,width,height,T);
                    status=am_gauss_jordan_elimination(T,4,a,1);

                    *Axx+=a[0][0];
                    *Ayx+=a[1][0];
                    *Ayy=*Axx;
                    *Axy=-(*Ayx);

                    dx=a[2][0];
                    dy=a[3][0];

                    break;
                case 2:
                    am_compute6by1ErrorVector(imgdiff,gradx,grady,width,height,a);
                    am_compute6by6GradientMatrix(gradx,grady,width,height,T);

                    status=am_gauss_jordan_elimination(T,6,a,1);

                    *Axx+=a[0][0];
                    *Ayx+=a[1][0];
                    *Axy+=a[2][0];
                    *Ayy+=a[3][0];

                    dx=a[4][0];
                    dy=a[5][0];
                    break;
            }
            *x2+=dx;
            *y2+=dy;

            /* old upper left corner - new upper left corner */
            ul_x-=*Axx*(-hw)+*Axy*hh+*x2;
            ul_y-=*Ayx*(-hw)+*Ayy*hh+*y2;

            /* old lower left corner - new lower left corner */
            ll_x-=*Axx*(-hw)+*Axy*(-hh)+*x2;
            ll_y-=*Ayx*(-hw)+*Ayy*(-hh)+*y2;

            /* old upper right corner - new upper right corner */
            ur_x-=*Axx*hw+*Axy*hh+*x2;
            ur_y-=*Ayx*hw+*Ayy*hh+*y2;

            /* old lower right corner - new lower right corner */
            lr_x-=*Axx*hw+*Axy*(-hh)+*x2;
            lr_y-=*Ayx*hw+*Ayy*(-hh)+*y2;

            convergence=(fabs(dx)<th&&fabs(dy)<th&&
                         fabs(ul_x)<th_aff&&fabs(ul_y)<th_aff&&
                         fabs(ll_x)<th_aff&&fabs(ll_y)<th_aff&&
                         fabs(ur_x)<th_aff&&fabs(ur_y)<th_aff&&
                         fabs(lr_x)<th_aff&&fabs(lr_y)<th_aff);
        }
        if (status==KLT_SMALL_DET) break;
        iteration++;

    } while(!convergence&&iteration<max_iterations);
    /*}  while ((fabs(dx)>=th||fabs(dy)>=th||(affine_map&&iteration<8))&&iteration<max_iterations); */
    am_free_matrix(T);
    am_free_matrix(a);

    /* Check whether window is out of bounds */
    if (*x2-hw<0.0f||nc2-(*x2+hw)<one_plus_eps ||
        *y2-hh<0.0f||nr2-(*y2+hh)<one_plus_eps)
        status=KLT_OOB;

    /* Check whether feature point has moved to much during iteration*/
    if ((*x2-old_x2)>mdd||(*y2-old_y2)>mdd )
        status = KLT_OOB;

    /* Check whether residue is too large */
    if (status==KLT_TRACKED)  {
        if (!affine_map){
            computeIntensityDifference(img1,img2,x1,y1,*x2,*y2,
                                       width,height,imgdiff);
        } else {
            am_computeIntensityDifferenceAffine(img1,img2,x1,y1,*x2,*y2,*Axx,*Ayx,*Axy,*Ayy,
                                                width,height,imgdiff);
        }
        if (sumAbsFloatWindow(imgdiff,width,height)/(width*height)>max_residue) {
            status=KLT_LARGE_RESIDUE;
        }
    }
    /* Free memory */
    free(imgdiff); free(gradx); free(grady);

    /* Return appropriate value */
    return status;
}
/*----------------------------------------------------------------------------
 * Tracks a feature point from one image to the next.
 *
 * RETURNS
 * KLT_SMALL_DET if feature is lost,
 * KLT_MAX_ITERATIONS if tracking stopped because iterations timed out,
 * KLT_TRACKED otherwise.
 *----------------------------------------------------------------------------*/
static int trackFeature(
        float x1,           /* location of window in first image */
        float y1,
        float *x2,          /* starting location of search in second image */
        float *y2,
        pfloatimg_t img1,
        pfloatimg_t gradx1,
        pfloatimg_t grady1,
        pfloatimg_t img2,
        pfloatimg_t gradx2,
        pfloatimg_t grady2,
        int width,           /* size of window */
        int height,
        float step_factor,   /* 2.0 comes from equations, 1.0 seems to avoid overshooting */
        int max_iterations,
        float small,         /* determinant threshold for declaring KLT_SMALL_DET */
        float th,            /* displacement threshold for stopping               */
        float max_residue,   /* residue threshold for declaring KLT_LARGE_RESIDUE */
        int lighting_insensitive)  /* whether to normalize for gain and bias */
{
    float* imgdiff,*gradx,*grady;
    float gxx,gxy,gyy,ex,ey,dx,dy;
    int iteration=0;
    int status;
    int hw=width/2;
    int hh=height/2;
    int nc=img1->ncols;
    int nr=img1->nrows;
    float one_plus_eps=1.001f;   /* To prevent rounding errors */


    /* Allocate memory for windows */
    imgdiff=allocateFloatWindow(width,height);
    gradx  =allocateFloatWindow(width,height);
    grady  =allocateFloatWindow(width,height);

    /* Iteratively update the window position */
    do  {

        /* If out of bounds, exit loop */
        if (x1-hw<0.0f||nc-( x1+hw)<one_plus_eps||
           *x2-hw<0.0f||nc-(*x2+hw)<one_plus_eps||
            y1-hh<0.0f||nr-( y1+hh)<one_plus_eps||
           *y2-hh<0.0f||nr-(*y2+hh)<one_plus_eps) {
            status=KLT_OOB;
            break;
        }
        /* Compute gradient and difference windows */
        if (lighting_insensitive) {
            computeIntensityDifferenceLightingInsensitive(img1,img2,x1,y1,*x2,*y2,
                                                           width,height,imgdiff);
            computeGradientSumLightingInsensitive(gradx1,grady1,gradx2,grady2,
                                                   img1,img2,x1,y1,*x2,*y2,width,height,gradx,grady);
        } else {
            computeIntensityDifference(img1, img2,x1,y1,*x2,*y2,
                                       width,height,imgdiff);
            computeGradientSum(gradx1,grady1,gradx2,grady2,
                               x1,y1,*x2,*y2,width,height,gradx,grady);
        }
        /* Use these windows to construct matrices */
        compute2by2GradientMatrix(gradx,grady,width,height,
                                  &gxx,&gxy,&gyy);
        compute2by1ErrorVector(imgdiff,gradx,grady,width,height,step_factor,
                               &ex,&ey);

        /* Using matrices, solve equation for new displacement */
        status=solveEquation(gxx,gxy,gyy,ex,ey,small,&dx,&dy);
        if (status==KLT_SMALL_DET) break;

        *x2+=dx;
        *y2+=dy;
        iteration++;

    } while((fabs(dx)>=th||fabs(dy)>=th)&&iteration<max_iterations);

    /* Check whether window is out of bounds */
    if (*x2-hw<0.0f||nc-(*x2+hw)<one_plus_eps||
        *y2-hh<0.0f||nr-(*y2+hh)<one_plus_eps)
        status=KLT_OOB;

    /* Check whether residue is too large */
    if (status == KLT_TRACKED)  {
        if (lighting_insensitive)
            computeIntensityDifferenceLightingInsensitive(img1,img2,x1,y1,*x2,*y2,
                                                          width,height,imgdiff);
        else
            computeIntensityDifference(img1,img2,x1,y1,*x2,*y2,
                                       width,height,imgdiff);
        if (sumAbsFloatWindow(imgdiff,width,height)/(width*height)>max_residue)
            status=KLT_LARGE_RESIDUE;
    }
    /* Free memory */
    free(imgdiff); free(gradx); free(grady);

    /* Return appropriate value */
    if (status==KLT_SMALL_DET) {
        return KLT_SMALL_DET;
    }
    else if (status==KLT_OOB          ) return KLT_OOB;
    else if (status==KLT_LARGE_RESIDUE) return KLT_LARGE_RESIDUE;
    else if (iteration>=max_iterations) return KLT_MAX_ITERATIONS;
    else return KLT_TRACKED;
}
/* KLT tracking feature points-------------------------------------------------
 * tracks feature points from one image to the next.
 * args:    tracking_context_t *tc  IO  tracking context
 *          unsigned char* img1     I   first image data
 *          unsigned char* img2     I   second image data
 *          int ncols,nrows         I   size if image
 *          pfeaturelist_t fl       O   tracked feature points list
 * return: none
 * ---------------------------------------------------------------------------*/
extern void klt_track_features(tracking_context_t* tc,
                               unsigned char *img1,unsigned char *img2,
                               int ncols,int nrows,
                               pfeaturelist_t featurelist)
{
    pfloatimg_t tmpimg, floatimg1, floatimg2;
    ppyramid_t pyramid1,pyramid1_gradx,pyramid1_grady,pyramid2,pyramid2_gradx,pyramid2_grady;
    float subsampling=(float)tc->subsampling;
    float xloc,yloc,xlocout,ylocout;
    int val=-1;
    int indx,r;
    KLT_BOOL floatimg1_created=false;
    int i;

    trace(3,"(KLT) Tracking %d features in a %d by %d image...  ",
          klt_count_remaining_features(featurelist),ncols,nrows);

    /* Check window size (and correct if necessary) */
    if (tc->window_width%2!=1) {
        tc->window_width=tc->window_width+1;
        trace(2,"Tracking context's window width must be odd.  "
                           "Changing to %d.\n",tc->window_width);
    }
    if (tc->window_height%2!=1) {
        tc->window_height=tc->window_height+1;
        trace(2,"Tracking context's window height must be odd.  "
                           "Changing to %d.\n",
              tc->window_height);
    }
    if (tc->window_width<3) {
        tc->window_width=3;
        trace(2,"Tracking context's window width must be at least three.\n"
                           "Changing to %d.\n",tc->window_width);
    }
    if (tc->window_height<3) {
        tc->window_height=3;
        trace(2,"Tracking context's window height must be at least three.\n"
                           "Changing to %d.\n",tc->window_height);
    }
    /* Create temporary image */
    tmpimg=create_floatimg(ncols,nrows);

    /* Process first image by converting to float, smoothing, computing */
    /* pyramid, and computing gradient pyramids
     * */
    if (tc->sequentialMode&&tc->pyramid_last!=NULL) {
        pyramid1=(ppyramid_t)tc->pyramid_last;
        pyramid1_gradx=(ppyramid_t)tc->pyramid_last_gradx;
        pyramid1_grady=(ppyramid_t)tc->pyramid_last_grady;
        if (pyramid1->ncols[0]!=ncols||pyramid1->nrows[0]!=nrows)
            trace(2,"(KLTTrackFeatures) Size of incoming image (%d by %d) "
                    "is different from size of previous image (%d by %d)\n",
                     ncols,nrows,pyramid1->ncols[0],pyramid1->nrows[0]);
        assert(pyramid1_gradx!=NULL);
        assert(pyramid1_grady!=NULL);
    }
    else  {
        floatimg1_created=true;
        floatimg1=create_floatimg(ncols,nrows);
        tofloatImage(img1,ncols,nrows,tmpimg);

        compute_smoothedimg(tmpimg,ComputeSmoothSigma(tc),floatimg1);
        pyramid1=klt_create_pyramid(ncols,nrows,(int)subsampling,tc->nPyramidLevels);

        klt_compute_pyramid(floatimg1,pyramid1,tc->pyramid_sigma_fact);
        pyramid1_gradx=klt_create_pyramid(ncols,nrows,(int)subsampling,tc->nPyramidLevels);
        pyramid1_grady=klt_create_pyramid(ncols,nrows,(int)subsampling,tc->nPyramidLevels);

        for (i=0;i<tc->nPyramidLevels;i++)
            compute_gradients(pyramid1->img[i],tc->grad_sigma,
                              pyramid1_gradx->img[i],
                              pyramid1_grady->img[i]);
    }
    /* Do the same thing with second image */
    floatimg2=create_floatimg(ncols,nrows);
    tofloatImage(img2,ncols,nrows,tmpimg);

    compute_smoothedimg(tmpimg,ComputeSmoothSigma(tc),floatimg2);
    pyramid2=klt_create_pyramid(ncols,nrows,(int)subsampling,
                                tc->nPyramidLevels);

    klt_compute_pyramid(floatimg2,pyramid2,tc->pyramid_sigma_fact);

    pyramid2_gradx=klt_create_pyramid(ncols,nrows,(int)subsampling,tc->nPyramidLevels);
    pyramid2_grady=klt_create_pyramid(ncols,nrows,(int)subsampling,tc->nPyramidLevels);

    for (i=0;i<tc->nPyramidLevels;i++) {
        compute_gradients(pyramid2->img[i],tc->grad_sigma,
                          pyramid2_gradx->img[i],
                          pyramid2_grady->img[i]);
    }
    /* Write internal images */
    if (tc->writeInternalImages) {
        char fname[80];
        for (i=0;i<tc->nPyramidLevels;i++) {
            sprintf(fname,"kltimg_tf_i%d.pgm",i);
            writeFloatImageToPGM(pyramid1->img[i],fname);

            sprintf(fname,"kltimg_tf_i%d_gx.pgm",i);
            writeFloatImageToPGM(pyramid1_gradx->img[i],fname);

            sprintf(fname,"kltimg_tf_i%d_gy.pgm",i);
            writeFloatImageToPGM(pyramid1_grady->img[i],fname);

            sprintf(fname,"kltimg_tf_j%d.pgm",i);
            writeFloatImageToPGM(pyramid2->img[i],fname);

            sprintf(fname,"kltimg_tf_j%d_gx.pgm",i);
            writeFloatImageToPGM(pyramid2_gradx->img[i],fname);

            sprintf(fname,"kltimg_tf_j%d_gy.pgm",i);
            writeFloatImageToPGM(pyramid2_grady->img[i],fname);
        }
    }
    /* For each feature, do ... */
    for (indx=0;indx<featurelist->n;indx++) {

        /* Only track features that are not lost */
        if (featurelist->feature[indx]->valid>=0) {

            xloc=featurelist->feature[indx]->x;
            yloc=featurelist->feature[indx]->y;

            /* Transform location to coarsest resolution */
            for (r=tc->nPyramidLevels-1;r>=0;r--) {
                xloc/=subsampling;
                yloc /= subsampling;
            }
            xlocout=xloc; ylocout=yloc;

            /* Beginning with coarsest resolution, do ... */
            for (r=tc->nPyramidLevels-1;r>=0;r--) {

                /* Track feature at current resolution */
                xloc*=subsampling;
                yloc*=subsampling;
                xlocout*=subsampling; ylocout*=subsampling;

                val=trackFeature(xloc,yloc,
                                 &xlocout,&ylocout,
                                 pyramid1->img[r],
                                 pyramid1_gradx->img[r],pyramid1_grady->img[r],
                                 pyramid2->img[r],
                                 pyramid2_gradx->img[r],pyramid2_grady->img[r],
                                 tc->window_width,tc->window_height,
                                 tc->step_factor,
                                 tc->max_iterations,
                                 tc->min_determinant,
                                 tc->min_displacement,
                                 tc->max_residue,
                                 tc->lighting_insensitive);

                if (val==KLT_SMALL_DET||val==KLT_OOB) {
                    break;
                }
            }
            /* Record feature */
            if (val==KLT_OOB) {
                featurelist->feature[indx]->x=-1.0f;
                featurelist->feature[indx]->y=-1.0f;
                featurelist->feature[indx]->valid=KLT_OOB;
                if (featurelist->feature[indx]->aff_img) free_floatimg(featurelist->feature[indx]->aff_img);
                if (featurelist->feature[indx]->aff_img_gradx) free_floatimg(featurelist->feature[indx]->aff_img_gradx);
                if (featurelist->feature[indx]->aff_img_grady) free_floatimg(featurelist->feature[indx]->aff_img_grady);

                featurelist->feature[indx]->aff_img=NULL;
                featurelist->feature[indx]->aff_img_gradx=NULL;
                featurelist->feature[indx]->aff_img_grady=NULL;

            }
            else if (outOfBounds(xlocout,ylocout,ncols,nrows,tc->borderx,tc->bordery)) {
                featurelist->feature[indx]->x    =-1.0f;
                featurelist->feature[indx]->y    =-1.0f;
                featurelist->feature[indx]->valid=KLT_OOB;
                if (featurelist->feature[indx]->aff_img) free_floatimg(featurelist->feature[indx]->aff_img);
                if (featurelist->feature[indx]->aff_img_gradx) free_floatimg(featurelist->feature[indx]->aff_img_gradx);
                if (featurelist->feature[indx]->aff_img_grady) free_floatimg(featurelist->feature[indx]->aff_img_grady);

                featurelist->feature[indx]->aff_img=NULL;
                featurelist->feature[indx]->aff_img_gradx=NULL;
                featurelist->feature[indx]->aff_img_grady=NULL;
            }
            else if (val==KLT_SMALL_DET)  {
                featurelist->feature[indx]->x=-1.0f;
                featurelist->feature[indx]->y=-1.0f;
                featurelist->feature[indx]->valid=KLT_SMALL_DET;
                if (featurelist->feature[indx]->aff_img) free_floatimg(featurelist->feature[indx]->aff_img);
                if (featurelist->feature[indx]->aff_img_gradx) free_floatimg(featurelist->feature[indx]->aff_img_gradx);
                if (featurelist->feature[indx]->aff_img_grady) free_floatimg(featurelist->feature[indx]->aff_img_grady);
                featurelist->feature[indx]->aff_img=NULL;
                featurelist->feature[indx]->aff_img_gradx=NULL;
                featurelist->feature[indx]->aff_img_grady=NULL;
            }
            else if (val==KLT_LARGE_RESIDUE) {
                featurelist->feature[indx]->x=-1.0f;
                featurelist->feature[indx]->y=-1.0f;
                featurelist->feature[indx]->valid=KLT_LARGE_RESIDUE;

                if (featurelist->feature[indx]->aff_img) free_floatimg(featurelist->feature[indx]->aff_img);
                if (featurelist->feature[indx]->aff_img_gradx) free_floatimg(featurelist->feature[indx]->aff_img_gradx);
                if (featurelist->feature[indx]->aff_img_grady) free_floatimg(featurelist->feature[indx]->aff_img_grady);

                featurelist->feature[indx]->aff_img=NULL;
                featurelist->feature[indx]->aff_img_gradx=NULL;
                featurelist->feature[indx]->aff_img_grady=NULL;
            }
            else if (val==KLT_MAX_ITERATIONS) {
                featurelist->feature[indx]->x=-1.0f;
                featurelist->feature[indx]->y=-1.0f;
                featurelist->feature[indx]->valid=KLT_MAX_ITERATIONS;
                if (featurelist->feature[indx]->aff_img) free_floatimg(featurelist->feature[indx]->aff_img);
                if (featurelist->feature[indx]->aff_img_gradx) free_floatimg(featurelist->feature[indx]->aff_img_gradx);
                if (featurelist->feature[indx]->aff_img_grady) free_floatimg(featurelist->feature[indx]->aff_img_grady);
                featurelist->feature[indx]->aff_img=NULL;
                featurelist->feature[indx]->aff_img_gradx=NULL;
                featurelist->feature[indx]->aff_img_grady=NULL;
            }
            else  {
                featurelist->feature[indx]->x=xlocout;
                featurelist->feature[indx]->y=ylocout;
                featurelist->feature[indx]->valid=KLT_TRACKED;
                if (tc->affineConsistencyCheck>=0&&val==KLT_TRACKED)  { /*for affine mapping*/
                    int border=2; /* add border for interpolation */

#ifdef DEBUG_AFFINE_MAPPING
                    glob_index = indx;
#endif
                    if(!featurelist->feature[indx]->aff_img){
                        /* save image and gradient for each feature at finest resolution after first successful track */
                        featurelist->feature[indx]->aff_img=create_floatimg((tc->affine_window_width+border),(tc->affine_window_height+border));
                        featurelist->feature[indx]->aff_img_gradx=create_floatimg((tc->affine_window_width+border),(tc->affine_window_height+border));
                        featurelist->feature[indx]->aff_img_grady=create_floatimg((tc->affine_window_width+border),(tc->affine_window_height+border));

                        am_getSubFloatImage(pyramid1->img[0],xloc,yloc,featurelist->feature[indx]->aff_img);
                        am_getSubFloatImage(pyramid1_gradx->img[0],xloc,yloc,featurelist->feature[indx]->aff_img_gradx);
                        am_getSubFloatImage(pyramid1_grady->img[0],xloc,yloc,featurelist->feature[indx]->aff_img_grady);

                        featurelist->feature[indx]->aff_x=xloc-(int)xloc+(tc->affine_window_width+border)/2;
                        featurelist->feature[indx]->aff_y=yloc-(int)yloc+(tc->affine_window_height+border)/2;;
                    }
                    else {
                        /* affine tracking */
                        val=am_trackFeatureAffine(featurelist->feature[indx]->aff_x, featurelist->feature[indx]->aff_y,
                                                  &xlocout, &ylocout,
                                                  featurelist->feature[indx]->aff_img,
                                                  featurelist->feature[indx]->aff_img_gradx,
                                                  featurelist->feature[indx]->aff_img_grady,
                                                  pyramid2->img[0],
                                                  pyramid2_gradx->img[0], pyramid2_grady->img[0],
                                                  tc->affine_window_width, tc->affine_window_height,
                                                  tc->step_factor,
                                                  tc->affine_max_iterations,
                                                  tc->min_determinant,
                                                  tc->min_displacement,
                                                  tc->affine_min_displacement,
                                                  tc->affine_max_residue,
                                                  tc->lighting_insensitive,
                                                  tc->affineConsistencyCheck,
                                                  tc->affine_max_displacement_differ,
                                                  &featurelist->feature[indx]->aff_Axx,
                                                  &featurelist->feature[indx]->aff_Ayx,
                                                  &featurelist->feature[indx]->aff_Axy,
                                                  &featurelist->feature[indx]->aff_Ayy
                        );
                        featurelist->feature[indx]->valid=val;
                        if (val!=KLT_TRACKED) {
                            featurelist->feature[indx]->x=-1.0f;
                            featurelist->feature[indx]->y=-1.0f;
                            featurelist->feature[indx]->aff_x=-1.0f;
                            featurelist->feature[indx]->aff_y=-1.0f;
                            /* free image and gradient for lost feature */
                            free_floatimg(featurelist->feature[indx]->aff_img);
                            free_floatimg(featurelist->feature[indx]->aff_img_gradx);
                            free_floatimg(featurelist->feature[indx]->aff_img_grady);
                            featurelist->feature[indx]->aff_img = NULL;
                            featurelist->feature[indx]->aff_img_gradx = NULL;
                            featurelist->feature[indx]->aff_img_grady = NULL;
                        } else {
                            /*featurelist->feature[indx]->x = xlocout;*/
                            /*featurelist->feature[indx]->y = ylocout;*/
                        }
                    }
                }

            }
        }
    }
    if (tc->sequentialMode) {
        tc->pyramid_last=pyramid2;
        tc->pyramid_last_gradx=pyramid2_gradx;
        tc->pyramid_last_grady=pyramid2_grady;
    }
    else  {
        free_pyramid(pyramid2);
        free_pyramid(pyramid2_gradx);
        free_pyramid(pyramid2_grady);
    }
    /* Free memory */
    free_floatimg(tmpimg);
    if (floatimg1_created) free_floatimg(floatimg1);
    free_floatimg(floatimg2);
    free_pyramid(pyramid1);
    free_pyramid(pyramid1_gradx);
    free_pyramid(pyramid1_grady);

    trace(3,"\n\t%d features successfully tracked.\n",
          klt_count_remaining_features(featurelist));
    return;
}
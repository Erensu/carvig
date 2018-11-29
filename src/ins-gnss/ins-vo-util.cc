/*------------------------------------------------------------------------------
* ins-vo-util.cc : visual odometry util functions
*
* reference :
*    [1] P.D.Groves, Principles of GNSS, Intertial, and Multisensor Integrated
*        Navigation System, Artech House, 2008
*    [2] Tedaldi D, Pretto A, Menegatti E. A robust and easy to implement method
*        for IMU calibration without external equipments,2014.
*    [3] Skog I, Handel P. Effects of time synchronization errors in GNSS-aided
*        INS 2008.
*    [4] Shin E H. Accuracy Improvement of Low Cost INS/GPS for Land Applications
*    [5] Weiss,Real-Time Metric State Estimation for Modular Vision-Inertial
*        Systems
*    [6] Kitt B , Geiger A , Lategahn H . Visual odometry based on stereo image
*        sequences with RANSAC-based outlier rejection scheme[C]// Intelligent
*        Vehicles Symposium. IEEE, 2010.
*
* version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
* history : 2018/05/09 1.0 new
*-----------------------------------------------------------------------------*/
#include <carvig.h>

/* constants------------------------------------------------------------------*/
#define COLINEAR_EPS    0.5       /* assume co-linear if within this threshold */

/* global variables-----------------------------------------------------------*/
static long int id_seed=1;        /* generate a new feature id */

/* check whether 3 points are co-linear---------------------------------------
 * args:    double *p1,*p2,*p3  I  input 3 points
 *          int dim             I  points in 2D or 3D
 *          int flag            I  0: indicating that p1, p2, p3 are homogneeous
 *                                    coordinates with arbitrary scale
 *                                 1: they are homogeneous with equal scale
 * return: 1: co-linear, 0: non-colinear
 * ---------------------------------------------------------------------------*/
extern int iscolinear(const double *p1,const double *p2,const double *p3,
                      int dim, int flag)
{
    double p1h[3],p2h[3],p3h[3],t[3],t1[3],t2[3];
    int i;

    trace(3,"iscolinear:\n");

    /* if data is 2D, assume they are 2D inhomogeneous coords.
     * make them homogeneous with scale 1.
     * */
    for (i=0;i<dim==2?2:3;i++) {
        p1h[i]=p1[i]; p2h[i]=p2[i]; p3h[i]=p3[i];
    }
    if (dim==2) p1h[2]=p2h[2]=p3h[2]=1.0;
    if (flag==0) {
        /* apply test that allows for homogeneous coords with arbitrary
         * scale.  p1 X p2 generates a normal vector to plane defined by
         * origin, p1 and p2. if the dot product of this normal with p3
         * is zero then p3 also lies in the plane, hence co-linear.
         * */
        cross3(p1h,p2h,t);
        return fabs(dot(t,p3h,3))<COLINEAR_EPS;
    }
    else if (flag==1) {
        /* assume inhomogeneous coords, or homogeneous coords
         * with equal scale
         * */
        for (i=0;i<3;i++) {
            t1[i]=p2h[i]-p1h[i]; t2[i]=p3h[i]-p1h[i];
        }
        cross3(t1,t2,t);
        return norm(t,3)<COLINEAR_EPS;
    }
}




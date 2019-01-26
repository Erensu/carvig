/*------------------------------------------------------------------------------
 * ins-vo-opencv-interface.cc : openCV interface common functions
 *
 * reference :
 *    [1] https://www.opencv.org/
 *
 * version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
 * history : 2018/12/02 1.0 new
 *-----------------------------------------------------------------------------*/
#include "carvig.h"
#include "relative-pose.h"

/* solve R|t using opencv------------------------------------------------------
 * args   :  voopt_t *opt      I  visual odometry options
 *           match_set_t *mf   I  feature list
 *           double *Tr        O  rotation and translation
 *           double *ratio     O  ratio of inliers
 * return : 1 (ok) or 0 (fail)
 * ----------------------------------------------------------------------------*/
extern int solveRt(const voopt_t *opt,const match_set_t *mf,double *Tr,double *ratio)
{
    trace(3,"solveRt:\n");
    

}
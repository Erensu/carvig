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

/* initial track---------------------------------------------------------------
 * args:    trackd_t *data  IO  track set data
 *          voopt_t *opt    I   track options
 * return: status (1: ok,0: fail)
 * ----------------------------------------------------------------------------*/
extern int inittrack(trackd_t *data,const voopt_t *opt)
{
    trace(3,"inittrack:\n");
}


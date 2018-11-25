/*------------------------------------------------------------------------------
* ins-features.cc : common features point function
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
*
* version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
* history : 2018/05/09 1.0 new
*-----------------------------------------------------------------------------*/
#include <carvig.h>

/* global variables-----------------------------------------------------------*/
static long int id_seed=1;        /* generate a new feature id */




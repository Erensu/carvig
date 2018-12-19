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
#include <cv.h>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

/* image buffer convert to Image struct--------------------------------------
 * args  :  img_t *img  I  image gray level measurement data
 *          Image **Img O  image rgb format measurement data
 * return: none
 * --------------------------------------------------------------------------*/
static int img2Img(const img_t *img,cv::Mat &Img)
{
    int i,j;
    for (i=0;i<Img.rows;i++) for (j=0;j<Img.cols;j++) {
            Img.at<cv::Vec3b>(i,j)[0]=img->data[i*img->w+j];
            Img.at<cv::Vec3b>(i,j)[1]=img->data[i*img->w+j];
            Img.at<cv::Vec3b>(i,j)[2]=img->data[i*img->w+j];
        }
    return 1;
}
/* draw tracks-----------------------------------------------------------------
 * args:    track_t *track  I  tracking data
 *          voopt_t *opt    I  options
 * return: none
 * ---------------------------------------------------------------------------*/
extern void drawtrack(const track_t *track,const voopt_t *opt)
{
    char allName[126][126];
    int i,j;

    trace(3,"drawtrack: n=%d\n",track->n);

    for (i=0;i<track->n;i++) { /* do for all tracks */
                                                                                                          
        /* for each track. */
        for (j=0;j<track->data[i].n;j++) {

            cv::Mat outImage(track->data[i].I[j].h,track->data[i].I[j].w,CV_8UC3);
            img2Img(&track->data[i].I[j],outImage);

            cv::Point point;
            point.x=(int)track->data[i].data[j].u;
            point.y=(int)track->data[i].data[j].v;

            cv::circle(outImage,point,2,cv::Scalar(0,0,255),1);
            cv::circle(outImage,point,8,cv::Scalar(0,0,255),1);

            sprintf(allName[j],"image_%ld_%d",track->data[i].uid,track->data[i].first_frame+j);
            cv::namedWindow(allName[j],CV_WINDOW_AUTOSIZE);
            cv::imshow(allName[j],outImage);                                                                                                                                  
            cv::waitKey(0);                                                                                                                                                     
        }
        for (j=0;j<track->data[i].n;j++) {
            cv::destroyWindow(allName[j]);
        }
    }
}
/* display image-------------------------------------------------------------*/
extern void dipsplyimg(const img_t *img)
{
    static int first=1;
    static char text[126];

    if (first) {
        cv::namedWindow("Display-Image",CV_WINDOW_AUTOSIZE);
        first=0;
    }
    cv::Mat outImage(img->h,img->w,CV_8UC3);
    img2Img(img,outImage);
    sprintf(text,"id=%d",img->id);

    cv::putText(outImage,text,cv::Point(8,15),cv::FONT_HERSHEY_PLAIN,1.0,cv::Scalar(0,255,255));
    cv::imshow("Display-Image",outImage);
    cv::waitKey(5);
}


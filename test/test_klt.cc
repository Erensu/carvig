/**********************************************************************
Finds the 100 best features in an image, and tracks these
features to the next image.  Saves the feature
locations (before and after tracking) to text files and to PPM files,
and prints the features to the screen.
**********************************************************************/
#include <vision.h>
#include <carvig.h>
int main()
{
    unsigned char *img1, *img2;
    ptracking_context_t tc;
    pfeaturelist_t fl;
    int nFeatures=100;
    int ncols, nrows;
    int i;

    tc=klt_create_tracking_context();
    fl=klt_create_featurelist(nFeatures);

    img1=pgmReadFile("/home/sujinglan/carvig/test/data/klt/img0.pgm",NULL,&ncols,&nrows);
    img2=pgmReadFile("/home/sujinglan/carvig/test/data/klt/img1.pgm",NULL,&ncols,&nrows);

    klt_select_goodfeatures(tc,img1,ncols,nrows,fl);

    printf("\nIn first image:\n");
    for (i = 0 ; i < fl->n ; i++)  {
        printf("Feature #%d:  (%f,%f) with value of %d\n",
               i, fl->feature[i]->x, fl->feature[i]->y,
               fl->feature[i]->valid);
    }
    writeFeatureListToPPM(fl, img1, ncols, nrows, "./feat1.ppm");

    klt_track_features(tc, img1, img2, ncols, nrows, fl);

    printf("\nIn second image:\n");
    for (i = 0 ; i < fl->n ; i++)  {
        printf("Feature #%d:  (%f,%f) with value of %d\n",
               i, fl->feature[i]->x, fl->feature[i]->y,
               fl->feature[i]->valid);
    }
    writeFeatureListToPPM(fl, img2, ncols, nrows, "./feat2.ppm");
    return 0;
}

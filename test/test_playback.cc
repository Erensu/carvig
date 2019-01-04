#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <unistd.h>

#include <carvig.h>

#include <cv.h>
#include <opencv2/opencv.hpp>
using namespace cv;

int main(int argc, char **argv)
{
    DIR *dir;
    struct dirent *ptr;
    const char basePath[1024]="/media/sujinglan/Files/mixed-sensor-data/m39/2018-11-22/camera";
    char imgath[1024];
    static char all_img_path[100000][156];
    int i=0;

    if ((dir=opendir(basePath)) == NULL) {
        return 0;
    }
    while ((ptr=readdir(dir)) != NULL) {
        if (strcmp(ptr->d_name,".")==0 || strcmp(ptr->d_name,"..")==0) continue;
        else if (ptr->d_type == 8) {

            if (strstr(ptr->d_name,"dat")) continue;

            sprintf(imgath,"%s/%s",basePath,ptr->d_name);

            int id;
            long int stamp,stamp1;
            sscanf(ptr->d_name,"%d_%ld_%ld.jpg",&id,&stamp,&stamp1);

            strcpy(all_img_path[id-1],imgath);

            i++;
        }

    }
    cv::namedWindow("image",CV_WINDOW_AUTOSIZE);
    for (int j=0;j<i;j+=5) {
        cv::Mat img=cv::imread(all_img_path[j],IMREAD_COLOR);

        if (img.empty())
        {
            fprintf(stderr, "Can not load image %s\n", ptr->d_name);
            continue;
        }
        cv::imshow("image",img);
        cv::waitKey(1);
    }
}
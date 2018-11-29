#include <carvig.h>
#include <iostream>
#include <png++/png.hpp>

using namespace std;

int main(int argc, char **argv)
{
    matchopt_t matchopt={0};
    bucketopt_t bucketopt;
    calib_t calib={0};

    voopt_t voopt;

    match_t match={{0}};

    matchopt.fu=645.24;
    matchopt.fv=645.24;
    matchopt.f =645.24; matchopt.cu=635.96; matchopt.cv=194.13;

    matchopt.nms_n=3;
    matchopt.nms_tau=50;
    matchopt.match_binsize=50;
    matchopt.match_radius =200;
    matchopt.match_disp_tol=2;
    matchopt.outlier_disp_tol=5;
    matchopt.outlier_flow_tol=5;

    matchopt.multi_stage=1;
    matchopt.half_res=1;
    matchopt.refine=1;

    matchopt.img_w=1344; matchopt.img_h=372;

    bucketopt.nmax=5;
    bucketopt.w=50;
    bucketopt.h=50;

    matchopt.bucket=bucketopt;

    calib.cu=635.96; calib.cv=194.13; calib.f=645.24;
    calib.fu=calib.fv=645.24;

    voopt.inlier_thres=0.00001;
    voopt.motion_thres=100.0;
    voopt.ransac_iters=2000;

    voopt.calib=calib;
    voopt.match=matchopt;

    match.opt=matchopt;
    match.opt.bucket=bucketopt;

    gtime_t t0={0};
    string dir="/home/sujinglan/libviso2/test/2010_03_09_drive_0019";
    img_t img_input={0};
    double Tr[19];

    init_match(&match,&matchopt);

    for (int i=0;i<370;i++) {
        /* input file names */
        char base_name[256]; sprintf(base_name,"%06d.png",i);
        string img_file_name =dir+"/I1_"+base_name;

        /* load  single input image */
        png::image<png::gray_pixel>img(img_file_name);

        /* image dimensions */
        int32_t width =img.get_width();
        int32_t height=img.get_height();

        /* convert input images to uint8_t buffer */
        initimg(&img_input,width,height,t0);
        int32_t k=0;
        for (int32_t v=0;v<height;v++) {
            for (int32_t u=0;u<width;u++) {
                img_input.data[k]=img.get_pixel(u,v);
                k++;
            }
        }
        img_input.w=width; img_input.h=height;
        
        /* match feature points */
        matchfeats(&match,&img_input);

        //for (int j=0;j<match.mp_bucket.n;j++) {
        //    match_point point;
        //
        //    point.up=match.mp_bucket.data[j].up; point.vp=match.mp_bucket.data[j].vp;
        //    point.uc=match.mp_bucket.data[j].uc; point.vc=match.mp_bucket.data[j].vc;
        //
        //    kltstatus(&point,&match.Ip,&match.Ic,NULL);
        //}

        /* vo estimate */
        //estmonort(&voopt,&match.mp_bucket,Tr);

        //trace(3,"transform matrix:\n");
        //tracemat(3,Tr,4,4,12,6);

        freeimg(&img_input);
    }
    return 0;
}
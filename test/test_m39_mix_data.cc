#include <carvig.h>
int main(int argc, char **argv)
{
    FILE *fp=fopen("/media/sujinglan/Files/mixed-sensor-data/m39/imu-cam-calib/181123065729.dat","rb") ;
    FILE *fpo=fopen("./m39-mix.out","w");
    char imgfile[1024];
    const char *output_dir="/media/sujinglan/Files/mixed-sensor-data/m39/imu-cam-calib/output";

    int flag=0,i=0;
    raw_t raw={{0}};

    init_raw(&raw,1);

    while (1) {
        flag=input_m39_mixf(&raw,fp);
        if (flag==-2) break;
        if (flag==11) {
            if (get_m39_img("/media/sujinglan/Files/mixed-sensor-data/m39/imu-cam-calib",
                            raw.m39.fts.tv_sec,raw.m39.fts.tv_nsec,imgfile)) {
                fprintf(fpo,"%.5lf  %.5lf  %ld-%ld  %ld-%ld  %s\n",
                        time2gpst(raw.m39.time,NULL),time2gpst(raw.m39.zda,NULL),
                        raw.m39.pps.tv_sec,raw.m39.pps.tv_nsec,
                        raw.m39.fts.tv_sec,raw.m39.fts.tv_nsec,imgfile);
            }
            else {
                exit(0);
            }
        }
    }
    fclose(fpo);
}
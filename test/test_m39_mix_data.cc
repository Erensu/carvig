#include <carvig.h>
#include <include/carvig.h>

int main(int argc, char **argv)
{
    FILE *fp=fopen("/media/sujinglan/Files/mixed-sensor-data/m39/2018-11-22/camera/181122074715.dat","rb") ;
    FILE *fpo=fopen("./m39-mix.out","w");
    char imgfile[1024];

    int flag=0,i=0;
    raw_t raw={{0}};

    init_raw(&raw,1);

    prcopt_t opt={0};
    strcpy(opt.monodir,"/media/sujinglan/Files/mixed-sensor-data/m39/2018-11-22/camera");

    raw.optp=&opt;

    while (1) {
        flag=input_m39_mixf(&raw,fp);
        if (flag==-2) break;
        if (flag==11) {
            if (get_m39_img("/media/sujinglan/Files/mixed-sensor-data/m39/2018-11-22/camera",
                            raw.m39.fts.tv_sec,raw.m39.fts.tv_nsec,imgfile)) {
                fprintf(fpo,"%.5lf  %.5lf  %ld-%ld  %ld-%ld  %s\n",
                        time2gpst(raw.m39.time,NULL),time2gpst(raw.m39.zda,NULL),
                        raw.m39.pps.tv_sec,raw.m39.pps.tv_nsec,
                        raw.m39.fts.tv_sec,raw.m39.fts.tv_nsec,imgfile);
                fflush(fpo);
            }
            else {
                exit(0);
            }
        }
    }
    fclose(fpo);
}
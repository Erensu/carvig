#include <carvig.h>

int main(int argc, char **argv)
{
    const char *file="/media/sujinglan/Files/mixed-sensor-data/m39/2018-11-22/gps-ins/181122154542.gps";

    int flag=0,i=0;
    raw_t raw={{0}};

    FILE *fp=fopen(file,"rb");

    init_raw(&raw,1);

    while (true) {
        flag=input_oem6f_pose(&raw,fp);
        if (flag==-1) break;
        if (flag==34) {
            fprintf(stderr,"%ld: %10.6lf  %10.6lf  %10.6lf  %10.6lf  %10.6lf  %10.6lf  %10.6lf  %d\n",
                    raw.pose.time.time,
                    raw.pose.rpy[0]*R2D,
                    raw.pose.rpy[1]*R2D,
                    raw.pose.rpy[2]*R2D,
                    SQRT(raw.pose.var[0]*R2D),
                    SQRT(raw.pose.var[1]*R2D),
                    SQRT(raw.pose.var[2]*R2D),raw.pose.len,raw.pose.stat);
        }
    }
    return 0;
}
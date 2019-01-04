#include <carvig.h>
int main(int argc, char **argv)
{
    const char *imufile="/media/sujinglan/Files/mixed-sensor-data/m39/imu-cam-calib/181123145613.imu";
    const char *imgdir="/media/sujinglan/Files/mixed-sensor-data/m39/imu-cam-calib";
    const char *datfile="/media/sujinglan/Files/mixed-sensor-data/m39/imu-cam-calib/181123065729.dat";

    kalibrrosbag(datfile,imufile,imgdir,imgdir);

    return 0;
}
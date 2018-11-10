/*----------------------------------------------------------------------------
* m39-mix.cc : decode m39/image time tag raw data
*
* version : $Revision:$ $Date:$
* history : 2017/03/08  1.0  new
*----------------------------------------------------------------------------*/
#include <carvig.h>

/* read m39 mix raw data (included image/imu)--------------------------------
 * args  :  raw_t *raw         IO  raw struct
 *          unsigned char data I   raw data
 * return >1:ok, 0:fail
 * --------------------------------------------------------------------------*/
extern int input_m39_mix(raw_t *raw, unsigned char data)
{
    trace(5,"input_m39_mix: data=%02x\n",data);
    return 0;
}
/* read m39-mix raw data from file-------------------------------------------
 * args   : raw_t  *raw   IO     receiver raw data control struct
 *          FILE   *fp    I      file pointer
 * return >1:ok, 0:fail
 * --------------------------------------------------------------------------*/
extern int input_m39_mixf(raw_t *raw,FILE *fp)
{
    int i,data,ret;

    trace(3,"input_m39_mixf:\n");

    for (i=0;i<4096;i++) {
        if ((data=fgetc(fp))==EOF) return -2;
        if ((ret=input_m39_mix(raw,(unsigned char)data))) return ret;
    }
    return 0; /* return at every 4k bytes */
}

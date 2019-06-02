/*------------------------------------------------------------------------------
* rnx2rtkp.cc : read rinex obs/nav files and compute receiver positions
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:55:16 $
*-----------------------------------------------------------------------------*/
#include <stdarg.h>
#include "carvig.h"

#define MAXFILE     16                  /* max number of input files */

/* receiver options table ----------------------------------------------------*/
#define TIMOPT  "0:gpst,1:utc,2:jst,3:tow"
#define CONOPT  "0:dms,1:deg,2:xyz,3:enu,4:pyl"
#define FLGOPT  "0:off,1:std+2:age/ratio/ns"
#define ISTOPT  "0:off,1:serial,2:file,3:tcpsvr,4:tcpcli,7:ntripcli,8:ftp,9:http"
#define OSTOPT  "0:off,1:serial,2:file,3:tcpsvr,4:tcpcli,6:ntripsvr,11:ntripc_c"
#define FMTOPT  "0:rtcm2,1:rtcm3,2:oem4,3:oem3,4:ubx,5:ss2,6:hemis,7:skytraq,8:gw10,9:javad,10:nvs,11:binex,12:rt17,13:sbf,14:cmr,15:tersus,18:sp3,19:rnxclk,20:sbas,21:nmea,22:gsof,23:ublox-evk-m8u,24:ublox-sol,25:m39,26:rinex,27:m39-mix,28:euroc-imu,29:euroc-img,30:karl-img,31:malaga-gnss,32:malaga-imu,33:malaga-img,34:oem6-sol,35:oem6-pose,36:oem6-raw"
#define NMEOPT  "0:off,1:latlon,2:single"
#define SOLOPT  "0:llh,1:xyz,2:enu,3:nmea,4:stat,5:gsif,6:ins"
#define MSGOPT  "0:all,1:rover,2:base,3:corr"

#define OPTSDIR     "."                 /* default config directory */
#define OPTSFILE    "rtkrcv.conf"       /* default config file */
#define MAXSTR      1024                /* max length of a stream */

/* help text -----------------------------------------------------------------*/
static const char *help[]={
        "",
        " usage: rnx2rtkp [option]... file file [...]",
        "",
        " Read RINEX OBS/NAV/GNAV/HNAV/CLK, SP3, SBAS message log files and ccompute ",
        " receiver (rover) positions and output position solutions.",
        " The first RINEX OBS file shall contain receiver (rover) observations. For the",
        " relative mode, the second RINEX OBS file shall contain reference",
        " (base station) receiver observations. At least one RINEX NAV/GNAV/HNAV",
        " file shall be included in input files. To use SP3 precise ephemeris, specify",
        " the path in the files. The extension of the SP3 file shall be .sp3 or .eph.",
        " All of the input file paths can include wild-cards (*). To avoid command",
        " line deployment of wild-cards, use \"...\" for paths with wild-cards.",
        " Command line options are as follows ([]:default). With -k option, the",
        " processing options are input from the configuration file. In this case,",
        " command line options precede options in the configuration file.",
        "",
        " -?        print help",
        " -p file   input options from configuration file [off]",
        " -o file   set output file [stdout]",
};
static int strtype[]={                  /* stream types */
        STR_SERIAL,STR_NONE,STR_NONE,STR_NONE,STR_NONE,STR_NONE,STR_NONE,STR_NONE,STR_NONE,
        STR_NONE,STR_NONE,STR_NONE,STR_NONE,STR_NONE
};
static char strpath[14][MAXSTR]={"","","","","","","","","","","","","",""}; /* stream paths */
static int strfmt[]={                   /* stream formats */
        STRFMT_UBX,STRFMT_RTCM3,STRFMT_SP3,STRFMT_UBXM8,STRFMT_UBXSOL,SOLF_LLH,SOLF_NMEA,
        SOLF_LLH,SOLF_LLH
};
static opt_t rcvopts[]={
        {"inpstr1-type",    3,  (void *)&strtype[0],         ISTOPT },
        {"inpstr2-type",    3,  (void *)&strtype[1],         ISTOPT },
        {"inpstr3-type",    3,  (void *)&strtype[2],         ISTOPT },
        {"inpstr4-type",    3,  (void *)&strtype[3],         ISTOPT },
        {"inpstr5-type",    3,  (void *)&strtype[4],         ISTOPT },
        {"inpstr6-type",    3,  (void *)&strtype[5],         ISTOPT },
        {"inpstr7-type",    3,  (void *)&strtype[6],         ISTOPT },
        {"inpstr1-path",    2,  (void *)strpath [0],         ""     },
        {"inpstr2-path",    2,  (void *)strpath [1],         ""     },
        {"inpstr3-path",    2,  (void *)strpath [2],         ""     },
        {"inpstr4-path",    2,  (void *)strpath [3],         ""     },
        {"inpstr5-path",    2,  (void *)strpath [4],         ""     },
        {"inpstr6-path",    2,  (void *)strpath [5],         ""     },
        {"inpstr7-path",    2,  (void *)strpath [6],         ""     },
        {"inpstr1-format",  3,  (void *)&strfmt [0],         FMTOPT },
        {"inpstr2-format",  3,  (void *)&strfmt [1],         FMTOPT },
        {"inpstr3-format",  3,  (void *)&strfmt [2],         FMTOPT },
        {"inpstr4-format",  3,  (void *)&strfmt [3],         FMTOPT },
        {"inpstr5-format",  3,  (void *)&strfmt [4],         FMTOPT },
        {"inpstr6-format",  3,  (void *)&strfmt [5],         FMTOPT },
        {"inpstr7-format",  3,  (void *)&strfmt [6],         FMTOPT },
        {"outstr1-type",    3,  (void *)&strtype[7],         OSTOPT },
        {"outstr2-type",    3,  (void *)&strtype[8],         OSTOPT },
        {"outstr1-path",    2,  (void *)strpath [7],         ""     },
        {"outstr2-path",    2,  (void *)strpath [8],         ""     },
        {"outstr1-format",  3,  (void *)&strfmt [7],         SOLOPT },
        {"outstr2-format",  3,  (void *)&strfmt [8],         SOLOPT },
        {"logstr1-type",    3,  (void *)&strtype[9],         OSTOPT },
        {"logstr2-type",    3,  (void *)&strtype[10],        OSTOPT },
        {"logstr3-type",    3,  (void *)&strtype[11],        OSTOPT },
        {"logstr4-type",    3,  (void *)&strtype[12],        OSTOPT },
        {"logstr5-type",    3,  (void *)&strtype[13],        OSTOPT },
        {"logstr1-path",    2,  (void *)strpath [9],         ""     },
        {"logstr2-path",    2,  (void *)strpath [10],        ""     },
        {"logstr3-path",    2,  (void *)strpath [11],        ""     },
        {"logstr4-path",    2,  (void *)strpath [12],        ""     },
        {"logstr5-path",    2,  (void *)strpath [13],        ""     },

        {"",0,NULL,""}
};
/* print help ----------------------------------------------------------------*/
static void printhelp(void)
{
    int i;
    for (i=0;i<(int)(sizeof(help)/sizeof(*help));i++) fprintf(stderr,"%s\n",help[i]);
    exit(0);
}
/* rnx2rtkp main -------------------------------------------------------------*/
int main(int argc, char **argv)
{
    prcopt_t prcopt=prcopt_default;
    solopt_t solopt=solopt_default;
    filopt_t filopt={""};
    gtime_t ts={0},te={0};
    double tint=0.0;
    int i,n,ret,port=0;
    char *infile[MAXFILE],*outfile="",file[1024];

    for (i=1;i<argc;i++) {
        if      (!strcmp(argv[i],"-p")&&i+1<argc) strcpy(file,argv[++i]);
        else if (!strcmp(argv[i],"-o")&&i+1<argc) outfile=argv[++i];
        else if (!strcmp(argv[i],"-m")&&i+1<argc) port=atoi(argv[++i]);
        else printhelp();
    }
    /* load options file */
    if (!*file) sprintf(file,"%s/%s",OPTSDIR,OPTSFILE);
    resetsysopts();
    if (!loadopts(file,rcvopts)||!loadopts(file,sysopts)||
        !loadopts(file,insopts)) {
        fprintf(stderr,"no options file: %s. defaults used\n",file);
    }
    getsysopts(&prcopt,&solopt,&filopt);
    for (i=0;i<MAXFILE;i++) {
        if (!(infile[i]=(char *)malloc(1024))) {
            for (;i>=0;i--) free(infile[i]); return -1;
        }
    }
    /* set input files */
    strcpy(infile[0],strpath[0]);
    strcpy(infile[1],strpath[1]);
    strcpy(infile[2],filopt.navfile);
    strcpy(infile[3],filopt.bdsfile);
    strcpy(infile[4],filopt.glofile);
    strcpy(infile[5],filopt.mixfile);
    n=6;

    if (n<=0) {
        showmsg("error : no input file");
        return -2;
    }
    ret=postpos(ts,te,tint,0.0,&prcopt,&solopt,&filopt,infile,n,outfile,"","",port);

    if (!ret) fprintf(stderr,"%40s\r","");
    for (i=0;i<MAXFILE;i++) {
        if (infile[i]) free(infile[i]);
    }
    return ret;
}
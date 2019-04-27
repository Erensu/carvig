/*------------------------------------------------------------------------------
* rnx2rtkp.cc : read rinex obs/nav files and compute receiver positions
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:55:16 $
*-----------------------------------------------------------------------------*/
#include <stdarg.h>
#include "carvig.h"

#define MAXFILE     16                  /* max number of input files */

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
        " -k file   input options from configuration file [off]",
        " -o file   set output file [stdout]",
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
    double tint=0.0,es[]={2000,1,1,0,0,0},ee[]={2000,12,31,23,59,59},pos[3];
    int i,j,n,ret;
    char *infile[MAXFILE],*outfile="",*p;

    prcopt.mode  =PMODE_KINEMA;
    prcopt.navsys=0;
    prcopt.refpos=1;
    prcopt.glomodear=1;
    solopt.timef=0;

    /* load options from configuration file */
    for (i=1;i<argc;i++) {
        if (!strcmp(argv[i],"-k")&&i+1<argc) {
            resetsysopts();
            if (!loadopts(argv[++i],sysopts)) return -1;
            getsysopts(&prcopt,&solopt,&filopt);
        }
    }
    for (i=1,n=0;i<argc;i++) {
        if      (!strcmp(argv[i],"-o")&&i+1<argc) outfile=argv[++i];
        else if (*argv[i]=='-?') printhelp();
        else if (n<MAXFILE) infile[n++]=argv[i];
    }
    if (n<=0) {
        showmsg("error : no input file");
        return -2;
    }
    ret=postpos(ts,te,tint,0.0,&prcopt,&solopt,&filopt,infile,n,outfile,"","");

    if (!ret) fprintf(stderr,"%40s\r","");
    return ret;
}
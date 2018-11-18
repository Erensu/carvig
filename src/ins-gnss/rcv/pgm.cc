/*-----------------------------------------------------------------------------
 * pgm.cc
 *
 * Various routines to manipulate PGM files.
 *----------------------------------------------------------------------------*/
#include <vision.h>
#include <carvig.h>

/* global variable -----------------------------------------------------------*/
#define LENGTH 80

/*----------------------------------------------------------------------------*/
static void _getNextString(FILE *fp,char *line)
{
    int i;

    line[0]='\0';
    while (line[0]=='\0') {
        fscanf(fp,"%s",line);
        i=-1;
        do {
            i++;
            if (line[i]=='#') {
                line[i]='\0';
                while (fgetc(fp)!='\n');
            }
        }  while (line[i]!='\0');
    }
}
/*----------------------------------------------------------------------------*/
static void pnmReadHeader(FILE *fp,int *magic,int *ncols, int *nrows,int *maxval)
{
    char line[LENGTH];

    /* Read magic number */
    _getNextString(fp, line);
    if (line[0] != 'P')
        trace(2,"Magic number does not begin with 'P', "
                "but with a '%c'",line[0]);
    sscanf(line,"P%d",magic);

    /* Read size, skipping comments */
    _getNextString(fp,line);
    *ncols=atoi(line);
    _getNextString(fp,line);
    *nrows=atoi(line);
    if (*ncols<0||*nrows<0||*ncols>10000||*nrows>10000) {
        trace(2,"The dimensions %d x %d are unacceptable",*ncols,*nrows);
    }
    /* Read maxval, skipping comments */
    _getNextString(fp,line);
    *maxval=atoi(line);
    fread(line,1,1,fp); /* Read newline which follows maxval */

    if (*maxval!=255) {
        trace(2,"Maxval is not 255, but %d",*maxval);
    }
}
/*----------------------------------------------------------------------------*/
static void pgmReadHeader(FILE *fp, int *magic, int *ncols, int *nrows,
                          int *maxval)
{
    pnmReadHeader(fp,magic,ncols,nrows,maxval);
    if (*magic!=5) {
        trace(2,"Magic number is not 'P5', but 'P%d'",*magic);
    }
}
/*----------------------------------------------------------------------------*/
static void ppmReadHeader(FILE *fp,int *magic,int *ncols, int *nrows,int *maxval)
{
    pnmReadHeader(fp,magic,ncols,nrows,maxval);
    if (*magic!=6) {
        trace(2,"Magic number is not 'P6', but 'P%d'",*magic);
    }
}
/*----------------------------------------------------------------------------*/
static void pgmReadHeaderFile(char *fname, int *magic, int *ncols, int *nrows,
                              int *maxval)
{
    FILE *fp;

    /* Open file */
    if ((fp=fopen(fname,"rb"))==NULL) {
        trace(2,"can't open file named '%s' for reading\n",fname);
    }
    /* Read header */
    pgmReadHeader(fp,magic,ncols,nrows,maxval);

    /* Close file */
    fclose(fp);
}
/*----------------------------------------------------------------------------*/
static void ppmReadHeaderFile(char *fname, int *magic, int *ncols, int *nrows,
                              int *maxval)
{
    FILE *fp;

    /* Open file */
    if ((fp = fopen(fname,"rb"))==NULL) {
        trace(2,"Can't open file named '%s' for reading\n",fname);
    }
    /* Read header */
    ppmReadHeader(fp,magic,ncols,nrows,maxval);

    /* Close file */
    fclose(fp);
}
/*----------------------------------------------------------------------------
 * pgmRead
 *
 * NOTE:  If img is NULL, memory is allocated.
 *--------------------------------------------------------------------------*/
extern unsigned char* pgmRead(FILE *fp,unsigned char *img,
                              int *ncols, int *nrows)
{
    unsigned char *ptr;
    int magic, maxval;
    int i;

    /* Read header */
    pgmReadHeader(fp,&magic,ncols,nrows,&maxval);

    /* Allocate memory, if necessary, and set pointer */
    if (img==NULL) {
        ptr=(unsigned char *)malloc(*ncols**nrows*sizeof(char));
        if (ptr==NULL)
            trace(2,"(pgmRead) Memory not allocated");
    }
    else
        ptr=img;

    /* Read binary image data */
    {
        unsigned char *tmpptr=ptr;
        for (i=0;i<*nrows;i++) {
            fread(tmpptr,*ncols,1,fp);
            tmpptr+=*ncols;
        }
    }
    return ptr;
}
/*----------------------------------------------------------------------------
 * pgmReadFile
 *
 * NOTE:  If img is NULL, memory is allocated.
 *----------------------------------------------------------------------------*/
extern unsigned char* pgmReadFile(char *fname,unsigned char *img,
                                  int *ncols, int *nrows)
{
    unsigned char *ptr;
    FILE *fp;

    /* Open file */
    if ((fp=fopen(fname,"rb")) == NULL)
        trace(2,"Can't open file named '%s' for reading\n",fname);

    /* Read file */
    ptr=pgmRead(fp,img,ncols,nrows);

    /* Close file */
    fclose(fp);
    return ptr;
}
/*----------------------------------------------------------------------------
 * pgmWrite
 *----------------------------------------------------------------------------*/
extern void pgmWrite(FILE *fp,unsigned char *img, int ncols, int nrows)
{
    int i;

    /* Write header */
    fprintf(fp,"P5\n");
    fprintf(fp,"%d %d\n",ncols,nrows);
    fprintf(fp,"255\n");

    /* Write binary data */
    for (i=0;i<nrows;i++) {
        fwrite(img,ncols,1,fp);
        img+=ncols;
    }
}
/*----------------------------------------------------------------------------
 * pgmWriteFile
 *---------------------------------------------------------------------------*/
extern void pgmWriteFile(char *fname, unsigned char *img, int ncols, int nrows)
{
    FILE *fp;

    /* Open file */
    if ((fp=fopen(fname,"wb")) == NULL)
        trace(2,"(pgmWriteFile) Can't open file named '%s' for writing\n",fname);

    /* Write to file */
    pgmWrite(fp,img,ncols,nrows);

    /* Close file */
    fclose(fp);
}
/*----------------------------------------------------------------------------
 * ppmWrite
 *---------------------------------------------------------------------------*/
extern void ppmWrite(FILE *fp,unsigned char *redimg,unsigned char *greenimg,
                     unsigned char *blueimg,
                     int ncols, int nrows)
{
    int i, j;

    /* Write header */
    fprintf(fp,"P6\n");
    fprintf(fp,"%d %d\n",ncols,nrows);
    fprintf(fp,"255\n");

    /* Write binary data */
    for (j=0;j<nrows;j++)  {
        for (i=0;i<ncols;i++)  {
            fwrite(redimg,1,1,fp);
            fwrite(greenimg,1,1,fp);
            fwrite(blueimg, 1,1,fp);
            redimg++; greenimg++; blueimg++;
        }
    }
}
/*----------------------------------------------------------------------------
 * ppmWriteFileRGB
 *---------------------------------------------------------------------------*/
extern void ppmWriteFileRGB(char *fname, unsigned char *redimg,
                            unsigned char *greenimg,unsigned char *blueimg,
                            int ncols, int nrows)
{
    FILE *fp;

    /* Open file */
    if ((fp=fopen(fname,"wb"))==NULL)
        trace(2,"(ppmWriteFileRGB) Can't open file named '%s' for writing\n",fname);

    /* Write to file */
    ppmWrite(fp,redimg,greenimg,blueimg,ncols,nrows);

    /* Close file */
    fclose(fp);
}
/*write float image to PGM ---------------------------------------------------
 * args:    pfloatimg_t *img  IO  float image data
 *          char *filename    O   write image data file path
 * return: none
 * ---------------------------------------------------------------------------*/
extern void writeFloatImageToPGM(pfloatimg_t img,char *filename)
{
    int npixs =img->ncols*img->nrows;
    float mmax=-999999.9f,mmin=999999.9f;
    float fact;
    float *ptr;
    unsigned char *byteimg,*ptrout;
    int i;

    /* calculate minimum and maximum values of float image */
    ptr=img->data;
    for (i=0;i<npixs;i++)  {
        mmax=MAX(mmax,*ptr);
        mmin=MIN(mmin,*ptr);
        ptr++;
    }
    /* allocate memory to hold converted image */
    byteimg=(unsigned char *)malloc(npixs*sizeof(unsigned char));

    /* Convert image from float to uchar */
    fact=255.0f/(mmax-mmin);
    ptr =img->data;
    ptrout=byteimg;
    for (i=0;i<npixs;i++) {
        *ptrout++=(unsigned char)((*ptr++-mmin)*fact);
    }
    /* write uchar image to PGM */
    pgmWriteFile(filename,byteimg,img->ncols,img->nrows);

    /* free memory */
    free(byteimg);
}
extern void writeFeatureListToPPM(pfeaturelist_t featurelist,
                                  unsigned char *greyimg,int ncols,int nrows,
                                  char *filename)
{
    int nbytes = ncols * nrows * sizeof(char);
    unsigned char *redimg, *grnimg, *bluimg;
    int offset;
    int x, y, xx, yy;
    int i;

    /* Allocate memory for component images */
    redimg = (unsigned char *)  malloc(nbytes);
    grnimg = (unsigned char *)  malloc(nbytes);
    bluimg = (unsigned char *)  malloc(nbytes);
    if (redimg == NULL || grnimg == NULL || bluimg == NULL)
        return;

    /* Copy grey image to component images */
    if (sizeof(unsigned char) != 1) {
        trace(2,"(KLTWriteFeaturesToPPM)  KLT_PixelType is not uchar");
    }
    memcpy(redimg, greyimg, nbytes);
    memcpy(grnimg, greyimg, nbytes);
    memcpy(bluimg, greyimg, nbytes);

    /* Overlay features in red */
    for (i = 0 ; i < featurelist->n ; i++)
        if (featurelist->feature[i]->valid >= 0)  {
            x = (int) (featurelist->feature[i]->x + 0.5);
            y = (int) (featurelist->feature[i]->y + 0.5);
            for (yy = y - 1 ; yy <= y + 1 ; yy++)
                for (xx = x - 1 ; xx <= x + 1 ; xx++)
                    if (xx >= 0 && yy >= 0 && xx < ncols && yy < nrows)  {
                        offset = yy * ncols + xx;
                        *(redimg + offset) = 255;
                        *(grnimg + offset) = 0;
                        *(bluimg + offset) = 0;
                    }
        }

    /* Write to PPM file */
    ppmWriteFileRGB(filename, redimg, grnimg, bluimg, ncols, nrows);

    /* Free memory */
    free(redimg);
    free(grnimg);
    free(bluimg);
}



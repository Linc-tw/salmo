

  //------------------------------------------------------//
  //--  FITSFunctions.h                                 --//
  //--  Version 2020.07.09                              --//
  //--                                                  --//
  //--  Copyright (C) 2020 - Chieh-An Lin               --//
  //--  GNU GPLv3 - https://www.gnu.org/licenses/       --//
  //------------------------------------------------------//


#include "commonHeader.h"
  
#ifndef __MFP_FITS_FUNCTIONS__
#define __MFP_FITS_FUNCTIONS__

#include <fitsio.h>


typedef struct {
  fitsfile *file;
  int status;
  int nbHDU;
  int HDUType;
  long nbRows;
  int nbColumns;
  int colType;
  long colRepeat;
  int bitpix;
  long resol[3];
} FITS_t;

typedef struct {
  long long nsidePix;
  long long nbPix;
  float *map;
} HPMap_t;


//-- Functions related to initialization
FITS_t *initializeTableReader_FITS_t(char *name);
FITS_t *initializeImageReader_FITS_t(char *name);
FITS_t *initializeTableWriter_FITS_t(char *name);
FITS_t *initializeImageWriter_FITS_t(char *name);
void free_FITS_t(FITS_t *fits);
void report_FITS_t(FITS_t *fits);

//-- Functions related to reading HDU
void printNbHDU(FITS_t *fits);
void chooseHDU(FITS_t *fits, int nbHDU);
void printHDU(FITS_t *fits);

//-- Functions related to writing TABLE and BINTABLE HDU
void getTFORM(int format, int nb, char *string);
char *getTFORMComment(int format);
void addColumn(FITS_t *fits, char colName[], int format, char unit[]);
void addWideColumn(FITS_t *fits, char colName[], int format, int nb, char unit[]);

//-- Functions related to reading keywords
int readKeyword_int(FITS_t *fits, char key[]);
double readKeyword_double(FITS_t *fits, char key[]);
unsigned long readKeyword_ulong(FITS_t *fits, char key[]);

//-- Functions related to writing keywords
void addKeyword(FITS_t *fits, int format, char key[], void *value, char comment[]);
void updateComment(FITS_t *fits, char key[], char comment[]);
void addComment(FITS_t *fits, char comment[]);
void addLineSpread(FITS_t *fits);

//-- Functions related to reading TABLE and BINTABLE extensions
void readTableInfo(FITS_t *fits);
void readTableColumnType(FITS_t *fits, int colInd);
void readTableColumn(FITS_t *fits, int colInd, void *stockArr);

//-- Functions related to reading IMAGE extensions
void readImageInfo(FITS_t *fits);
void read2DImage(FITS_t *fits, void *stockArr);
void read2DSubImage(FITS_t *fits, void *stockArr, int limit[4]);

//-- Functions related to writing TABLE and BINTABLE extensions
void writeTableColumn(FITS_t *fits, int colInd, int nbData, void *data);
void writeTableColumn_value(FITS_t *fits, int colInd, int start, int nbData, void *value);
void nextRow(FITS_t *fits);
void updateNbRows(FITS_t *fits, int nbData);
void copyTableRow(FITS_t *from, FITS_t *to, long start, long nbRows);

//-- Functions related to writing IMAGE extension
void write2DImage_double(FITS_t *fits, int inputFormat, int outputFormat, int N1, int N2, void *data);

//-- Functions related to HPMap_t
HPMap_t *initialize_HPMap_t(long long nsidePix);
void free_HPMap_t(HPMap_t *full);
void read_HPMap_t(char name[], HPMap_t *full, int verbose);
void outFits_HPMap_t(char name[], HPMap_t *full, int verbose);

//-- Functions related to fast manipulation
void outFits2DImage(char name[], int N1, int N2, double *image, int verbose);

#endif


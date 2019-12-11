

  //------------------------------------------------------//
  //--  FITSFunctions.c					--//
  //--  Version 2018.09.25				--//
  //--  						--//
  //--  Copyright (C) 2018 - Chieh-An Lin		--//
  //--  GNU GPLv3 - https://www.gnu.org/licenses/	--//
  //------------------------------------------------------//

  
#include "FITSFunctions.h"


//----------------------------------------------------------------------
//-- Functions related to initialization

FITS_t *initializeTableReader_FITS_t(char *name)
{
  FITS_t *fits = malloc(sizeof(FITS_t));
  fits->nbRows = 0;
  fits->status = 0;
  fits_open_table(&fits->file, name, READONLY, &fits->status);
  
  readTableInfo(fits);
  report_FITS_t(fits);
  return fits;
}

FITS_t *initializeImageReader_FITS_t(char *name)
{
  FITS_t *fits = malloc(sizeof(FITS_t));
  fits->bitpix = 0;
  fits->status = 0;
  fits_open_image(&fits->file, name, READONLY, &fits->status);
  
  readImageInfo(fits);
  report_FITS_t(fits);
  return fits;
}

FITS_t *initializeTableWriter_FITS_t(char *name)
{
  FITS_t *fits = malloc(sizeof(FITS_t));
  fits->nbColumns = 0;
  fits->nbRows    = 0;
  fits->status    = 0;
  
  char name2[STRING_LENGTH_MAX];
  sprintf(name2, "!%s", name);
  fits_create_file(&fits->file, name2, &fits->status);
  fits_create_tbl(fits->file, BINARY_TBL, 0, 0, NULL, NULL, NULL, NULL, &fits->status);
  report_FITS_t(fits);
  return fits;
}

FITS_t *initializeImageWriter_FITS_t(char *name)
{
  FITS_t *fits = malloc(sizeof(FITS_t));
  fits->bitpix = 0;
  fits->status = 0;
  
  char name2[STRING_LENGTH_MAX];
  sprintf(name2, "!%s", name);
  fits_create_file(&fits->file, name2, &fits->status);
  fits_create_img(fits->file, BYTE_IMG, 0, NULL, &fits->status); //-- This is a dummy image; naxis = 0, naxes = NULL
  report_FITS_t(fits);
  return fits;
}

void free_FITS_t(FITS_t *fits)
{
  if (fits) {
    fits_close_file(fits->file, &fits->status);
    report_FITS_t(fits);
    free(fits); fits = NULL;
  }
  return;
}

void report_FITS_t(FITS_t *fits)
{
  if (fits->status) fits_report_error(stderr, fits->status);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to reading HDU

void printNbHDU(FITS_t *fits)
{
  fits_get_num_hdus(fits->file, &fits->nbHDU, &fits->status);
  printf("Number of HDU = %d\n", fits->nbHDU);
  report_FITS_t(fits);
  return;
}

void chooseHDU(FITS_t *fits, int nbHDU)
{
  fits_movabs_hdu(fits->file, nbHDU+1, &fits->HDUType, &fits->status);
  printf("HDU No.%d is chosen.\n", nbHDU);
  report_FITS_t(fits);
  return;
}

void printHDU(FITS_t *fits)
{
  char key[STRING_LENGTH_MAX];
  int i, nbKeys;
  fits_get_hdrspace(fits->file, &nbKeys, NULL, &fits->status);
  
  for (i=1; i<=nbKeys; i++)  { 
    fits_read_record(fits->file, i, key, &fits->status); 
    printf("%s\n", key);
  }
  printf("END\n");
  
  report_FITS_t(fits);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to writing TABLE and BINTABLE HDU

void getTFORM(int format, int nb, char *string)
{
  if (nb > 1) sprintf(string, "%d", nb);
  else        sprintf(string, "%s", "");
  
  if      (format == TSHORT)    sprintf(string, "%s%c", string, 'I');
  else if (format == TINT)      sprintf(string, "%s%c", string, 'J');
  else if (format == TINT32BIT) sprintf(string, "%s%c", string, 'J');
  else if (format == TLONGLONG) sprintf(string, "%s%c", string, 'K');
  else if (format == TFLOAT)    sprintf(string, "%s%c", string, 'E');
  else if (format == TDOUBLE)   sprintf(string, "%s%c", string, 'D');
  else if (format == TSTRING)   sprintf(string, "%s%c", string, 'A');
  else                          sprintf(string, "%s%c", string, 'A');
  return;
}

char *getTFORMComment(int format)
{
  if (format == TSHORT)    return "2-byte integer";
  if (format == TINT)      return "4-byte integer";
  if (format == TINT32BIT) return "4-byte integer";
  if (format == TLONGLONG) return "8-byte integer";
  if (format == TFLOAT)    return "4-byte real number";
  if (format == TDOUBLE)   return "8-byte real number";
  if (format == TSTRING)   return "string";
  return "string";
}

void addColumn(FITS_t *fits, char colName[], int format, char unit[])
{
  char buffer[STRING_LENGTH_MAX];
  getTFORM(format, 1, buffer); //-- nb = 1
  fits_insert_col(fits->file, fits->nbColumns+1, colName, buffer, &fits->status);
  
  sprintf(buffer, "TTYPE%d", fits->nbColumns+1);
  updateComment(fits, buffer, NULL);
  
  sprintf(buffer, "TFORM%d", fits->nbColumns+1);
  updateComment(fits, buffer, getTFORMComment(format));
  
  sprintf(buffer, "TUNIT%d", fits->nbColumns+1);
  addKeyword(fits, TSTRING, buffer, unit, NULL);
  fits->nbColumns++;
  return;
}

void addWideColumn(FITS_t *fits, char colName[], int format, int nb, char unit[])
{
  char buffer[STRING_LENGTH_MAX];
  getTFORM(format, nb, buffer);
  fits_insert_col(fits->file, fits->nbColumns+1, colName, buffer, &fits->status);
  
  sprintf(buffer, "TTYPE%d", fits->nbColumns+1);
  updateComment(fits, buffer, NULL);
  
  sprintf(buffer, "TFORM%d", fits->nbColumns+1);
  updateComment(fits, buffer, getTFORMComment(format));
  
  sprintf(buffer, "TUNIT%d", fits->nbColumns+1);
  addKeyword(fits, TSTRING, buffer, unit, NULL);
  fits->nbColumns++;
  return;
}

//----------------------------------------------------------------------
//-- Functions related to reading keywords

int readKeyword_int(FITS_t *fits, char key[])
{
  int value;
  fits_read_key(fits->file, TINT, key, &value, NULL, &fits->status);
  report_FITS_t(fits);
  return value;
}

double readKeyword_double(FITS_t *fits, char key[])
{
  double value;
  fits_read_key(fits->file, TDOUBLE, key, &value, NULL, &fits->status);
  report_FITS_t(fits);
  return value;
}

unsigned long readKeyword_ulong(FITS_t *fits, char key[])
{
  unsigned long value;
  fits_read_key(fits->file, TULONG, key, &value, NULL, &fits->status);
  report_FITS_t(fits);
  return value;
}

//----------------------------------------------------------------------
//-- Functions related to writing keywords

void addKeyword(FITS_t *fits, int format, char key[], void *value, char comment[])
{
  fits_update_key(fits->file, format, key, value, comment, &fits->status);
  report_FITS_t(fits);
  return;
}

void updateComment(FITS_t *fits, char key[], char comment[])
{
  //-- Use "TTYPEn" for key[]
  fits_modify_comment(fits->file, key, comment, &fits->status);
  report_FITS_t(fits);
  return;
}

void addComment(FITS_t *fits, char comment[])
{
  fits_write_comment(fits->file, comment, &fits->status);
  report_FITS_t(fits);
  return;
}

void addLineSpread(FITS_t *fits)
{
  addComment(fits, " ");
  report_FITS_t(fits);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to reading TABLE and BINTABLE extensions

void readTableInfo(FITS_t *fits)
{
  fits_get_hdu_type(fits->file, &fits->HDUType, &fits->status);
  fits_get_num_cols(fits->file, &fits->nbColumns, &fits->status);
  fits_get_num_rows(fits->file, &fits->nbRows, &fits->status);
  report_FITS_t(fits);
  return;
}

void readTableColumnType(FITS_t *fits, int colInd)
{
  long width;
  fits_get_coltype(fits->file, colInd+1, &fits->colType, &fits->colRepeat, &width, &fits->status);
  report_FITS_t(fits);
  return;
}

void readTableColumn(FITS_t *fits, int colInd, void *stockArr)
{
  if (fits->nbRows == 0) readTableInfo(fits);
  readTableColumnType(fits, colInd);
  if (fits->colType == 41) fits->colType = 31; //-- Correction for int
  fits_read_col(fits->file, fits->colType, colInd+1 , 1, 1, fits->nbRows*fits->colRepeat, NULL, stockArr, NULL, &fits->status); //-- firstrow = 1, firstelem = 1
  report_FITS_t(fits);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to reading IMAGE extensions

void readImageInfo(FITS_t *fits)
{
  int dim;
  fits_get_hdu_type(fits->file, &fits->HDUType, &fits->status);
  fits_get_img_param(fits->file, 2, &fits->bitpix, &dim, fits->resol, &fits->status);
  if      (fits->bitpix == FLOAT_IMG)    fits->colType = TFLOAT;
  else if (fits->bitpix == DOUBLE_IMG)   fits->colType = TDOUBLE;
  else if (fits->bitpix == SHORT_IMG)    fits->colType = TSHORT;
  else if (fits->bitpix == LONG_IMG)     fits->colType = TINT32BIT;
  else if (fits->bitpix == LONGLONG_IMG) fits->colType = TLONGLONG;
  else {
    printf("FITS image type error\n");
    exit(1);
  }
  report_FITS_t(fits);
  return;
}

void read2DImage(FITS_t *fits, void *stockArr)
{
  if (fits->bitpix == 0) readImageInfo(fits);
  long length = fits->resol[0] * fits->resol[1];
  int anyNull;
  fits_read_img(fits->file, fits->colType, 1, length, NULL, stockArr, &anyNull, &fits->status);
  report_FITS_t(fits);
  return;
}

void read2DSubImage(FITS_t *fits, void *stockArr, int limit[4])
{
  if (fits->bitpix == 0) readImageInfo(fits);
  
  long begin[2], end[2], inc[2];
  begin[0] = (long)limit[0] + 1;
  begin[1] = (long)limit[2] + 1;
  end[0]   = (long)limit[1];
  end[1]   = (long)limit[3];
  inc[0]   = 1;
  inc[1]   = 1;
  
  int anyNull;
  fits_read_subset(fits->file, fits->colType, begin, end, inc, NULL, stockArr, &anyNull, &fits->status);
  report_FITS_t(fits);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to writing TABLE and BINTABLE extensions

void writeTableColumn(FITS_t *fits, int colInd, int nbData, void *data)
{
  //-- Write value into column
  readTableColumnType(fits, colInd);
  fits_write_col(fits->file, fits->colType, colInd+1, fits->nbRows+1, 1, nbData, data, &fits->status);
  report_FITS_t(fits);
  return;
}

void writeTableColumn_value(FITS_t *fits, int colInd, int start, int nbData, void *value)
{
  //-- Write a constant "value" into column from start to start + nbData
  readTableColumnType(fits, colInd);
  int i;
  for (i=0; i<nbData; i++) {
    fits_write_col(fits->file, fits->colType, colInd+1, start+i+1, 1, 1, value, &fits->status);
    report_FITS_t(fits);
  }
  return;
}

void nextRow(FITS_t *fits)
{
  fits->nbRows += 1;
  return;
}

void updateNbRows(FITS_t *fits, int nbData)
{
  fits->nbRows += nbData;
  return;
}

void copyTableRow(FITS_t *from, FITS_t *to, long start, long nbRows)
{
  fits_copy_rows(from->file, to->file, start+1, nbRows, &from->status);
  to->nbRows += nbRows;
  report_FITS_t(from);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to writing IMAGE extension

void write2DImage_double(FITS_t *fits, int inputFormat, int outputFormat, int N1, int N2, void *data)
{
  long resol[2] = {(long)N1, (long)N2};
  fits_create_img(fits->file, outputFormat, 2, resol, &fits->status); //-- naxis = 2
  fits_write_img(fits->file, inputFormat, 1, N1*N2, data, &fits->status); //-- firstelem = 1
  report_FITS_t(fits);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to HPMap_t

HPMap_t *initialize_HPMap_t(long long nsidePix)
{
  HPMap_t *full  = (HPMap_t*)malloc(sizeof(HPMap_t));
  full->nsidePix = nsidePix;
  full->nbPix    = 12 * nsidePix * nsidePix;
  full->map      = (float*)malloc(full->nbPix * sizeof(float));
  return full;
}

void free_HPMap_t(HPMap_t *full)
{
  if (full) {
    if (full->map) {free(full->map); full->map = NULL;}
    free(full); full = NULL;
  }
  return;
}

void read_HPMap_t(char name[], HPMap_t *full, int verbose)
{
  FITS_t *fits = initializeTableReader_FITS_t(name);
  
  if (fits->nbRows * 1024 != full->nbPix) {
    printf("Wrong nside\n");
    return;
  }
  
  readTableColumn(fits, 0, (void*)full->map); //-- colInd = 0
  free_FITS_t(fits);
  if (verbose) printf("Read \"%s\"\n", name);
  return;
}

void outFits_HPMap_t(char name[], HPMap_t *full, int verbose)
{
  FITS_t *fits     = initializeTableWriter_FITS_t(name);
  long long nbRows = full->nbPix / 1024;
  float *map       = full->map;
  addWideColumn(fits, "VALUE", TFLOAT, 1024, "-       "); //-- nb = 1024
  
  int i;
  for (i=0; i<nbRows; i++) {
    writeTableColumn(fits, 0, 1024, map+i*1024); //-- colInd = 0, nbData = 1024
    nextRow(fits);
  }
  
  free_FITS_t(fits);
  if (verbose) printf("Outputed \"%s\"\n", name);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to fast manipulation

void outFits2DImage(char name[], int N1, int N2, double *image, int verbose)
{
  FITS_t *fits = initializeImageWriter_FITS_t(name);
  write2DImage_double(fits, TDOUBLE, FLOAT_IMG, N1, N2, (void*)image);
  free_FITS_t(fits);
  if (verbose) printf("Outputed \"%s\"\n", name);
  return;
}

//----------------------------------------------------------------------


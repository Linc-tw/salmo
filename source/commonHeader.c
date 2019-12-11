

  //------------------------------------------------------//
  //--  commonHeader.c					--//
  //--  Version 2019.04.15				--//
  //--  						--//
  //--  Copyright (C) 2018 - Chieh-An Lin		--//
  //--  GNU GPLv3 - https://www.gnu.org/licenses/	--//
  //------------------------------------------------------//


#include "commonHeader.h"


//----------------------------------------------------------------------
//-- Functions related to array

void reset_double(double *lfArr, int length)
{
  int i;
  for (i=0; i<length; i++) lfArr[i] = 0.0;
  return;
}

void rescale_double(double *lfArr, int length, double factor)
{
  int i;
  for (i=0; i<length; i++) lfArr[i] *= factor;
  return;
}

void reset_fftw_complex(fftw_complex *table, int length)
{
  double *c;
  int i;
  for (i=0; i<length; i++) {
    c = table[i];
    c[0] = 0.0;
    c[1] = 0.0;
  }
  return;
}

void rescaleReal_fftw_complex(fftw_complex *table, int length, double factor)
{
  int i;
  for (i=0; i<length; i++) table[i][0] *= factor;
  return;
}

void rescale_fftw_complex(fftw_complex *table, int length, double factor)
{
  double *c;
  int i;
  for (i=0; i<length; i++) {
    c = table[i];
    c[0] *= factor;
    c[1] *= factor;
  }
  return;
}

void multiplication_fftw_complex(fftw_complex *table1, fftw_complex *table2, fftw_complex *product, int length)
{
  double value, *c1, *c2;
  int i;
  for (i=0; i<length; i++) {
    c1 = table1[i];
    c2 = table2[i];
    value         = c1[0] * c2[0] - c1[1] * c2[1];
    product[i][1] = c1[0] * c2[1] + c1[1] * c2[0];
    product[i][0] = value;
  }
  return;
}

void copy_fftw_complex(fftw_complex *from, fftw_complex *to, int length)
{
  double *c1, *c2;
  int i;
  for (i=0; i<length; i++) {
    c1 = from[i];
    c2 = to[i];
    c2[0] = c1[0];
    c2[1] = c1[1];
  }
  return;
}

void print_fftw_complex(fftw_complex *table, int N1)
{
  int i, j;
  for (j=0; j<10; j++) {
    for (i=0; i<10; i++) {
      printf(" % .3e ", table[i+j*N1][0]);
    }
    printf("\n");
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to int_arr

int_arr *initialize_int_arr(int length)
{
  int_arr *intArr = (int_arr*)malloc(sizeof(int_arr));
  intArr->length  = length;
  intArr->array   = (int*)calloc(length, sizeof(int));
  return intArr;
}

void free_int_arr(int_arr *intArr)
{
  if (intArr) {
    if (intArr->array) {free(intArr->array); intArr->array = NULL;}
    free(intArr); intArr = NULL;
  }
  return;
}

void print_int_arr(int_arr *intArr)
{
  printf("# int_arr\n");
  
  int L = (int)fmin(intArr->length, 20);
  int i;
  for (i=0; i<L; i++) printf(" %4d ", intArr->array[i]);
  printf("\n");
  return;
}

//----------------------------------------------------------------------
//-- Functions related to int_mat

int_mat *initialize_int_mat(int N1, int N2)
{
  int_mat *intMat = (int_mat*)malloc(sizeof(int_mat));
  intMat->N1         = N1;
  intMat->N2         = N2;
  intMat->length     = N1 * N2;
  intMat->matrix     = (int*)calloc(intMat->length, sizeof(int));
  return intMat;
}

void free_int_mat(int_mat *intMat)
{
  if (intMat) {
    if (intMat->matrix) {free(intMat->matrix); intMat->matrix = NULL;}
    free(intMat); intMat = NULL;
  }
  return;
}

void print_int_mat(int_mat *intMat)
{
  printf("# int_mat\n");
  
  int L1 = (int)fmin(intMat->N1, 8);
  int L2 = (int)fmin(intMat->N2, 8);
  
  int i, j;
  for (j=0; j<L2; j++) {
    for (i=0; i<L1; i++) {
      printf(" %4d ", intMat->matrix[i+j*intMat->N1]);
    }
    printf("\n");
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to int_ten3

int_ten3 *initialize_int_ten3(int N1, int N2, int N3)
{
  int_ten3 *intTen = (int_ten3*)malloc(sizeof(int_ten3));
  intTen->N1          = N1;
  intTen->N2          = N2;
  intTen->N3          = N3;
  intTen->length      = N1 * N2 * N3;
  intTen->tensor      = (int*)calloc(intTen->length, sizeof(int));
  return intTen;
}

void free_int_ten3(int_ten3 *intTen)
{
  if (intTen) {
    if (intTen->tensor) {free(intTen->tensor); intTen->tensor = NULL;}
    free(intTen); intTen = NULL;
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to short_mat

short_mat *initialize_short_mat(int N1, int N2)
{
  short_mat *shrtMat = (short_mat*)malloc(sizeof(short_mat));
  shrtMat->N1        = N1;
  shrtMat->N2        = N2;
  shrtMat->length    = N1 * N2;
  shrtMat->matrix    = (short*)calloc(N1 * N2, sizeof(short));
  return shrtMat;
}

void free_short_mat(short_mat *shrtMat)
{
  if (shrtMat) {
    if (shrtMat->matrix) {free(shrtMat->matrix); shrtMat->matrix = NULL;}
    free(shrtMat); shrtMat = NULL;
  }
  return;
}

void print_short_mat(short_mat *shrtMat)
{
  printf("# short_mat\n");
  
  int L1 = (int)fmin(shrtMat->N1, 8);
  int L2 = (int)fmin(shrtMat->N2, 8);
  
  int i, j;
  for (j=0; j<L2; j++) {
    for (i=0; i<L1; i++) {
      printf(" %4d ", (int)shrtMat->matrix[i+j*shrtMat->N1]);
    }
    printf("\n");
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to double_arr

double_arr *initialize_double_arr(int length)
{
  double_arr *dblArr = (double_arr*)malloc(sizeof(double_arr));
  dblArr->length     = length;
  dblArr->array      = (double*)calloc(length, sizeof(double));
  return dblArr;
}

void free_double_arr(double_arr *dblArr)
{
  if (dblArr) {
    if (dblArr->array) {free(dblArr->array); dblArr->array = NULL;}
    free(dblArr); dblArr = NULL;
  }
  return;
}

void print_double_arr(double_arr *dblArr)
{
  printf("# double_arr\n");
  
  int L = (int)fmin(dblArr->length, 20);
  int i;
  for (i=0; i<L; i++) printf(" % .3e ", dblArr->array[i]);
  printf("\n");
  return;
}

//----------------------------------------------------------------------
//-- Functions related to double_mat

double_mat *initialize_double_mat(int N1, int N2)
{
  double_mat *dblMat = (double_mat*)malloc(sizeof(double_mat));
  dblMat->N1         = N1;
  dblMat->N2         = N2;
  dblMat->length     = N1 * N2;
  dblMat->matrix     = (double*)calloc(dblMat->length, sizeof(double));
  return dblMat;
}

void free_double_mat(double_mat *dblMat)
{
  if (dblMat) {
    if (dblMat->matrix) {free(dblMat->matrix); dblMat->matrix = NULL;}
    free(dblMat); dblMat = NULL;
  }
  return;
}

void print_double_mat(double_mat *dblMat)
{
  printf("# double_mat\n");
  
  int L1 = (int)fmin(dblMat->N1, 8);
  int L2 = (int)fmin(dblMat->N2, 8);
  
  int i, j;
  for (j=0; j<L2; j++) {
    for (i=0; i<L1; i++) {
      printf(" % .3e ", dblMat->matrix[i+j*dblMat->N1]);
    }
    printf("\n");
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to double_ten3

double_ten3 *initialize_double_ten3(int N1, int N2, int N3)
{
  double_ten3 *dblTen = (double_ten3*)malloc(sizeof(double_ten3));
  dblTen->N1          = N1;
  dblTen->N2          = N2;
  dblTen->N3          = N3;
  dblTen->length      = N1 * N2 * N3;
  dblTen->tensor      = (double*)calloc(dblTen->length, sizeof(double));
  return dblTen;
}

void free_double_ten3(double_ten3 *dblTen)
{
  if (dblTen) {
    if (dblTen->tensor) {free(dblTen->tensor); dblTen->tensor = NULL;}
    free(dblTen); dblTen = NULL;
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to float_ten5

float_ten5 *initialize_float_ten5(int N1, int N2, int N3, int N4, int N5)
{
  float_ten5 *fltTen = (float_ten5*)malloc(sizeof(float_ten5));
  fltTen->N1         = N1;
  fltTen->N2         = N2;
  fltTen->N3         = N3;
  fltTen->N4         = N4;
  fltTen->N5         = N5;
  fltTen->length     = N1 * N2 * N3 * N4 * N5;
  fltTen->tensor     = (float*)calloc(fltTen->length, sizeof(float));
  return fltTen;
}

void free_float_ten5(float_ten5 *fltTen)
{
  if (fltTen) {
    if (fltTen->tensor) {free(fltTen->tensor); fltTen->tensor = NULL;}
    free(fltTen); fltTen = NULL;
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to interpolator_t

interpolator_t *initialize_interpolator_t(int length)
{
  //-- length = number of points
  //-- dx     = interval width
  //-- *x     = coordinate
  //-- *value = value
  interpolator_t *inter = (interpolator_t*)malloc(sizeof(interpolator_t));
  inter->length         = length;
  inter->dx             = 0.0;
  inter->x              = (double*)calloc(length, sizeof(double));
  inter->value          = (double*)calloc(length, sizeof(double));
  return inter;
}

void free_interpolator_t(interpolator_t *inter)
{
  if (inter) {
    if (inter->x)     {free(inter->x);     inter->x     = NULL;}
    if (inter->value) {free(inter->value); inter->value = NULL;}
    free(inter); inter = NULL;
  }
  return;
}

void print_interpolator_t(interpolator_t *inter)
{
  printf("# interpolator_t\n");
  printf("# Number of points = %d\n", inter->length);
  printf("# Interval = %g\n", inter->dx);
  printf("#\n");
  printf("#        x      value\n");
  
  int i;
  for (i=0; i<inter->length; i++) printf(" % .3e  % .3e\n", inter->x[i], inter->value[i]);
  return;
}

double execute_interpolator_t(interpolator_t *inter, double x, int border)
{
  int length    = inter->length;
  double *coor  = inter->x;
  double *value = inter->value;
  double r;
  
  //-- 0-padding border
  if (border == 0) {
    if (x < coor[0] || x > coor[length-1]) return 0.0;
  }
  
  //-- Constant border
  else if (border == 1) {
    if (x < coor[0])        return value[0];
    if (x > coor[length-1]) return value[length-1];
  }
  
  //-- Linear extrapolation border
  else if (border == 2) {
    if (x < coor[0]) {
      r = (x - coor[0]) / (coor[1] - coor[0]); 
      return (1-r)*value[0] + r*value[1];
    }
    if (x > coor[length-1]) {
      r = (x - coor[length-2]) / (coor[length-1] - coor[length-2]); 
      return (1-r)*value[length-2] + r*value[length-1];
    }
  }
  
  //-- Exit with error
  else {
    if (x < coor[0] || x > coor[length-1]) {
      printf("Value out of range of interpolator\n");
      exit(1);
    }
  }
    
  int i;
  for (i=0; i<length-1; i++) {
    if (x < coor[i+1]) break;
  }
  r = (x - coor[i]) / (coor[i+1] - coor[i]); 
  return (1-r)*value[i] + r*value[i+1];
}

//----------------------------------------------------------------------
//-- Functions related to sampler_t

sampler_t *initialize_sampler_t(int length)
{
  //-- length  = number of points
  //-- dx      = interval width
  //-- *x      = bin values
  //-- *pdf    = normalized pdf
  //-- *cdf    = normalized cdf
  //-- totPdf  = integration over pdf, before or after normalization (see set_sampler_t)
  //-- x_mean  = integration over x*p(x), before or after normalization (see set_sampler_t)
  //-- 
  //-- Usage:
  //--   - initialize it with length
  //--   - fill dx, x[0], rest of x
  //--   - fill pdf (not necessarily normalized)
  //--   - use set_sampler_t to set
  
  sampler_t *samp = (sampler_t*)malloc(sizeof(sampler_t));
  samp->length    = length;
  samp->ind       = 0;
  samp->x         = (double*)calloc(length, sizeof(double));
  samp->pdf       = (double*)calloc(length, sizeof(double));
  samp->cdf       = (double*)calloc(length, sizeof(double));
  samp->x_mean    = 0.0;
  return samp;
}

void free_sampler_t(sampler_t *samp)
{
  if (samp) {
    if (samp->x)   {free(samp->x);   samp->x   = NULL;}
    if (samp->pdf) {free(samp->pdf); samp->pdf = NULL;}
    if (samp->cdf) {free(samp->cdf); samp->cdf = NULL;}
    free(samp); samp = NULL;
  }
  return;
}

void print_sampler_t(sampler_t *samp)
{
  printf("# sampler_t\n");
  int L = (int)fmin(samp->length, 20);
  int i;
  
  printf("# length = %d\n", samp->length);
  printf("# dx     = %g\n", samp->dx);
  printf("# x      =");
  for (i=0; i<L; i++) printf(" % .3e ", samp->x[i]);
  printf("\n");
  printf("# pdf    =");
  for (i=0; i<L; i++) printf(" % .3e ", samp->pdf[i]);
  printf("\n");
  printf("# cdf    =");
  for (i=0; i<L; i++) printf(" % .3e ", samp->cdf[i]);
  printf("\n");
  printf("# totPdf = %g\n", samp->totPdf);
  printf("# x_mean = %g\n", samp->x_mean);
  return;
}

void set_sampler_t(sampler_t *samp, int mode, int setTotalToOne)
{
  //-- mode:
  //--   0 = Continuous & linear-linear, interpolations are linear between bin edges.
  //--   1 = Continuous & linear-log, interpolations are linear between log value of each edges.
  //--   2 = Continuous & log-linear, bins are considered as log of true values, x_mean is calculated with true x values.
  //--   3 = Continuous & log-log, interpolations are linear between log value of log edges.
  //--   4 = Continuous histogram & linear, inputs are histogram not pdf, bins are lower edges, interpolations are linear.
  //--   5 = Discrete histogram, no interpolation
  //-- When doLog = 1, the integration over [x_1, x_2] is (x_2 - x_1)(y_2 - y_1)/ln(y_2/y_1).
  //--
  //-- Allow non-equal binwidth
  //-- Let p(x) = input pdf, q(x) = normalized pdf.
  //--
  //-- If setTotalToOne = 0, totPdf and x_mean will continue to contain physical information:
  //--   totPdf = integration over p(x)
  //--   x_mean = integration over x*p(x)
  //-- If setTotalToOne = 1, totPdf and x_mean will be normalized, so:
  //--   totPdf = integration over q(x) = 1.0
  //--   x_mean = integration over x*q(x)
  //-- In both cases, samp->pdf is always normalized: q(x).
  
  double *x     = samp->x;
  double *pdf   = samp->pdf;
  double *cdf   = samp->cdf;
  int length    = samp->length;
  double x_mean = 0.0;
  double trapeze, xPdf1, xPdf2, totPdf;
  int i;
  
  //-- Compute cdf, total pdf and <x>
  cdf[0] = 0.0;
  
  //-- 0 = Continuous & linear-linear
  if (mode == 0) {
    for (i=1; i<length; i++) {
      trapeze = 0.5 * fabs(x[i] - x[i-1]) * (pdf[i-1] + pdf[i]); //-- Trapeze in linear space
      cdf[i]  = cdf[i-1] + trapeze;
      x_mean += 0.5 * fabs(x[i] - x[i-1]) * (pdf[i-1] * x[i-1] + pdf[i] * x[i]);
    }
    totPdf = cdf[length-1];
  }
  
  //-- 1 = Continuous & linear-log
  else if (mode == 1) {
    for (i=1; i<length; i++) {
      trapeze = (pdf[i-1] == 0 || pdf[i] == 0) ? 0.0 : (fabs(x[i] - x[i-1]) * (pdf[i] - pdf[i-1]) / log(pdf[i] / pdf[i-1])); //-- Trapeze in log space
      cdf[i]  = cdf[i-1] + trapeze;
      xPdf1   = pdf[i-1] * x[i-1];
      xPdf2   = pdf[i] * x[i];
      x_mean += (xPdf1 == 0 || xPdf2 == 0) ? 0.0 : (fabs(x[i] - x[i-1]) * (xPdf2 - xPdf1) / log(xPdf2 / xPdf1));
    }
    totPdf = cdf[length-1];
  }
  
  //-- 2 = Continuous & log-linear
  else if (mode == 2) {
    for (i=1; i<length; i++) {
      trapeze = 0.5 * fabs(x[i] - x[i-1]) * (pdf[i-1] + pdf[i]); //-- Trapeze in linear space
      cdf[i]  = cdf[i-1] + trapeze;
      x_mean += 0.5 * fabs(x[i] - x[i-1]) * (pdf[i-1] * pow(10, x[i-1]) + pdf[i] * pow(10, x[i]));
    }
    totPdf = cdf[length-1];
  }
  
  //-- 3 = Continuous & log-log (mass function)
  else if (mode == 3) {
    for (i=1; i<length; i++) {
      trapeze = (pdf[i-1] == 0 || pdf[i] == 0) ? 0.0 : (fabs(x[i] - x[i-1]) * (pdf[i] - pdf[i-1]) / log(pdf[i] / pdf[i-1])); //-- Trapeze in log space
      cdf[i]  = cdf[i-1] + trapeze;
      xPdf1   = pdf[i-1] * pow(10, x[i-1]);
      xPdf2   = pdf[i] * pow(10, x[i]);
      x_mean += (xPdf1 == 0 || xPdf2 == 0) ? 0.0 : (fabs(x[i] - x[i-1]) * (xPdf2 - xPdf1) / log(xPdf2 / xPdf1));
    }
    totPdf = cdf[length-1];
  }
  
  //-- 4 = Continuous histogram & linear
  else if (mode == 4) {
    x_mean += pdf[0] * 0.5 * (x[0] + x[1]);
    for (i=1; i<length; i++) {
      cdf[i] += cdf[i-1] + pdf[i-1];
      x_mean += pdf[i] * 0.5 * (x[i] + x[i+1]);
    }
    totPdf = cdf[length-1] + pdf[length-1];
  }
  
  //-- 5 = Discrete histogram
  else {
    x_mean += pdf[0] * x[0];
    for (i=1; i<length; i++) {
      cdf[i] += cdf[i-1] + pdf[i-1];
      x_mean += pdf[i] * x[i];
    }
    totPdf = cdf[length-1] + pdf[length-1];
  }
  
  samp->totPdf = (setTotalToOne == 1) ? 1.0 : totPdf;
  samp->x_mean = (setTotalToOne == 1) ? x_mean / totPdf : x_mean;
  
  //-- Normalization
  for (i=0; i<length; i++) {
    pdf[i] /= totPdf;
    cdf[i] /= totPdf;
  }
  return;
}

double execute_sampler_t(sampler_t *samp, double x, int mode)
{
  //-- mode:
  //--   0 = Continuous & linear-linear, interpolations are linear between bin edges.
  //--   1 = Continuous & linear-log, interpolations are linear between log value of each edges.
  //--   2 = Continuous & log-linear, bins are considered as log of true values, x_mean is calculated with true x values.
  //--   3 = Continuous & log-log, interpolations are linear between log value of log edges.
  //--   4 = Continuous histogram & linear, inputs are histogram not pdf, bins are lower edges, interpolations are linear.
  //--   5 = Discrete histogram, no interpolation
  //-- When doLog = 1, the integration over [x_1, x_2] is (x_2 - x_1)(y_2 - y_1)/ln(y_2/y_1).
  
  int length    = samp->length;
  double *cdf   = samp->cdf;
  double *value = samp->x;
    
  int i;
  for (i=0; i<=length-2; i++) {
    if (x < cdf[i+1]) break;
  }
  
  double r;
  
  //-- 0 = Continuous & linear-linear
  //-- 2 = Continuous & log-linear
  //-- 4 = Continuous histogram & linear
  if (mode == 0 || mode == 2 || mode == 4) {
    if (i == length-1) i -= 1;
    r = (x - cdf[i]) / (cdf[i+1] - cdf[i]); 
    return (1-r)*value[i] + r*value[i+1];
  }
  
  //-- 1 = Continuous & linear-log
  //-- 3 = Continuous & log-log (mass function)
  if (mode == 1 || mode == 3) {
    if (i == length-1) i -= 1;
    r = (x - cdf[i]) / (cdf[i+1] - cdf[i]); 
    return pow(value[i], 1-r) * pow(value[i+1], r);
  }
  
  //-- 5 = Discrete histogram
  return samp->x[i];
}

//----------------------------------------------------------------------
//-- Functions related to sampler_arr

sampler_arr *initialize_sampler_arr(int N_type, int N_array)
{
  sampler_arr *sampArr = (sampler_arr*)malloc(sizeof(sampler_arr));
  sampArr->length      = N_array;
  sampArr->array       = (sampler_t**)malloc(N_array * sizeof(sampler_t*));
  int i;
  for (i=0; i<N_array; i++) {
    sampArr->array[i]      = initialize_sampler_t(N_type);
    sampArr->array[i]->ind = i;
  }
  return sampArr;
}

void free_sampler_arr(sampler_arr *sampArr)
{
  int i;
  if (sampArr) {
    if (sampArr->array) {
      for (i=0; i<sampArr->length; i++) {free_sampler_t(sampArr->array[i]); sampArr->array[i] = NULL;}
      free(sampArr->array); sampArr->array = NULL;
    }
    free(sampArr); sampArr = NULL;
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to sampler_mat

sampler_mat *initialize_sampler_mat(int N_type, int N1_array, int N2_array)
{
  sampler_mat *sampMat = (sampler_mat*)malloc(sizeof(sampler_mat));
  sampMat->N1          = N1_array;
  sampMat->N2          = N2_array;
  sampMat->length      = N1_array * N2_array;
  sampMat->matrix      = (sampler_t**)malloc(sampMat->length * sizeof(sampler_t*));
  int i;
  for (i=0; i<sampMat->length; i++) {
    sampMat->matrix[i]      = initialize_sampler_t(N_type);
    sampMat->matrix[i]->ind = i;
  }
  return sampMat;
}

void free_sampler_mat(sampler_mat *sampMat)
{
  int i;
  if (sampMat) {
    if (sampMat->matrix) {
      for (i=0; i<sampMat->length; i++) {free_sampler_t(sampMat->matrix[i]); sampMat->matrix[i] = NULL;}
      free(sampMat->matrix); sampMat->matrix = NULL;
    }
    free(sampMat); sampMat = NULL;
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to FFT_t

FFT_t *initialize_FFT_t(int N)
{
  //-- 'before' and 'kernel' should be set to 0 at the beginning.
  
  FFT_t *transformer      = (FFT_t*)malloc(sizeof(FFT_t));
  transformer->N          = N;
  transformer->length     = N * N;
  transformer->normFactor = 1.0 / (double)transformer->length;
  transformer->ind        = 0;
  
  //-- Allocate fftw_complex elements
  transformer->before     = (fftw_complex*)fftw_malloc(transformer->length * sizeof(fftw_complex));
  transformer->kernel     = (fftw_complex*)fftw_malloc(transformer->length * sizeof(fftw_complex));
  transformer->after      = (fftw_complex*)fftw_malloc(transformer->length * sizeof(fftw_complex));
  
  //-- Allocate fftw_plan elements
  transformer->before_f   = fftw_plan_dft_2d(N, N, transformer->before, transformer->before, FFTW_FORWARD,  FFTW_ESTIMATE);
  transformer->kernel_f   = fftw_plan_dft_2d(N, N, transformer->kernel, transformer->kernel, FFTW_FORWARD,  FFTW_ESTIMATE);
  transformer->after_b    = fftw_plan_dft_2d(N, N, transformer->after,  transformer->after,  FFTW_BACKWARD, FFTW_ESTIMATE);
  transformer->before_b   = fftw_plan_dft_2d(N, N, transformer->before, transformer->before, FFTW_BACKWARD, FFTW_ESTIMATE);
  
  return transformer;
}

void free_FFT_t(FFT_t *transformer)
{
  if (transformer) {
    if (transformer->before)   {fftw_free(transformer->before);           transformer->before   = NULL;}
    if (transformer->kernel)   {fftw_free(transformer->kernel);           transformer->kernel   = NULL;}
    if (transformer->after)    {fftw_free(transformer->after);            transformer->after    = NULL;}
    
    if (transformer->before_f) {fftw_destroy_plan(transformer->before_f); transformer->before_f = NULL;}
    if (transformer->kernel_f) {fftw_destroy_plan(transformer->kernel_f); transformer->kernel_f = NULL;}
    if (transformer->after_b)  {fftw_destroy_plan(transformer->after_b);  transformer->after_b  = NULL;}
    if (transformer->before_b) {fftw_destroy_plan(transformer->before_b); transformer->before_b = NULL;}
    
    free(transformer); transformer = NULL;
  }
  return;
}

void reset_FFT_t(FFT_t *transformer)
{
  reset_fftw_complex(transformer->before, transformer->length);
  return;
}

void execute_FFT_t(FFT_t *transformer)
{
  int length    = transformer->length;
  double factor = 1.0 / (double)length;
  fftw_complex *after = transformer->after;
  
  //-- Go to Fourier space
  fftw_execute(transformer->before_f);
  
  //-- Multiplication
  multiplication_fftw_complex(transformer->before, transformer->kernel, after, length);
  
  //-- Go to direct space
  fftw_execute(transformer->after_b);
  
  //-- Rescale
  rescale_fftw_complex(after, length, factor);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to FFT_arr

FFT_arr *initialize_FFT_arr(int N_type, int N_array)
{
  FFT_arr *transArr = (FFT_arr*)malloc(sizeof(FFT_arr));
  transArr->length  = N_array;
  transArr->array   = (FFT_t**)malloc(N_array * sizeof(FFT_t*));
  int i;
  for (i=0; i<N_array; i++) {
    transArr->array[i]      = initialize_FFT_t(N_type);
    transArr->array[i]->ind = i;
  }
  return transArr;
}

void free_FFT_arr(FFT_arr *transArr)
{
  int i;
  if (transArr) {
    if (transArr->array) {
      for (i=0; i<transArr->length; i++) {free_FFT_t(transArr->array[i]); transArr->array[i] = NULL;}
      free(transArr->array); transArr->array = NULL;
    }
    free(transArr); transArr = NULL;
  }
  return;
}

void reset_FFT_arr(FFT_arr *transArr)
{ 
  int i;
  for (i=0; i<transArr->length; i++) reset_FFT_t(transArr->array[i]);
  return;
}

void execute_FFT_arr(FFT_arr *transArr)
{
  int i;
  for (i=0; i<transArr->length; i++) execute_FFT_t(transArr->array[i]);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to hist_t

hist_t *initialize_hist_t(int length)
{
  //-- length   = number of bins
  //-- dx       = bin width
  //-- n_tot    = total counts
  //-- *x_lower = lower limit of each bin
  //-- x_max;   = maximal limit
  //-- *n       = histogram
  hist_t *hist  = (hist_t*)malloc(sizeof(hist_t));
  hist->length  = length;
  hist->dx      = 0.0;
  hist->x_lower = (double*)calloc(length, sizeof(double));
  hist->x_max   = 0.0;
  hist->n       = (int*)calloc(length, sizeof(int));
  hist->n_tot   = 0;
  return hist;
}

void free_hist_t(hist_t *hist)
{
  if (hist) {
    if (hist->x_lower) {free(hist->x_lower); hist->x_lower = NULL;}
    if (hist->n)       {free(hist->n);       hist->n       = NULL;}
    free(hist); hist = NULL;
  }
  return;
}

void print_hist_t(hist_t *hist)
{
  printf("# hist_t\n");
  printf("# Number of bins = %d\n", hist->length);
  printf("# Bin width = %g\n", hist->dx);
  printf("# Total counts = %d\n", hist->n_tot);
  printf("#\n");
  printf("#                   bin        n\n");
  
  int i;
  for (i=0; i<hist->length-1; i++) printf("  [% .3e, % .3e]  %8d\n", hist->x_lower[i], hist->x_lower[i+1], hist->n[i]);
  printf("  [% .3e, % .3e]  %8d\n", hist->x_lower[hist->length-1], hist->x_max, hist->n[hist->length-1]);
  return;
}

void set_hist_t(hist_t *hist, double x_min, double x_max)
{
  int length  = hist->length;
  double dx   = (x_max - x_min) / (double)length;
  hist->dx    = dx;
  hist->n_tot = 0;
  int i;
  for (i=0; i<length; i++) {
    hist->x_lower[i] = x_min + i * dx;
    hist->n[i]       = 0;
  }
  hist->x_max = hist->x_lower[length-1] + dx;
  return;
}

void reset_hist_t(hist_t *hist)
{
  int i;
  for (i=0; i<hist->length; i++) hist->n[i] = 0;
  hist->n_tot = 0;
  return;
}

void push_hist_t(hist_t *hist, double x, int verbose)
{
  double *x_lower = hist->x_lower;
  if (x < x_lower[0] || x >= hist->x_max) {
    if (verbose == 1) printf("Value out of range, not pushed into histogram\n");
    return;
  }
  
  int i;
  for (i=0; i<hist->length-1; i++) {
    if (x_lower[i+1] > x) break;
  }
  hist->n[i]  += 1;
  hist->n_tot += 1;
  return;
}

//----------------------------------------------------------------------
//-- Functions related to KDE_t

KDE_t *initialize_KDE_t(int length)
{
  //-- length    = number of samples
  //-- mean      = mean
  //-- std       = standard deviation
  //-- two_h_sq  = 2 * bandwidth^2
  //-- amplitude = common amplitude factor
  //-- *sample   = sample points
  KDE_t *estimator  = (KDE_t*)malloc(sizeof(KDE_t));
  estimator->length = length;
  estimator->sample = (double*)malloc(length * sizeof(double));
  return estimator;
}

void free_KDE_t(KDE_t *estimator)
{
  if (estimator) {
    if (estimator->sample) {free(estimator->sample); estimator->sample = NULL;}
    free(estimator); estimator = NULL;
  }
  return;
}

void set_KDE_t(KDE_t *estimator, double h)
{
  estimator->mean     = gsl_stats_mean(estimator->sample, 1, estimator->length);
  estimator->variance = gsl_stats_variance_m(estimator->sample, 1, estimator->length, estimator->mean);
  
  if (h > 0) estimator->two_h_sq = 2 * h * h;
  else {
    //-- Use Silverman's rule
    double factor       = 4.0 / (3.0 * estimator->length);
    estimator->two_h_sq = 2 * pow(factor, 0.4) * estimator->variance;
  }
  estimator->amplitude = 1.0 / (estimator->length * sqrt(PI * estimator->two_h_sq));
  
  return;
}

double execute_KDE_t(KDE_t *estimator, double x)
{
  double dist_sq, sum = 0.0;
  int i;
  for (i=0; i<estimator->length; i++) {
    dist_sq = pow(x - estimator->sample[i], 2);
    sum += exp(-dist_sq / estimator->two_h_sq);
  }
  return estimator->amplitude * sum;
}

double integrate_KDE_t(KDE_t *estimator, int N, double x)
{
  //-- Fixed path at 10-sigma / N
  double dx  = 10 * sqrt(estimator->variance) / (double)N;
  double sum = 0.0;
  
  int i;
  if (x > estimator->mean) {
    sum += 0.5 * (execute_KDE_t(estimator, x + N*dx) + execute_KDE_t(estimator, x));
    for (i=1; i<N; i++) sum += execute_KDE_t(estimator, x + i*dx);
    sum *= dx;
    sum = 1 - sum;
  }
  else {
    sum += 0.5 * (execute_KDE_t(estimator, x - N*dx) + execute_KDE_t(estimator, x));
    for (i=1; i<N; i++) sum += execute_KDE_t(estimator, x - i*dx);
    sum *= dx;
  }
  
  return sum;
}

//----------------------------------------------------------------------
//-- Functions related to RNG

u_int32_t renewSeed()
{
  FILE *file = fopen("/dev/urandom","r");
  u_int32_t seed;
  if (file == NULL) seed = 0;
  else {
    fread(&seed, sizeof(u_int32_t), 1, file);
    fclose(file);
  }
  return seed;
}

gsl_rng *initializeGenerator(u_int32_t seed)
{
  gsl_rng *generator = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(generator, seed);
  return generator;
}

void printGenerator(gsl_rng *generator)
{
  printf("Generator type = %s\n", gsl_rng_name(generator));
  printf("Seed range = [%lu, %lu]\n", gsl_rng_min(generator), gsl_rng_max(generator));
  return;
}

//-- Use the following functions:
//--   double gsl_ran_flat(const gsl_rng *r, double a, double b)
//--   double gsl_ran_gaussian(const gsl_rng *r, double sigma)

//----------------------------------------------------------------------
//-- Functions related to rotation & projection

void rotate(double oldPos[2], double rotAng, double newPos[2])
{
  //-- rotAng in [rad]
  double sinTheta = sin(rotAng);
  double cosTheta = cos(rotAng);
  double buffer = cosTheta * oldPos[0] - sinTheta * oldPos[1];
  newPos[1] = sinTheta * oldPos[0] + cosTheta * oldPos[1];
  newPos[0] = buffer;
  return;
}

void project(double RADEC[2], double center[4], double thetaXY[2])
{
  //-- RADEC is (RA, DEC) of the point to be projected, in [rad]
  //-- center is (RA, DEC, sin(DEC), cos(DEC)) of the projection center, in [rad]
  //-- thetaXY is (theta_x, theta_y) after projection, in [rad]
  
  double sin_DEC     = sin(RADEC[1]);
  double cos_DEC     = cos(RADEC[1]);
  double dRA         = RADEC[0] - center[0];
  double cos_dRA     = cos(dRA);
  double denominator = sin_DEC * center[2] + cos_DEC * center[3] * cos_dRA;
  thetaXY[0] = cos_DEC * sin(dRA) / denominator;
  thetaXY[1] = (center[3] * sin_DEC - center[2] * cos_DEC * cos_dRA) / denominator;
  return;
}

void deproject(double thetaXY[2], double center[4], double RADEC[2])
{
  //-- thetaXY is (theta_x, theta_y) of the point to be deprojected, in [rad]
  //-- center is (RA, DEC, sin(DEC), cos(DEC)) of the projection center, in [rad]
  //-- RADEC is (RA, DEC) after deprojection, in [rad]
  
  double sin_ctrRA = sin(center[0]);
  double cos_ctrRA = cos(center[0]);
  double value1 = center[3] - thetaXY[1] * center[2];
  double value2 = (center[2] + thetaXY[1] * center[3]) / sqrt(1.0 + SQ(thetaXY[0]) + SQ(thetaXY[1]));
  RADEC[0] = atan2(sin_ctrRA * value1 + thetaXY[0] * cos_ctrRA, cos_ctrRA * value1 - thetaXY[0] * sin_ctrRA);
  RADEC[1] = HALF_PI - acos(value2);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to comparison

int compare_float(const void *a, const void *b)
{
  float a2 = *(float*)a;
  float b2 = *(float*)b;
  return (a2 > b2) - (a2 < b2);
}

int compare_longLong(const void *a, const void *b)
{
  long long a2 = *(long long*)a;
  long long b2 = *(long long*)b;
  return (a2 > b2) - (a2 < b2);
}

//----------------------------------------------------------------------
//-- Functions related to stopwatch

//-- Use the following functions:
//--   clock_t start = clock();
//--   clock_t finish = clock();
//--   double duration = (double)(finish - start) / CLOCKS_PER_SEC;

void printTime(clock_t start, clock_t stop)
{
  double duration = (double)(stop - start) / CLOCKS_PER_SEC;
  int hours       = (int)(duration / 3600);
  double remain   = fmod(duration, 3600);
  int minutes     = (int)(remain / 60);
  double seconds  = fmod(remain, 60);
  
  if (hours == 0) {
    if (minutes == 0) printf("Computation time = %.2f secs\n", duration);
    else              printf("Computation time = %d m %d s\n", minutes, (int)seconds);
  }
  else                printf("Computation time = %d h %d m %d s\n", hours, minutes, (int)seconds);
  return;
}

void printTime2(clock_t start, clock_t stop, char *string)
{
  double duration = (double)(stop - start) / CLOCKS_PER_SEC;
  int hours       = (int)(duration / 3600);
  double remain   = fmod(duration, 3600);
  int minutes     = (int)(remain / 60);
  double seconds  = fmod(remain, 60);
  
  if (hours == 0) {
    if (minutes == 0) printf("%s time = %.2f secs\n", string, duration);
    else              printf("%s time = %d m %d s\n", string, minutes, (int)seconds);
  }
  else                printf("%s time = %d h %d m %d s\n", string, hours, minutes, (int)seconds);
  return;
}

//----------------------------------------------------------------------

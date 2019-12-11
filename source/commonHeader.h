

  //------------------------------------------------------//
  //--  commonHeader.h					--//
  //--  Version 2018.09.26				--//
  //--  						--//
  //--  Copyright (C) 2018 - Chieh-An Lin		--//
  //--  GNU GPLv3 - https://www.gnu.org/licenses/	--//
  //------------------------------------------------------//


#ifndef __MFP_COMMON_HEADER__
#define __MFP_COMMON_HEADER__

#ifndef __MFP_USE_HEALPIX_CXX__
  #define __MFP_USE_HEALPIX_CXX__
#endif
// #ifndef __MFP_USE_MPI__
//   #define __MFP_USE_MPI__
// #endif

#ifndef _GNU_SOURCE
  #define _GNU_SOURCE 1
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <fftw3.h>

//-- Constants
//-- The precision of 64 bits double is 52 bits, and 2^-52 is roughly 2e-16, so we only keep 17 digits in base 10.

//-- Mathematical constants
#define PI                 3.1415926535897932
#define HALF_PI            1.5707963267948966
#define TWO_PI             6.2831853071795865
#define FOUR_PI            12.566370614359172
#define FOUR_PI_OVER_THREE 4.1887902047863910 
#define PI_SQ              9.8696044010893586
#define PI_INV             0.31830988618379067
#define ONE_OVER_TWO_PI    0.15915494309189534

#define SQRT_2             1.4142135623730950
#define SQRT_2_INV         0.70710678118654752

#define EXP_1              2.7182818284590452
#define LN_10              2.3025850929940457
#define LOG10_E            0.43429448190325183
#define LOG2_E             1.4426950408889634

//-- Unit conversion
#define DEGREE_TO_RADIAN       0.017453292519943296
#define ARCMIN_TO_RADIAN       2.9088820866572160e-04
#define ARCSEC_TO_RADIAN       4.8481368110953599e-06
#define RADIAN_TO_DEGREE       57.295779513082321
#define RADIAN_TO_ARCMIN       3437.7467707849393
#define RADIAN_TO_ARCSEC       206264.80624709636
#define DEGREE_SQ_TO_RADIAN_SQ 3.0461741978670860e-04
#define ARCMIN_SQ_TO_RADIAN_SQ 8.4615949940752389e-08
#define ARCSEC_SQ_TO_RADIAN_SQ 2.3504430539097886e-11
#define RADIAN_SQ_TO_DEGREE_SQ 3.2828063500117438e+03
#define RADIAN_SQ_TO_ARCMIN_SQ 1.1818102860042278e+07
#define RADIAN_SQ_TO_ARCSEC_SQ 4.2545170296152200e+10
#define MEGA_PARSEC_TO_METER   3.0856775814671916e+22

//-- Physical constants
#define LIGHT_SPEED          299792458              //-- [m/s]
#define CRITICAL_DENSITY     2.7753619786618317e+11 //-- [M_sol h^2 / Mpc^3], mass density, 3 H^2 / 8 pi G
#define HUBBLE_DISTANCE      2997.92458             //-- [Mpc/h]; LIGHT SPEED is in [m/s], so we should write H_0 in 100000 h [m/s/Mpc]
#define FOUR_PI_G_OVER_C2    6.0135402045290702e-19 //-- [Mpc / M_sol]
#define FULL_SKY             41252.961249419271     //-- [deg^2]

//-- Computational constants
#define STRING_LENGTH_MAX 1024                   //-- Maximal length of a string
#define EPS_MIN           4.4408920985006262e-16 //-- 2^-51
#define EPS_NUM           1e-12 //-- Numerical tolerance


typedef struct {
  int length;
  int *array;
} int_arr;

typedef struct {
  int N1, N2, length;
  int *matrix;
} int_mat;

typedef struct {
  int N1, N2, N3, length;
  int *tensor;
} int_ten3;

typedef struct {
  int N1, N2, length;
  short *matrix;
} short_mat;

typedef struct {
  int length;
  double *array;
} double_arr;

typedef struct {
  int N1, N2, length;
  double *matrix;
} double_mat;

typedef struct {
  int N1, N2, N3, length;
  double *tensor;
} double_ten3;

typedef struct {
  int N1, N2, N3, N4, N5, length;
  float *tensor;
} float_ten5;

typedef struct {
  int length;    //-- Number of points
  double dx;     //-- Interval width
  double *x;     //-- Coordinates
  double *value; //-- Values
} interpolator_t;

typedef struct {
  int length;    //-- Number of points
  int ind;       //-- Index of the sampler in an array
  double dx;     //-- Interval width
  double *x;     //-- Coordinates
  double *pdf;   //-- Normalized pdf
  double *cdf;   //-- Normalized cdf
  double totPdf; //-- Integration over pdf, before or after normalization (see set_sampler_t)
  double x_mean; //-- Integration over x*p(x), before or after normalization (see set_sampler_t)
} sampler_t;

typedef struct {
  int length;
  sampler_t **array;
} sampler_arr;

typedef struct {
  int N1, N2, length;
  sampler_t **matrix;
} sampler_mat;

typedef struct {
  int N;             //-- Resolution, should be a square
  int length;        //-- Number of pixels
  int ind;           //-- Index of the transformer in an array
  double normFactor; //-- Normalization factor
  fftw_complex *before, *kernel, *after; //-- fftw_complex elements, transformations are in-place
  fftw_plan before_f, kernel_f, after_b; //-- fftw_plan elements, _f = forward, _b = backward
  fftw_plan before_b;                    //-- Only used for KS inversion
} FFT_t;

typedef struct {
  int length;
  FFT_t **array;
} FFT_arr;

typedef struct {
  int length;      //-- Number of bins
  double dx;       //-- Bin width
  int n_tot;       //-- Total counts
  double *x_lower; //-- Lower limit of each bin
  double x_max;    //-- Maximal limit
  int *n;          //-- Histogram
} hist_t;

typedef struct {
  int length;       //-- Number of samples
  double mean;      //-- Mean
  double variance;  //-- Standard deviation
  double two_h_sq;  //-- 2 * bandwidth^2
  double amplitude; //-- Common amplitude factor
  double *sample;   //-- Sample points
} KDE_t;


//-- Math functions
#define SQ(a)             (pow((a), 2))
#define CB(a)             (pow((a), 3))
#define SUM_SQ_2(a, b)    (pow((a), 2) + pow((b), 2))
#define SUM_SQ_3(a, b, c) (pow((a), 2) + pow((b), 2) + pow((c), 2))
#define POS_MOD(N, i)     ((i) >= 0 ? (i)%(N) : ((i)%(N) + N)%(N))                                           //-- Return value in range [0, N-1]
#define CEN_MOD(N, i)     ((i)%(N) > (N)/2 ? (i)%(N) - (N) : (i)%(N) <= (N)/2 - (N) ? (i)%(N) + N : (i)%(N)) //-- Return value in range [N/2 - N + 1, N/2], [-1, 0, 1] if N = 3 and [-1, 0, 1, 2] if N = 4

//-- Distance functions
#define DIST_2D_SQ(a, b)  (pow(a[0]-b[0], 2) + pow(a[1]-b[1], 2))
#define DIST_3D_SQ(a, b)  (pow(a[0]-b[0], 2) + pow(a[1]-b[1], 2) + pow(a[2]-b[2], 2))
#define NORM_2D_SQ(a)     (pow(a[0], 2) + pow(a[1], 2))
#define NORM_3D_SQ(a)     (pow(a[0], 2) + pow(a[1], 2) + pow(a[2], 2))
#define SPHE_DIST(a, b)   (2.0 * asin(sqrt(pow(sin(0.5*(a[1]-b[1])), 2) + cos(a[1]) * cos(b[1]) * pow(sin(0.5*(a[0]-b[0])), 2))))

//-- Bitwise operations
#define SET_BIT(c, r)     ((c) |  (1 << (r)))
#define CLEAN_BIT(c, r)   ((c) & ~(1 << (r)))
#define TOGGLE_BIT(c, r)  ((c) ^  (1 << (r)))
#define CHECK_BIT(c, r)   (((c) >> (r)) & 1)
#define SET_BIT_64(c, r)     ((c) |  (1UL << (r)))
#define CLEAN_BIT_64(c, r)   ((c) & ~(1UL << (r)))
#define TOGGLE_BIT_64(c, r)  ((c) ^  (1UL << (r)))
#define CHECK_BIT_64(c, r)   (((c) >> (r)) & 1UL)

//-- Functions related to array
void reset_double(double *lfArr, int length);
void rescale_double(double *lfArr, int length, double factor);
void reset_fftw_complex(fftw_complex *table, int length);
void rescaleReal_fftw_complex(fftw_complex *table, int length, double factor);
void rescale_fftw_complex(fftw_complex *table, int length, double factor);
void multiplication_fftw_complex(fftw_complex *table1, fftw_complex *table2, fftw_complex *product, int length);
void copy_fftw_complex(fftw_complex *from, fftw_complex *to, int length);
void print_fftw_complex(fftw_complex *table, int N1);

//-- Functions related to int_arr
int_arr *initialize_int_arr(int length);
void free_int_arr(int_arr *intArr);
void print_int_arr(int_arr *intArr);

//-- Functions related to int_mat
int_mat *initialize_int_mat(int N1, int N2);
void free_int_mat(int_mat *intMat);
void print_int_mat(int_mat *intMat);

//-- Functions related to int_ten3
int_ten3 *initialize_int_ten3(int N1, int N2, int N3);
void free_int_ten3(int_ten3 *intTen);
  
//-- Functions related to short_mat
short_mat *initialize_short_mat(int N1, int N2);
void free_short_mat(short_mat *shrtMat);
void print_short_mat(short_mat *shrtMat);

//-- Functions related to double_arr
double_arr *initialize_double_arr(int length);
void free_double_arr(double_arr *dblArr);
void print_double_arr(double_arr *dblArr);

//-- Functions related to double_mat
double_mat *initialize_double_mat(int N1, int N2);
void free_double_mat(double_mat *dblMat);
void print_double_mat(double_mat *dblMat);

//-- Functions related to double_ten3
double_ten3 *initialize_double_ten3(int N1, int N2, int N3);
void free_double_ten3(double_ten3 *dblTen);

//-- Functions related to float_ten5
float_ten5 *initialize_float_ten5(int N1, int N2, int N3, int N4, int N5);
void free_float_ten5(float_ten5 *fltTen);

//-- Functions related to interpolator_t
interpolator_t *initialize_interpolator_t(int length);
void free_interpolator_t(interpolator_t *inter);
void print_interpolator_t(interpolator_t *inter);
double execute_interpolator_t(interpolator_t *inter, double x, int border);

//-- Functions related to sampler_t
sampler_t *initialize_sampler_t(int length);
void free_sampler_t(sampler_t *samp);
void print_sampler_t(sampler_t *samp);
void set_sampler_t(sampler_t *samp, int mode, int setTotalToOne);
double execute_sampler_t(sampler_t *samp, double x, int mode);

//-- Functions related to sampler_arr
sampler_arr *initialize_sampler_arr(int N_type, int N_array);
void free_sampler_arr(sampler_arr *sampArr);

//-- Functions related to sampler_mat
sampler_mat *initialize_sampler_mat(int N_type, int N1_array, int N2_array);
void free_sampler_mat(sampler_mat *sampMat);

//-- Functions related to FFT_t
FFT_t *initialize_FFT_t(int N);
void free_FFT_t(FFT_t *transformer);
void reset_FFT_t(FFT_t *transformer);
void execute_FFT_t(FFT_t *transformer);

//-- Functions related to FFT_arr
FFT_arr *initialize_FFT_arr(int N_type, int N_array);
void free_FFT_arr(FFT_arr *transArr);
void reset_FFT_arr(FFT_arr *transArr);
void execute_FFT_arr(FFT_arr *transArr);

//-- Functions related to hist_t
hist_t *initialize_hist_t(int length);
void free_hist_t(hist_t *hist);
void print_hist_t(hist_t *hist);
void set_hist_t(hist_t *hist, double x_min, double x_max);
void reset_hist_t(hist_t *hist);
void push_hist_t(hist_t *hist, double x, int verbose);

//-- Functions related to KDE_t
KDE_t *initialize_KDE_t(int length);
void free_KDE_t(KDE_t *estimator);
void set_KDE_t(KDE_t *estimator, double h);
double execute_KDE_t(KDE_t *estimator, double x);
double integrate_KDE_t(KDE_t *estimator, int N, double x);

//-- Functions related to RNG
u_int32_t renewSeed();
gsl_rng *initializeGenerator(u_int32_t seed);
void printGenerator(gsl_rng *generator);

//-- Functions related to rotation & projection
void rotate(double oldPos[2], double rotAng, double newPos[2]);
void project(double RADEC[2], double center[4], double thetaXY[2]);
void deproject(double thetaXY[2], double center[4], double RADEC[2]);

//-- Functions related to comparison
int compare_float(const void *a, const void *b);
int compare_longLong(const void *a, const void *b);

//-- Functions related to stopwatch
void printTime(clock_t start, clock_t stop);
void printTime2(clock_t start, clock_t stop, char *string);

//----------------------------------------------------------------------

#endif


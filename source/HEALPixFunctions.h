

  //------------------------------------------------------//
  //--  HEALPixFunctions.h                              --//
  //--  Version 2020.07.09                              --//
  //--                                                  --//
  //--  Copyright (C) 2020 - Chieh-An Lin               --//
  //--  GNU GPLv3 - https://www.gnu.org/licenses/       --//
  //------------------------------------------------------//


#include "commonHeader.h"

#ifndef __SALMO_HEALPIX_FUNCTIONS__
#define __SALMO_HEALPIX_FUNCTIONS__

#include "chealpix.h"


//-- Functions related to base 2 calculations
int log2_int(int a);
int log2_long(long a);
int log2_longLong(long long a);
int isPowerOf2(int a);

//-- Functions related to HEALPix conversion
void pixelToPatch(long long nsidePat, long long nsidePix, long long pix, long long *patch, int doNest);
void patchToPixels(long long nsidePat, long long nsidePix, long long patch, long long *pix, int doNest, int doSort);
void nsideToLevels(long long nside, long long *lenArr, long long *cumLenArr);
void patchToCap(long long nside, long long patch, int doNest, int *cap);

//-- Functions related to Cartesian ordering
void ijPixToLocalNest(int resol, int i_pix, int j_pix, int *localNest);
void localNestToIJPix(int resol, int localNest, int *i_pix, int *j_pix);
void CartesianToLocalNest(int resol, int carte, int *localNest);
void localNestToCartesian(int resol, int localNest, int *carte);
void resolToIJPix(int resol, int *i_pix, int *j_pix);
void resolToLocalNest(int resol, int *localNest);
void CartesianToRingInPatch(long long nside, long long patch, int resol, long long *pix);

//-- Functions related to neighbors
void decompose(long long nside, long long patch, int doNest, int *cap, int *level, int *length, int *off, int *j);
void decompose2(long long nside, long long patch, int doNest, long long *lenArr, long long *cumLenArr, long long v1, long long v2, int *cap, int *level, int *length, int *off, int *j);
void recombine(long long nside, int level, int off, int j, long long *patch);
void findNeighborsForCap(long long nside, int level, int off, int j, long long neiInfo[3][8]);
void findNeighborsForBelt(long long nside, int level, int off, int j, long long neiInfo[3][8]);
void findNeighbors(long long nside, long long patch, int doNest, long long neighbors[8]);

//-- Functions related to boundary
void getLimits(long long nside, long long patch, int doNest, double limits[4]);
void adjustRA(double *RA, int N, double half, int sign);
double HEALPixCapLine(double psi, double l);
double HEALPixBeltLine(double psi, double l);
double HEALPixCapArea(double psi_0, double psi, double l, double z_0);
double HEALPixCapInverseArea_plus(double A, double cst1, double cst2, double cst3);
double HEALPixCapInverseArea_minus(double A, double cst1, double cst2, double cst3);
double HEALPixBeltInverseArea(double A, double cst4, double cst5);
void capSampling(gsl_rng *generator, long long nside, int N, int level, int length, int j, double z_0, double *pos1, double *pos2);
void beltSampling(gsl_rng *generator, long long nside, int N, double *pos1, double *pos2);
void patchSampling1(gsl_rng *generator, long long nside, long long patch, int N, int doNest, double *pos1, double *pos2);
void patchSampling2(gsl_rng *generator, long long nside, long long patch, int N, int doNest, double *pos1, double *pos2);
void patchSampling3(gsl_rng *generator, long long nside, long long patch, int N, int doNest, long long *lenArr, long long *cumLenArr, long long v1, long long v2, double *pos1, double *pos2);
void patchSampling4(gsl_rng *generator, long long nside, int cap, int level, int length, int off, int j, double z0, double ctrPhi, double *pos1, double *pos2);
void comparePatchSampling();

#endif


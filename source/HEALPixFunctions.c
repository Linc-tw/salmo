

  //------------------------------------------------------//
  //--  HEALPixFunctions.c				--//
  //--  Version 2019.01.25				--//
  //--  						--//
  //--  Copyright (C) 2019 - Chieh-An Lin		--//
  //--  GNU GPLv3 - https://www.gnu.org/licenses/	--//
  //------------------------------------------------------//


#include "HEALPixFunctions.h"


//----------------------------------------------------------------------
//-- Functions related to base 2 calculations

int log2_int(int a)
{
  if (a < 1)  return -1;
  if (a == 1) return 0;
  return 1 + log2_int(a >> 1);
}

int log2_long(long a)
{
  if (a < 1)  return -1;
  if (a == 1) return 0;
  return 1 + log2_long(a >> 1);
}

int log2_longLong(long long a)
{
  if (a < 1)  return -1;
  if (a == 1) return 0;
  return 1 + log2_longLong(a >> 1);
}

int isPowerOf2(int a)
{
  if (a < 1)  return 0;
  if (a == 1) return 1;
  return (1 - a % 2) && isPowerOf2(a / 2);
}

//----------------------------------------------------------------------
//-- Functions related to HEALPix conversion

void pixelToPatch(long long nsidePat, long long nsidePix, long long pix, long long *patch, int doNest)
{
  long long pixNest = pix;
  if (doNest == 0) ring2nest64(nsidePix, pix, &pixNest); //-- Require nest
  
  long long resol  = nsidePix / nsidePat;
  long long length = resol * resol;
  
  patch[0] = pixNest / length;
  if (doNest == 0) nest2ring64(nsidePat, patch[0], patch);
  return;
}

void patchToPixels(long long nsidePat, long long nsidePix, long long patch, long long *pix, int doNest, int doSort)
{
  //-- *pix should be initialized to have length = (nsidePix / nsidePat)^2
  
  long long patNest = patch;
  if (doNest == 0) ring2nest64(nsidePat, patch, &patNest); //-- Require nest
  
  long long resol  = nsidePix / nsidePat;
  long long length = resol * resol;
  
  long long i;
  for (i=0; i<length; i++) pix[i] = patNest * length + i;
  
  if (doNest == 0) for (i=0; i<length; i++) nest2ring64(nsidePix, pix[i], &pix[i]); //-- Require nest
  if (doSort == 1) qsort(pix, length, sizeof(long long), compare_longLong);
  return;
}

void nsideToLevels(long long nside, long long *lenArr, long long *cumLenArr)
{
  long long beltLength     = 4 * nside;
  long long nbLevelsInBelt = 2 * nside - 1;
  
  cumLenArr[0] = 0;
  long long i, j = 1;
  
  for (i=1; i<=nside; i++, j++) {
    lenArr[j-1]  = 4 * i;
    cumLenArr[j] = cumLenArr[j-1] + lenArr[j-1];
  }
  
  for (i=0; i<nbLevelsInBelt; i++, j++) {
    lenArr[j-1]  = beltLength;
    cumLenArr[j] = cumLenArr[j-1] + lenArr[j-1];
  }
  
  for (i=nside; i>0; i--, j++) {
    lenArr[j-1]  = 4 * i;
    cumLenArr[j] = cumLenArr[j-1] + lenArr[j-1];
  }
  return;
}

void patchToCap(long long nside, long long patch, int doNest, int *cap)
{
  long long patRing = patch;
  if (doNest == 1) nest2ring64(nside, patch, &patRing); //-- Require ring
  
  long long v1 = 2 * nside * (nside + 1);
  long long v2 = 2 * nside * (5 * nside - 1);
  cap[0] = -1 + (int)(patRing < v2) + (int)(patRing < v1);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to Cartesian ordering

void ijPixToLocalNest(int resol, int i_pix, int j_pix, int *localNest)
{
  localNest[0] = 0;
  
  int halfNbBits = log2_int(resol);
  int i, k;
  
  for (i=0; i<halfNbBits; i++) {
    k = 1 << i;
    localNest[0] += (i_pix & k) * k;     //-- Bitwise and
    localNest[0] += (j_pix & k) * k * 2; //-- Bitwise and
  }
  return;
}

void localNestToIJPix(int resol, int localNest, int *i_pix, int *j_pix)
{
  i_pix[0] = 0;
  j_pix[0] = 0;
  
  int halfNbBits = log2_int(resol);
  int i, k;
  
  for (i=0; i<halfNbBits; i++) {
    k = 1 << i;
    i_pix[0] += (localNest & (k * k)) / k;
    j_pix[0] += (localNest & (2 * k * k)) / k;
  }
  j_pix[0] /= 2;
  return;
}

void CartesianToLocalNest(int resol, int carte, int *localNest)
{
  int i_pix = carte % resol;
  int j_pix = carte / resol;
  ijPixToLocalNest(resol, i_pix, j_pix, localNest);
  return;
}

void localNestToCartesian(int resol, int localNest, int *carte)
{
  int i_pix, j_pix;
  localNestToIJPix(resol, localNest, &i_pix, &j_pix);
  carte[0] = i_pix + j_pix * resol;
  return;
}

void resolToIJPix(int resol, int *i_pix, int *j_pix)
{
  int length = resol * resol;
  int i;
  for (i=0; i<length; i++) localNestToIJPix(resol, i, &i_pix[i], &j_pix[i]);
  return;
}

void resolToLocalNest(int resol, int *localNest)
{
  int length = resol * resol;
  int i;
  for (i=0; i<length; i++) CartesianToLocalNest(resol, i, &localNest[i]);
  return;
}

void CartesianToRingInPatch(long long nside, long long patch, int resol, long long *pix)
{
  int length = resol * resol;
  long long first;
  ring2nest64(nside, patch, &first);
  first *= length;
  
  long long nsidePix = nside * resol;
  int i, localNest;
  
  for (i=0; i<length; i++) {
    CartesianToLocalNest(resol, i, &localNest);
    pix[i] = first + localNest;
    nest2ring64(nsidePix, pix[i], &pix[i]);
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to neighbors

void decompose(long long nside, long long patch, int doNest, int *cap, int *level, int *length, int *off, int *j)
{
  long long *lenArr    = (long long*)malloc((4*nside-1) * sizeof(long long));
  long long *cumLenArr = (long long*)malloc((4*nside) * sizeof(long long));
  
  nsideToLevels(nside, lenArr, cumLenArr);
  patchToCap(nside, patch, doNest, cap);
  
  long long patRing = patch;
  if (doNest == 1) nest2ring64(nside, patch, &patRing); //-- Require ring
  
  //-- Determine level
  level[0] = 0;
  while (patRing >= cumLenArr[level[0]+1]) level[0]++;
  
  //-- Determine length, off, j
  long long rest = patRing - cumLenArr[level[0]];
  length[0] = lenArr[level[0]] / 4;
  off[0]    = (rest / length[0] + 2) % 4 - 2;
  j[0]      = rest % length[0];
  
  free(lenArr);
  free(cumLenArr);
  return;
}

void decompose2(long long nside, long long patch, int doNest, long long *lenArr, long long *cumLenArr, long long v1, long long v2, int *cap, int *level, int *length, int *off, int *j)
{
  long long patRing = patch;
  if (doNest == 1) nest2ring64(nside, patch, &patRing); //-- Require ring
  
  //-- Determine cap
  cap[0] = -1 + (int)(patRing < v2) + (int)(patRing < v1);
  
  //-- Determine level
  level[0] = 0;
  while (patRing >= cumLenArr[level[0]+1]) level[0]++;
  
  //-- Determine length, off, j
  long long rest = patRing - cumLenArr[level[0]];
  length[0] = lenArr[level[0]] / 4;
  off[0]    = (rest / length[0] + 2) % 4 - 2;
  j[0]      = rest % length[0];
  return;
}

void recombine(long long nside, int level, int off, int j, long long *patch)
{
  if (level == -1) {
    patch[0] = -1;
    return;
  }
  
  long long *lenArr    = (long long*)malloc((4*nside-1) * sizeof(long long));
  long long *cumLenArr = (long long*)malloc((4*nside) * sizeof(long long));
  
  nsideToLevels(nside, lenArr, cumLenArr);
  off      = (off + 4) % 4;
  patch[0] = cumLenArr[level] + off * lenArr[level] / 4 + j;
  
  free(lenArr);
  free(cumLenArr);
  return;
}

void findNeighborsForCap(long long nside, int level, int off, int j, long long neiInfo[3][8])
{
  //-- neiInfor contains level, length, off, j of 8 neighbors
  //-- in the order of SW, W, NW, N, NE, E, SE, S
  
      neiInfo[0][0] = level+1; neiInfo[1][0] = off;   neiInfo[2][0] = j;
      neiInfo[0][1] = level;   neiInfo[1][1] = off;   neiInfo[2][1] = j-1;
      neiInfo[0][2] = level-1; neiInfo[1][2] = off;   neiInfo[2][2] = j-1;
      neiInfo[0][3] = level-2; neiInfo[1][3] = off;   neiInfo[2][3] = j-1;
      neiInfo[0][4] = level-1; neiInfo[1][4] = off;   neiInfo[2][4] = j;
      neiInfo[0][5] = level;   neiInfo[1][5] = off;   neiInfo[2][5] = j+1;
      neiInfo[0][6] = level+1; neiInfo[1][6] = off;   neiInfo[2][6] = j+1;
      neiInfo[0][7] = level+2; neiInfo[1][7] = off;   neiInfo[2][7] = j+1;
  
  if (j == 0) {
      neiInfo[0][1] = level+1; neiInfo[1][1] = off-1; neiInfo[2][1] = level+1;
      neiInfo[0][2] = level;   neiInfo[1][2] = off-1; neiInfo[2][2] = level;
      neiInfo[0][3] = level-1; neiInfo[1][3] = off-1; neiInfo[2][3] = level-1;
  }
  
  if (j == level) {
      neiInfo[0][3] = level-1; neiInfo[1][3] = off+1; neiInfo[2][3] = 0;
      neiInfo[0][4] = level;   neiInfo[1][4] = off+1; neiInfo[2][4] = 0;
      neiInfo[0][5] = level+1; neiInfo[1][5] = off+1; neiInfo[2][5] = 0;
  }
  
  if (level == nside - 1) {
      neiInfo[0][7] = level+2; neiInfo[1][7] = off;   neiInfo[2][7] = j;
    if (j == 0) {
      neiInfo[0][1] = -1;
    }
    if (j == level) {
      neiInfo[0][5] = -1;
      neiInfo[0][6] = nside;   neiInfo[1][6] = off+1; neiInfo[2][6] = 0;
    }
  }
  
  if (level == 0) {
      neiInfo[0][0] = 1;       neiInfo[1][0] = off;   neiInfo[2][0] = 0;
      neiInfo[0][1] = 1;       neiInfo[1][1] = off-1; neiInfo[2][1] = 1;
      neiInfo[0][2] = 0;       neiInfo[1][2] = off-1; neiInfo[2][2] = 0;
      neiInfo[0][3] = 0;       neiInfo[1][3] = off-2; neiInfo[2][3] = 0;
      neiInfo[0][4] = 0;       neiInfo[1][4] = off+1; neiInfo[2][4] = 0;
      neiInfo[0][5] = 1;       neiInfo[1][5] = off+1; neiInfo[2][5] = 0;
      neiInfo[0][6] = 1;       neiInfo[1][6] = off;   neiInfo[2][6] = 1;
      neiInfo[0][7] = 2;       neiInfo[1][7] = off;   neiInfo[2][7] = 1;
  }
  
  int i;
  for (i=0; i<8; i++) neiInfo[1][i] = (neiInfo[1][i] + 4) % 4;
  return; 
}

void findNeighborsForBelt(long long nside, int level, int off, int j, long long neiInfo[3][8])
{
  //-- neiInfor contains level, off, j of 8 neighbors
  //-- in the order of SW, W, NW, N, NE, E, SE, S
  
  if (level % 2 == 0) {
      neiInfo[0][0] = level+1; neiInfo[1][0] = off;   neiInfo[2][0] = j-1;
      neiInfo[0][1] = level;   neiInfo[1][1] = off;   neiInfo[2][1] = j-1;
      neiInfo[0][2] = level-1; neiInfo[1][2] = off;   neiInfo[2][2] = j-1;
      neiInfo[0][3] = level-2; neiInfo[1][3] = off;   neiInfo[2][3] = j;
      neiInfo[0][4] = level-1; neiInfo[1][4] = off;   neiInfo[2][4] = j;
      neiInfo[0][5] = level;   neiInfo[1][5] = off;   neiInfo[2][5] = j+1;
      neiInfo[0][6] = level+1; neiInfo[1][6] = off;   neiInfo[2][6] = j;
      neiInfo[0][7] = level+2; neiInfo[1][7] = off;   neiInfo[2][7] = j;
  
    if (j == 0) {
      neiInfo[0][0] = level+1; neiInfo[1][0] = off-1; neiInfo[2][0] = nside-1;
      neiInfo[0][1] = level;   neiInfo[1][1] = off-1; neiInfo[2][1] = nside-1;
      neiInfo[0][2] = level-1; neiInfo[1][2] = off-1; neiInfo[2][2] = nside-1;
    }
    
    if (j == nside - 1) {
      neiInfo[0][5] = level;   neiInfo[1][5] = off+1; neiInfo[2][5] = 0;
    }
    
    if (level == nside) {
      neiInfo[0][3] = level-2; neiInfo[1][3] = off;   neiInfo[2][3] = j-1;
      if (j == 0) {
      neiInfo[0][3] = -1;
      }
    }
    if (level == 3*nside-2) {
      neiInfo[0][7] = level+2; neiInfo[1][7] = off;   neiInfo[2][7] = j-1;
      if (j == 0) {
      neiInfo[0][7] = -1;
      }
    }
  }
  
  else {
      neiInfo[0][0] = level+1; neiInfo[1][0] = off;   neiInfo[2][0] = j;
      neiInfo[0][1] = level;   neiInfo[1][1] = off;   neiInfo[2][1] = j-1;
      neiInfo[0][2] = level-1; neiInfo[1][2] = off;   neiInfo[2][2] = j;
      neiInfo[0][3] = level-2; neiInfo[1][3] = off;   neiInfo[2][3] = j;
      neiInfo[0][4] = level-1; neiInfo[1][4] = off;   neiInfo[2][4] = j+1;
      neiInfo[0][5] = level;   neiInfo[1][5] = off;   neiInfo[2][5] = j+1;
      neiInfo[0][6] = level+1; neiInfo[1][6] = off;   neiInfo[2][6] = j+1;
      neiInfo[0][7] = level+2; neiInfo[1][7] = off;   neiInfo[2][7] = j;
      
    if (j == 0) {
      neiInfo[0][1] = level;   neiInfo[1][1] = off-1; neiInfo[2][1] = nside-1;
    }
      
    if (j == nside - 1) {
      neiInfo[0][4] = level-1; neiInfo[1][4] = off+1; neiInfo[2][4] = 0;
      neiInfo[0][5] = level;   neiInfo[1][5] = off+1; neiInfo[2][5] = 0;
      neiInfo[0][6] = level+1; neiInfo[1][6] = off+1; neiInfo[2][6] = 0;
    }
  }
  
  int i;
  for (i=0; i<8; i++) neiInfo[1][i] = (neiInfo[1][i] + 4) % 4;
  return; 
}

void findNeighbors(long long nside, long long patch, int doNest, long long neighbors[8])
{
  long long neiInfo[3][8];
  int cap, level, length, off, j, i;
  
  long long patRing = patch;
  if (doNest == 1) nest2ring64(nside, patch, &patRing); //-- Require ring
  decompose(nside, patRing, doNest, &cap, &level, &length, &off, &j);
  
  if      (cap > 0)  findNeighborsForCap(nside, level, off, j, neiInfo);
  else if (cap == 0) findNeighborsForBelt(nside, level, off, j, neiInfo);
  else  {
    level = 4*nside-2 - level;
    findNeighborsForCap(nside, level, off, j, neiInfo);
    for (i=0; i<8; i++) neiInfo[0][i] = (neiInfo[0][i] == -1) ? -1 : 4*nside-2 - neiInfo[0][i];
    for (i=0; i<3; i++) {
      neighbors[0]  = neiInfo[i][0];
      neiInfo[i][0] = neiInfo[i][2];
      neiInfo[i][2] = neighbors[0];
      neighbors[0]  = neiInfo[i][3];
      neiInfo[i][3] = neiInfo[i][7];
      neiInfo[i][7] = neighbors[0];
      neighbors[0]  = neiInfo[i][4];
      neiInfo[i][4] = neiInfo[i][6];
      neiInfo[i][6] = neighbors[0];
    }
  }
  
  for (i=0; i<8; i++) recombine(nside, neiInfo[0][i], neiInfo[1][i], neiInfo[2][i], &neighbors[i]);
  return; 
}

//----------------------------------------------------------------------
//-- Functions related to boundary

void getLimits(long long nside, long long patch, int doNest, double limits[4])
{
  //-- Return z_min, z_max, phi_min, phi_max
  
  int cap, level, length, off, j;
  
  long long patRing = patch;
  if (doNest == 1) nest2ring64(nside, patch, &patRing); //-- Require ring
  decompose(nside, patRing, doNest, &cap, &level, &length, &off, &j);
  
  if (cap == 1) {
    if (level == nside-1) limits[0] = 2.0/3.0 * (2.0 - (double)(level+2) / (double)nside);
    else                  limits[0] = 1.0 - pow((double)(level+2) / (double)nside, 2) / 3.0;
    limits[1] = 1.0 - pow((double)level / (double)nside, 2) / 3.0;
    limits[2] = HALF_PI * ((double)j / (double)length + off);
    limits[3] = HALF_PI * ((double)(j+1) / (double)length + off);
  }
  
  else if (cap == 0) {
    limits[0] = 2.0/3.0 * (2.0 - (double)(level+2) / (double)nside);
    limits[1] = 2.0/3.0 * (2.0 - (double)level / (double)nside);
    limits[2] = HALF_PI * (((double)j - 0.5 * (double)((level+1)%2)) / (double)nside + off);
    limits[3] = HALF_PI * (((double)j - 0.5 * (double)((level+1)%2) + 1) / (double)nside + off);
  }
  
  else {
    limits[0] = -(1.0 - pow((double)(4*nside-2 - level) / (double)nside, 2) / 3.0);
    if (level == 3*nside-1) limits[1] = 2.0/3.0 * (2.0 - (double)level / (double)nside);
    else                    limits[1] = -(1.0 - pow((double)(4*nside - level) / (double)nside, 2) / 3.0);
    limits[2] = HALF_PI * ((double)j / (double)length + off);
    limits[3] = HALF_PI * ((double)(j+1) / (double)length + off);
  }
  return;
}

void adjustRA(double *RA, int N, double half, int sign)
{
  int i;
  if (sign > 0) {
    for (i=0; i<N; i++) RA[i] = fmod(RA[i], 2.0*half);
  }
  else if (sign < 0) {
    for (i=0; i<N; i++) RA[i] = fmod(RA[i], 2.0*half) - 2.0*half;
  }
  else {
    for (i=0; i<N; i++) RA[i] = fmod(RA[i] + half, 2.0*half) - half;
  }
  return;
}

double HEALPixCapLine(double psi, double l)
{
  return 1.0 - pow(l / psi, 2) / 3.0;
}

double HEALPixBeltLine(double psi, double l)
{
  return 4.0/3.0 * (0.5 - l + psi);
}

double HEALPixCapArea(double psi_0, double psi, double l, double z_0)
{
  return (1.0-z_0) * (psi-psi_0) + (l*l)/3.0 * (1.0/psi-1.0/psi_0);
}

double HEALPixCapInverseArea_plus(double A, double cst1, double cst2, double cst3)
{
  return (A + cst1 + sqrt(pow(A + cst1, 2) - cst2)) * cst3;
}

double HEALPixCapInverseArea_minus(double A, double cst1, double cst2, double cst3)
{
  return (A + cst1 - sqrt(pow(A + cst1, 2) - cst2)) * cst3;
}

double HEALPixBeltInverseArea(double A, double cst4, double cst5)
{
  return 0.75 * (-cst4 - sqrt(cst5 + 8.0/3.0*A));
}

void capSampling(gsl_rng *generator, long long nside, int N, int level, int length, int j, double z_0, double *pos1, double *pos2)
{
  //-- Determine corners
  double dummy        = -1.0;
  double j2           = (double)j;
  double length2      = (double)length;
  double cornerPsi[4] = {dummy, j2/length2, (j2+1.0)/(length2+1.0), (j2+1.0)/length2};
  if (level > 0)         cornerPsi[0] = j2 / (length2-1.0);
  if (level == nside -1) cornerPsi[2] = 0.5 * (cornerPsi[1] + cornerPsi[3]);
  
  //-- Determine psi
  double psi0Arr[4]   = {cornerPsi[1], 1.0-cornerPsi[3], 1.0-cornerPsi[2], cornerPsi[2]};
  if (j == 0)     psi0Arr[0] = dummy;
  if (j == level) psi0Arr[1] = dummy;
  
  //-- Determine l
  double dl   = 1.0 / (double)nside;
  double lArr[4];
  lArr[0] = j2 * dl;
  lArr[1] = (length2-1.0-j2) * dl;
  lArr[2] = lArr[1] + dl;
  lArr[3] = lArr[0] + dl;
  
  double halfDenom = 3.0 * pow((double)nside, 2);
  double A_0       = 1.0 / (2.0 * halfDenom);
  double AArr[4]   = {A_0, A_0, A_0, A_0};
  
  //-- Determine areas
  if (level == 0 && nside > 1) {
    AArr[2] = -HEALPixCapArea(1.0-cornerPsi[2], 1.0-cornerPsi[1], lArr[2], z_0);
    AArr[3] = -HEALPixCapArea(    cornerPsi[2],     cornerPsi[3], lArr[3], z_0);
  }
  else {
    if (level < nside - 1) {
      AArr[2] = -HEALPixCapArea(1.0-cornerPsi[2], 1.0-cornerPsi[1], lArr[2], z_0);
      AArr[3] = -HEALPixCapArea(    cornerPsi[2],     cornerPsi[3], lArr[3], z_0);
    }
    
    if (j == 0) {
      AArr[0] = 0.0;
      AArr[1] = HEALPixCapArea(1.0-cornerPsi[3], 1.0-cornerPsi[0], lArr[1], z_0);
    }
    else if (j == level) {
      AArr[0] = HEALPixCapArea(    cornerPsi[1],     cornerPsi[0], lArr[0], z_0);
      AArr[1] = 0.0;
    }
    else {
      AArr[0] = HEALPixCapArea(    cornerPsi[1],     cornerPsi[0], lArr[0], z_0);
      AArr[1] = HEALPixCapArea(1.0-cornerPsi[3], 1.0-cornerPsi[0], lArr[1], z_0);
    }
  }
  
  double cst3 = 0.5 / (1.0 - z_0);
  double cst1[4], cst2[4], cst4[4], cst5[4];
  double ACum[5];
  int i;
  
  //-- Original sampling functions:
  //--   cSampP = lambda A, l, psi_0: (A +(1.0-z_0)*psi_0+l**2/(3.0*psi_0) + np.sqrt((A +(1.0-z_0)*psi_0+l**2/(3.0*psi_0))**2 - 4.0*(1-z_0)*l**2/3.0)) * 0.5 / (1.0 - z_0)
  //--   cSampM = lambda A, l, psi_0: (A +(1.0-z_0)*psi_0+l**2/(3.0*psi_0) - np.sqrt((A +(1.0-z_0)*psi_0+l**2/(3.0*psi_0))**2 - 4.0*(1-z_0)*l**2/3.0)) * 0.5 / (1.0 - z_0)
  //-- Belows are simplification.
  ACum[0] = 0.0;
  for (i=0; i<4; i++) {
    ACum[i+1] = ACum[i] + AArr[i];
    cst1[i] = (1.0 - z_0) * psi0Arr[i] + pow(lArr[i], 2) / (3.0 * psi0Arr[i]);
    cst2[i] = 4.0/3.0 * (1 - z_0) * pow(lArr[i], 2);
    cst4[i] = 4.0/3.0 * (0.5 - lArr[i]) - z_0;
    cst5[i] = pow(cst4[i] + 4.0/3.0*psi0Arr[i], 2);
  }
  
  int count = 0;
  double p, q, rest, psi, z;
  
  while (count < N) {
    //-- Sample and determine which zone
    p = gsl_ran_flat(generator, 0.0, ACum[4]);
    q = gsl_ran_flat(generator, 0.0, 1.0);
    for (i=0; i<4; i++) {
      if (ACum[i+1] > p) break;
    }
    rest = p - ACum[i];
    
    //-- Transform into psi & z
    if (i < 2) {
      if (level == 0) {
	psi = rest * halfDenom;
	z   = q / halfDenom;
      }
      else {
	psi = HEALPixCapInverseArea_plus(rest, cst1[i], cst2[i], cst3);
	z   = (HEALPixCapLine(psi, lArr[i]) - z_0) * q;
      }
    }
    else {
      if (level == nside - 1) {
	psi = HEALPixBeltInverseArea(-rest, cst4[i], cst5[i]);
	z   = (HEALPixBeltLine(psi, lArr[i]) - z_0) * q;
      }
      else {
	psi = HEALPixCapInverseArea_minus(-rest, cst1[i], cst2[i], cst3);
	z   = (HEALPixCapLine(psi, lArr[i]) - z_0) * q;
      }
    }
    
    //-- Flip to the original referential
    pos1[count] = (i == 1 || i == 2) ? 1.0 - psi: psi;
    pos2[count] = z;
    count++;
  }
  return;
}

void beltSampling(gsl_rng *generator, long long nside, int N, double *pos1, double *pos2)
{
  int count = 0;
  double psi_half = 0.5;
  double z_half   = 2.0 / 3.0;
  double factor   = 1.0 / (double)nside;
  double psi, z;
  
  while (count < N) {
    psi = gsl_ran_flat(generator, -psi_half, psi_half);
    z   = gsl_ran_flat(generator, 0.0, z_half);
    
    if (4.0 * fabs(psi) > 3.0 * z) {
      if (psi > 0.0) psi -= psi_half;
      else           psi += psi_half;
    }
    else z -= z_half;
    
    pos1[count] = psi * factor;
    pos2[count] = z * factor;
    count++;
  }
  return;
}

void patchSampling1(gsl_rng *generator, long long nside, long long patch, int N, int doNest, double *pos1, double *pos2)
{
  int cap, level, length, off, j;
  double ctrTheta, ctrPhi, z0;
  
  //-- Decompose & get z0
  decompose(nside, patch, doNest, &cap, &level, &length, &off, &j);
  if (doNest == 1) pix2ang_nest(nside, patch, &ctrTheta, &ctrPhi);
  else             pix2ang_ring(nside, patch, &ctrTheta, &ctrPhi);
  z0 = cos(ctrTheta);
  
  int i;
  
  //-- Sample
  if (cap == 1) {
    capSampling(generator, nside, N, level, length, j, z0, pos1, pos2);
    
    for (i=0; i<N; i++) {
      pos1[i]  = (pos1[i] + off) * HALF_PI;
      pos2[i] += z0;
    }
  }
  else if (cap == 0) {
    beltSampling(generator, nside, N, pos1, pos2);
    adjustRA(&ctrPhi, 1, PI, 0); //-- N = 1, half = PI, sign = 0
    
    for (i=0; i<N; i++) {
      pos1[i]  = pos1[i] * HALF_PI + ctrPhi;
      pos2[i] += z0;
    }
  }
  else {
    capSampling(generator, nside, N, 4*nside-2-level, length, length-1-j, -z0, pos1, pos2);
    
    for (i=0; i<N; i++) {
      pos1[i] = (1.0 - pos1[i] + off) * HALF_PI;
      pos2[i] = z0 - pos2[i];
    }
  }
  
  //-- Conversion
  for (i=0; i<N; i++) pos2[i]  = HALF_PI - acos(pos2[i]);
  return;
}

void patchSampling2(gsl_rng *generator, long long nside, long long patch, int N, int doNest, double *pos1, double *pos2)
{
  int count = 0;
  double limits[4];
  double z, theta, phi;
  long long patch2;
  
  getLimits(nside, patch, doNest, limits);
  
  while (count < N) {
    z     = gsl_ran_flat(generator, limits[0], limits[1]);
    phi   = gsl_ran_flat(generator, limits[2], limits[3]);
    theta = acos(z);
    if (doNest == 1) ang2pix_nest64(nside, theta, phi, &patch2);
    else             ang2pix_ring64(nside, theta, phi, &patch2);
    
    if (patch2 == patch) {
      pos1[count] = phi;
      pos2[count] = HALF_PI - theta;
      count++;
    }
  }
  return;
}

void patchSampling3(gsl_rng *generator, long long nside, long long patch, int N, int doNest, long long *lenArr, long long *cumLenArr, long long v1, long long v2, double *pos1, double *pos2)
{
  int cap, level, length, off, j;
  double ctrTheta, ctrPhi, z0;
  
  //-- Decompose & get z0
  decompose2(nside, patch, doNest, lenArr, cumLenArr, v1, v2, &cap, &level, &length, &off, &j);
  if (doNest == 1) pix2ang_nest64(nside, patch, &ctrTheta, &ctrPhi);
  else             pix2ang_ring64(nside, patch, &ctrTheta, &ctrPhi);
  z0 = cos(ctrTheta);
  
  int i;
  
  //-- Sample
  if (cap == 1) {
    capSampling(generator, nside, N, level, length, j, z0, pos1, pos2);
    
    for (i=0; i<N; i++) {
      pos1[i]  = (pos1[i] + off) * HALF_PI;
      pos2[i] += z0;
    }
  }
  else if (cap == 0) {
    beltSampling(generator, nside, N, pos1, pos2);
    adjustRA(&ctrPhi, 1, PI, 0); //-- N = 1, half = PI, sign = 0
    
    for (i=0; i<N; i++) {
      pos1[i]  = pos1[i] * HALF_PI + ctrPhi;
      pos2[i] += z0;
    }
  }
  else {
    capSampling(generator, nside, N, 4*nside-2-level, length, length-1-j, -z0, pos1, pos2);
    
    for (i=0; i<N; i++) {
      pos1[i] = (1.0 - pos1[i] + off) * HALF_PI;
      pos2[i] = z0 - pos2[i];
    }
  }
  
  //-- Conversion
  for (i=0; i<N; i++) pos2[i] = HALF_PI - acos(pos2[i]);
  return;
}

void patchSampling4(gsl_rng *generator, long long nside, int cap, int level, int length, int off, int j, double z0, double ctrPhi, double *pos1, double *pos2)
{
  //-- Sample, N = 1
  if (cap == 1) {
    capSampling(generator, nside, 1, level, length, j, z0, pos1, pos2);
    
    pos1[0]  = (pos1[0] + off) * HALF_PI;
    pos2[0] += z0;
  }
  else if (cap == 0) {
    beltSampling(generator, nside, 1, pos1, pos2);
    adjustRA(&ctrPhi, 1, PI, 0); //-- N = 1, half = PI, sign = 0
    
    pos1[0]  = pos1[0] * HALF_PI + ctrPhi;
    pos2[0] += z0;
  }
  else {
    capSampling(generator, nside, 1, 4*nside-2-level, length, length-1-j, -z0, pos1, pos2);
    
    pos1[0] = (1.0 - pos1[0] + off) * HALF_PI;
    pos2[0] = z0 - pos2[0];
  }
  
  //-- Conversion
  pos2[0] = HALF_PI - acos(pos2[0]);
  return;
}

void comparePatchSampling()
{
  u_int32_t seed     = renewSeed();
  gsl_rng *generator = initializeGenerator(seed);
  
  int N = 10000000;
  long long nside = 2048;
  long long patch = 30331648;
  int doNest = 0;
  
  double *pos1 = (double*)malloc(N * sizeof(double));
  double *pos2 = (double*)malloc(N * sizeof(double));
  long long *lenArr    = (long long*)malloc((4*nside-1) * sizeof(long long));
  long long *cumLenArr = (long long*)malloc((4*nside) * sizeof(long long));
  long long v1 = 2 * nside * (nside + 1);
  long long v2 = 2 * nside * (5 * nside - 1);
  
  int cap, level, length, off, j;
  double ctrTheta, ctrPhi, z0;
  int i;
  
  clock_t stop0 = clock();
  
  patchSampling1(generator, nside, patch, N, doNest, pos1, pos2);
  clock_t stop1 = clock();
  
  patchSampling2(generator, nside, patch, N, doNest, pos1, pos2);
  clock_t stop2 = clock();
  
  nsideToLevels(nside, lenArr, cumLenArr);
  patchSampling3(generator, nside, patch, N, doNest, lenArr, cumLenArr, v1, v2, pos1, pos2);
  clock_t stop3 = clock();
  
  nsideToLevels(nside, lenArr, cumLenArr);
  decompose2(nside, patch, doNest, lenArr, cumLenArr, v1, v2, &cap, &level, &length, &off, &j);
  pix2ang_ring64(nside, patch, &ctrTheta, &ctrPhi);
  z0 = cos(ctrTheta);
  for (i=0; i<N; i++) patchSampling4(generator, nside, cap, level, length, off, j, z0, ctrPhi, &pos1[i], &pos2[i]);
  clock_t stop4 = clock();
  
  printf("HEALPix uniform sampling for %d identical patch:\n", N);
  printTime2(stop0, stop1, "patchSampling1");
  printTime2(stop1, stop2, "patchSampling2");
  printTime2(stop2, stop3, "patchSampling3");
  printTime2(stop3, stop4, "patchSampling4");
  
  N = 200000;
  long long *patchList = (long long*)malloc(N * sizeof(long long));
  for (i=0; i<N; i++) patchList[i] = gsl_rng_uniform_int(generator, 12*nside*nside);
  
  clock_t stop5 = clock();
  
  for (i=0; i<N; i++) patchSampling1(generator, nside, patchList[i], 1, doNest, &pos1[i], &pos2[i]);
  clock_t stop6 = clock();
  
  for (i=0; i<N; i++) patchSampling2(generator, nside, patchList[i], 1, doNest, &pos1[i], &pos2[i]);
  clock_t stop7 = clock();
  
  nsideToLevels(nside, lenArr, cumLenArr);
  for (i=0; i<N; i++) patchSampling3(generator, nside, patchList[i], 1, doNest, lenArr, cumLenArr, v1, v2, &pos1[i], &pos2[i]);
  clock_t stop8 = clock();
  
  nsideToLevels(nside, lenArr, cumLenArr);
  for (i=0; i<N; i++) {
    decompose2(nside, patchList[i], doNest, lenArr, cumLenArr, v1, v2, &cap, &level, &length, &off, &j);
    pix2ang_ring64(nside, patchList[i], &ctrTheta, &ctrPhi);
    z0 = cos(ctrTheta);
    patchSampling4(generator, nside, cap, level, length, off, j, z0, ctrPhi, &pos1[i], &pos2[i]);
  }
  clock_t stop9 = clock();
  
  printf("\n");
  printf("HEALPix uniform sampling for %d different patches:\n", N);
  printTime2(stop5, stop6, "patchSampling1");
  printTime2(stop6, stop7, "patchSampling2");
  printTime2(stop7, stop8, "patchSampling3");
  printTime2(stop8, stop9, "patchSampling4");
  printf("\n");
  printf("Always use patchSampling3\n");
  return;
}

//----------------------------------------------------------------------


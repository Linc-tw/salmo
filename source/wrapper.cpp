

  //------------------------------------------------------//
  //--  wrapper.cpp					--//
  //--  Version 2019.12.11				--//
  //--  						--//
  //--  Copyright (C) 2019 - Chieh-An Lin		--//
  //--  GNU GPLv3 - https://www.gnu.org/licenses/	--//
  //------------------------------------------------------//


#ifdef __MFP_USE_HEALPIX_CXX__
#include <stdio.h>
#include <unistd.h>

#ifdef __cplusplus
#include <Healpix_cxx/alm_healpix_tools.h>
#include <cxxsupport/arr.h>
#include <../src/cxx/cxxsupport/pointing.h>
#include <../src/cxx/Healpix_cxx/healpix_base.h>
#include <../src/cxx/Healpix_cxx/alm.h>
#include <../src/cxx/Healpix_cxx/alm_healpix_tools.h>
#include <../src/cxx/Healpix_cxx/alm_powspec_tools.h>
#include <../src/cxx/Healpix_cxx/healpix_map.h>
#include <../src/cxx/Healpix_cxx/healpix_map_fitsio.h>
#include <../src/cxx/cxxsupport/xcomplex.h>


extern "C" {

#include "FITSFunctions.h"
#include "HEALPixFunctions.h"
#include "wrapper.h"
#endif


//----------------------------------------------------------------------
//-- Functions related to a_lm

void readWeight(long long nside, double *weight)
{
  char name[1024];
  char name2[1024];
  gethostname(name, 1024);
  
  if (!strcmp(name, "cuillin") || strstr(name, "worker") != NULL) sprintf(name2, "/home/calin/03_Codes/Healpix_3.31/data/weight_ring_n%05lld.fits", nside);
  else                                                            sprintf(name2, "/usr/local/lib/python3.5/dist-packages/healpy/data/weight_ring_n%05lld.fits", nside);
  FITS_t *fits = initializeTableReader_FITS_t(name2);
  readTableColumn(fits, 0, weight);
  long long i;
  for (i=0; i<2*nside; i++) weight[i] += 1.0;
  return;
}

void mapToAlm(long long nside, float *map, int l_maxx, double *alm, double *weight)
{
  arr<float> fltArr_cxx(map, 12*nside*nside);
  Healpix_Map<float> map_cxx(fltArr_cxx, RING);
  Alm< xcomplex<float> > alm_cxx(l_maxx, l_maxx);
  
  double *weight2;
  if (weight == NULL) {
    weight2 = (double*)malloc(2*nside * sizeof(double));
    readWeight(nside, weight2);
  }
  else weight2 = weight;
  arr<double> weight_cxx(weight2, 2*nside);
  
  int l, m;
  
  //-- Set to 0
  for (l=0; l<=l_maxx; l++) {
    for (m=0; m<=l; m++) {
//       alm_cxx(l, m).Set(0.0, 0.0);
      alm_cxx(l, m).real(0.0);
      alm_cxx(l, m).imag(0.0);
    }
  }
  
  map2alm(map_cxx, alm_cxx, weight_cxx);
  
  int index;
  
  //-- Retrieve
  for (m=0; m<=l; m++) {
    for (l=m; l<=l_maxx; l++) {
      index = m*(2*l_maxx+1-m)/2+l;
      alm[0+2*index] = alm_cxx(l, m).real();
      alm[1+2*index] = alm_cxx(l, m).imag();
    }
  }
  
  if (weight == NULL) free(weight2);
  return;
}

void almToMap_spin2(int l_maxx, double *alm, long long nside, float *map1, float *map2)
{
  Alm< xcomplex<float> > alm1_cxx(l_maxx, l_maxx);
  Alm< xcomplex<float> > alm2_cxx(l_maxx, l_maxx);
  Healpix_Map<float> map1_cxx((int)nside, RING, SET_NSIDE);
  Healpix_Map<float> map2_cxx((int)nside, RING, SET_NSIDE);
  
  int l, m, index;
  
  //-- Fill
  for (m=0; m<=l; m++) {
    for (l=m; l<=l_maxx; l++) {
      index = m*(2*l_maxx+1-m)/2+l;
//       alm1_cxx(l, m).Set(alm[0+2*index], alm[1+2*index]);
//       alm2_cxx(l, m).Set(0.0, 0.0); //-- No B modes
      alm1_cxx(l, m).real(alm[0+2*index]);
      alm1_cxx(l, m).imag(alm[1+2*index]);
      alm2_cxx(l, m).real(0.0);
      alm2_cxx(l, m).imag(0.0);
    }
  }
  
  alm2map_spin(alm1_cxx, alm2_cxx, map1_cxx, map2_cxx, 2); //-- spin = 2
  
  long long i;
  //-- Retrieve
  for (i=0; i<12*nside*nside; i++) {
    map1[i] = (float)map1_cxx[i];
    map2[i] = (float)map2_cxx[i];
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to pixel interpolation

void getNgbAndWeight(double pos[2], long long nside, long long neighbor[4], double weight[4])
{
  pointing pos_cxx(pos[0], pos[1]);
  fix_arr<int64, 4> neighbor_cxx;
  fix_arr<double, 4> weight_cxx;
  Healpix_Base2 base((int64)nside, RING, SET_NSIDE);
  base.get_interpol(pos_cxx, neighbor_cxx, weight_cxx);
  
  int i;
  for (i=0; i<4; i++) {
    neighbor[i] = neighbor_cxx[i];
    weight[i]   = weight_cxx[i];
  }
  return;
}

//----------------------------------------------------------------------

#ifdef __cplusplus
}
#endif
#endif


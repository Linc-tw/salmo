

  //------------------------------------------------------//
  //--  sampling.c					--//
  //--  Version 2020.04.01				--//
  //--  						--//
  //--  Copyright (C) 2019 - Chieh-An Lin		--//
  //--  GNU GPLv3 - https://www.gnu.org/licenses/	--//
  //------------------------------------------------------//


#include "sampling.h"


//----------------------------------------------------------------------
//-- Functions related to initialization

type_map *initialize_type_map(long long nsidePix)
{
  type_map *tMap = (type_map*)malloc(sizeof(type_map));
  tMap->nsidePix = nsidePix;
  tMap->nbPix    = 12 * nsidePix * nsidePix;
  tMap->map      = (long long*)malloc(tMap->nbPix * sizeof(long long));
  return tMap;
}

void free_type_map(type_map *tMap)
{
  if (tMap) {
    if (tMap->map) {free(tMap->map); tMap->map = NULL;}
    free(tMap); tMap = NULL;
  }
  return;
}

void reset_type_map(type_map *tMap)
{
  long long i;
  for (i=0; i<tMap->nbPix; i++) tMap->map[i] = 0;
  return;
}

gal_node *initialize_gal_node()
{
  gal_node *gNode = (gal_node*)malloc(sizeof(gal_node));
  gNode->g        = (gal_t*)malloc(sizeof(gal_t));
  gNode->next     = NULL;
  return gNode;
}

gal_list *initialize_gal_list()
{
  gal_list *gList = (gal_list*)malloc(sizeof(gal_list));
  gList->length   = 0;
  gList->size     = 0;
  gList->first    = NULL;
  gList->last     = NULL;
  return gList;
}

void free_gal_list(gal_list *gList)
{
  gal_node *gNode;
  if (gList) {
    while (gList->first != NULL) {
      gNode        = gList->first;
      gList->first = gNode->next;
      if (gNode->g) {free(gNode->g); gNode->g = NULL;}
      free(gNode); gNode = NULL;
    }
    free(gList); gList = NULL;
  }
  return;
}

gal_list_mat *initialize_gal_list_mat(int N_z_map, int nbTypes)
{
  gal_list_mat *gListMat = (gal_list_mat*)malloc(sizeof(gal_list_mat));
  gListMat->N1           = N_z_map;
  gListMat->N2           = nbTypes;
  gListMat->length       = gListMat->N1 * gListMat->N2;
  gListMat->total        = 0;
  gListMat->matrix       = (gal_list**)malloc(gListMat->length * sizeof(gal_list*));
  int i;
  for (i=0; i<gListMat->length; i++) gListMat->matrix[i] = initialize_gal_list();
  return gListMat;
}

void free_gal_list_mat(gal_list_mat *gListMat)
{
  int i;
  if (gListMat) {
    if (gListMat->matrix) {
      for (i=0; i<gListMat->length; i++) {free_gal_list(gListMat->matrix[i]); gListMat->matrix[i] = NULL;}
      free(gListMat->matrix); gListMat->matrix = NULL;
    }
    free(gListMat); gListMat = NULL;
  }
  return;
}

HPMap_arr *initialize_HPMap_arr(long long nside, int length)
{
  HPMap_arr *fullArr = (HPMap_arr*)malloc(sizeof(HPMap_arr));
  fullArr->length    = length;
  fullArr->array     = (HPMap_t**)malloc(length * sizeof(HPMap_t*));
  int i;
  for (i=0; i<length; i++) fullArr->array[i] = initialize_HPMap_t(nside);
  return fullArr;
}

void free_HPMap_arr(HPMap_arr *fullArr)
{
  int i;
  if (fullArr) {
    if (fullArr->array) {
      for (i=0; i<fullArr->length; i++) {free_HPMap_t(fullArr->array[i]); fullArr->array[i] = NULL;}
      free(fullArr->array); fullArr->array = NULL;
    }
    free(fullArr); fullArr = NULL;
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to sampling

void setTypeMap(MFP_param *mPar, type_map *tMap, HPMap_t *mask)
{
  long long i;
  int k;
  for (k=0; k<mPar->nbTypes; k++) {
    read_HPMap_t(mPar->maskPath[k], mask, 1); //-- verbose = 1
    for (i=0; i<mask->nbPix; i++) {
      if (mask->map[i] > 0) tMap->map[i] = SET_BIT_64(tMap->map[i], k); //-- 0 = masked, 1 = activated
    }
  }
  if (mPar->verbose < 3) printf("Set up angular selections\n");
  return;
}

void setRatioMaps(MFP_param *mPar, HPMap_arr *rMapArr)
{
  int k;
  for (k=0; k<mPar->nbTypes; k++) read_HPMap_t(mPar->maskPath[k], rMapArr->array[k], 1); //-- verbose = 1
  if (mPar->verbose < 3) printf("Set up angular selections\n");
  return;
}

void readAsciiNOfZ(char name[], interpolator_t *inter, int verbose)
{
  //-- Reset
  int j;
  for (j=0; j<inter->length; j++) {
    inter->x[j]     = 0.0;
    inter->value[j] = 0.0;
  }
  
  //-- Open
  int zInd   = 0;
  int nInd   = 1;
  int count  = 0;
  FILE *file = fopen(name, "r");
  
  char kv[STRING_LENGTH_MAX][STRING_LENGTH_MAX], buffer[STRING_LENGTH_MAX];
  double pos[2], z, ww, M;
  
  //-- Read
  int c = fgetc(file);
  while (c != EOF) {
    if (c == (int)'#') fgets(buffer, STRING_LENGTH_MAX, file);
    else {
      ungetc(c, file);
      fgets(buffer, STRING_LENGTH_MAX, file);
      getKeyAndValues(buffer, kv);
      
      inter->x[count]       = atof(kv[zInd]);
      inter->value[count+1] = atof(kv[nInd]); //-- To make cdf
      count++;
    }
    c = fgetc(file);
  }
  fclose(file);
  
  //-- Fill the rest of x
  for (j=count; j<inter->length; j++) inter->x[j] = inter->x[count-1];
  
  //-- Make cdf, no need to normalize
  for (j=1; j<inter->length; j++) inter->value[j] += inter->value[j-1];
  
  //-- Check
  if (verbose < 3) printf("Read \"%s\"   \n", name);
  return;
}

void setNOfZArr(MFP_param *mPar, interpolator_t *inter, sampler_arr *nOfZArr)
{
  double *bin_z_map = mPar->bin_z_map;
  sampler_t *nOfZ;
  double z, cdf_lower, cdf_upper;
  int j, k;
  
  for (k=0; k<nOfZArr->length; k++) {
    readAsciiNOfZ(mPar->nOfZPath[k], inter, 1);             //-- verbose = 1, reset inside
    nOfZ       = nOfZArr->array[k];
    
    z          = bin_z_map[0];
    cdf_upper  = execute_interpolator_t(inter, z, 1);       //-- border = 1 (constant)
    nOfZ->x[0] = z;
    
    for (j=1; j<nOfZ->length; j++) {
      cdf_lower      = cdf_upper;
      z              = bin_z_map[j];
      cdf_upper      = execute_interpolator_t(inter, z, 1); //-- border = 1 (constant)
      nOfZ->x[j]     = z;
      nOfZ->pdf[j-1] = cdf_upper - cdf_lower;
    }
    nOfZ->pdf[nOfZ->length-1] = 0.0;
    
    set_sampler_t(nOfZ, 4, 1);                                                   //-- mode = 4, setTotalToOne = 1 (see commonHeader.c)
    for (j=0; j<nOfZ->length; j++) nOfZ->pdf[j] *= mPar->A_pix * mPar->n_gal[k]; //-- Rescale by A_pix * n_gal[k]
  }
  
  if (mPar->verbose < 3) printf("Set up redshift selections\n");
  return;
}

void setDepthMaps(MFP_param *mPar, HPMap_arr *vdMapArr)
{
  int k;
  for (k=0; k<mPar->nbDepthMaps; k++) read_HPMap_t(mPar->depthMapPath[k], vdMapArr->array[k], 1); //-- verbose = 1
  if (mPar->verbose < 3) printf("Set up depth maps\n");
  return;
}

void setVDNOfZArr(MFP_param *mPar, interpolator_t *inter, sampler_arr *VD_nOfZArr)
{
  double *bin_z_map = mPar->bin_z_map;
  sampler_t *nOfZ;
  double z, cdf_lower, cdf_upper;
  int j, k;
  
  for (k=0; k<VD_nOfZArr->length; k++) {
    readAsciiNOfZ(mPar->VD_nOfZPath[k], inter, 1);          //-- verbose = 1, reset inside
    nOfZ       = VD_nOfZArr->array[k];
    
    z          = bin_z_map[0];
    cdf_upper  = execute_interpolator_t(inter, z, 1);       //-- border = 1 (constant)
    nOfZ->x[0] = z;
    
    for (j=1; j<nOfZ->length; j++) {
      cdf_lower      = cdf_upper;
      z              = bin_z_map[j];
      cdf_upper      = execute_interpolator_t(inter, z, 1); //-- border = 1 (constant)
      nOfZ->x[j]     = z;
      nOfZ->pdf[j-1] = cdf_upper - cdf_lower;
    }
    nOfZ->pdf[nOfZ->length-1] = 0.0;
    
    set_sampler_t(nOfZ, 4, 1);                                  //-- mode = 4, setTotalToOne = 1 (see commonHeader.c)
    for (j=0; j<nOfZ->length; j++) nOfZ->pdf[j] *= mPar->A_pix; //-- Rescale by A_pix
  }
  
  if (mPar->verbose < 3) printf("Set up variable depth redshift selections\n");
  return;
}

void set_gal_t(gal_t *g, double RA, double DEC, double z, double p, long long pix, double sigma_eps)
{
  g->RA        = RA;
  g->DEC       = DEC;
  g->z         = z;
  g->p         = p;
  g->pix       = pix;
  g->sigma_eps = sigma_eps; //-- sigma_eps comes in here because in the variable depth case, sigma_eps varies with pixels.
//   g->kappa;
//   g->gamma_1;
//   g->gamma_2;
//   g->e_1;
//   g->e_2;
  return;
}

void append_gal_list(gal_list *gList, double RA, double DEC, double z, double p, long long pix, double sigma_eps)
{
  if (gList->length == 0) {
    gList->first = initialize_gal_node();
    gList->last  = gList->first;
    gList->length++;
  }
  else if (gList->length == gList->size) {
    gList->last->next = initialize_gal_node();
    gList->last       = gList->last->next;
    gList->length++;
  }
  else if (gList->size == 0) gList->last = gList->first;
  else                       gList->last = gList->last->next;
  
  set_gal_t(gList->last->g, RA, DEC, z, p, pix, sigma_eps);
  gList->size++;
  return;
}

void samplePos(gsl_rng *generator, long long nside, long long pix, int N_max, int N, long long *lenArr, long long *cumLenArr, long long v1, long long v2, double *RAArr, double *DECArr, 
	       double z1, double z2, gal_list *gList, double sigma_eps)
{
  int q = N / N_max;
  int r = N % N_max;
  double p, z;
  int l, l2;
  
  for (l2=0; l2<q; l2++) {
    patchSampling3(generator, nside, pix, N_max, 0, lenArr, cumLenArr, v1, v2, RAArr, DECArr); //-- doNest = 0
    for (l=0; l<N_max; l++) {
      p = gsl_ran_flat(generator, 0.0, 1.0);
      z = z1 + p * (z2 - z1);
      append_gal_list(gList, RAArr[l], DECArr[l], z, p, pix, sigma_eps);
    }
  }
  
  patchSampling3(generator, nside, pix, r, 0, lenArr, cumLenArr, v1, v2, RAArr, DECArr); //-- doNest = 0
  for (l=0; l<r; l++) {
    p = gsl_ran_flat(generator, 0.0, 1.0);
    z = z1 + p * (z2 - z1);
    append_gal_list(gList, RAArr[l], DECArr[l], z, p, pix, sigma_eps);
  }
  return;
}

void samplePos_projCL(gsl_rng *generator, long long nside, long long pix, int N_max, int N, long long *lenArr, long long *cumLenArr, long long v1, long long v2, double *RAArr, double *DECArr, 
		      sampler_t *nOfZ, gal_list *gList, double sigma_eps)
{
  int q = N / N_max;
  int r = N % N_max;
  double p, z;
  int l, l2;
  
  for (l2=0; l2<q; l2++) {
    patchSampling3(generator, nside, pix, N_max, 0, lenArr, cumLenArr, v1, v2, RAArr, DECArr); //-- doNest = 0
    for (l=0; l<N_max; l++) {
      p = gsl_ran_flat(generator, 0.0, 1.0);
      z = execute_sampler_t(nOfZ, p, 4); //-- mode = 4 (continuous histogram & linear)
      append_gal_list(gList, RAArr[l], DECArr[l], z, p, pix, sigma_eps);
    }
  }
  
  patchSampling3(generator, nside, pix, r, 0, lenArr, cumLenArr, v1, v2, RAArr, DECArr); //-- doNest = 0
  for (l=0; l<r; l++) {
    p = gsl_ran_flat(generator, 0.0, 1.0);
    z = execute_sampler_t(nOfZ, p, 4); //-- mode = 4 (continuous histogram & linear)
    append_gal_list(gList, RAArr[l], DECArr[l], z, p, pix, sigma_eps);
  }
  return;
}

void sampleGalaxies(MFP_param *mPar, type_map *tMap, HPMap_arr *rMapArr, sampler_arr *nOfZArr, HPMap_arr *vdMapArr, sampler_arr *VD_nOfZArr, gal_list_mat *gListMat, HPMap_t *delta, int N_max)
{
  int nbTypes         = mPar->nbTypes;
  int nbDepthMaps     = mPar->nbDepthMaps;
  int N_depth         = mPar->N_depth;
  int nbTomo          = mPar->nbTomo;
  long long nside     = mPar->nside;
  long long nbPix     = mPar->nbPix;
  double *a_n_gal     = mPar->a_n_gal;
  double *b_n_gal     = mPar->b_n_gal;
  double *a_sigma_eps = mPar->a_sigma_eps;
  double *b_sigma_eps = mPar->b_sigma_eps;
  gsl_rng *generator  = mPar->generator;
  
  long long *tArr     = (tMap == NULL) ? NULL : tMap->map;
  float *dArr         = delta->map;
  sampler_t **sampArr = nOfZArr->array;
  
  long long *lenArr    = (long long*)malloc((4*nside-1) * sizeof(long long));
  long long *cumLenArr = (long long*)malloc((4*nside) * sizeof(long long));
  long long v1         = 2 * nside * (nside + 1);
  long long v2         = 2 * nside * (5 * nside - 1);
  double *RAArr        = (double*)malloc(N_max * sizeof(double));
  double *DECArr       = (double*)malloc(N_max * sizeof(double));
  nsideToLevels(nside, lenArr, cumLenArr);
  
  gal_list *gList;
  char name[STRING_LENGTH_MAX];
  double z1, z2, value, sigma_eps;
  long long i;
  int j, k, l, N;
  
  //-- Only for variable depth
  sampler_t **VD_sampArr = (VD_nOfZArr == NULL) ? NULL : VD_nOfZArr->array;
  double depth, n_gal;
  int m, m2, n, k2;
  
  //-- Loop over redshift slices
  for (j=0; j<gListMat->N1; j++) {
    //-- Read delta map
    sprintf(name, "%s%s_f1z%d.fits", mPar->denPrefix, mPar->runTag, j+1);
    read_HPMap_t(name, delta, 1); //-- verbose = 1
    z1 = mPar->bin_z_map[j];
    z2 = mPar->bin_z_map[j+1];
    
    //-- Loop over pixels
    for (i=0; i<nbPix; i++) {
      dArr[i] += 1.0; //-- delta to 1+delta
      
      //-- Loop over types
      for (k=0; k<gListMat->N2; k++) {
	if (k >= nbTypes) continue;
	sigma_eps = mPar->sigma_eps[k];
	
	if (tArr == NULL)                  value = rMapArr->array[k]->map[i]; //-- Use ratio maps
	else if (CHECK_BIT_64(tArr[i], k)) value = 1.0;                       //-- Use type map, ratio = 1
	else continue;                                                        //-- Masked
	
	value *= dArr[i] * sampArr[k]->pdf[j];      //-- Compute <N_gal> = ratio * (1+delta) * N(z) = ratio * (1+delta) * (p(z) * n_gal * A_pix)
	                                            //-- n_gal = N_gal / (A_pix * sum_i ratio)
	N      = gsl_ran_poisson(generator, value); //-- Number of galaxies for pixel i, redshift slice j, and type k
	if (N == 0) continue;
	gList  = gListMat->matrix[j+k*gListMat->N1];
	
	samplePos(generator, nside, i, N_max, N, lenArr, cumLenArr, v1, v2, RAArr, DECArr, z1, z2, gList, sigma_eps);
      }
      
      //-- Loop over depth maps
      for (m=0; m<nbDepthMaps; m++) {
	depth = vdMapArr->array[m]->map[i]; //-- Use depth maps
	if (depth < 1) continue; //-- Pixel masked
	
	//-- Determine which depth bin
	for (m2=-1; m2<N_depth; m2++) {
	  if (mPar->bin_depth[m2+1] > depth) break;
	}
	if (m2 < 0 || m2 >= N_depth) continue;
	
	//-- Loop over tomographic bins
	for (n=0; n<mPar->nbTomo; n++) {
	  k2        = m2 + n * N_depth;
	  n_gal     = a_n_gal[n] * depth + b_n_gal[n];
	  n_gal    *= RADIAN_SQ_TO_ARCMIN_SQ;
	  sigma_eps = a_sigma_eps[n] * depth + b_sigma_eps[n];
	  value     = dArr[i] * n_gal * VD_sampArr[k2]->pdf[j]; //-- Compute <N_gal> = (1+delta) * n_gal * pdf(z) = (1+delta) * n_gal * (p(z) * A_pix)
	  
	  N     = gsl_ran_poisson(generator, value); //-- Number of galaxies for pixel i, redshift slice j, and type k
	  if (N == 0) continue;
	  k     = nbTypes + m + n * nbDepthMaps;
	  gList = gListMat->matrix[j+k*gListMat->N1];
	  
	  samplePos(generator, nside, i, N_max, N, lenArr, cumLenArr, v1, v2, RAArr, DECArr, z1, z2, gList, sigma_eps);
	}
      }
    }
  }
  
  free(lenArr);
  free(cumLenArr);
  free(RAArr);
  free(DECArr);
  if (mPar->verbose < 3) printf("Did galaxy sampling\n");
  return;
}

void smoothDeltaMap(MFP_param *mPar, HPMap_t *delta, int resol)
{
  int length      = resol*resol;
  long long nside = mPar->nside;
  long long nbPix = mPar->nbPix;
  float *dArr     = delta->map;
  
  double value;
  long long pix, first, pixNest;
  
  for (first=0; first<nbPix; first+=length) {
    value = 0.0;
    
    for (pixNest=first; pixNest<first+length; pixNest++) {
      nest2ring64(nside, pixNest, &pix);
      value += dArr[pix];
    }
    
    value /= length;
    
    for (pixNest=first; pixNest<first+length; pixNest++) {
      nest2ring64(nside, pixNest, &pix);
      dArr[pix] = value;
    }
  }
  return;
}

void sampleGalaxies_projCL(MFP_param *mPar, type_map *tMap, HPMap_arr *rMapArr, sampler_arr *nOfZArr, HPMap_arr *vdMapArr, sampler_arr *VD_nOfZArr, gal_list_mat *gListMat, 
			   HPMap_t *delta, int N_max, int resol, int nbSplits)
{
  int N_z_map         = gListMat->N1 / nbSplits; //WARNING Special settings
  long long nside     = mPar->nside;
  long long nbPix     = mPar->nbPix;
  gsl_rng *generator  = mPar->generator;
  
  long long *tArr     = (tMap == NULL) ? NULL : tMap->map;
  float *dArr         = delta->map;
  sampler_t **sampArr = nOfZArr->array;
  
  long long *lenArr    = (long long*)malloc((4*nside-1) * sizeof(long long));
  long long *cumLenArr = (long long*)malloc((4*nside) * sizeof(long long));
  long long v1         = 2 * nside * (nside + 1);
  long long v2         = 2 * nside * (5 * nside - 1);
  double *RAArr        = (double*)malloc(N_max * sizeof(double));
  double *DECArr       = (double*)malloc(N_max * sizeof(double));
  nsideToLevels(nside, lenArr, cumLenArr);
  
  gal_list *gList;
  char name[STRING_LENGTH_MAX];
  double z1, z2, value, sigma_eps, factor;
  long long i;
  int j, k, l, N, j2;
  
  //-- Loop over redshift slices
  for (j=0; j<N_z_map; j++) {
    j2 = j * nbSplits;
    
    //-- Read delta map
    sprintf(name, "%s%s_f1z%d.fits", mPar->denPrefix, mPar->runTag, j+1);
    read_HPMap_t(name, delta, 1); //-- verbose = 1
    smoothDeltaMap(mPar, delta, resol);
    factor = mPar->A_pix * mPar->n_gal[j2];
    
    //-- Loop over pixels
    for (i=0; i<nbPix; i++) {
      dArr[i] += 1.0; //-- delta to 1+delta
      
      //-- Loop over types
      for (k=j2; k<j2+nbSplits; k++) {
	sigma_eps = mPar->sigma_eps[k];
	
	if (tArr == NULL)                  value = rMapArr->array[k]->map[i]; //-- Use ratio maps
	else if (CHECK_BIT_64(tArr[i], k)) value = 1.0;                       //-- Use type map, ratio = 1
	else continue;                                                        //-- Masked
	
	value *= dArr[i] * factor;                  //-- Compute <N_gal> = ratio * (1+delta) * N(z) = ratio * (1+delta) * (n_gal * A_pix)
						    //-- n_gal = N_gal / (A_pix * sum_i ratio)
	N      = gsl_ran_poisson(generator, value); //-- Number of galaxies for pixel i, redshift slice j, and type k
	if (N == 0) continue;
	gList  = gListMat->matrix[j+k*gListMat->N1];
	
	samplePos_projCL(generator, nside, i, N_max, N, lenArr, cumLenArr, v1, v2, RAArr, DECArr, sampArr[k], gList, sigma_eps);
      }
    }
  }
  
  free(lenArr);
  free(cumLenArr);
  free(RAArr);
  free(DECArr);
  if (mPar->verbose < 3) printf("Did galaxy sampling\n");
  return;
}

//----------------------------------------------------------------------
//-- Functions related to lensing

void readLensingMaps(char name[], HPMap_t *kappa, HPMap_t *gamma_1, HPMap_t *gamma_2, int verbose)
{
  FITS_t *fits = initializeTableReader_FITS_t(name);
  
  if (fits->nbRows * 1024 != kappa->nbPix) {
    printf("Wrong nside\n");
    return;
  }
  
  readTableColumn(fits, 0, (void*)kappa->map);   //-- colInd = 0
  readTableColumn(fits, 1, (void*)gamma_1->map); //-- colInd = 1
  readTableColumn(fits, 2, (void*)gamma_2->map); //-- colInd = 2
  free_FITS_t(fits);
  if (verbose) printf("Read \"%s\"\n", name);
  return;
}

double interpolateBetweenPixels(HPMap_t *full, long long neighbor[4], double weight[4])
{
  double value = 0.0;
  long long pix;
  int i;
  
  for (i=0; i<4; i++) {
    pix = neighbor[i];
    value += full->map[pix] * weight[i];
  }
  return value;
}

void interpolateLensing(MFP_param *mPar, gal_list_mat *gListMat, HPMap_t *k_lower, HPMap_t *g1_lower, HPMap_t *g2_lower, HPMap_t *k_upper, HPMap_t *g1_upper, HPMap_t *g2_upper)
{
  long long i;
  
  //-- Lensing is null at z = 0
  for (i=0; i<k_upper->nbPix; i++) {
    k_upper->map[i]  = 0.0;
    g1_upper->map[i] = 0.0;
    g2_upper->map[i] = 0.0;
  }
  
  gal_list *gList;
  gal_node *gNode;
  gal_t *g;
  
  long long neighbor[4];
  double thetaPhi[2], weight[4];
  float *buffer;
  double k_lower_itp, g1_lower_itp, g2_lower_itp, k_upper_itp, g1_upper_itp, g2_upper_itp;
  
  char name[STRING_LENGTH_MAX];
  double p;
  int j, k, l;
  
  //-- Loop over redshift
  for (j=0; j<gListMat->N1; j++) {
    //-- Swap
    buffer        = k_lower->map;
    k_lower->map  = k_upper->map;
    k_upper->map  = buffer;
    
    buffer        = g1_lower->map;
    g1_lower->map = g1_upper->map;
    g1_upper->map = buffer;
    
    buffer        = g2_lower->map;
    g2_lower->map = g2_upper->map;
    g2_upper->map = buffer;
    
    //-- Read 3 new maps
    sprintf(name, "%s%s_f2z%d.fits", mPar->lenPrefix, mPar->runTag, j+1);
    readLensingMaps(name, k_upper, g1_upper, g2_upper, 1); //-- verbose = 1
    
    //-- Loop over types
    for (k=0; k<gListMat->N2; k++) {
      gList = gListMat->matrix[j+k*gListMat->N1];
      
      if (k >= mPar->nbTypes || mPar->doLensing[k]) {
	//-- Loop over galaxies
	for (l=0, gNode=gList->first; l<gList->size; l++, gNode=gNode->next) {
	  g = gNode->g;
	  p = g->p;
	  
#ifdef __MFP_USE_HEALPIX_CXX__
	  //-- Interpolate between pixels
	  thetaPhi[0] = HALF_PI - g->DEC; //-- [rad]
	  thetaPhi[1] = g->RA;            //-- [rad]
	  getNgbAndWeight(thetaPhi, mPar->nside, neighbor, weight);
	  k_lower_itp  = interpolateBetweenPixels(k_lower, neighbor, weight);
	  g1_lower_itp = interpolateBetweenPixels(g1_lower, neighbor, weight);
	  g2_lower_itp = interpolateBetweenPixels(g2_lower, neighbor, weight);
	  k_upper_itp  = interpolateBetweenPixels(k_upper, neighbor, weight);
	  g1_upper_itp = interpolateBetweenPixels(g1_upper, neighbor, weight);
	  g2_upper_itp = interpolateBetweenPixels(g2_upper, neighbor, weight);
	  
	  //-- Interpolate between redshift slices
	  g->kappa   = k_lower_itp  + p * (k_upper_itp  - k_lower_itp);
	  g->gamma_1 = g1_lower_itp + p * (g1_upper_itp - g1_lower_itp);
	  g->gamma_2 = g2_lower_itp + p * (g2_upper_itp - g2_lower_itp);
#else
	  //-- Interpolate directly between redshift slices
	  i = g->pix;
	  g->kappa   = k_lower->map[i]  + p * (k_upper->map[i]  - k_lower->map[i]);
	  g->gamma_1 = g1_lower->map[i] + p * (g1_upper->map[i] - g1_lower->map[i]);
	  g->gamma_2 = g2_lower->map[i] + p * (g2_upper->map[i] - g2_lower->map[i]);
#endif
	  g->gamma_1 = -g->gamma_1; //-- The minus sign comes from a Flask bug.
	  g->gamma_2 = -g->gamma_2; //-- The minus sign comes from a Flask bug.
	}
      }
      
      else {
	//-- Loop over galaxies, no lensing
	for (l=0, gNode=gList->first; l<gList->size; l++, gNode=gNode->next) {
	  g = gNode->g;
	  g->kappa   = 0.0;
	  g->gamma_1 = 0.0;
	  g->gamma_2 = 0.0;
	}
      }
    }
  }
  
  if (mPar->verbose < 3) {
#ifdef __MFP_USE_HEALPIX_CXX__
    printf("Has compiled with healpix_cxx\n");
    printf("Interpolated lensing signals between pixels and redshifts\n");
#else
    printf("Has compiled without healpix_cxx\n");
    printf("Interpolated lensing signals only between redshifts, not between pixels\n");
#endif
  }
  return;
}

void assignLensing(MFP_param *mPar, gal_list_mat *gListMat, HPMap_t *kappa, HPMap_t *gamma_1, HPMap_t *gamma_2)
{
  gal_list *gList;
  gal_node *gNode;
  gal_t *g;
  
  long long neighbor[4];
  double thetaPhi[2], weight[4];
  
  char name[STRING_LENGTH_MAX];
  long long i;
  int j, k, l;
  
  //-- Loop over redshift
  for (j=0; j<gListMat->N1; j++) {
    //-- Read 3 new maps
    sprintf(name, "%s%s_f2z%d.fits", mPar->lenPrefix, mPar->runTag, j+1);
    readLensingMaps(name, kappa, gamma_1, gamma_2, 1); //-- verbose = 1
    
    //-- Loop over types
    for (k=0; k<gListMat->N2; k++) {
      gList = gListMat->matrix[j+k*gListMat->N1];
      
      if (k >= mPar->nbTypes || mPar->doLensing[k]) {
	//-- Loop over galaxies
	for (l=0, gNode=gList->first; l<gList->size; l++, gNode=gNode->next) {
	  g = gNode->g;
	  
#ifdef __MFP_USE_HEALPIX_CXX__
	  //-- Interpolate between pixels
	  thetaPhi[0] = HALF_PI - g->DEC; //-- [rad]
	  thetaPhi[1] = g->RA;            //-- [rad]
	  getNgbAndWeight(thetaPhi, mPar->nside, neighbor, weight);
	  g->kappa   = interpolateBetweenPixels(kappa, neighbor, weight);
	  g->gamma_1 = interpolateBetweenPixels(gamma_1, neighbor, weight);
	  g->gamma_2 = interpolateBetweenPixels(gamma_2, neighbor, weight);
#else
	  //-- Assign directly values
	  i = g->pix;
	  g->kappa   = kappa->map[i];
	  g->gamma_1 = gamma_1->map[i];
	  g->gamma_2 = gamma_2->map[i];
#endif
	  g->gamma_1 = -g->gamma_1; //-- The minus sign comes from a Flask bug.
	  g->gamma_2 = -g->gamma_2; //-- The minus sign comes from a Flask bug.
	}
      }
      
      else {
	//-- Loop over galaxies, no lensing
	for (l=0, gNode=gList->first; l<gList->size; l++, gNode=gNode->next) {
	  g = gNode->g;
	  g->kappa   = 0.0;
	  g->gamma_1 = 0.0;
	  g->gamma_2 = 0.0;
	}
      }
    }
  }
  
  if (mPar->verbose < 3) {
#ifdef __MFP_USE_HEALPIX_CXX__
    printf("Has compiled with healpix_cxx\n");
    printf("Interpolated lensing signals between pixels\n");
#else
    printf("Has compiled without healpix_cxx\n");
    printf("Assigned lensing signals, no interpolation between pixels\n");
#endif
  }
  return;
}

void assignLensing_projCL(MFP_param *mPar, gal_list_mat *gListMat, HPMap_t *kappa, HPMap_t *gamma_1, HPMap_t *gamma_2, int nbSplits)
{
  int N_z_map = gListMat->N1 / nbSplits; //WARNING Special settings
  
  gal_list *gList;
  gal_node *gNode;
  gal_t *g;
  
  long long neighbor[4];
  double thetaPhi[2], weight[4];
  
  char name[STRING_LENGTH_MAX];
  long long i;
  int j, k, l, j2;
  
  //-- Loop over redshift
  for (j=0; j<N_z_map; j++) {
    j2 = j * nbSplits;
    
    if (mPar->doLensing[j2]) {
      //-- Read 3 new maps
      sprintf(name, "%s%s_f2z%d.fits", mPar->lenPrefix, mPar->runTag, j+1);
      readLensingMaps(name, kappa, gamma_1, gamma_2, 1); //-- verbose = 1
      
      //-- Loop over types
      for (k=j2; k<j2+nbSplits; k++) {
	gList = gListMat->matrix[j+k*gListMat->N1];
	
	//-- Loop over galaxies
	for (l=0, gNode=gList->first; l<gList->size; l++, gNode=gNode->next) {
	  g = gNode->g;
	  
#ifdef __MFP_USE_HEALPIX_CXX__
	  //-- Interpolate between pixels
	  thetaPhi[0] = HALF_PI - g->DEC; //-- [rad]
	  thetaPhi[1] = g->RA;            //-- [rad]
	  getNgbAndWeight(thetaPhi, mPar->nside, neighbor, weight);
	  g->kappa   = interpolateBetweenPixels(kappa, neighbor, weight);
	  g->gamma_1 = interpolateBetweenPixels(gamma_1, neighbor, weight);
	  g->gamma_2 = interpolateBetweenPixels(gamma_2, neighbor, weight);
#else
	  //-- Assign directly values
	  i = g->pix;
	  g->kappa   = kappa->map[i];
	  g->gamma_1 = gamma_1->map[i];
	  g->gamma_2 = gamma_2->map[i];
#endif
	  g->gamma_1 = -g->gamma_1; //-- The minus sign comes from a Flask bug.
	  g->gamma_2 = -g->gamma_2; //-- The minus sign comes from a Flask bug.
	}
      }
    }
    
    else {
      //-- Loop over types
      for (k=j2; k<j2+nbSplits; k++) {
	gList = gListMat->matrix[j+k*gListMat->N1];
	
	//-- Loop over galaxies, no lensing
	for (l=0, gNode=gList->first; l<gList->size; l++, gNode=gNode->next) {
	  g = gNode->g;
	  g->kappa   = 0.0;
	  g->gamma_1 = 0.0;
	  g->gamma_2 = 0.0;
	}
      }
    }
  }
  
  if (mPar->verbose < 3) {
#ifdef __MFP_USE_HEALPIX_CXX__
    printf("Has compiled with healpix_cxx\n");
    printf("Interpolated lensing signals between pixels\n");
#else
    printf("Has compiled without healpix_cxx\n");
    printf("Assigned lensing signals, no interpolation between pixels\n");
#endif
  }
  return;
}

void makeG(MFP_param *mPar, gal_list_mat *gListMat)
{
  gal_list *gList;
  gal_node *gNode;
  gal_t *g;
  double factor;
  int j, k, l;
  
  for (k=0; k<gListMat->N2; k++) {
    if (k < mPar->nbTypes && mPar->doLensing[k] == 0) continue;
    
    for (j=0; j<gListMat->N1; j++) {
      gList = gListMat->matrix[j+k*gListMat->N1];
      
      for (l=0, gNode=gList->first; l<gList->size; l++, gNode=gNode->next) {
	g = gNode->g;
	factor = 1.0 / (1.0 - g->kappa);
	g->gamma_1 *= factor;
	g->gamma_2 *= factor;
      }
    }
  }
  
  if (mPar->verbose < 3) printf("Computed reduced shear\n");
  return;
}

void addNoise(MFP_param *mPar, gal_list_mat *gListMat)
{
  gal_list *gList;
  gal_node *gNode;
  gal_t *g;
  double g1, g2, e1, e2, A, B, C, g1_sq, g2_sq, factor;
  int j, k, l;
  
  if (mPar->doNoise) {
    //-- Let e be source ellipticity, epsilon be observed ellipticity
    //-- epsilon = (e + g) / (1 + g^* e)
    //-- where
    //--   e   = e1 + i e2
    //--   g   = g1 + i g2
    //--   g^* = g1 - i g2
    //--
    //-- If x + iy = (a + ib) / (c + id)
    //-- then x = (ac + bd) / (c^2 + d^2)
    //--      y = (bc - ad) / (c^2 + d^2)
    //-- Here,
    //--   a = e1 + g1
    //--   b = e2 + g2
    //--   c = 1 + g1 e1 + g2 e2
    //--   d = g1 e2 - g2 e1
    //-- such that,
    //--   ac + bd   = e1 + g1 + g1(e1^2 + e2^2) + 2 g1 g2 e2 + e1(g1^2 - g2^2)
    //--   bc - ad   = e2 + g2 + g2(e1^2 + e2^2) + 2 g1 g2 e1 - e2(g1^2 - g2^2)
    //--   c^2 + d^2 = 1 + 2 g1 e1 + 2 g2 e2 + (g1^2 + g2^2)(e1^2 + e2^2)
    //--   
    //-- So,
    //--   epsilon1 = [e1 + g1 + g1 C + g1 B + e1(g1^2 - g2^2)] / [1 + A + B + C(g1^2 + g2^2)]
    //--   epsilon2 = [e2 + g2 + g2 C + g2 A - e2(g1^2 - g2^2)] / [1 + A + B + C(g1^2 + g2^2)]
    //-- where
    //--   A = 2 g1 e1
    //--   B = 2 g2 e2
    //--   C = e1^2 + e2^2
    //-- Finally,
    //--   epsilon1 = [e1(1 + g1^2 - g2^2) + g1(1 + C + B)] / [1 + A + B + C(g1^2 + g2^2)]
    //--   epsilon2 = [e2(1 - g1^2 + g2^2) + g2(1 + C + A)] / [1 + A + B + C(g1^2 + g2^2)]
    
    for (k=0; k<gListMat->N2; k++) {
      if (k < mPar->nbTypes && mPar->doLensing[k] == 0) continue;
      
      for (j=0; j<gListMat->N1; j++) {
	gList = gListMat->matrix[j+k*gListMat->N1];
	
	for (l=0, gNode=gList->first; l<gList->size; l++, gNode=gNode->next) {
	  g  = gNode->g;
	  g1 = g->gamma_1; //-- Already reduced shear
	  g2 = g->gamma_2;
	  
	  do e1 = gsl_ran_gaussian(mPar->generator, g->sigma_eps);
	  while (fabs(e1) >= 1.0);
	  do e2 = gsl_ran_gaussian(mPar->generator, g->sigma_eps);
	  while (fabs(e2) >= 1.0);
	  
	  A = 2.0 * g1 * e1;
	  B = 2.0 * g2 * e2;
	  C = e1 * e1 + e2 * e2;
	  g1_sq = g1 * g1;
	  g2_sq = g2 * g2;
	  factor = 1.0 / (1.0 + A + B + (g1_sq + g2_sq) * C);
	  
	  g->e_1 = factor * (e1 * (1.0 + g1_sq - g2_sq) + g1 * (1.0 + C + B));
	  g->e_2 = factor * (e2 * (1.0 - g1_sq + g2_sq) + g2 * (1.0 + C + A));
	}
      }
    }
  }
  
  if (mPar->verbose < 3) {
    if (mPar->doNoise) printf("Added noise\n");
    else               printf("No noise\n");
  }
  return;
}

void flipSign(MFP_param *mPar, gal_list_mat *gListMat)
{
  gal_list *gList;
  gal_node *gNode;
  gal_t *g;
  int j, k, l;
  
  //-- Athena convention, do nothing
  if (mPar->signConv == 0) ;
  
  //-- treecorr convention, flip g_2
  else if (mPar->signConv == 1) {
    for (k=0; k<gListMat->N2; k++) {
      if (k < mPar->nbTypes && mPar->doLensing[k] == 0) continue;
      
      for (j=0; j<gListMat->N1; j++) {
	gList = gListMat->matrix[j+k*gListMat->N1];
	
	for (l=0, gNode=gList->first; l<gList->size; l++, gNode=gNode->next) {
	  g = gNode->g;
	  g->gamma_2 = -g->gamma_2;
	  g->e_2     = -g->e_2;
	}
      }
    }
  }
  
  //-- Flask convention, flip g_1 & g_2
  else if (mPar->signConv == 1) {
    for (k=0; k<gListMat->N2; k++) {
      if (k < mPar->nbTypes && mPar->doLensing[k] == 0) continue;
      
      for (j=0; j<gListMat->N1; j++) {
	gList = gListMat->matrix[j+k*gListMat->N1];
	
	for (l=0, gNode=gList->first; l<gList->size; l++, gNode=gNode->next) {
	  g = gNode->g;
	  g->gamma_1 = -g->gamma_1;
	  g->gamma_2 = -g->gamma_2;
	  g->e_1     = -g->e_1;
	  g->e_2     = -g->e_2;
	}
      }
    }
  }
  
  if (mPar->verbose < 3) {
    if      (mPar->signConv == 0) printf("Flipped sign to the Athena convention\n");
    else if (mPar->signConv == 1) printf("Flipped sign to the treecorr convention\n");
    else if (mPar->signConv == 2) printf("Flipped sign to the Flask convention\n");
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to outputs

void outFits_gal_t(FITS_t *fits, gal_t *g, double factor, int doNoise, int doWgt)
{
  float fBuff;
  fBuff = g->RA * factor;  writeTableColumn(fits, 0, 1, &fBuff);
  fBuff = g->DEC * factor; writeTableColumn(fits, 1, 1, &fBuff);
  fBuff = g->z;            writeTableColumn(fits, 2, 1, &fBuff);
  
  if (doNoise == 0) {
    fBuff = g->gamma_1;    writeTableColumn(fits, 3, 1, &fBuff);
    fBuff = g->gamma_2;    writeTableColumn(fits, 4, 1, &fBuff);
  }
  else if (doNoise == 1) {
    fBuff = g->e_1;        writeTableColumn(fits, 3, 1, &fBuff);
    fBuff = g->e_2;        writeTableColumn(fits, 4, 1, &fBuff);
  }
  else if (doNoise == 2) {
    fBuff = g->gamma_1;    writeTableColumn(fits, 3, 1, &fBuff);
    fBuff = g->gamma_2;    writeTableColumn(fits, 4, 1, &fBuff);
    fBuff = g->e_1;        writeTableColumn(fits, 5, 1, &fBuff);
    fBuff = g->e_2;        writeTableColumn(fits, 6, 1, &fBuff);
  }
  
  if (doWgt == 1) {
    fBuff = 1.0;
    if (doNoise == 0 || doNoise == 1) writeTableColumn(fits, 5, 1, &fBuff);
    else if (doNoise == 2)            writeTableColumn(fits, 7, 1, &fBuff);
  }
  nextRow(fits);
  return;
}

void outFitsParam(FITS_t *fits, MFP_param *mPar, int k2)
{
  int outStyle   = mPar->outStyle;
  int *doLensing = mPar->doLensing;
  
  char sBuff[STRING_LENGTH_MAX];
  double lfBuff;
  int k, m, n, m2;
  
  addLineSpread(fits);
  addKeyword(fits, TSTRING,   "PARPATH",  mPar->parPath,       "[-] Path of the parameter file");
  addKeyword(fits, TSTRING,   "SEED",     mPar->seed,          "[-] Seed for generating random numbers");
  addKeyword(fits, TSTRING,   "RUNTAG",   mPar->runTag,        "[-] Tag for identifying different runs");
  
  addLineSpread(fits);
  addKeyword(fits, TLONGLONG, "NSIDE",    &mPar->nside,        "[-] N_side of input maps");
  addKeyword(fits, TLONGLONG, "NBPIX",    &mPar->nbPix,        "[-] Number of pixels");
  lfBuff = mPar->A_pix * RADIAN_SQ_TO_ARCMIN_SQ;
  addKeyword(fits, TDOUBLE,   "APIX",     &lfBuff,             "[arcmin^2] Area of pixels");
  
  addLineSpread(fits);
  addKeyword(fits, TINT,      "NZMAP",    &mPar->N_z_map,      "[-] Number of redshift slices of maps");
  addKeyword(fits, TDOUBLE,   "ZMAPMIN",  &mPar->zMapRange[0], "[-] Minimal map redshift");
  addKeyword(fits, TDOUBLE,   "ZMAPMAX",  &mPar->zMapRange[1], "[-] Maximal map redshift");
  addKeyword(fits, TDOUBLE,   "DZMAP",    &mPar->zMapRange[2], "[-] Map redshift binwidth (can be 0)");
  
  addLineSpread(fits);
  addKeyword(fits, TINT,      "NBTYPES",  &mPar->totNbTypes,   "[-] Number of types of galaxies");
  if (outStyle == 64) {
    if (k2 < mPar->nbTypes) {
      addKeyword(fits, TSTRING, "MASKPATH", mPar->maskPath[k2],  "[-] Path to mask for type k");
      addKeyword(fits, TSTRING, "NOFZPATH", mPar->nOfZPath[k2],  "[-] Path to n(z) for type k");
    }
  }
  for (k=0; k<mPar->nbTypes; k++) {
    if (outStyle == 1 || (outStyle == 2 && doLensing[k] == k2) || (outStyle == 64 && k == k2)) {
      sprintf(sBuff, "NGAL%d", k);
      lfBuff = mPar->n_gal[k] / RADIAN_SQ_TO_ARCMIN_SQ;
      addKeyword(fits, TDOUBLE, sBuff,    &lfBuff,             "[-] Galaxy number density for type k");
    }
  }
  
  addLineSpread(fits);
  addKeyword(fits, TINT,      "DONOISE",  &mPar->doNoise,      "[-] 0 = no, 1 = yes, 2 = output both");
  if (k2<mPar->nbTypes) {
    addKeyword(fits, TDOUBLE, "SIGMAEPS", &mPar->sigma_eps,    "[-] Ellipticity dispersion");
  }
  for (k=0; k<mPar->nbTypes; k++) {
    if (outStyle == 1 || (outStyle == 2 && doLensing[k] == k2) || (outStyle == 64 && k == k2)) {
      sprintf(sBuff, "DOLEN%d", k);
      addKeyword(fits, TINT,  sBuff,      &doLensing[k],       "[-] Do lensing for type k");
    }
  }
  addKeyword(fits, TINT,      "OUTSTYLE", &mPar->outStyle,     "[-] Output style");
  
  //-- Variable depth
  if (outStyle == 64) {
    if (k2 >= mPar->nbTypes) {
      m = (k2 - mPar->nbTypes) % mPar->nbDepthMaps;
      n = (k2 - mPar->nbTypes) / mPar->nbDepthMaps;
      addKeyword(fits, TINT,    "NDEPTH",   &mPar->N_depth,         "[-] Number of depth bins");
      addKeyword(fits, TSTRING, "DMAPPATH", mPar->depthMapPath[m],  "[-] Path to depth map for type k");
      addKeyword(fits, TDOUBLE, "ANGAL",    &mPar->a_n_gal[n],      "[-] Slope parameter for n_gal = a*depth + b");
      addKeyword(fits, TDOUBLE, "BNGAL",    &mPar->b_n_gal[n],      "[-] Intercept parameter for n_gal = a*depth + b");
      addKeyword(fits, TDOUBLE, "ASIGEPS",  &mPar->a_sigma_eps[n],  "[-] Slope parameter for sigma_eps = a*depth + b");
      addKeyword(fits, TDOUBLE, "BSIGEPS",  &mPar->b_sigma_eps[n],  "[-] Intercept parameter for sigma_eps = a*depth + b");
      for (m2=0; m2<mPar->N_depth; m2++) {
	k = m2 + n * mPar->N_depth;
	addKeyword(fits, TSTRING, "VDNZPATH", mPar->VD_nOfZPath[k], "[-] Path to variable depth n(z)");
      }
    }
  }
  return;
}

void outFitsGalListMat(MFP_param *mPar, gal_list_mat *gListMat, int verbose)
{
  if (mPar->verbose < 2) printf("Outputing...\r");
  
  int nbFiles             = (mPar->outStyle == 64) ? mPar->totNbTypes : mPar->outStyle;
  char *defaultFileTag[3] = {"_all", "_noLen", "_len"};
  long long *nbGal        = (long long*)malloc(nbFiles * sizeof(long long));
  FITS_t **fitsArr        = (FITS_t**)malloc(nbFiles * sizeof(FITS_t*));
  
  char name[STRING_LENGTH_MAX];
  char fileTag[STRING_LENGTH_MAX];
  FITS_t *fits;
  int k2;
  
  for (k2=0; k2<nbFiles; k2++) {
    nbGal[k2] = 0;
    if      (mPar->outStyle == 1)            sprintf(fileTag, "%s", defaultFileTag[0]);
    else if (mPar->outStyle == 2 && k2 == 0) sprintf(fileTag, "%s", defaultFileTag[1]);
    else if (mPar->outStyle == 2 && k2 == 1) sprintf(fileTag, "%s", defaultFileTag[2]);
    else                                     sprintf(fileTag, "_type%d", k2);
    
    sprintf(name, "%s%s%s.fits", mPar->outPrefix, mPar->runTag, fileTag);
    fitsArr[k2] = initializeTableWriter_FITS_t(name);
    fits = fitsArr[k2];
    
    addColumn(fits, "ALPHA_J2000", TFLOAT, "deg");
    addColumn(fits, "DELTA_J2000", TFLOAT, "deg");
    addColumn(fits, "z_spec",      TFLOAT, "-       ");
    if (mPar->doNoise == 0) {
      addColumn(fits, "g1",       TFLOAT, "-       ");
      addColumn(fits, "g2",       TFLOAT, "-       ");
    }  
    else if (mPar->doNoise == 1) {
      addColumn(fits, "e1",       TFLOAT, "-       ");
      addColumn(fits, "e2",       TFLOAT, "-       ");
    }  
    else if (mPar->doNoise == 2) {
      addColumn(fits, "g1",       TFLOAT, "-       ");
      addColumn(fits, "g2",       TFLOAT, "-       ");
      addColumn(fits, "e1",       TFLOAT, "-       ");
      addColumn(fits, "e2",       TFLOAT, "-       ");
    }
    if (mPar->doWgt == 1) {
      addColumn(fits, "weight",   TFLOAT, "-       ");
    }
  }
  
  gal_list *gList;
  gal_node *gNode;
  int j, k, l;
  
  //-- Loop over types
  for (k=0; k<gListMat->N2; k++) {
    //-- Determine which file to output
    if      (mPar->outStyle == 1) k2 = 0;
    else if (mPar->outStyle == 2) k2 = (k >= mPar->nbTypes || mPar->doLensing[k]);
    else                          k2 = k;
    
    for (j=0; j<gListMat->N1; j++) {
      gList      = gListMat->matrix[j+k*gListMat->N1];
      nbGal[k2] += gList->size;
      for (l=0, gNode=gList->first; l<gList->size; l++, gNode=gNode->next) outFits_gal_t(fitsArr[k2], gNode->g, RADIAN_TO_DEGREE, mPar->doNoise, mPar->doWgt); //-- factor = RADIAN_TO_DEGREE
    }
  }
  
  for (k2=0; k2<nbFiles; k2++) {
    fits = fitsArr[k2];
    addLineSpread(fits);
    addKeyword(fits, TLONGLONG, "NBGAL", &nbGal[k2], "[-] Number of galaxies in this file");
    addLineSpread(fits);
    outFitsParam(fits, mPar, k2);
    
    free_FITS_t(fits);
    if      (mPar->outStyle == 1)            sprintf(fileTag, "%s", defaultFileTag[0]);
    else if (mPar->outStyle == 2 && k2 == 0) sprintf(fileTag, "%s", defaultFileTag[1]);
    else if (mPar->outStyle == 2 && k2 == 1) sprintf(fileTag, "%s", defaultFileTag[2]);
    else                                     sprintf(fileTag, "_type%d", k2);
    sprintf(name, "%s%s%s.fits", mPar->outPrefix, mPar->runTag, fileTag);
    if (verbose) printf("Outputed \"%s\"\n", name);
  }
  
  free(nbGal);
  free(fitsArr);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to lensing maps

double E_sq_of_z(double z, double Omega_m)
{
  return (1.0-Omega_m) + Omega_m * pow(1+z, 3);
}

void fillWInterpolator(interpolator_t *w_inter, double dz, double Omega_m)
{
  w_inter->x[0] = 0.0;;
  w_inter->value[0] = 0.0;
  w_inter->dx = dz;
  
  double E_sq;
  int i;
  
  for (i=1; i<w_inter->length; i++) {
    w_inter->x[i]     = w_inter->x[i-1] + dz;
    E_sq = E_sq_of_z(w_inter->x[i-1]+0.5*dz, Omega_m);
    w_inter->value[i] = w_inter->value[i-1] + 1.0 / sqrt(E_sq);
  }
  
  double factor = HUBBLE_DISTANCE * w_inter->dx;
  for (i=1; i<w_inter->length; i++) w_inter->value[i] *= factor;
  return;
}

double comovDist(interpolator_t *w_inter, double z)
{
  return execute_interpolator_t(w_inter, z, 2); //-- border = 2 (linear extrapolation)
}

void makeKappaMap(MFP_param *mPar, HPMap_t *delta, HPMap_t *kappa, interpolator_t *w_inter, double Omega_m, int zInd)
{
  //-- Reset
  long long i;
  for (i=0; i<kappa->nbPix; i++) kappa->map[i] = 0.0;
  
  double z_s    = mPar->bin_z_map[zInd+1];
  double w_s    = comovDist(w_inter, z_s);
  double factor = FOUR_PI_G_OVER_C2 * CRITICAL_DENSITY * HUBBLE_DISTANCE * Omega_m / w_s;
  
  char name[STRING_LENGTH_MAX];
  double z_l, dz_l, w_l, E_sq, kernel;
  int j;
  
  for (j=0; j<=zInd; j++) {
    sprintf(name, "%s%s_f1z%d.fits", mPar->denPrefix, mPar->runTag, j+1);
    read_HPMap_t(name, delta, mPar->verbose<3);
    
    z_l    = mPar->z_map[j];
    dz_l   = 2.0 * mPar->half_dz_map[j];
    w_l    = comovDist(w_inter, z_l);
    E_sq   = E_sq_of_z(z_l, Omega_m);
    kernel = fmax(0.0, w_s-w_l) * w_l * (1.0 + z_l) * dz_l / sqrt(E_sq);
    
    for (i=0; i<kappa->nbPix; i++) kappa->map[i] += kernel * delta->map[i];
  }
  
  for (i=0; i<kappa->nbPix; i++) kappa->map[i] *= factor;
  return;
}

void kappaToGamma(MFP_param *mPar, HPMap_t *kappa, HPMap_t *gamma1, HPMap_t *gamma2, double_mat *kAlm, double_mat *gAlm, double_arr *weight, int l_maxx)
{
#ifdef __MFP_USE_HEALPIX_CXX__
  mapToAlm(mPar->nside, kappa->map, l_maxx, kAlm->matrix, weight->array);
  
  double factor;
  int l, m, index;
  
  //-- Hu (2000) - PRD, 62, 043007
  for (m=0; m<=l; m++) {
    for (l=m; l<=l_maxx; l++) {
      index = m*(2*l_maxx+1-m)/2+l;
      if (l < 2) { //-- No monopole & dipole for shear
	gAlm->matrix[0+2*index] = 0.0; //-- Real part
	gAlm->matrix[1+2*index] = 0.0; //-- Imaginary part
      }
      else {
	factor = -sqrt((l + 2.0) * (l - 1.0) / (l * (l + 1.0))); //-- Be careful here for the sign
	gAlm->matrix[0+2*index] = kAlm->matrix[0+2*index] * factor;
	gAlm->matrix[1+2*index] = kAlm->matrix[1+2*index] * factor;
      }
    }
  }
  
  almToMap_spin2(l_maxx, gAlm->matrix, mPar->nside, gamma1->map, gamma2->map);
  if (mPar->verbose < 3) printf("Transformed kappa to gamma\n");
#endif
  return;
}

void flipSignForFlask(HPMap_t *gamma1, HPMap_t *gamma2)
{
  long long i;
  for (i=0; i<gamma1->nbPix; i++) {
    gamma1->map[i] *= -1;
    gamma2->map[i] *= -1;
  }
  return;
}

void outFitsLensingMaps(MFP_param *mPar, HPMap_t *kappa, HPMap_t *gamma1, HPMap_t *gamma2, int zInd, int verbose)
{
  char name[STRING_LENGTH_MAX];
  sprintf(name, "%s%s_f2z%d.fits", mPar->lenPrefix, mPar->runTag, zInd+1);
  
  FITS_t *fits     = initializeTableWriter_FITS_t(name);
  long long nbRows = kappa->nbPix / 1024;
  addWideColumn(fits, "KAPPA",  TFLOAT, 1024, "-       "); //-- nb = 1024
  addWideColumn(fits, "GAMMA1", TFLOAT, 1024, "-       "); //-- nb = 1024
  addWideColumn(fits, "GAMMA2", TFLOAT, 1024, "-       "); //-- nb = 1024
  
  int i;
  for (i=0; i<nbRows; i++) {
    writeTableColumn(fits, 0, 1024, kappa->map+i*1024); //-- colInd = 0, nbData = 1024
    writeTableColumn(fits, 1, 1024, gamma1->map+i*1024); //-- colInd = 0, nbData = 1024
    writeTableColumn(fits, 2, 1024, gamma2->map+i*1024); //-- colInd = 0, nbData = 1024
    nextRow(fits);
  }
  
  free_FITS_t(fits);
  if (verbose) printf("Outputed \"%s\"\n\n", name);
  return;
}

//----------------------------------------------------------------------
//-- Main function

#define N_MAX_1 4000
#define N_MAX_2 10000

void processLensingMaps(MFP_param *mPar)
{
#ifdef __MFP_USE_HEALPIX_CXX__
  //-- Initialize interpolator
  double Omega_m = 0.2905;
  double dz = 0.0002;
  interpolator_t *w_inter = initialize_interpolator_t(10001);
  fillWInterpolator(w_inter, dz, Omega_m);
  
  //-- Initialize multiples
  int l_maxx         = 3 * mPar->nside - 1;
  int N_l            = (l_maxx+1) * (l_maxx+2) / 2;
  double_mat *kAlm   = initialize_double_mat(2, N_l);
  double_mat *gAlm   = initialize_double_mat(2, N_l);
  double_arr *weight = initialize_double_arr(2*mPar->nside);
  readWeight(mPar->nside, weight->array);
  if (mPar->verbose < 3) printf("Read HEALPix ring weights\n");
  
  HPMap_arr *bufferMap = initialize_HPMap_arr(mPar->nside, 3);
  int j;
  
  for (j=0; j<mPar->N_z_map; j++) {
    makeKappaMap(mPar, bufferMap->array[1], bufferMap->array[0], w_inter, Omega_m, j);
    kappaToGamma(mPar, bufferMap->array[0], bufferMap->array[1], bufferMap->array[2], kAlm, gAlm, weight, l_maxx);
    flipSignForFlask(bufferMap->array[1], bufferMap->array[2]); //-- WARNING Go to Flask sign convention
    outFitsLensingMaps(mPar, bufferMap->array[0], bufferMap->array[1], bufferMap->array[2], j, mPar->verbose<3);
  }
  
  free_interpolator_t(w_inter);
  free_double_mat(kAlm);
  free_double_mat(gAlm);
  free_double_arr(weight);
  free_HPMap_arr(bufferMap);
#else
  printf("healpix_cxx not found; nothing is done.");
#endif
  printf("------------------------------------------------------------------------\n");
  return;
}

void processMock_typeMap_LOS(MFP_param *mPar)
{
  type_map *tMap         = initialize_type_map(mPar->nside);
  interpolator_t *inter  = initialize_interpolator_t(N_MAX_1);
  sampler_arr *nOfZArr   = initialize_sampler_arr(mPar->N_z_map+1, mPar->nbTypes); //-- +1 is necessary
  gal_list_mat *gListMat = initialize_gal_list_mat(mPar->N_z_map, mPar->nbTypes);
  HPMap_arr *bufferMap   = initialize_HPMap_arr(mPar->nside, 6);
  
  reset_type_map(tMap);
  setTypeMap(mPar, tMap, bufferMap->array[0]);
  setNOfZArr(mPar, inter, nOfZArr);
  sampleGalaxies(mPar, tMap, NULL, nOfZArr, NULL, NULL, gListMat, bufferMap->array[0], N_MAX_2);
  
  if (mPar->skipLensing == 1) {
    if (mPar->verbose < 3) printf("Skipped lensing module\n");
  }
  else {
    interpolateLensing(mPar, gListMat, bufferMap->array[0], bufferMap->array[1], bufferMap->array[2], bufferMap->array[3], bufferMap->array[4], bufferMap->array[5]);
    makeG(mPar, gListMat);
    addNoise(mPar, gListMat);
    flipSign(mPar, gListMat);
  }
  
  outFitsGalListMat(mPar, gListMat, mPar->verbose<3);
  
  free_type_map(tMap);
  free_interpolator_t(inter);
  free_sampler_arr(nOfZArr);
  free_gal_list_mat(gListMat);
  free_HPMap_arr(bufferMap);
  printf("------------------------------------------------------------------------\n");
  return;
}

void processMock_ratioMap_LOS(MFP_param *mPar)
{
  HPMap_arr *rMapArr     = initialize_HPMap_arr(mPar->nside, mPar->nbTypes);
  interpolator_t *inter  = initialize_interpolator_t(N_MAX_1);
  sampler_arr *nOfZArr   = initialize_sampler_arr(mPar->N_z_map+1, mPar->nbTypes); //-- +1 is necessary
  gal_list_mat *gListMat = initialize_gal_list_mat(mPar->N_z_map, mPar->nbTypes);
  HPMap_arr *bufferMap   = initialize_HPMap_arr(mPar->nside, 6);
  
  setRatioMaps(mPar, rMapArr);
  setNOfZArr(mPar, inter, nOfZArr);
  sampleGalaxies(mPar, NULL, rMapArr, nOfZArr, NULL, NULL, gListMat, bufferMap->array[0], N_MAX_2);
  
  if (mPar->skipLensing == 1) {
    if (mPar->verbose < 3) printf("Skipped lensing module\n");
  }
  else {
    interpolateLensing(mPar, gListMat, bufferMap->array[0], bufferMap->array[1], bufferMap->array[2], bufferMap->array[3], bufferMap->array[4], bufferMap->array[5]);
    makeG(mPar, gListMat);
    addNoise(mPar, gListMat);
    flipSign(mPar, gListMat);
  }
  
  outFitsGalListMat(mPar, gListMat, mPar->verbose<3);
  
  free_HPMap_arr(rMapArr);
  free_interpolator_t(inter);
  free_sampler_arr(nOfZArr);
  free_gal_list_mat(gListMat);
  free_HPMap_arr(bufferMap);
  printf("------------------------------------------------------------------------\n");
  return;
}

void processMock_ratioMap_CL(MFP_param *mPar)
{
  HPMap_arr *rMapArr     = initialize_HPMap_arr(mPar->nside, mPar->nbTypes);
  interpolator_t *inter  = initialize_interpolator_t(N_MAX_1);
  sampler_arr *nOfZArr   = initialize_sampler_arr(mPar->N_z_map+1, mPar->nbTypes); //-- +1 is necessary
  gal_list_mat *gListMat = initialize_gal_list_mat(mPar->N_z_map, mPar->nbTypes);
  HPMap_arr *bufferMap   = initialize_HPMap_arr(mPar->nside, 3);
  
  setRatioMaps(mPar, rMapArr);
  setNOfZArr(mPar, inter, nOfZArr);
  sampleGalaxies(mPar, NULL, rMapArr, nOfZArr, NULL, NULL, gListMat, bufferMap->array[0], N_MAX_2);
  
  if (mPar->skipLensing == 1) {
    if (mPar->verbose < 3) printf("Skipped lensing module\n");
  }
  else {
    assignLensing(mPar, gListMat, bufferMap->array[0], bufferMap->array[1], bufferMap->array[2]);
    makeG(mPar, gListMat);
    addNoise(mPar, gListMat);
    flipSign(mPar, gListMat);
  }
  
  outFitsGalListMat(mPar, gListMat, mPar->verbose<3);
  
  free_HPMap_arr(rMapArr);
  free_interpolator_t(inter);
  free_sampler_arr(nOfZArr);
  free_gal_list_mat(gListMat);
  free_HPMap_arr(bufferMap);
  printf("------------------------------------------------------------------------\n");
  return;
}

void processMock_depthMap_LOS(MFP_param *mPar)
{
  HPMap_arr *rMapArr      = initialize_HPMap_arr(mPar->nside, mPar->nbTypes);
  interpolator_t *inter   = initialize_interpolator_t(N_MAX_1);
  sampler_arr *nOfZArr    = initialize_sampler_arr(mPar->N_z_map+1, mPar->nbTypes); //-- +1 is necessary
  HPMap_arr *vdMapArr     = initialize_HPMap_arr(mPar->nside, mPar->nbDepthMaps);
  sampler_arr *VD_nOfZArr = initialize_sampler_arr(mPar->N_z_map+1, mPar->VD_nbNOfZ); //-- +1 is necessary
  gal_list_mat *gListMat  = initialize_gal_list_mat(mPar->N_z_map, mPar->totNbTypes);
  HPMap_arr *bufferMap    = initialize_HPMap_arr(mPar->nside, 6);
  
  setRatioMaps(mPar, rMapArr);
  setNOfZArr(mPar, inter, nOfZArr);
  setDepthMaps(mPar, vdMapArr);
  setVDNOfZArr(mPar, inter, VD_nOfZArr);
  sampleGalaxies(mPar, NULL, rMapArr, nOfZArr, vdMapArr, VD_nOfZArr, gListMat, bufferMap->array[0], N_MAX_2);
  
  if (mPar->skipLensing == 1) {
    if (mPar->verbose < 3) printf("Skipped lensing module\n");
  }
  else {
    interpolateLensing(mPar, gListMat, bufferMap->array[0], bufferMap->array[1], bufferMap->array[2], bufferMap->array[3], bufferMap->array[4], bufferMap->array[5]);
    makeG(mPar, gListMat);
    addNoise(mPar, gListMat);
    flipSign(mPar, gListMat);
  }
  
  outFitsGalListMat(mPar, gListMat, mPar->verbose<3);
  
  free_HPMap_arr(rMapArr);
  free_interpolator_t(inter);
  free_sampler_arr(nOfZArr);
  free_HPMap_arr(vdMapArr);
  free_sampler_arr(VD_nOfZArr);
  free_gal_list_mat(gListMat);
  free_HPMap_arr(bufferMap);
  printf("------------------------------------------------------------------------\n");
  return;
}

void processMock_ratioMap_projCL(MFP_param *mPar, int resol, int nbSplits)
{
  HPMap_arr *rMapArr     = initialize_HPMap_arr(mPar->nside, mPar->nbTypes);
  interpolator_t *inter  = initialize_interpolator_t(N_MAX_1);
  sampler_arr *nOfZArr   = initialize_sampler_arr(mPar->N_z_map+1, mPar->nbTypes); //-- +1 is necessary
  gal_list_mat *gListMat = initialize_gal_list_mat(mPar->nbTypes, mPar->nbTypes);
  HPMap_arr *bufferMap   = initialize_HPMap_arr(mPar->nside, 3);
  
  setRatioMaps(mPar, rMapArr);
  setNOfZArr(mPar, inter, nOfZArr);
  sampleGalaxies_projCL(mPar, NULL, rMapArr, nOfZArr, NULL, NULL, gListMat, bufferMap->array[0], N_MAX_2, resol, nbSplits);
  
  if (mPar->skipLensing == 1) {
    if (mPar->verbose < 3) printf("Skipped lensing module\n");
  }
  else {
    assignLensing_projCL(mPar, gListMat, bufferMap->array[0], bufferMap->array[1], bufferMap->array[2], nbSplits);
    makeG(mPar, gListMat);
    addNoise(mPar, gListMat);
    flipSign(mPar, gListMat);
  }
  
  outFitsGalListMat(mPar, gListMat, mPar->verbose<3);
  
  free_HPMap_arr(rMapArr);
  free_interpolator_t(inter);
  free_sampler_arr(nOfZArr);
  free_gal_list_mat(gListMat);
  free_HPMap_arr(bufferMap);
  printf("------------------------------------------------------------------------\n");
  return;
}
#undef N_MAX_1
#undef N_MAX_2

//----------------------------------------------------------------------


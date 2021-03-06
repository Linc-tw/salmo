

  //------------------------------------------------------//
  //--  sampling.h                                      --//
  //--  Version 2020.07.09                              --//
  //--                                                  --//
  //--  Copyright (C) 2020 - Chieh-An Lin               --//
  //--  GNU GPLv3 - https://www.gnu.org/licenses/       --//
  //------------------------------------------------------//


#include "commonHeader.h"

#ifndef __SALMO_SAMPLING__
#define __SALMO_SAMPLING__

#ifdef __SALMO_USE_HEALPIX_CXX__
#include "wrapper.h"
#endif

#include "FITSFunctions.h"
#include "HEALPixFunctions.h"
#include "parameters.h"


typedef struct {
  long long nsidePix; //-- nside of the map
  long long nbPix;    //-- Number of pixels
  long long *map;     //-- Map, allow maximum 64
} type_map;

typedef struct {
  int length;
  HPMap_t **array;
} HPMap_arr;

typedef struct {
  double RA;        //-- [rad]
  double DEC;       //-- [rad]
  double z;         //-- [-]
  double p;         //-- [-] For interpolation
  long long pix;    //-- [-] Pixel index
  double sigma_eps; //-- [-] sigma_eps for variable depth
  double kappa;     //-- [-] Convergence
  double gamma_1;   //-- [-] Shear or reduced shear
  double gamma_2;   //-- [-] Shear or reduced shear
  double e_1;       //-- [-] Noisy reduced shear
  double e_2;       //-- [-] Noisy reduced shear
} gal_t;

typedef struct gal_node {
  gal_t *g;
  struct gal_node *next;
} gal_node;

typedef struct {
  int length;      //-- Number of nodes allocated
  int size;        //-- Number of nodes containing information
  gal_node *first; //-- Begin of the list
  gal_node *last;  //-- End of the list
} gal_list;

typedef struct {
  int N1, N2, length; //-- Lengths
  int total;          //-- Number of galaxies
  gal_list **matrix;  //-- Array of galaxy lists
} gal_list_mat;


//-- Functions related to initialization
type_map *initialize_type_map(long long nsidePix);
void free_type_map(type_map *tMap);
void reset_type_map(type_map *tMap);
gal_node *initialize_gal_node();
gal_list *initialize_gal_list();
void free_gal_list(gal_list *gList);
gal_list_mat *initialize_gal_list_mat(int N_z_map, int nbTypes);
void free_gal_list_mat(gal_list_mat *gListMat);

//-- Functions related to sampling
void setTypeMap(Salmo_param *sPar, type_map *tMap, HPMap_t *mask);
void setRatioMaps(Salmo_param *sPar, HPMap_arr *rMapArr);
void readAsciiNOfZ(char name[], interpolator_t *inter, int verbose);
void setNOfZArr(Salmo_param *sPar, interpolator_t *inter, sampler_arr *nOfZArr);
void setDepthMaps(Salmo_param *sPar, HPMap_arr *vdMapArr);
void setVDNOfZArr(Salmo_param *sPar, interpolator_t *inter, sampler_arr *VD_nOfZArr);
void set_gal_t(gal_t *g, double RA, double DEC, double z, double p, long long pix, double sigma_eps);
void append_gal_list(gal_list *gList, double RA, double DEC, double z, double p, long long pix, double sigma_eps);
void samplePos(gsl_rng *generator, long long nside, long long pix, int N_max, int N, long long *lenArr, long long *cumLenArr, long long v1, long long v2, double *RAArr, double *DECArr, 
	       double z1, double z2, gal_list *gList, double sigma_eps);
void samplePos_projCL(gsl_rng *generator, long long nside, long long pix, int N_max, int N, long long *lenArr, long long *cumLenArr, long long v1, long long v2, double *RAArr, double *DECArr, 
		      sampler_t *nOfZ, gal_list *gList, double sigma_eps);
void sampleGalaxies(Salmo_param *sPar, type_map *tMap, HPMap_arr *rMapArr, sampler_arr *nOfZArr, HPMap_arr *vdMapArr, sampler_arr *VD_nOfZArr, gal_list_mat *gListMat, HPMap_t *delta, int N_max);
void smoothDeltaMap(Salmo_param *sPar, HPMap_t *delta, int resol);
void sampleGalaxies_projCL(Salmo_param *sPar, type_map *tMap, HPMap_arr *rMapArr, sampler_arr *nOfZArr, HPMap_arr *vdMapArr, sampler_arr *VD_nOfZArr, gal_list_mat *gListMat, 
			   HPMap_t *delta, int N_max, int resol, int nbSplits);

//-- Functions related to lensing
void readLensingMaps(char name[], HPMap_t *kappa, HPMap_t *gamma_1, HPMap_t *gamma_2, int verbose);
double interpolateBetweenPixels(HPMap_t *full, long long neighbor[4], double weight[4]);
void interpolateLensing(Salmo_param *sPar, gal_list_mat *gListMat, HPMap_t *k_lower, HPMap_t *g1_lower, HPMap_t *g2_lower, HPMap_t *k_upper, HPMap_t *g1_upper, HPMap_t *g2_upper);
void assignLensing(Salmo_param *sPar, gal_list_mat *gListMat, HPMap_t *kappa, HPMap_t *gamma_1, HPMap_t *gamma_2);
void assignLensing_projCL(Salmo_param *sPar, gal_list_mat *gListMat, HPMap_t *kappa, HPMap_t *gamma_1, HPMap_t *gamma_2, int nbSplits);
void makeG(Salmo_param *sPar, gal_list_mat *gListMat);
void addNoise(Salmo_param *sPar, gal_list_mat *gListMat);
void flipSign(Salmo_param *sPar, gal_list_mat *gListMat);

//-- Functions related to outputs
void outFits_gal_t(FITS_t *fits, gal_t *g, double factor, int doNoise, int doWgt);
void outFitsParam(FITS_t *fits, Salmo_param *sPar, int k2);
void outFitsGalListMat(Salmo_param *sPar, gal_list_mat *gListMat, int verbose);

//-- Functions related to lensing maps
double E_sq_of_z(double z, double Omega_m);
void fillWInterpolator(interpolator_t *w_inter, double dz, double Omega_m);
double comovDist(interpolator_t *w_inter, double z);
void makeKappaMap(Salmo_param *sPar, HPMap_t *delta, HPMap_t *kappa, interpolator_t *w_inter, double Omega_m, int zInd);
void kappaToGamma(Salmo_param *sPar, HPMap_t *kappa, HPMap_t *gamma1, HPMap_t *gamma2, double_mat *kAlm, double_mat *gAlm, double_arr *weight, int l_maxx);
void flipSignForFlask(HPMap_t *gamma1, HPMap_t *gamma2);
void outFitsLensingMaps(Salmo_param *sPar, HPMap_t *kappa, HPMap_t *gamma1, HPMap_t *gamma2, int zInd, int verbose);

//-- Main function
void processLensingMaps(Salmo_param *sPar);
void processMock_typeMap_LOS(Salmo_param *sPar);
void processMock_ratioMap_LOS(Salmo_param *sPar);
void processMock_ratioMap_CL(Salmo_param *sPar);
void processMock_depthMap_LOS(Salmo_param *sPar);
void processMock_ratioMap_projCL(Salmo_param *sPar, int resol, int nbSplits);

#endif


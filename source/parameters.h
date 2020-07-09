

  //------------------------------------------------------//
  //--  parameters.h                                    --//
  //--  Version 2020.07.09                              --//
  //--                                                  --//
  //--  Copyright (C) 2020 - Chieh-An Lin               --//
  //--  GNU GPLv3 - https://www.gnu.org/licenses/       --//
  //------------------------------------------------------//


#include "commonHeader.h"

#ifndef __SALMO_PARAMETERS__
#define __SALMO_PARAMETERS__

#include "FITSFunctions.h"
#include "HEALPixFunctions.h"


//-- Parameters
typedef struct {
  //-- When adding new keys:
  //-- - Check salmoParam.par
  //-- - Check parameters.h
  //-- - Check find Key
  //-- - Check print
  //-- - Check set
  //-- - Check free
  
  //----------------------------------------------------------------------
  //-- Customized part (see salmoParam.par for documentations)
  
  //-- 1. Generality
  char parPath[STRING_LENGTH_MAX];
  char seed[64];
  int verbose;
  
  //-- 2. Input maps
  long long nside;
  int N_z_map;
  double *bin_z_map;
  char denPrefix[STRING_LENGTH_MAX];
  char lenPrefix[STRING_LENGTH_MAX];
  char runTag[STRING_LENGTH_MAX];
  
  //-- 3. Selection functions
  int nbTypes;
  char **maskPath;
  char **nOfZPath;
  double *n_gal;
  
  //-- 4. Lensing & outputs
  int doNoise;
  int doWgt;
  int signConv;
  double *sigma_eps;
  int *doLensing;
  char outPrefix[STRING_LENGTH_MAX];
  int outStyle;
  
  //-- 5. Variable depth
  int doVariableDepth;
  int nbDepthMaps;
  char **depthMapPath;
  int N_depth;
  double *bin_depth;
  int nbTomo;
  double *a_n_gal;
  double *b_n_gal;
  double *a_sigma_eps;
  double *b_sigma_eps;
  char **VD_nOfZPath;
  
  //----------------------------------------------------------------------
  //-- Precomputed part
  
  //-- Generality
  gsl_rng *generator;		//-- [-] Random number generator
  long long nbPix;		//-- [int] Number of pixels
  double A_pix;			//-- [float | rad^2] Area of pixels
  double *z_map;		//-- [float array] Centers of map redshift slices
  double *half_dz_map;		//-- [float array] Half width of map redshift slices
  double zMapRange[3];		//-- [3 float] z_map_min, z_map_max, dz_map
  int skipLensing;		//-- [int] 0 = any doLensing, 1 = skip
  
  int VD_nbNOfZ;
  int VD_nbTypes;
  int totNbTypes;
  
  //----------------------------------------------------------------------
  //-- Running part
  
  int MPISize;			//-- [int] Number of MPI processors
  int MPIInd;			//-- [int] Index for MPI processors
  
  //----------------------------------------------------------------------
} Salmo_param;


//-- Functions related to file reading
void ignoreComments(char line[]);
int getKeyAndValues(char line[], char kv[][STRING_LENGTH_MAX]);
int *makeIntArray(char kv[][STRING_LENGTH_MAX], int count);
double *makeDoubleArray(char kv[][STRING_LENGTH_MAX], int count);
void setPathWhichCanBeBlank(char path[], char kv[][STRING_LENGTH_MAX], int count);

//-- Functions related to Salmo_param
Salmo_param *initialize_Salmo_param();
void free_Salmo_param(Salmo_param *sPar);
int findParameterKey(Salmo_param *sPar, char kv[][STRING_LENGTH_MAX], int count);
void readParameters(char name[], Salmo_param *sPar, int verbose);
int updateFromCommandLine(int argc, char *argv[], Salmo_param *sPar);
void setParameters(Salmo_param *sPar);

//-- Functions related to printing
void printIntArray(int *iArr, int length);
void printDoubleArray(double *lfArr, int length, double factor, int digit);
void printParam(Salmo_param *sPar);
void printCompleteParam(Salmo_param *sPar);

#endif


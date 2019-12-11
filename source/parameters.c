

  //------------------------------------------------------//
  //--  parameters.c					--//
  //--  Version 2019.10.02				--//
  //--  						--//
  //--  Copyright (C) 2018 - Chieh-An Lin		--//
  //--  GNU GPLv3 - https://www.gnu.org/licenses/	--//
  //------------------------------------------------------//


#include "parameters.h"


//----------------------------------------------------------------------
//-- Functions related to file reading

void ignoreComments(char line[])
{
  char *end = strchr(line, '#');
  if (end == line) sprintf(end, "%s", "");
  else if (end != NULL) sprintf(end, "\n");
  return;
}

int getKeyAndValues(char line[], char kv[][STRING_LENGTH_MAX])
{
  if (!strcmp(line, "-h") || !strcmp(line, "-H") || !strcmp(line, "--help")) return -1; //-- This is for parameters updated from the command line.
  
  char buffer[STRING_LENGTH_MAX];
  sprintf(buffer, "%s", line);
  
  int count = 0;
  char *token = strtok(buffer, " ,=\"\t\n");
  while (token != NULL) {
    sprintf(kv[count], "%s", token);
    count++;
    token = strtok(NULL, " ,=\"\t\n");
  }
  return count;
}

int *makeIntArray(char kv[][STRING_LENGTH_MAX], int count)
{
  int *array = (count == 1) ? NULL : (int*)malloc((count-1) * sizeof(int));
  int i;
  for (i=0; i<count-1; i++) array[i] = atoi(kv[i+1]);
  return array;
}

double *makeDoubleArray(char kv[][STRING_LENGTH_MAX], int count)
{
  double *array = (count == 1) ? NULL : (double*)malloc((count-1) * sizeof(double));
  int i;
  for (i=0; i<count-1; i++) array[i] = atof(kv[i+1]);
  return array;
}

void setPathWhichCanBeBlank(char path[], char kv[][STRING_LENGTH_MAX], int count)
{
  if (count == 1) sprintf(path, "%s", "");
  else            sprintf(path, "%s", kv[1]);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to MFP_param

MFP_param *initialize_MFP_param()
{
  MFP_param *mPar = (MFP_param*)malloc(sizeof(MFP_param));
  //-- Essential
  mPar->N_depth    = 0;
  mPar->nbTomo     = 0;
  return mPar;
}

void free_MFP_param(MFP_param *mPar)
{
  int i;
  if (mPar) {
    if (mPar->bin_z_map)                  {free(mPar->bin_z_map);         mPar->bin_z_map       = NULL;}
    if (mPar->maskPath) {
      for (i=0; i<mPar->nbTypes; i++)     {free(mPar->maskPath[i]);       mPar->maskPath[i]     = NULL;}
                                           free(mPar->maskPath);          mPar->maskPath        = NULL;
    }
    if (mPar->nOfZPath) {
      for (i=0; i<mPar->nbTypes; i++)     {free(mPar->nOfZPath[i]);       mPar->nOfZPath[i]     = NULL;}
                                           free(mPar->nOfZPath);          mPar->nOfZPath        = NULL;
    }
    if (mPar->n_gal)                      {free(mPar->n_gal);             mPar->n_gal           = NULL;}
    if (mPar->doLensing)                  {free(mPar->doLensing);         mPar->doLensing       = NULL;}
    if (mPar->depthMapPath) {
      for (i=0; i<mPar->nbDepthMaps; i++) {free(mPar->depthMapPath[i]);   mPar->depthMapPath[i] = NULL;}
                                           free(mPar->depthMapPath);      mPar->depthMapPath    = NULL;
    }
    if (mPar->bin_depth)                  {free(mPar->bin_depth);         mPar->bin_depth       = NULL;}
    if (mPar->a_n_gal)                    {free(mPar->a_n_gal);           mPar->a_n_gal         = NULL;}
    if (mPar->b_n_gal)                    {free(mPar->b_n_gal);           mPar->b_n_gal         = NULL;}
    if (mPar->a_sigma_eps)                {free(mPar->a_sigma_eps);       mPar->a_sigma_eps     = NULL;}
    if (mPar->b_sigma_eps)                {free(mPar->b_sigma_eps);       mPar->b_sigma_eps     = NULL;}
    if (mPar->VD_nOfZPath) {
      for (i=0; i<mPar->VD_nbNOfZ; i++)   {free(mPar->VD_nOfZPath[i]);    mPar->VD_nOfZPath[i]  = NULL;}
                                           free(mPar->VD_nOfZPath);       mPar->VD_nOfZPath     = NULL;
    }
    
    if (mPar->generator)                  {gsl_rng_free(mPar->generator); mPar->generator       = NULL;}
    if (mPar->z_map)                      {free(mPar->z_map);             mPar->z_map           = NULL;}
    if (mPar->half_dz_map)                {free(mPar->half_dz_map);       mPar->half_dz_map     = NULL;}
    free(mPar); mPar = NULL;
  }
  return;
}

int findParameterKey(MFP_param *mPar, char kv[][STRING_LENGTH_MAX], int count)
{
  int i;
  
  if (count == 0) return 0; //-- unknown = 0
  
  //-- 1. Generality
  if (!strcmp(kv[0], "seed"))                 setPathWhichCanBeBlank(mPar->seed, kv, count);
  else if (!strcmp(kv[0], "verbose"))         mPar->verbose = atoi(kv[1]);
  
  //-- 2. Input maps
  else if (!strcmp(kv[0], "nside"))           mPar->nside     = strtol(kv[1], NULL, 10);
  else if (!strcmp(kv[0], "N_z_map"))         mPar->N_z_map   = atoi(kv[1]);
  else if (!strcmp(kv[0], "bin_z_map")) {
                                              if (mPar->bin_z_map) free(mPar->bin_z_map);
                                              mPar->bin_z_map = makeDoubleArray(kv, count);
  }
  else if (!strcmp(kv[0], "denPrefix"))       setPathWhichCanBeBlank(mPar->denPrefix, kv, count);
  else if (!strcmp(kv[0], "lenPrefix"))       setPathWhichCanBeBlank(mPar->lenPrefix, kv, count);
  else if (!strcmp(kv[0], "runTag"))          setPathWhichCanBeBlank(mPar->runTag, kv, count);
  
  //-- 3. Selection functions
  else if (!strcmp(kv[0], "nbTypes")) {
                                              if (mPar->maskPath) {
                                                for (i=0; i<mPar->nbTypes; i++) {free(mPar->maskPath[i]); mPar->maskPath[i] = NULL;}
                                                free(mPar->maskPath); mPar->maskPath = NULL;
                                              }
                                              if (mPar->nOfZPath) {
                                                for (i=0; i<mPar->nbTypes; i++) {free(mPar->nOfZPath[i]); mPar->nOfZPath[i] = NULL;}
                                                free(mPar->nOfZPath); mPar->nOfZPath = NULL;
                                              }
                                              mPar->nbTypes  = atoi(kv[1]);
                                              mPar->maskPath = (char**)malloc(mPar->nbTypes * sizeof(char*));
                                              mPar->nOfZPath = (char**)malloc(mPar->nbTypes * sizeof(char*));
                                              for (i=0; i<mPar->nbTypes; i++) {
                                                mPar->maskPath[i] = (char*)malloc(STRING_LENGTH_MAX * sizeof(char));
                                                mPar->nOfZPath[i] = (char*)malloc(STRING_LENGTH_MAX * sizeof(char));
                                              }
  }
  else if (!strcmp(kv[0], "maskPath"))        sprintf(mPar->maskPath[atoi(kv[1])], "%s", kv[2]);
  else if (!strcmp(kv[0], "nOfZPath"))        sprintf(mPar->nOfZPath[atoi(kv[1])], "%s", kv[2]);
  else if (!strcmp(kv[0], "n_gal")) {
                                              if (mPar->n_gal) free(mPar->n_gal);
                                              mPar->n_gal = makeDoubleArray(kv, count);
                                              for (i=0; i<mPar->nbTypes; i++) mPar->n_gal[i] /= ARCMIN_SQ_TO_RADIAN_SQ;
  }
  
  //-- 4. Lensing & outputs
  else if (!strcmp(kv[0], "doNoise"))         mPar->doNoise   = atoi(kv[1]);
  else if (!strcmp(kv[0], "doWgt"))           mPar->doWgt     = atoi(kv[1]);
  else if (!strcmp(kv[0], "signConv"))        mPar->signConv  = atoi(kv[1]);
  else if (!strcmp(kv[0], "sigma_eps"))       mPar->sigma_eps = atof(kv[1]);
  else if (!strcmp(kv[0], "doLensing")) {
                                              if (mPar->doLensing) free(mPar->doLensing);
                                              mPar->doLensing = makeIntArray(kv, count);
  }
  else if (!strcmp(kv[0], "outPrefix"))       setPathWhichCanBeBlank(mPar->outPrefix, kv, count);
  else if (!strcmp(kv[0], "outStyle"))        mPar->outStyle = atoi(kv[1]);
  
  //-- 5. Variable depth
  else if (!strcmp(kv[0], "doVariableDepth")) mPar->doVariableDepth = atoi(kv[1]);
  else if (!strcmp(kv[0], "nbDepthMaps")) {
                                              if (mPar->depthMapPath) {
                                                for (i=0; i<mPar->nbDepthMaps; i++) {free(mPar->depthMapPath[i]); mPar->depthMapPath[i] = NULL;}
                                                free(mPar->depthMapPath); mPar->depthMapPath = NULL;
                                              }
                                              mPar->nbDepthMaps  = atoi(kv[1]);
                                              mPar->depthMapPath = (char**)malloc(mPar->nbDepthMaps * sizeof(char*));
                                              for (i=0; i<mPar->nbDepthMaps; i++) mPar->depthMapPath[i] = (char*)malloc(STRING_LENGTH_MAX * sizeof(char));
  }
  else if (!strcmp(kv[0], "depthMapPath"))    sprintf(mPar->depthMapPath[atoi(kv[1])], "%s", kv[2]);
  else if (!strcmp(kv[0], "N_depth")) {
                                              if (mPar->VD_nOfZPath) {
                                                for (i=0; i<mPar->VD_nbNOfZ; i++) {free(mPar->VD_nOfZPath[i]); mPar->VD_nOfZPath[i] = NULL;}
                                                free(mPar->VD_nOfZPath); mPar->VD_nOfZPath = NULL;
                                              }
                                              mPar->N_depth   = atoi(kv[1]);
                                              mPar->VD_nbNOfZ = mPar->N_depth * mPar->nbTomo;
                                              if (mPar->VD_nbNOfZ > 0) mPar->VD_nOfZPath = (char**)malloc(mPar->VD_nbNOfZ * sizeof(char*));
                                              for (i=0; i<mPar->VD_nbNOfZ; i++) mPar->VD_nOfZPath[i] = (char*)malloc(STRING_LENGTH_MAX * sizeof(char));
  }
  else if (!strcmp(kv[0], "bin_depth")) {
                                              if (mPar->bin_depth) free(mPar->bin_depth);
                                              mPar->bin_depth = makeDoubleArray(kv, count);
  }
  else if (!strcmp(kv[0], "nbTomo")) {
                                              if (mPar->VD_nOfZPath) {
                                                for (i=0; i<mPar->VD_nbNOfZ; i++) {free(mPar->VD_nOfZPath[i]); mPar->VD_nOfZPath[i] = NULL;}
                                                free(mPar->VD_nOfZPath); mPar->VD_nOfZPath = NULL;
                                              }
                                              mPar->nbTomo    = atoi(kv[1]);
                                              mPar->VD_nbNOfZ = mPar->N_depth * mPar->nbTomo;
                                              if (mPar->VD_nbNOfZ > 0) mPar->VD_nOfZPath = (char**)malloc(mPar->VD_nbNOfZ * sizeof(char*));
                                              for (i=0; i<mPar->VD_nbNOfZ; i++) mPar->VD_nOfZPath[i] = (char*)malloc(STRING_LENGTH_MAX * sizeof(char));
  }
  else if (!strcmp(kv[0], "a_n_gal")) {
                                              if (mPar->a_n_gal) free(mPar->a_n_gal);
                                              mPar->a_n_gal = makeDoubleArray(kv, count);
  }
  else if (!strcmp(kv[0], "b_n_gal")) {
                                              if (mPar->b_n_gal) free(mPar->b_n_gal);
                                              mPar->b_n_gal = makeDoubleArray(kv, count);
  }
  else if (!strcmp(kv[0], "a_sigma_eps")) {
                                              if (mPar->a_sigma_eps) free(mPar->a_sigma_eps);
                                              mPar->a_sigma_eps = makeDoubleArray(kv, count);
  }
  else if (!strcmp(kv[0], "b_sigma_eps")) {
                                              if (mPar->b_sigma_eps) free(mPar->b_sigma_eps);
                                              mPar->b_sigma_eps = makeDoubleArray(kv, count);
  }
  else if (!strcmp(kv[0], "VD_nOfZPath"))     sprintf(mPar->VD_nOfZPath[atoi(kv[1])], "%s", kv[2]);
  
  else return 1; //-- unknown = 1
  
  return 0;
}

void readParameters(char name[], MFP_param *mPar, int verbose)
{
  sprintf(mPar->parPath, "%s", name);
  FILE *file = fopen(name, "r");
  int verbose2 = (mPar->MPIInd == 0) && (mPar->verbose < 99) && verbose;
  
  char line[STRING_LENGTH_MAX], kv[STRING_LENGTH_MAX][STRING_LENGTH_MAX], *buffer;
  int count, unknown;
  
  //-- Read
  buffer = fgets(line, STRING_LENGTH_MAX, file);
  while (buffer != NULL) {
    ignoreComments(line);
    count   = getKeyAndValues(line, kv);
    unknown = findParameterKey(mPar, kv, count);
    buffer  = fgets(line, STRING_LENGTH_MAX, file);
    
    if (verbose == 0 && !strcmp(kv[0], "verbose")) {
      fclose(file);
      return;
    }
    if (unknown == 1 && verbose2) printf("Detected the keyword \"%s\" but unable to interpret, skipped\n", kv[0]);
  }
  
  fclose(file);
  return;
}

int updateFromCommandLine(int argc, char *argv[], MFP_param *mPar)
{
  //-- Update from the command line
  
  int verbose = (mPar->MPIInd == 0) && (mPar->verbose < 99);
  char kv[STRING_LENGTH_MAX][STRING_LENGTH_MAX];
  int i, count, unknown;
  
  for (i=3; i<argc; i++) {
    count = getKeyAndValues(argv[i], kv);
    if (count == -1) break; //-- Help detected
    unknown = findParameterKey(mPar, kv, count); 
    if (unknown == 1 && count > 1 && verbose) printf("Detected the keyword \"%s\" but unable to interpret, skipped\n", kv[0]);
  }
  
  return (int)(count == -1);
}

void setParameters(MFP_param *mPar)
{
  int verbose = mPar->MPIInd == 0 && mPar->verbose < 99;
  u_int32_t seed;
  int i;
  
  //-- Precomputed part - Generality
  if (strchr(mPar->seed, 'r') == NULL) seed = strtoul(mPar->seed, NULL, 10);
  else {
    seed = renewSeed();
    sprintf(mPar->seed, "%u (random)", seed);
  }
  mPar->generator = initializeGenerator(seed);
  mPar->nbPix     = 12 * mPar->nside * mPar->nside;
  mPar->A_pix     = FULL_SKY * DEGREE_SQ_TO_RADIAN_SQ / mPar->nbPix;
  if (mPar->N_z_map < 0) {
    mPar->zMapRange[0] = mPar->bin_z_map[0];
    mPar->zMapRange[1] = mPar->bin_z_map[1];
    mPar->zMapRange[2] = mPar->bin_z_map[2];
    mPar->N_z_map      = (int)round((mPar->zMapRange[1] - mPar->zMapRange[0]) / mPar->zMapRange[2]);
    if (mPar->bin_z_map) free(mPar->bin_z_map);
    mPar->bin_z_map    = (double*)malloc((mPar->N_z_map+1) * sizeof(double));
    for (i=0; i<=mPar->N_z_map; i++) mPar->bin_z_map[i] = mPar->zMapRange[0] + i * mPar->zMapRange[2];
  }
  else {
    mPar->zMapRange[0] = mPar->bin_z_map[0];
    mPar->zMapRange[1] = mPar->bin_z_map[mPar->N_z_map];
    mPar->zMapRange[2] = 0.0; //-- dz_map is not defined.
  }
  if (mPar->z_map)       free(mPar->z_map);
  if (mPar->half_dz_map) free(mPar->half_dz_map);
  mPar->z_map       = (double*)malloc(mPar->N_z_map * sizeof(double));
  mPar->half_dz_map = (double*)malloc(mPar->N_z_map * sizeof(double));
  for (i=0; i<mPar->N_z_map; i++) {
    mPar->z_map[i]       = 0.5 * (mPar->bin_z_map[i] + mPar->bin_z_map[i+1]);
    mPar->half_dz_map[i] = 0.5 * (mPar->bin_z_map[i+1] - mPar->bin_z_map[i]);
  }
  mPar->skipLensing = 1;
  for (i=0; i<mPar->nbTypes; i++) {
    if (mPar->doLensing[i] == 1) {
      mPar->skipLensing = 0;
      break;
    }
  }
  
  if (mPar->doVariableDepth == 0) {
    mPar->nbDepthMaps = 0;
    mPar->N_depth     = 0;
    mPar->nbTomo      = 0;
    mPar->VD_nbNOfZ   = 0;
    mPar->VD_nbTypes  = 0;
  }
  else {
    mPar->skipLensing = 0;
//     mPar->VD_nbNOfZ; //-- Already initialized
    mPar->VD_nbTypes = mPar->nbDepthMaps * mPar->nbTomo;
  }
  mPar->totNbTypes = mPar->nbTypes + mPar->VD_nbTypes;
  
  
  //-- Running part
  //mPar->MPISize;  //-- Should not be initialized here
  //mPar->MPIInd;   //-- Should not be initialized here
  return;
}

//----------------------------------------------------------------------
//-- Functions related to printing

void printIntArray(int *iArr, int length)
{
  if (length == 0) {
    printf("NULL");
    return;
  }
  
  int i;
  printf("%d", iArr[0]);
  for (i=1; i<length; i++) printf(", %d", iArr[i]);
  return;
}

void printDoubleArray(double *lfArr, int length, double factor, int digit)
{
  if (lfArr == NULL || length == 0) {
    printf("NULL");
    return;
  }
  
  int i;
  if (digit == 0) {
    printf("%f", lfArr[0]*factor);
    for (i=1; i<length; i++) printf(", %f", lfArr[i]*factor);
  }
  else if (digit == 1) {
    printf("%.1f", lfArr[0]*factor);
    for (i=1; i<length; i++) printf(", %.1f", lfArr[i]*factor);
  }
  else if (digit == 2) {
    printf("%.2f", lfArr[0]*factor);
    for (i=1; i<length; i++) printf(", %.2f", lfArr[i]*factor);
  }
  else if (digit == 3) {
    printf("%.3f", lfArr[0]*factor);
    for (i=1; i<length; i++) printf(", %.3f", lfArr[i]*factor);
  }
  else if (digit == 4) {
    printf("%.4f", lfArr[0]*factor);
    for (i=1; i<length; i++) printf(", %.4f", lfArr[i]*factor);
  }
  else if (digit == 5) {
    printf("%.5f", lfArr[0]*factor);
    for (i=1; i<length; i++) printf(", %.5f", lfArr[i]*factor);
  }
  return;
}

void printParam(MFP_param *mPar)
{
  int i;
  printf("seed         = %s\n", mPar->seed);
  printf("verbose      = %d\n", mPar->verbose);
  
  printf("\n");
  printf("nside        = %lld\n", mPar->nside);
  printf("N_z_map      = %d\n", mPar->N_z_map);
  printf("bin_z_map    = "); printDoubleArray(mPar->bin_z_map, mPar->N_z_map+1, 1.0, 2); printf(" [-]\n");
  printf("denPrefix    = \"%s\"\n", mPar->denPrefix);
  printf("lenPrefix    = \"%s\"\n", mPar->lenPrefix);
  printf("runTag       = \"%s\"\n", mPar->runTag);
  
  printf("\n");
  printf("nbTypes      = %d\n", mPar->nbTypes);
  for (i=0; i<mPar->nbTypes;i++) printf("maskPath[%2d] = \"%s\"\n", i, mPar->maskPath[i]);
  for (i=0; i<mPar->nbTypes;i++) printf("nOfZPath[%2d] = \"%s\"\n", i, mPar->nOfZPath[i]);
  printf("n_gal        = "); printDoubleArray(mPar->n_gal, mPar->nbTypes, 1.0/RADIAN_SQ_TO_ARCMIN_SQ, 2); printf(" [arcmin^-2]\n");
  
  printf("\n");
  printf("doNoise      = %d\n", mPar->doNoise);
  printf("doWgt        = %d\n", mPar->doWgt);
  printf("signConv     = %d\n", mPar->signConv);
  printf("sigma_eps    = %f\n", mPar->sigma_eps);
  printf("doLensing    = "); printIntArray(mPar->doLensing, mPar->nbTypes); printf("\n");
  printf("outPrefix    = \"%s\"\n", mPar->outPrefix);
  printf("outStyle     = %d\n", mPar->outStyle);
  
  printf("\n");
  printf("doVariableDepth  = %d\n", mPar->doVariableDepth);
  
  if (mPar->doVariableDepth == 1) {
    printf("nbDepthMaps      = %d\n", mPar->nbDepthMaps);
    for (i=0; i<mPar->nbDepthMaps;i++) printf("depthMapPath[%2d] = \"%s\"\n", i, mPar->depthMapPath[i]);
    printf("N_depth          = %d\n", mPar->N_depth);
    printf("bin_depth        = "); printDoubleArray(mPar->bin_depth, mPar->N_depth+1, 1.0, 2); printf(" [-]\n");
    printf("nbTomo           = %d\n", mPar->nbTomo);
    printf("a_n_gal          = "); printDoubleArray(mPar->a_n_gal, mPar->nbTomo, 1.0, 3); printf("\n");
    printf("b_n_gal          = "); printDoubleArray(mPar->b_n_gal, mPar->nbTomo, 1.0, 3); printf("\n");
    printf("a_sigma_eps      = "); printDoubleArray(mPar->a_sigma_eps, mPar->nbTomo, 1.0, 3); printf("\n");
    printf("b_sigma_eps      = "); printDoubleArray(mPar->b_sigma_eps, mPar->nbTomo, 1.0, 3); printf("\n");
    for (i=0; i<mPar->VD_nbNOfZ;i++) printf("VD_nOfZPath[%2d]  = \"%s\"\n", i, mPar->VD_nOfZPath[i]);
  }
  return;
}

void printCompleteParam(MFP_param *mPar)
{
  int i;
  printf("------ 1. Generality ------\n");
  printf("parPath      = \"%s\"\n", mPar->parPath);
  printf("seed         = %s\n", mPar->seed);
  printf("verbose      = %d\n", mPar->verbose);
  
  printf("------ 2. Input maps ------\n");
  printf("nside        = %lld\n", mPar->nside);
  printf("N_z_map      = %d\n", mPar->N_z_map);
  printf("bin_z_map    = "); printDoubleArray(mPar->bin_z_map, mPar->N_z_map+1, 1.0, 2); printf(" [-]\n");
  printf("denPrefix    = \"%s\"\n", mPar->denPrefix);
  printf("lenPrefix    = \"%s\"\n", mPar->lenPrefix);
  printf("runTag       = \"%s\"\n", mPar->runTag);
  
  printf("------ 3. Selection functions ------\n");
  printf("nbTypes      = %d\n", mPar->nbTypes);
  for (i=0; i<mPar->nbTypes;i++) printf("maskPath[%2d] = \"%s\"\n", i, mPar->maskPath[i]);
  for (i=0; i<mPar->nbTypes;i++) printf("nOfZPath[%2d] = \"%s\"\n", i, mPar->nOfZPath[i]);
  printf("n_gal        = "); printDoubleArray(mPar->n_gal, mPar->nbTypes, 1.0/RADIAN_SQ_TO_ARCMIN_SQ, 2); printf(" [arcmin^-2]\n");
  
  printf("------ 4. Lensing & outputs ------\n");
  printf("doNoise      = %d\n", mPar->doNoise);
  printf("doWgt        = %d\n", mPar->doWgt);
  printf("signConv     = %d\n", mPar->signConv);
  printf("sigma_eps    = %f\n", mPar->sigma_eps);
  printf("doLensing    = "); printIntArray(mPar->doLensing, mPar->nbTypes); printf("\n");
  printf("outPrefix    = \"%s\"\n", mPar->outPrefix);
  printf("outStyle     = %d\n", mPar->outStyle);
  
  printf("------ 5. Variable depth ------\n");
  printf("doVariableDepth  = %d\n", mPar->doVariableDepth);
  printf("nbDepthMaps      = %d\n", mPar->nbDepthMaps);
  for (i=0; i<mPar->nbDepthMaps;i++) printf("depthMapPath[%2d] = \"%s\"\n", i, mPar->depthMapPath[i]);
  printf("N_depth          = %d\n", mPar->N_depth);
  printf("bin_depth        = "); printDoubleArray(mPar->bin_depth, mPar->N_depth+1, 1.0, 2); printf(" [-]\n");
  printf("nbTomo           = %d\n", mPar->nbTomo);
  printf("a_n_gal          = "); printDoubleArray(mPar->a_n_gal, mPar->nbTomo, 1.0, 3); printf("\n");
  printf("b_n_gal          = "); printDoubleArray(mPar->b_n_gal, mPar->nbTomo, 1.0, 3); printf("\n");
  printf("a_sigma_eps      = "); printDoubleArray(mPar->a_sigma_eps, mPar->nbTomo, 1.0, 3); printf("\n");
  printf("b_sigma_eps      = "); printDoubleArray(mPar->b_sigma_eps, mPar->nbTomo, 1.0, 3); printf("\n");
  for (i=0; i<mPar->VD_nbNOfZ;i++) printf("VD_nOfZPath[%2d]  = \"%s\"\n", i, mPar->VD_nOfZPath[i]);
  
  printf("------ Precomputed part ------\n");
  printf("generator    = What do you expect from this?\n");
  printf("nbPix        = %lld\n", mPar->nbPix);
  printf("A_pix        = %f [arcmin^2]\n", mPar->A_pix*RADIAN_SQ_TO_ARCMIN_SQ);
  printf("z_map        = "); printDoubleArray(mPar->z_map, mPar->N_z_map, 1.0, 2); printf(" [-]\n");
  printf("half_dz_map  = "); printDoubleArray(mPar->half_dz_map, mPar->N_z_map, 1.0, 2); printf(" [-]\n");
  printf("zMapRange    = %.3f, %.3f, %.3f [-]\n", mPar->zMapRange[0], mPar->zMapRange[1], mPar->zMapRange[2]);
  printf("skipLensing  = %d\n", mPar->skipLensing);
  printf("VD_nbNOfZ    = %d\n", mPar->VD_nbNOfZ);
  printf("VD_nbTypes   = %d\n", mPar->VD_nbTypes);
  printf("totNbTypes   = %d\n", mPar->totNbTypes);
  return;
}

//----------------------------------------------------------------------


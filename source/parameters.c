

  //------------------------------------------------------//
  //--  parameters.c                                    --//
  //--  Version 2020.07.09                              --//
  //--                                                  --//
  //--  Copyright (C) 2020 - Chieh-An Lin               --//
  //--  GNU GPLv3 - https://www.gnu.org/licenses/       --//
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
//-- Functions related to Salmo_param

Salmo_param *initialize_Salmo_param()
{
  Salmo_param *sPar = (Salmo_param*)malloc(sizeof(Salmo_param));
  //-- Essential
  sPar->N_depth    = 0;
  sPar->nbTomo     = 0;
  return sPar;
}

void free_Salmo_param(Salmo_param *sPar)
{
  int i;
  if (sPar) {
    if (sPar->bin_z_map)                  {free(sPar->bin_z_map);         sPar->bin_z_map       = NULL;}
    if (sPar->maskPath) {
      for (i=0; i<sPar->nbTypes; i++)     {free(sPar->maskPath[i]);       sPar->maskPath[i]     = NULL;}
                                           free(sPar->maskPath);          sPar->maskPath        = NULL;
    }
    if (sPar->nOfZPath) {
      for (i=0; i<sPar->nbTypes; i++)     {free(sPar->nOfZPath[i]);       sPar->nOfZPath[i]     = NULL;}
                                           free(sPar->nOfZPath);          sPar->nOfZPath        = NULL;
    }
    if (sPar->n_gal)                      {free(sPar->n_gal);             sPar->n_gal           = NULL;}
    
    if (sPar->sigma_eps)                  {free(sPar->sigma_eps);         sPar->sigma_eps       = NULL;}
    if (sPar->doLensing)                  {free(sPar->doLensing);         sPar->doLensing       = NULL;}
    if (sPar->depthMapPath) {
      for (i=0; i<sPar->nbDepthMaps; i++) {free(sPar->depthMapPath[i]);   sPar->depthMapPath[i] = NULL;}
                                           free(sPar->depthMapPath);      sPar->depthMapPath    = NULL;
    }
    if (sPar->bin_depth)                  {free(sPar->bin_depth);         sPar->bin_depth       = NULL;}
    if (sPar->a_n_gal)                    {free(sPar->a_n_gal);           sPar->a_n_gal         = NULL;}
    if (sPar->b_n_gal)                    {free(sPar->b_n_gal);           sPar->b_n_gal         = NULL;}
    if (sPar->a_sigma_eps)                {free(sPar->a_sigma_eps);       sPar->a_sigma_eps     = NULL;}
    if (sPar->b_sigma_eps)                {free(sPar->b_sigma_eps);       sPar->b_sigma_eps     = NULL;}
    if (sPar->VD_nOfZPath) {
      for (i=0; i<sPar->VD_nbNOfZ; i++)   {free(sPar->VD_nOfZPath[i]);    sPar->VD_nOfZPath[i]  = NULL;}
                                           free(sPar->VD_nOfZPath);       sPar->VD_nOfZPath     = NULL;
    }
    
    if (sPar->generator)                  {gsl_rng_free(sPar->generator); sPar->generator       = NULL;}
    if (sPar->z_map)                      {free(sPar->z_map);             sPar->z_map           = NULL;}
    if (sPar->half_dz_map)                {free(sPar->half_dz_map);       sPar->half_dz_map     = NULL;}
    free(sPar); sPar = NULL;
  }
  return;
}

int findParameterKey(Salmo_param *sPar, char kv[][STRING_LENGTH_MAX], int count)
{
  int i;
  
  if (count == 0) return 0; //-- unknown = 0
  
  //-- 1. Generality
  if (!strcmp(kv[0], "seed"))                 setPathWhichCanBeBlank(sPar->seed, kv, count);
  else if (!strcmp(kv[0], "verbose"))         sPar->verbose = atoi(kv[1]);
  
  //-- 2. Input maps
  else if (!strcmp(kv[0], "nside"))           sPar->nside     = strtol(kv[1], NULL, 10);
  else if (!strcmp(kv[0], "N_z_map"))         sPar->N_z_map   = atoi(kv[1]);
  else if (!strcmp(kv[0], "bin_z_map")) {
                                              if (sPar->bin_z_map) free(sPar->bin_z_map);
                                              sPar->bin_z_map = makeDoubleArray(kv, count);
  }
  else if (!strcmp(kv[0], "denPrefix"))       setPathWhichCanBeBlank(sPar->denPrefix, kv, count);
  else if (!strcmp(kv[0], "lenPrefix"))       setPathWhichCanBeBlank(sPar->lenPrefix, kv, count);
  else if (!strcmp(kv[0], "runTag"))          setPathWhichCanBeBlank(sPar->runTag, kv, count);
  
  //-- 3. Selection functions
  else if (!strcmp(kv[0], "nbTypes")) {
                                              if (sPar->maskPath) {
                                                for (i=0; i<sPar->nbTypes; i++) {free(sPar->maskPath[i]); sPar->maskPath[i] = NULL;}
                                                free(sPar->maskPath); sPar->maskPath = NULL;
                                              }
                                              if (sPar->nOfZPath) {
                                                for (i=0; i<sPar->nbTypes; i++) {free(sPar->nOfZPath[i]); sPar->nOfZPath[i] = NULL;}
                                                free(sPar->nOfZPath); sPar->nOfZPath = NULL;
                                              }
                                              sPar->nbTypes  = atoi(kv[1]);
                                              sPar->maskPath = (char**)malloc(sPar->nbTypes * sizeof(char*));
                                              sPar->nOfZPath = (char**)malloc(sPar->nbTypes * sizeof(char*));
                                              for (i=0; i<sPar->nbTypes; i++) {
                                                sPar->maskPath[i] = (char*)malloc(STRING_LENGTH_MAX * sizeof(char));
                                                sPar->nOfZPath[i] = (char*)malloc(STRING_LENGTH_MAX * sizeof(char));
                                              }
  }
  else if (!strcmp(kv[0], "maskPath"))        sprintf(sPar->maskPath[atoi(kv[1])], "%s", kv[2]);
  else if (!strcmp(kv[0], "nOfZPath"))        sprintf(sPar->nOfZPath[atoi(kv[1])], "%s", kv[2]);
  else if (!strcmp(kv[0], "n_gal")) {
                                              if (sPar->n_gal) free(sPar->n_gal);
                                              sPar->n_gal = makeDoubleArray(kv, count);
                                              for (i=0; i<sPar->nbTypes; i++) sPar->n_gal[i] /= ARCMIN_SQ_TO_RADIAN_SQ;
  }
  
  //-- 4. Lensing & outputs
  else if (!strcmp(kv[0], "doNoise"))         sPar->doNoise   = atoi(kv[1]);
  else if (!strcmp(kv[0], "doWgt"))           sPar->doWgt     = atoi(kv[1]);
  else if (!strcmp(kv[0], "signConv"))        sPar->signConv  = atoi(kv[1]);
  else if (!strcmp(kv[0], "sigma_eps")) {
                                              if (sPar->sigma_eps) free(sPar->sigma_eps);
                                              sPar->sigma_eps = makeDoubleArray(kv, count);
  }
  else if (!strcmp(kv[0], "doLensing")) {
                                              if (sPar->doLensing) free(sPar->doLensing);
                                              sPar->doLensing = makeIntArray(kv, count);
  }
  else if (!strcmp(kv[0], "outPrefix"))       setPathWhichCanBeBlank(sPar->outPrefix, kv, count);
  else if (!strcmp(kv[0], "outStyle"))        sPar->outStyle = atoi(kv[1]);
  
  //-- 5. Variable depth
  else if (!strcmp(kv[0], "doVariableDepth")) sPar->doVariableDepth = atoi(kv[1]);
  else if (!strcmp(kv[0], "nbDepthMaps")) {
                                              if (sPar->depthMapPath) {
                                                for (i=0; i<sPar->nbDepthMaps; i++) {free(sPar->depthMapPath[i]); sPar->depthMapPath[i] = NULL;}
                                                free(sPar->depthMapPath); sPar->depthMapPath = NULL;
                                              }
                                              sPar->nbDepthMaps  = atoi(kv[1]);
                                              sPar->depthMapPath = (char**)malloc(sPar->nbDepthMaps * sizeof(char*));
                                              for (i=0; i<sPar->nbDepthMaps; i++) sPar->depthMapPath[i] = (char*)malloc(STRING_LENGTH_MAX * sizeof(char));
  }
  else if (!strcmp(kv[0], "depthMapPath"))    sprintf(sPar->depthMapPath[atoi(kv[1])], "%s", kv[2]);
  else if (!strcmp(kv[0], "N_depth")) {
                                              if (sPar->VD_nOfZPath) {
                                                for (i=0; i<sPar->VD_nbNOfZ; i++) {free(sPar->VD_nOfZPath[i]); sPar->VD_nOfZPath[i] = NULL;}
                                                free(sPar->VD_nOfZPath); sPar->VD_nOfZPath = NULL;
                                              }
                                              sPar->N_depth   = atoi(kv[1]);
                                              sPar->VD_nbNOfZ = sPar->N_depth * sPar->nbTomo;
                                              if (sPar->VD_nbNOfZ > 0) sPar->VD_nOfZPath = (char**)malloc(sPar->VD_nbNOfZ * sizeof(char*));
                                              for (i=0; i<sPar->VD_nbNOfZ; i++) sPar->VD_nOfZPath[i] = (char*)malloc(STRING_LENGTH_MAX * sizeof(char));
  }
  else if (!strcmp(kv[0], "bin_depth")) {
                                              if (sPar->bin_depth) free(sPar->bin_depth);
                                              sPar->bin_depth = makeDoubleArray(kv, count);
  }
  else if (!strcmp(kv[0], "nbTomo")) {
                                              if (sPar->VD_nOfZPath) {
                                                for (i=0; i<sPar->VD_nbNOfZ; i++) {free(sPar->VD_nOfZPath[i]); sPar->VD_nOfZPath[i] = NULL;}
                                                free(sPar->VD_nOfZPath); sPar->VD_nOfZPath = NULL;
                                              }
                                              sPar->nbTomo    = atoi(kv[1]);
                                              sPar->VD_nbNOfZ = sPar->N_depth * sPar->nbTomo;
                                              if (sPar->VD_nbNOfZ > 0) sPar->VD_nOfZPath = (char**)malloc(sPar->VD_nbNOfZ * sizeof(char*));
                                              for (i=0; i<sPar->VD_nbNOfZ; i++) sPar->VD_nOfZPath[i] = (char*)malloc(STRING_LENGTH_MAX * sizeof(char));
  }
  else if (!strcmp(kv[0], "a_n_gal")) {
                                              if (sPar->a_n_gal) free(sPar->a_n_gal);
                                              sPar->a_n_gal = makeDoubleArray(kv, count);
  }
  else if (!strcmp(kv[0], "b_n_gal")) {
                                              if (sPar->b_n_gal) free(sPar->b_n_gal);
                                              sPar->b_n_gal = makeDoubleArray(kv, count);
  }
  else if (!strcmp(kv[0], "a_sigma_eps")) {
                                              if (sPar->a_sigma_eps) free(sPar->a_sigma_eps);
                                              sPar->a_sigma_eps = makeDoubleArray(kv, count);
  }
  else if (!strcmp(kv[0], "b_sigma_eps")) {
                                              if (sPar->b_sigma_eps) free(sPar->b_sigma_eps);
                                              sPar->b_sigma_eps = makeDoubleArray(kv, count);
  }
  else if (!strcmp(kv[0], "VD_nOfZPath"))     sprintf(sPar->VD_nOfZPath[atoi(kv[1])], "%s", kv[2]);
  
  else return 1; //-- unknown = 1
  
  return 0;
}

void readParameters(char name[], Salmo_param *sPar, int verbose)
{
  sprintf(sPar->parPath, "%s", name);
  FILE *file = fopen(name, "r");
  int verbose2 = (sPar->MPIInd == 0) && (sPar->verbose < 99) && verbose;
  
  char line[STRING_LENGTH_MAX], kv[STRING_LENGTH_MAX][STRING_LENGTH_MAX], *buffer;
  int count, unknown;
  
  //-- Read
  buffer = fgets(line, STRING_LENGTH_MAX, file);
  while (buffer != NULL) {
    ignoreComments(line);
    count   = getKeyAndValues(line, kv);
    unknown = findParameterKey(sPar, kv, count);
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

int updateFromCommandLine(int argc, char *argv[], Salmo_param *sPar)
{
  //-- Update from the command line
  
  int verbose = (sPar->MPIInd == 0) && (sPar->verbose < 99);
  char kv[STRING_LENGTH_MAX][STRING_LENGTH_MAX];
  int i, count, unknown;
  
  for (i=3; i<argc; i++) {
    count = getKeyAndValues(argv[i], kv);
    if (count == -1) break; //-- Help detected
    unknown = findParameterKey(sPar, kv, count); 
    if (unknown == 1 && count > 1 && verbose) printf("Detected the keyword \"%s\" but unable to interpret, skipped\n", kv[0]);
  }
  
  return (int)(count == -1);
}

void setParameters(Salmo_param *sPar)
{
  int verbose = sPar->MPIInd == 0 && sPar->verbose < 99;
  u_int32_t seed;
  int i;
  
  //-- Precomputed part - Generality
  if (strchr(sPar->seed, 'r') == NULL) seed = strtoul(sPar->seed, NULL, 10);
  else {
    seed = renewSeed();
    sprintf(sPar->seed, "%u (random)", seed);
  }
  sPar->generator = initializeGenerator(seed);
  sPar->nbPix     = 12 * sPar->nside * sPar->nside;
  sPar->A_pix     = FULL_SKY * DEGREE_SQ_TO_RADIAN_SQ / sPar->nbPix;
  if (sPar->N_z_map < 0) {
    sPar->zMapRange[0] = sPar->bin_z_map[0];
    sPar->zMapRange[1] = sPar->bin_z_map[1];
    sPar->zMapRange[2] = sPar->bin_z_map[2];
    sPar->N_z_map      = (int)round((sPar->zMapRange[1] - sPar->zMapRange[0]) / sPar->zMapRange[2]);
    if (sPar->bin_z_map) free(sPar->bin_z_map);
    sPar->bin_z_map    = (double*)malloc((sPar->N_z_map+1) * sizeof(double));
    for (i=0; i<=sPar->N_z_map; i++) sPar->bin_z_map[i] = sPar->zMapRange[0] + i * sPar->zMapRange[2];
  }
  else {
    sPar->zMapRange[0] = sPar->bin_z_map[0];
    sPar->zMapRange[1] = sPar->bin_z_map[sPar->N_z_map];
    sPar->zMapRange[2] = 0.0; //-- dz_map is not defined.
  }
  if (sPar->z_map)       free(sPar->z_map);
  if (sPar->half_dz_map) free(sPar->half_dz_map);
  sPar->z_map       = (double*)malloc(sPar->N_z_map * sizeof(double));
  sPar->half_dz_map = (double*)malloc(sPar->N_z_map * sizeof(double));
  for (i=0; i<sPar->N_z_map; i++) {
    sPar->z_map[i]       = 0.5 * (sPar->bin_z_map[i] + sPar->bin_z_map[i+1]);
    sPar->half_dz_map[i] = 0.5 * (sPar->bin_z_map[i+1] - sPar->bin_z_map[i]);
  }
  sPar->skipLensing = 1;
  for (i=0; i<sPar->nbTypes; i++) {
    if (sPar->doLensing[i] == 1) {
      sPar->skipLensing = 0;
      break;
    }
  }
  
  if (sPar->doVariableDepth == 0) {
    sPar->nbDepthMaps = 0;
    sPar->N_depth     = 0;
    sPar->nbTomo      = 0;
    sPar->VD_nbNOfZ   = 0;
    sPar->VD_nbTypes  = 0;
  }
  else {
    sPar->skipLensing = 0;
//     sPar->VD_nbNOfZ; //-- Already initialized
    sPar->VD_nbTypes = sPar->nbDepthMaps * sPar->nbTomo;
  }
  sPar->totNbTypes = sPar->nbTypes + sPar->VD_nbTypes;
  
  
  //-- Running part
  //sPar->MPISize;  //-- Should not be initialized here
  //sPar->MPIInd;   //-- Should not be initialized here
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

void printParam(Salmo_param *sPar)
{
  int i;
  printf("seed         = %s\n", sPar->seed);
  printf("verbose      = %d\n", sPar->verbose);
  
  printf("\n");
  printf("nside        = %lld\n", sPar->nside);
  printf("N_z_map      = %d\n", sPar->N_z_map);
  printf("bin_z_map    = "); printDoubleArray(sPar->bin_z_map, sPar->N_z_map+1, 1.0, 2); printf(" [-]\n");
  printf("denPrefix    = \"%s\"\n", sPar->denPrefix);
  printf("lenPrefix    = \"%s\"\n", sPar->lenPrefix);
  printf("runTag       = \"%s\"\n", sPar->runTag);
  
  printf("\n");
  printf("nbTypes      = %d\n", sPar->nbTypes);
  for (i=0; i<sPar->nbTypes;i++) printf("maskPath[%2d] = \"%s\"\n", i, sPar->maskPath[i]);
  for (i=0; i<sPar->nbTypes;i++) printf("nOfZPath[%2d] = \"%s\"\n", i, sPar->nOfZPath[i]);
  printf("n_gal        = "); printDoubleArray(sPar->n_gal, sPar->nbTypes, 1.0/RADIAN_SQ_TO_ARCMIN_SQ, 2); printf(" [arcmin^-2]\n");
  
  printf("\n");
  printf("doNoise      = %d\n", sPar->doNoise);
  printf("doWgt        = %d\n", sPar->doWgt);
  printf("signConv     = %d\n", sPar->signConv);
  printf("sigma_eps    = "); printDoubleArray(sPar->sigma_eps, sPar->nbTypes, 1.0, 5); printf("\n");
  printf("doLensing    = "); printIntArray(sPar->doLensing, sPar->nbTypes); printf("\n");
  printf("outPrefix    = \"%s\"\n", sPar->outPrefix);
  printf("outStyle     = %d\n", sPar->outStyle);
  
  printf("\n");
  printf("doVariableDepth  = %d\n", sPar->doVariableDepth);
  
  if (sPar->doVariableDepth == 1) {
    printf("nbDepthMaps      = %d\n", sPar->nbDepthMaps);
    for (i=0; i<sPar->nbDepthMaps;i++) printf("depthMapPath[%2d] = \"%s\"\n", i, sPar->depthMapPath[i]);
    printf("N_depth          = %d\n", sPar->N_depth);
    printf("bin_depth        = "); printDoubleArray(sPar->bin_depth, sPar->N_depth+1, 1.0, 2); printf(" [-]\n");
    printf("nbTomo           = %d\n", sPar->nbTomo);
    printf("a_n_gal          = "); printDoubleArray(sPar->a_n_gal, sPar->nbTomo, 1.0, 3); printf("\n");
    printf("b_n_gal          = "); printDoubleArray(sPar->b_n_gal, sPar->nbTomo, 1.0, 3); printf("\n");
    printf("a_sigma_eps      = "); printDoubleArray(sPar->a_sigma_eps, sPar->nbTomo, 1.0, 3); printf("\n");
    printf("b_sigma_eps      = "); printDoubleArray(sPar->b_sigma_eps, sPar->nbTomo, 1.0, 3); printf("\n");
    for (i=0; i<sPar->VD_nbNOfZ;i++) printf("VD_nOfZPath[%2d]  = \"%s\"\n", i, sPar->VD_nOfZPath[i]);
  }
  return;
}

void printCompleteParam(Salmo_param *sPar)
{
  int i;
  printf("\n");
  
  printf("------ 1. Generality ------\n");
  printf("parPath      = \"%s\"\n", sPar->parPath);
  printf("seed         = %s\n", sPar->seed);
  printf("verbose      = %d\n", sPar->verbose);
  printf("\n");
  
  printf("------ 2. Input maps ------\n");
  printf("nside        = %lld\n", sPar->nside);
  printf("N_z_map      = %d\n", sPar->N_z_map);
  printf("bin_z_map    = "); printDoubleArray(sPar->bin_z_map, sPar->N_z_map+1, 1.0, 2); printf(" [-]\n");
  printf("denPrefix    = \"%s\"\n", sPar->denPrefix);
  printf("lenPrefix    = \"%s\"\n", sPar->lenPrefix);
  printf("runTag       = \"%s\"\n", sPar->runTag);
  printf("\n");
  
  printf("------ 3. Selection functions ------\n");
  printf("nbTypes      = %d\n", sPar->nbTypes);
  for (i=0; i<sPar->nbTypes;i++) printf("maskPath[%2d] = \"%s\"\n", i, sPar->maskPath[i]);
  for (i=0; i<sPar->nbTypes;i++) printf("nOfZPath[%2d] = \"%s\"\n", i, sPar->nOfZPath[i]);
  printf("n_gal        = "); printDoubleArray(sPar->n_gal, sPar->nbTypes, 1.0/RADIAN_SQ_TO_ARCMIN_SQ, 2); printf(" [arcmin^-2]\n");
  printf("\n");
  
  printf("------ 4. Lensing & outputs ------\n");
  printf("doNoise      = %d\n", sPar->doNoise);
  printf("doWgt        = %d\n", sPar->doWgt);
  printf("signConv     = %d\n", sPar->signConv);
  printf("sigma_eps    = "); printDoubleArray(sPar->sigma_eps, sPar->nbTypes, 1.0, 5); printf("\n");
  printf("doLensing    = "); printIntArray(sPar->doLensing, sPar->nbTypes); printf("\n");
  printf("outPrefix    = \"%s\"\n", sPar->outPrefix);
  printf("outStyle     = %d\n", sPar->outStyle);
  printf("\n");
  
  printf("------ 5. Variable depth ------\n");
  printf("doVariableDepth  = %d\n", sPar->doVariableDepth);
  printf("nbDepthMaps      = %d\n", sPar->nbDepthMaps);
  for (i=0; i<sPar->nbDepthMaps;i++) printf("depthMapPath[%2d] = \"%s\"\n", i, sPar->depthMapPath[i]);
  printf("N_depth          = %d\n", sPar->N_depth);
  printf("bin_depth        = "); printDoubleArray(sPar->bin_depth, sPar->N_depth+1, 1.0, 2); printf(" [-]\n");
  printf("nbTomo           = %d\n", sPar->nbTomo);
  printf("a_n_gal          = "); printDoubleArray(sPar->a_n_gal, sPar->nbTomo, 1.0, 3); printf("\n");
  printf("b_n_gal          = "); printDoubleArray(sPar->b_n_gal, sPar->nbTomo, 1.0, 3); printf("\n");
  printf("a_sigma_eps      = "); printDoubleArray(sPar->a_sigma_eps, sPar->nbTomo, 1.0, 3); printf("\n");
  printf("b_sigma_eps      = "); printDoubleArray(sPar->b_sigma_eps, sPar->nbTomo, 1.0, 3); printf("\n");
  for (i=0; i<sPar->VD_nbNOfZ;i++) printf("VD_nOfZPath[%2d]  = \"%s\"\n", i, sPar->VD_nOfZPath[i]);
  printf("\n");
  
  printf("------ Precomputed part ------\n");
  printf("generator    = %s (default)\n", gsl_rng_name(sPar->generator));
  printf("nbPix        = %lld\n", sPar->nbPix);
  printf("A_pix        = %f [arcmin^2]\n", sPar->A_pix*RADIAN_SQ_TO_ARCMIN_SQ);
  printf("z_map        = "); printDoubleArray(sPar->z_map, sPar->N_z_map, 1.0, 2); printf(" [-]\n");
  printf("half_dz_map  = "); printDoubleArray(sPar->half_dz_map, sPar->N_z_map, 1.0, 2); printf(" [-]\n");
  printf("zMapRange    = %.3f, %.3f, %.3f [-]\n", sPar->zMapRange[0], sPar->zMapRange[1], sPar->zMapRange[2]);
  printf("skipLensing  = %d\n", sPar->skipLensing);
  printf("VD_nbNOfZ    = %d\n", sPar->VD_nbNOfZ);
  printf("VD_nbTypes   = %d\n", sPar->VD_nbTypes);
  printf("totNbTypes   = %d\n", sPar->totNbTypes);
  printf("\n");
  return;
}

//----------------------------------------------------------------------


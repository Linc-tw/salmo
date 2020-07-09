

  //------------------------------------------------------//
  //--  main.c                                          --//
  //--  Version 2020.07.09                              --//
  //--                                                  --//
  //--  Copyright (C) 2020 - Chieh-An Lin               --//
  //--  GNU GPLv3 - https://www.gnu.org/licenses/       --//
  //------------------------------------------------------//


#include "main.h"


int main(int argc, char *argv[])
{
  //-- This code has not been implemented with MPI.
  //-- Here the related variables & functions are only kept for future developments.
  int MPISize = 1;
  int MPIInd  = 0;
  
#ifdef __SALMO_USE_MPI__
  //-- MPI initialization
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &MPISize);
  MPI_Comm_rank(MPI_COMM_WORLD, &MPIInd);
#endif
  
  //-- Start stopwatch
  clock_t start = clock();
  
  //-- Load inputs
  //-- Below are some strings that will be automatically recognized as keys for parameter files.
  char *parPath = (argc >= 2) ? argv[1] : "../param/salmoParam.par";
  if (!strcmp(parPath, "default"))        parPath = "../param/salmoParam.par";
  else if (!strcmp(parPath, "anser"))     parPath = "../param/MFPParam_anser.par";
  else if (!strcmp(parPath, "buceros"))   parPath = "../param/MFPParam_buceros.par";
  else if (!strcmp(parPath, "corvus"))    parPath = "../param/MFPParam_corvus.par";
  else if (!strcmp(parPath, "cygnus"))    parPath = "../param/MFPParam_cygnus.par";
  else if (!strcmp(parPath, "cuculus"))   parPath = "../param/MFPParam_cuculus.par";
  else if (!strcmp(parPath, "diomedea"))  parPath = "../param/MFPParam_diomedea.par";
  else if (!strcmp(parPath, "egretta"))   parPath = "../param/MFPParam_egretta.par";
  else if (!strcmp(parPath, "falco"))     parPath = "../param/MFPParam_falco.par";
  else if (!strcmp(parPath, "gallus"))    parPath = "../param/MFPParam_gallus.par";
  else if (!strcmp(parPath, "grus"))      parPath = "../param/MFPParam_grus.par";
  else if (!strcmp(parPath, "hirundo"))   parPath = "../param/MFPParam_hirundo.par";
  else if (!strcmp(parPath, "icterus"))   parPath = "../param/MFPParam_icterus.par";
  else if (!strcmp(parPath, "jacana"))    parPath = "../param/MFPParam_jacana.par";
  else if (!strcmp(parPath, "larus"))     parPath = "../param/MFPParam_larus.par";
  else if (!strcmp(parPath, "milvus"))    parPath = "../param/MFPParam_milvus.par";
  else if (!strcmp(parPath, "numida"))    parPath = "../param/MFPParam_numida.par";
  else if (!strcmp(parPath, "otis"))      parPath = "../param/MFPParam_otis.par";
  else if (!strcmp(parPath, "pavo"))      parPath = "../param/MFPParam_pavo.par";
  else if (!strcmp(parPath, "raphus"))    parPath = "../param/MFPParam_raphus.par";
  else if (!strcmp(parPath, "strix"))     parPath = "../param/MFPParam_strix.par";
  else if (!strcmp(parPath, "turdus"))    parPath = "../param/MFPParam_turdus.par";
  else if (!strcmp(parPath, "upupa"))     parPath = "../param/MFPParam_upupa.par";
  else if (!strcmp(parPath, "vultur"))    parPath = "../param/MFPParam_vultur.par";
  else if (!strcmp(parPath, "zosterops")) parPath = "../param/MFPParam_zosterops.par";
  
  //-- Read the task index. 
  //-- If not provided, set to -1 (print instructions).
  int task = (argc >= 3) ? atoi(argv[2]) : -1;
  
  //-- Initialization
  Salmo_param *sPar = initialize_Salmo_param();
  readParameters(parPath, sPar, 0); //-- Get only verbose
  sPar->MPISize = MPISize;
  sPar->MPIInd  = MPIInd;
  
  //-- Print the header after read pkPar->verbose
  if (MPIInd == 0 && sPar->verbose < 99) {
    printf("\n");
    printf("        #################################################\n");
    printf("        ##                    Salmo                    ##\n");
    printf("        ##                                             ##\n");
    printf("        ##  Copyright (C) 2020 - Chieh-An Lin          ##\n");
    printf("        ##  GNU GPLv3 - https://www.gnu.org/licenses/  ##\n");
    printf("        #################################################\n");
    printf("\n");
  }
  
  //-- Initialization (continued)
  readParameters(parPath, sPar, 1);                   //-- Read again
  int help = updateFromCommandLine(argc, argv, sPar); //-- Update parameters from the command line
  setParameters(sPar);                                //-- Precalculate some parameters
  
  //-- Print parameters
  if (MPIInd == 0 && sPar->verbose < 99) {
    printf("Initialized parameters\n");
    printf("------------------------------------------------------------------------\n");
    
    if (task == -1 || task == -2) ;
    else {
      printParam(sPar);
      printf("------------------------------------------------------------------------\n");
    }
  }
  
#ifdef __SALMO_USE_MPI__
  //-- Stopping lock for no MPI task
  if (MPISize > 1 && !(task == 1)) {
    printf("Cannot use multiple processors for task %d\n", task);
    MPI_Finalize();
    return 1;
  }
  
  //-- Synchronize all MPI processes
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
  //-- Sandbox
  if (task == 0) {
    sandbox(sPar);
  }
  
  //-- Print instructions
  else if (task == -1) printInstructions(-1, 1);
  
  //-- Print complete parameters
  else if (task == -2) {
    printCompleteParam(sPar);
    printf("------------------------------------------------------------------------\n");
  }
  
  //-- WARNING This function has been deactivated in the public version.
  //-- Make kappa+gamma maps (3 extensions in a FITS file) given a kappa map
  else if (task == 1) {
    if (argc < 3 || help) printInstructions(task, 1);
    else {
      processLensingMaps(sPar);
    }
  }
  
  //-- Create mocks using a type map.
  //-- A type map is an int64 array.
  //-- Each int64 of the array represents a pixel.
  //-- The bit of each int64 indicates whether the selection function of a certain tracer is activated at the corresponding pixel or not.
  //-- Since int64 has 64 bits, this approach allows up to 64 different types of tracers.
  else if (task == 2) {
    if (argc < 3 || help) printInstructions(task, 1);
    else {
      processMock_typeMap_LOS(sPar);
    }
  }
  
  //-- Create mocks using a ratio map in LOS mode.
  //-- LOS mode means that the lensing maps provided are considered as lensing signals at the outer edge of the shell (not the center).
  //-- A ratio map is a float32 array.
  //-- Each float32 of the array represents the filling factor of a certain tracer in a pixel.
  //-- The filling factor can be any positive real number or zero.
  //-- It is the ratio of the expected tracer density in that pixel with regard to the reference value (usually the global mean density).
  //-- If filling factors are either 0 or 1, we find the case of the classical binary mask.
  //-- A ratio map can only provide information for one tracer.
  //-- Therefore, users need to specify paths of all ratio maps if multiple tracers are used.
  else if (task == 3) {
    if (argc < 3 || help) printInstructions(task, 1);
    else {
      if (sPar->doVariableDepth) processMock_depthMap_LOS(sPar);
      else                       processMock_ratioMap_LOS(sPar);
    }
  }
  
  //-- Create mocks using a ratio map in C_ell mode.
  //-- C_ell mode means that the lensing maps provided are considered as lensing signals at the center of the shell (not the outer edge).
  else if (task == 4) {
    if (argc < 3 || help) printInstructions(task, 1);
    else {
      processMock_ratioMap_CL(sPar);
    }
  }
  
  //-- You won't want to use this, but I will still explain.
  //-- Given a tracer, given its n(z), one can calculate its projected C_ell.
  //-- If we use this C_ell to generate a matter shell from Flask, then this matter shell becomes a projected density map (i.e. all redshift included).
  //-- In practice, we may have N tracers.
  //-- We can then calculate N sets of projected C_ell & generate N projected density maps.
  //-- Salmo will take that & sample redshifts beside to create catalogues.
  else if (task == 5) {
    if (argc < 5 || help) printInstructions(task, 1);
    else {
      int resol    = atoi(argv[3]);
      int nbSplits = atoi(argv[4]);
      printf("task = 5, resol = %d, nb of splits = %d\n", resol, nbSplits);
      processMock_ratioMap_projCL(sPar, resol, nbSplits);
    }
  }
  
  //-- Stop stopwatch
  clock_t finish = clock();
  if (MPIInd == 0 && sPar->verbose < 99) {
    printTime(start, finish);
    printf("\n");
  }
  
  free_Salmo_param(sPar);
#ifdef __SALMO_USE_MPI__
  MPI_Finalize();
#endif
  return 0;
}

void printInstructions(int task, int doHelp)
{
  if (task >= 1 && task <= 30) {
    printHeader(1, 0);
    printDetails(task, doHelp);
  }
  
  else { //-- Complete menu
    printHeader(1, doHelp);
    printDetails(-2, 0);
    printDetails(-1, 0);
    printDetails(0, 0); printDetails(1, 0); 
    printDetails(2, 0); printDetails(3, 0); printDetails(4, 0); 
  }
  
  printf("------------------------------------------------------------------------\n");
  return;
}

//-- Print header of instructions
//-- Provide examples if doHelp is on
void printHeader(int task, int doHelp)
{
  if (task >= 0 && task <= 30) {
    if (doHelp) {
      printf("Commands:\n");
      printf("  ./salmo PATH TASK\n");
      printf("  ./salmo PATH TASK -h\n");
      printf("  ./salmo PATH TASK KEY=VALUE\n");
      printf("\n");
      printf("Notes:\n");
      printf("  - PATH = path to the salmoParam.par file\n");
      printf("  - Can replace PATH with the string \"default\", equivalent to \"../param/salmoParam.par\"\n");
      printf("  - TASK = number of the task to do\n");
      printf("  - Add -h after TASK for more detailed instructions\n");
      printf("  - Non-array parameters can be updated by KEY=VALUE.\n");
      printf("\n");
      printf("Examples:\n");
      printf("  ./salmo default 0\n");
      printf("  ./salmo ../param/salmoParam.par 1\n");
      printf("  ./salmo default 1 runTag=\"_hello\"\n");
      printf("\n");
    }
      printf("Tasks:\n");
  }
  return; 
}

//-- Print detaileds of each option
//-- Provide variable format if doHelp is on
void printDetails(int task, int doHelp)
{
  if (task == 0) {
      printf("  0 = Sandbox\n");
    if (doHelp) {
      printf("\n");
      printf("Commands:\n");
      printf("  ./salmo PATH 0    # Sandbox\n");
    }
  }
  else if (task == -1) {
      printf(" -1 = Print help messages\n");
    if (doHelp) {
      printf("\n");
      printf("Commands:\n");
      printf("  ./salmo PATH -1   # Print help messages\n");
    }
  }
  else if (task == -2) {
      printf(" -2 = Print all parameters\n");
    if (doHelp) {
      printf("\n");
      printf("Commands:\n");
      printf("  ./salmo PATH -2   # Print all parameters\n");
    }
  }
  else if (task == 1) {
//       printf("  1 = Compute lensing maps\n");
      printf("  1 = (deactivated)\n");
    if (doHelp) {
      printf("\n");
      printf("Commands:\n");
//       printf("  ./salmo PATH 1    # delta to lensing\n");
      printf("  ./salmo PATH 1    # (deactivated)\n");
    }
  }
  else if (task == 2) {
      printf("  2 = Generate a mock, type map, LOS\n");
    if (doHelp) {
      printf("\n");
      printf("Commands:\n");
      printf("  ./salmo PATH 2    # Generate a mock using type map, lenMap interpreted as at the upper bin edges\n");
    }
  }
  else if (task == 3) {
      printf("  3 = Generate a mock, LOS\n");
    if (doHelp) {
      printf("\n");
      printf("Commands:\n");
      printf("  ./salmo PATH 3    # Generate a mock, lenMap interpreted as at the upper bin edges\n");
    }
  }
  else if (task == 4) {
      printf("  4 = Generate a mock, C_l\n");
    if (doHelp) {
      printf("\n");
      printf("Commands:\n");
      printf("  ./salmo PATH 4    # Generate a mock, lenMap interpreted as at the middle of bins\n");
    }
  }
  else if (task == 5) {
      printf("  5 = Generate galaxies from projected maps\n");
    if (doHelp) {
      printf("\n");
      printf("Commands:\n");
      printf("  ./salmo PATH 5 resol nbSplits  # Generate galaxies from projected maps\n");
      printf("                                 # Matter fields smoothed with a top hat of resol^2\n");
      printf("                                 # nbSplits tells how many types of galaxies should a map accounts for\n");
    }
  }
  return;
}

void sandbox(Salmo_param *sPar)
{
  char name[STRING_LENGTH_MAX];
  char name2[STRING_LENGTH_MAX];
  printf("\n");
  printf("This is a sandbox in Salmo.\n");
  printf("If you are a developer, please feel free to use this function for quick tests.\n");
  printf("\n");
  
//   printCompleteParam(sPar);
  
  //mpirun -n 3 ./salmo default 35
  printf("------------------------------------------------------------------------\n");
  return;
}


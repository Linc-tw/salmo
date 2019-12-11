

  //------------------------------------------------------//
  //--  main.c						--//
  //--  Version 2019.10.02				--//
  //--  						--//
  //--  Copyright (C) 2019 - Chieh-An Lin		--//
  //--  GNU GPLv3 - https://www.gnu.org/licenses/	--//
  //------------------------------------------------------//


#include "main.h"


int main(int argc, char *argv[])
{
  int MPISize = 1;
  int MPIInd  = 0;
  
#ifdef __MFP_USE_MPI__
  //-- MPI initialization
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &MPISize);
  MPI_Comm_rank(MPI_COMM_WORLD, &MPIInd);
#endif
  
  //-- Start stopwatch
  clock_t start = clock();
  
  //-- Read inputs
  //-- anser, buceros, corvus, diomedea, egretta, 
  //-- falco, gallus, hirundo, icterus, jacana, 
  //-- larus, milvus, numida, otis, pavo, 
  //-- raphus, strix, turdus, upupa, vultur, 
  //-- zosterops
  char *parPath = (argc >= 2) ? argv[1] : "../param/MFPParam.par";
  if (!strcmp(parPath, "default"))        parPath = "../param/MFPParam.par";
  else if (!strcmp(parPath, "anser"))     parPath = "../param/MFPParam_anser.par";
  else if (!strcmp(parPath, "buceros"))   parPath = "../param/MFPParam_buceros.par";
  else if (!strcmp(parPath, "corvus"))    parPath = "../param/MFPParam_corvus.par";
  else if (!strcmp(parPath, "diomedea"))  parPath = "../param/MFPParam_diomedea.par";
  else if (!strcmp(parPath, "egretta"))   parPath = "../param/MFPParam_egretta.par";
  else if (!strcmp(parPath, "falco"))     parPath = "../param/MFPParam_falco.par";
  else if (!strcmp(parPath, "gallus"))    parPath = "../param/MFPParam_gallus.par";
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
  int task = (argc >= 3) ? atoi(argv[2]) : -1;
  
  //-- Initialization
  MFP_param *mPar = initialize_MFP_param();
  readParameters(parPath, mPar, 0); //-- Get only verbose
  mPar->MPISize = MPISize;
  mPar->MPIInd  = MPIInd;
  
  //-- Print the header after read pkPar->verbose
  if (MPIInd == 0 && mPar->verbose < 99) {
    printf("\n");
    printf("        #################################################\n");
    printf("        ##               mockFootprint                 ##\n");
    printf("        ##                                             ##\n");
    printf("        ##  Copyright (C) 2018 - Chieh-An Lin          ##\n");
    printf("        ##  GNU GPLv3 - https://www.gnu.org/licenses/  ##\n");
    printf("        #################################################\n");
    printf("\n");
  }
  
  //-- Initialization (continued)
  readParameters(parPath, mPar, 1);                   //-- Read again
  int help = updateFromCommandLine(argc, argv, mPar); //-- Update parameters from the command line
  setParameters(mPar);                                //-- Precalculate some parameters
  
  //-- Print parameters
  if (MPIInd == 0 && mPar->verbose < 99) {
    printf("Initialized parameters\n");
    printf("------------------------------------------------------------------------\n");
    
    if (task == -3) ;
    else {
      printParam(mPar);
      printf("------------------------------------------------------------------------\n");
    }
  }
  
#ifdef __MFP_USE_MPI__
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
    sandbox(mPar);
  }
  
  else if (task == -3) {
    printCompleteParam(mPar);
    printf("------------------------------------------------------------------------\n");
  }
  
  else if (task == 1) {
    if (argc < 3 || help) printInstructions(task, 1);
    else {
      processLensingMaps(mPar);
    }
  }
  
  else if (task == 2) {
    if (argc < 3 || help) printInstructions(task, 1);
    else {
      processMock_typeMap_LOS(mPar);
    }
  }
  
  else if (task == 3) {
    if (argc < 3 || help) printInstructions(task, 1);
    else {
      if (mPar->doVariableDepth) processMock_depthMap_LOS(mPar);
      else                       processMock_ratioMap_LOS(mPar);
    }
  }
  
  else if (task == 4) {
    if (argc < 3 || help) printInstructions(task, 1);
    else {
      processMock_ratioMap_CL(mPar);
    }
  }
  
  else if (task == 5) {
    if (argc < 5 || help) printInstructions(task, 1);
    else {
      int resol    = atoi(argv[3]);
      int nbSplits = atoi(argv[4]);
      printf("task = 5, resol = %d, nb of splits = %d\n", resol, nbSplits);
      processMock_ratioMap_projCL(mPar, resol, nbSplits);
    }
  }
  
  else if (task == -1) printInstructions(-1, 1);
  else if (task == -2) printInstructions(-2, 0);
  else                 printInstructions(-1, 0);
  
  //-- Stop stopwatch
  clock_t finish = clock();
  if (MPIInd == 0 && mPar->verbose < 99) {
    printTime(start, finish);
    printf("\n");
  }
  
  free_MFP_param(mPar);
#ifdef __MFP_USE_MPI__
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

  else if (task == -2) { //-- Simple menu
    printHeader(1, doHelp);
    printDetails(0, 0); printDetails(1, 0); printDetails(2, 0); printDetails(3, 0); printDetails(4, 0); 
  }
  
  else { //-- Complete menu
    printHeader(1, doHelp);
    printDetails(0, 0); printDetails(1, 0); printDetails(2, 0); printDetails(3, 0); printDetails(4, 0); 
  }
  
  printf("------------------------------------------------------------------------\n");
  return;
}

void printHeader(int task, int doHelp)
{
  if (task >= 0 && task <= 20) {
    if (doHelp) {
      printf("Commands:\n");
      printf("  ./mockFootprint PATH TASK\n");
      printf("  ./mockFootprint PATH TASK -h\n");
      printf("  ./mockFootprint PATH TASK KEY=VALUE\n");
      printf("\n");
      printf("Notes:\n");
      printf("  - PATH = path to the camelusParam.par file\n");
      printf("  - Can replace PATH with the string \"default\", equivalent to \"../param/camelusParam.par\"\n");
      printf("  - TASK = number of the task to do\n");
      printf("  - Add -h after TASK for more detailed instructions\n");
      printf("  - Some parameters can be updated by KEY=VALUE.\n");
      printf("\n");
      printf("Examples:\n");
      printf("  ./mockFootprint default 0\n");
      printf("  ./mockFootprint ../param/camelusParam.par 1\n");
      printf("  ./mockFootprint default 1 runTag=\"_hello\"\n");
      printf("\n");
    }
      printf("Tasks:\n");
  }
  return; 
}

void printDetails(int task, int doHelp)
{
  if (task == 0) {
      printf("  0 = Sandbox\n");
    if (doHelp) {
      printf("\n");
      printf("Commands:\n");
      printf("  ./mockFootprint PATH 0    # Sandbox\n");
    }
  }
  else if (task == 1) {
      printf("  1 = Compute lensing maps\n");
    if (doHelp) {
      printf("\n");
      printf("Commands:\n");
      printf("  ./mockFootprint PATH 1    # delta to lensing\n");
    }
  }
  else if (task == 2) {
      printf("  2 = Generate a mock, type map, LOS\n");
    if (doHelp) {
      printf("\n");
      printf("Commands:\n");
      printf("  ./mockFootprint PATH 2    # Generate a mock using type map, lenMap interpreted as at the upper bin edges\n");
    }
  }
  else if (task == 3) {
      printf("  3 = Generate a mock, LOS\n");
    if (doHelp) {
      printf("\n");
      printf("Commands:\n");
      printf("  ./mockFootprint PATH 3    # Generate a mock, lenMap interpreted as at the upper bin edges\n");
    }
  }
  else if (task == 4) {
      printf("  4 = Generate a mock, C_l\n");
    if (doHelp) {
      printf("\n");
      printf("Commands:\n");
      printf("  ./mockFootprint PATH 4    # Generate a mock, lenMap interpreted as at the middle of bins\n");
    }
  }
  else if (task == 5) {
      printf("  5 = Generate galaxies from projected maps\n");
    if (doHelp) {
      printf("\n");
      printf("Commands:\n");
      printf("  ./mockFootprint PATH 5 resol nbSplits  # Generate galaxies from projected maps\n");
      printf("                                         # Matter fields smoothed with a top hat of resol^2\n");
      printf("                                         # nbSplits tells how many types of galaxies should a map accounts for\n");
    }
  }
  return;
}

void sandbox(MFP_param *mPar)
{
  char name[STRING_LENGTH_MAX];
  char name2[STRING_LENGTH_MAX];
  printf("This is a sandbox of MFP.\n\n");
  
  printCompleteParam(mPar);
  
  //mpirun -n 3 ./camelus default 35
  //-- Hold root processor if others have not finished
  printf("------------------------------------------------------------------------\n");
  return;
}


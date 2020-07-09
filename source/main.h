

  //------------------------------------------------------//
  //--  main.h                                          --//
  //--  Version 2020.07.09                              --//
  //--                                                  --//
  //--  Copyright (C) 2020 - Chieh-An Lin               --//
  //--  GNU GPLv3 - https://www.gnu.org/licenses/       --//
  //------------------------------------------------------//


#include "commonHeader.h"

#ifndef __SALMO_MAIN__
#define __SALMO_MAIN__

#ifdef __SALMO_USE_HEALPIX_CXX__
#include "wrapper.h"
#endif

#ifdef __SALMO_USE_MPI__
#include <mpi/mpi.h>
#endif

#include "FITSFunctions.h"
#include "HEALPixFunctions.h"
#include "parameters.h"
#include "sampling.h"


int main(int argc, char *argv[]);
void printInstructions(int task, int doHelp);
void printHeader(int task, int doHelp);
void printDetails(int task, int doHelp);
void sandbox(Salmo_param *mPar);

#endif


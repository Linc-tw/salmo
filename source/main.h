

  //------------------------------------------------------//
  //--  main.h						--//
  //--  Version 2019.01.25				--//
  //--  						--//
  //--  Copyright (C) 2019 - Chieh-An Lin		--//
  //--  GNU GPLv3 - https://www.gnu.org/licenses/	--//
  //------------------------------------------------------//


#include "commonHeader.h"

#ifndef __MFP_MAIN__
#define __MFP_MAIN__

#ifdef __MFP_USE_HEALPIX_CXX__
#include "wrapper.h"
#endif

#ifdef __MFP_USE_MPI__
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
void sandbox(MFP_param *mPar);

#endif


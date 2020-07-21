

    ##############################
    ##  flask_salmo_wrapper.py  ##
    ##  Chieh-An Lin            ##
    ##  Version 2020.07.21      ##
    ##############################


import time
import subprocess as spc

import numpy as np


################################################################################
## Global variables

GALAXY_BIAS = 2.0
FLASK_EXEC  = '../../flask/bin/flask' #WARNING Change this!
FLASK_PARAM = 'flask_input/flaskParam.config'
DEMO_PATH   = './'
DEMO_PATH_FROM_BUILD = '../demo/'
BUILD_PATH  = '../build/'

################################################################################
## Call functions

def printTime(start, stop):
  """
  Print elapsed computation time.
  """
  duration = stop - start
  hours    = int(duration / 3600.0)
  remain   = duration % 3600.0
  minutes  = int(remain / 60.0)
  seconds  = remain % 60.0
  
  if hours == 0:
    if minutes == 0:
      print('Computation time = %.2f secs' % duration)
    else:
      print('Computation time = %d m %d s' % (minutes, int(seconds)))
  else:
    print('Computation time = %d h %d m %d s' % (hours, minutes, int(seconds)))
  return

def callFlask(runInd):
  """
  Call Flask using the inputs in demo.
  """
  
  start = time.time()
  
  ## If runInd is 0, we will use a specific seed, so that the unitary test will work.
  if runInd == 0:
    seed = 2718281
  else:
    seed = np.random.randint(0, 10000000)
  
  RNDSEED           = '%d' % seed
  SCALE_CLS         = '%g' % GALAXY_BIAS**2
  FIELDS_INFO       = '%sflask_input/fieldInfo.dat' % DEMO_PATH
  CL_PREFIX         = '%sflask_input/CL_' % DEMO_PATH
  MAPFITS_PREFIX    = '%sinput/denMap_run%d_' % (DEMO_PATH, runInd)
  SHEAR_FITS_PREFIX = '%sinput/lenMap_run%d_' % (DEMO_PATH, runInd)
  
  ## We will call Flask twice.
  
  print('################################################################################')
  print('## Flask - galaxy density maps - bias = %f' % GALAXY_BIAS)
  print()
  
  ## In the 1st call, we scale the matter C_ell by bias^2, so that 
  ## it becomes galaxy C_ell. 
  ## We don't calculate lensing by asking Flask to exit at MAPFITS_PREFIX.
  ## The resulting maps are galaxy number density maps.
  spc.run([
    FLASK_EXEC, FLASK_PARAM,
    'RNDSEED:', RNDSEED, 
    'FIELDS_INFO:', FIELDS_INFO,
    'CL_PREFIX:', CL_PREFIX,
    'SCALE_CLS:', SCALE_CLS, 
    'EXIT_AT:', 'MAPFITS_PREFIX', 
    'DENS2KAPPA:', '0',
    'MAPFITS_PREFIX:', MAPFITS_PREFIX, 
    'SHEAR_FITS_PREFIX:', '0'
  ])
  
  start1 = time.time()
  printTime(start, start1)
  print('################################################################################')
  print('## Flask - lensing maps with the same seed')
  print()
  
  ## In the 2nd call, we do not scale C_ell.
  ## We calculate lensing maps with the LOS approach, and do not output density maps.
  ## Also, we use the same seed as the 1st call.
  ## This will results in lensing maps which are in phase with the previously-
  ## created galaxy density maps.
  spc.run([
    FLASK_EXEC, FLASK_PARAM, 
    'RNDSEED:', RNDSEED, 
    'FIELDS_INFO:', FIELDS_INFO,
    'CL_PREFIX:', CL_PREFIX,
    'SCALE_CLS:', '1.0', 
    'EXIT_AT:', 'SHEAR_FITS_PREFIX', 
    'DENS2KAPPA:', '1',
    'MAPFITS_PREFIX:', '0', 
    'SHEAR_FITS_PREFIX:', SHEAR_FITS_PREFIX
  ])
  
  printTime(start1, time.time())
  print('################################################################################')
  return

def callSalmo(runInd):
  """
  Call Salmo using the inputs in demo.
  """
  
  start = time.time()
  print('## Salmo')
  
  ## If runInd is 0, we will use a specific seed, so that the unitary test will work.
  if runInd == 0:
    seed = str(2718281828)
  else:
    seed = 'random'
  
  spc.run([
    './salmo', 'default', '3', 
    'seed=%s' % seed,
    'verbose=2', 
    'runTag=_run%d' % runInd, 
    'denPrefix=%sinput/denMap' % DEMO_PATH_FROM_BUILD, 
    'lenPrefix=%sinput/lenMap' % DEMO_PATH_FROM_BUILD, 
    'outPrefix=%soutput/galCat' % DEMO_PATH_FROM_BUILD
  ], cwd=BUILD_PATH)
  
  printTime(start, time.time())
  print('################################################################################')
  return

def process(runInd):
  """
  Call Flask and Salmo together using the inputs in demo.
  """
  
  start = time.time()
  print()
  
  callFlask(runInd)
  callSalmo(runInd)
  
  print('## Total runtime')
  print()
  printTime(start, time.time())
  print('################################################################################')
  print()
  return

################################################################################
## Main

if __name__ == '__main__':
  process(0)

###############################################################################


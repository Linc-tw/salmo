

    ##############################
    ##  unitary_test.py         ##
    ##  Chieh-An Lin            ##
    ##  Version 2020.07.21      ##
    ##############################


#import unittest

import numpy as np
import astropy.io.fits as fits

################################################################################
## Global variables

EPS_NUM   = 1e-12
DEMO_PATH = './'

################################################################################
## Unitary test

def loadFitsGalCat(name, verbose=True):
  data = fits.getdata(name, 1)
  if verbose:
    print('Loaded \"%s\"' % name)
  return data

def compare(name1, name2):
  hdr1  = fits.getheader(name1, 1)
  data1 = loadFitsGalCat(name1, verbose=True)
  data2 = loadFitsGalCat(name2, verbose=True)
  nbCols = 7
  print()
  
  for i in range(nbCols):
    key = hdr1['TTYPE%d' % (i+1)]
    diff = np.fabs(data1.field(i) - data2.field(i))
    boolean = (diff < EPS_NUM).all()
    print('column = %12s, match = %s' % (key, boolean))
  return

def unitaryTest():
  """
  Call Flask and Salmo together using the inputs in demo.
  """
  name1 = '%soutput/galCat_ref_type0.fits' % DEMO_PATH
  name2 = '%soutput/galCat_ref_type1.fits' % DEMO_PATH
  name3 = '%soutput/galCat_run0_type0.fits' % DEMO_PATH
  name4 = '%soutput/galCat_run0_type1.fits' % DEMO_PATH
  
  print()
  print('## Unitary test - type 0')
  print()
  compare(name1, name3)
  
  print()
  print('## Unitary test - type 1')
  print()
  compare(name2, name4)
  print()
  return

#class TestGalCat(unittest.TestCase):

  #def test_type0(self):
    #print()
    #name1  = '%soutput/galCat_ref_type0.fits' % DEMO_PATH
    #name2  = '%soutput/galCat_run0_type0.fits' % DEMO_PATH
    #data1  = loadFitsGalCat(name1, verbose=True)
    #data2  = loadFitsGalCat(name2, verbose=True)
    #nbCols = len(data1[0])
      
    #for i in range(nbCols):
      #self.assertListEqual(data1.field(i).tolist(), data2.field(i).tolist())
  
  #def test_type1(self):
    #print()
    #name1  = '%soutput/galCat_ref_type1.fits' % DEMO_PATH
    #name2  = '%soutput/galCat_run0_type1.fits' % DEMO_PATH
    #data1  = loadFitsGalCat(name1, verbose=True)
    #data2  = loadFitsGalCat(name2, verbose=True)
    #nbCols = len(data1[0])
      
    #for i in range(nbCols):
      #self.assertListEqual(data1.field(i).tolist(), data2.field(i).tolist())

################################################################################
## Main

if __name__ == '__main__':
  unitaryTest()
  #unittest.main()

###############################################################################


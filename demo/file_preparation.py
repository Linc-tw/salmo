

    ##############################
    ##  file_preparation.py     ##
    ##  Chieh-An Lin            ##
    ##  Version 2020.07.21      ##
    ##############################


import numpy as np
import astropy.io.fits as fits
import healpy as hp


################################################################################
## Global variables

DEMO_PATH = './'

################################################################################
## Prepare masks

def patchToPixels(nsidePat, nsidePix, patch, nest=False, sort=False):
  """
  Return all pixel ring or nest numbers related to a given patch.
  
  Parameters
  ----------
  nsidePat : int
    nside of the patch
  nsidePix : int
    nside of the child patch
  patch : int
    ring or nest number of the patch
  nest : bool, optional
    consider nest
  sort : bool, optional
    sort the returned array
  
  Returns
  -------
  pix : int array
    ring or nest numbers of all pixels
  """
  patch  = patch if nest == True else hp.ring2nest(nsidePat, patch)
  length = (nsidePix // nsidePat)**2
  pix    = np.arange(patch*length, (patch+1)*length)
  pix    = pix if nest == True else hp.nest2ring(nsidePix, pix)
  if sort == True:
    pix.sort()
  return pix

def saveFitsFullMap(name, full, verbose=True):
  """
  Save a HEALPix map as a FITS file under a specific convention.
  
  Parameters
  ----------
  name : string
    full name of the file
  full : numpy array
    the HEALPix map to save, will be saved as a float32 array
  verbose : bool, optional
    print verbose message
  
  Returns
  -------
  No returns
  """
  full   = full.astype(np.float32)
  nbRows = full.size // 1024
  full   = full.reshape(nbRows, 1024)
  nside  = hp.npix2nside(full.size)
  
  HDU1 = fits.PrimaryHDU()
  HDU2 = fits.BinTableHDU.from_columns([
    fits.Column(name='VALUE', format='1024E', unit='-       ', array=full)
  ])
  
  hdr = HDU2.header
  hdr.append(('COMMENT',  'HEALPIX pixelisation'),                                        bottom=True)
  hdr.append(('ORDERING', 'RING    ',    'Pixel ordering scheme'),                        bottom=True)
  hdr.append(('COORDSYS', 'C       ',    'Ecliptic, Galactic or Celestial (equatorial)'), bottom=True)
  hdr.append(('NSIDE',    nside,         'nside of the pixel'),                           bottom=True)
  hdr.append(('FIRSTPIX', 0,             'First pixel # (0 based)'),                      bottom=True)
  hdr.append(('LASTPIX',  12*nside**2-1, 'Last pixel # (0 based)'),                       bottom=True)
  hdr.append(('INDXSCHM', 'IMPLICIT',    'Indexing: IMPLICIT or EXPLICIT'),               bottom=True)
  
  fits.HDUList([HDU1, HDU2]).writeto(name, overwrite=True)
  if verbose == True:
    print('Saved \"%s\"' % name)
  return

def saveFitsMask_demo():
  """
  Save two masks used in demo.
  """
  
  ## Create masks of N_side = 256
  nside = 256
  nbPix = 12 * nside * nside
  mask  = np.zeros(nbPix, dtype=np.float32)
  
  ## For type 1, we define the mask as a patch of N_side = 4.
  nsidePat = 4
  patch    = 88
  pix  = patchToPixels(nsidePat, nside, patch)
  mask[pix] = 1
  
  ## Save
  name = '%sinput/mask_type1.fits' % DEMO_PATH
  saveFitsFullMap(name, mask, verbose=True)
  
  ## Reset
  mask[pix] = 0
  
  ## For type 0, we define the mask as 2 patches of N_side = 8
  nsidePat = 8
  patch    = 305
  pix  = patchToPixels(nsidePat, nside, patch)
  mask[pix] = 1
  
  patch    = 337
  pix  = patchToPixels(nsidePat, nside, patch)
  mask[pix] = 1
  
  ## Save
  name = '%sinput/mask_type0.fits' % DEMO_PATH
  saveFitsFullMap(name, mask, verbose=True)
  return

################################################################################
## Prepare n(z)

def cdfToHist(xArr, cdf, bins):
  """
  Return a histogram-like distribution from a cumulative distribution 
  given a set of bin edges.
  
  Here, a "histogram-like distribution" does not know what the pdf value
  is at any point, but the integrated pdf over a given interval. It assumes
  that the distribution is uniform within each interval.
  
  Parameters
  ----------
  xArr : float array
    sample points of the cumulative distribution
  cdf : float array
    values of the cumulative distribution
  bins : float array
    bin edges of the histogram
  
  Returns
  -------
  nArr : float array
    histogram-like distribution, with length of len(bins)-1
  """
  import scipy.interpolate as itp
  inter = itp.interp1d(xArr, cdf, bounds_error=False, fill_value=(0.0, 1.0))
  nArr  = inter(bins)
  nArr  = nArr[1:] - nArr[:-1]
  return nArr

def saveAsciiNOfZ(name, nArr, bins, verbose=True):
  """
  Save n(z) files under ASCII format.
  
  Parameters
  ----------
  name : string
    full name of the file
  nArr : float array
    histogram-like distribution
  bins : float array
    bin edges of the histogram
  verbose : bool, optional
    print verbose message
    
  Returns
  -------
  No returns
  """
  
  f = open(name, 'w')
  f.write('# hist\n')
  
  for z, n in zip(bins, nArr):
    f.write('%f  %e\n' % (z, n))
  f.close()
  
  if verbose:
    print('Saved \"%s\"' % name)
  return

def saveAsciiNOfZ_fromCdf(name, x_cdf, cdf):
  """
  Save n(z) file from a cdf.
  
  Parameters
  ----------
  name : string
    full name of the file
  x_cdf : float array
    sample points of the cumulative distribution
  cdf : float array
    values of the cumulative distribution
    
  Returns
  -------
  No returns
  """
  
  ## Define bin range & width
  bins = np.arange(0.0, 1.000001, 0.02)
  
  ## Convert cdf to histogram
  n    = cdfToHist(x_cdf, cdf, bins)
  
  ## Add a dummy zero at the end
  n    = np.insert(n, len(n), 0.0)
  
  ## Save
  saveAsciiNOfZ(name, n, bins, verbose=True)
  return

def saveAsciiNOfZ_demo():
  """
  Save two n(z) used in demo.
  """
  
  x_cdf = [0.0, 0.2, 0.3, 0.4, 0.5, 0.6]
  cdf0  = [0.0, 0.0, 0.5, 1.0, 1.0, 1.0]
  cdf1  = [0.0, 0.0, 0.0, 0.2, 0.6, 1.0]
  
  ## Save for type 0
  name = '%sinput/nOfZ_hist_type0.dat' % DEMO_PATH
  saveAsciiNOfZ_fromCdf(name, x_cdf, cdf0)
  
  ## Save for type 1
  name = '%sinput/nOfZ_hist_type1.dat' % DEMO_PATH
  saveAsciiNOfZ_fromCdf(name, x_cdf, cdf1)
  return

################################################################################
## Main

if __name__ == '__main__':
  saveFitsMask_demo()
  saveAsciiNOfZ_demo()

###############################################################################


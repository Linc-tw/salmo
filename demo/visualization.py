

    ##############################
    ##  visualization.py        ##
    ##  Chieh-An Lin            ##
    ##  Version 2020.07.21      ##
    ##############################


import numpy as np
import astropy.io.fits as fits
import healpy as hp
import matplotlib.pyplot as plt


################################################################################
## Global variables

DEMO_PATH = './'

###############################################################################
## Functions related to mask plotting

def loadFitsFullMap(name, verbose=True):
  full = fits.getdata(name, 1).field(0).flatten()
  if verbose == True:
    print('Loaded \"%s\"' % name)
  return full

def makeFullMapAxes(flip='geo', ctrRA=0.0):
  """
  Create an Axes object for a mollview plot
  
  Parameters
  ----------
    flip : 'geo' or 'astro'
    
    ctrRA : float
      RA value in degree on which the map is centered
      
  Returns
  -------
    ax : matplotlib Axes instance
  """
  fig = plt.gcf()
  ax  = hp.projaxes.HpxMollweideAxes(fig, [0, 0, 1, 1], coord='E', rot=(ctrRA, 0, 0), format='%g', flipconv=flip)
  fig.add_axes(ax)
  return ax

def plotFullMap(ax, full, nest=False, flip='geo', resol=1600, cmap='magma', vmin=None, vmax=None, alpha=1):
  transparent = (1, 1, 1, 0)
  full = hp.pixelfunc.ma_to_array(full)
  ax.projmap(full, nest=nest, xsize=resol, coord=None, cmap=cmap, badcolor=transparent, vmin=vmin, vmax=vmax, norm=None, alpha=alpha)
  return

def loadFitsRADEC(name, verbose=True):
  data  = fits.getdata(name, 1)
  RADEC = [data.field(0), data.field(1)]
  if verbose:
    print('Loaded \"%s\"' % name)
  return RADEC

def RADECToPatch(nside, RA, DEC=None, nest=False, lonlat=True, invRotMat=np.diag(np.ones(3))):
  if DEC is None:
    DEC = RA[1]
    RA  = RA[0]
  uXYZ  = hp.dir2vec(RA, DEC, lonlat=lonlat)
  uXYZ  = invRotMat.dot(uXYZ)
  patch = hp.vec2pix(nside, uXYZ[0], uXYZ[1], uXYZ[2], nest=nest)
  return patch

def plotFullCat(ax, RADEC, nside=512, nest=False, flip='geo', resol=1600, log=False, cmap=None, vmin=None, vmax=None):
  nbPix = 12 * nside * nside
  full  = np.zeros(nbPix, dtype=float)
  pix   = RADECToPatch(nside, RADEC[0], RADEC[1], nest=nest)
  for i in pix:
    full[i] += 1
  if log != False or log != 0:
    full = np.log10(full+float(log))
  plotFullMap(ax, full, nest=nest, flip=flip, resol=resol, cmap=cmap, vmin=vmin, vmax=vmax)
  return

def showFullComp(name1, name2, saveName, nside=512, nest=False, flip='geo', ctrRA=0.0, resol=1600, log=False, cmap='magma', vmin=None, vmax=None, verbose=True):
  fig = plt.gcf()
  fig.clf()
  ax1 = makeFullMapAxes(flip=flip, ctrRA=ctrRA)
  ax2 = makeFullMapAxes(flip=flip, ctrRA=ctrRA)
  
  full  = loadFitsFullMap(name1, verbose=verbose)
  RADEC = loadFitsRADEC(name2, verbose=verbose)
  
  ## Plot
  plotFullMap(ax1, full, nest=nest, flip=flip, resol=resol, cmap=cmap, vmin=vmin, vmax=vmax)
  plotFullCat(ax2, RADEC, nside=nside, nest=nest, flip=flip, resol=resol, log=log, cmap=cmap, vmin=vmin, vmax=vmax)
  
  ## Save
  fig.set_size_inches(6, 6)
  ax1.set_position([0.01, 0.5, 0.98, 0.5])
  ax2.set_position([0.01, 0, 0.98, 0.5])
  fig.patch.set_alpha(1.0)
  
  fig.savefig(saveName)
  if verbose:
    print('Saved \"%s\"' % saveName)
  return

###############################################################################
## Functions related to n(z) plotting

def loadAscii(name, sep=None, cmt='#', verbose=True):
  data = np.loadtxt(name, comments=cmt, delimiter=sep)
  if verbose:
    print('Loaded \"%s\"' % name)
  return data.T

def plotNOfZ_hist(ax, data):
  zArr = data[0]
  nArr = data[1]
  ax.step(zArr, nArr, where='post', color='#222288', ls='-')
  return

def loadFitsRedshift(name, verbose=True):
  data = fits.getdata(name, 1)
  if verbose:
    print('Loaded \"%s\"' % name)
  return data.field(2)

def centerOfBins(bins, area=False):
  bins  = np.array(bins, dtype=float)
  left  = bins[:-1]
  right = bins[1:]
  if area is True:
    return np.sqrt(0.5 * (left**2 + right**2))
  return 0.5 * (left + right)

def makeHist(data, bins, wgt=None, factor=1.0, pdf=False):
  """
  Make the histogram such that the output can be plotted directly
  
  Parameters
  ----------
  data : array-like
  bins : (1, N) float array
    bin edges
  factor : float, optional
    rescaling factor for the histogram
  pdf : bool, optional
    make the output a pdf, i.e. normalized by the binwidth & the total counts
  
  Returns
  -------
  n : (1, N) float array
    number counts, could be rescaled
  ctrBin : (1, N) float array
    center of the bins
    n & ctrBin have the same size.
  """
  nArr, bins = np.histogram(data, bins, weights=wgt)
  ctrBin     = centerOfBins(bins)
  if pdf == True:
    nArr = nArr.astype(float) / (float(sum(nArr)) * (bins[1:] - bins[:-1]))
  else:
    nArr = nArr.astype(float) * factor
  return nArr, ctrBin

def plotNOfZ_cat(ax, data):
  bins = np.arange(0.0, 1.000001, 0.02)
  nArr, ctrBins = makeHist(data, bins)
  nArr /= nArr.sum()
  ax.plot(ctrBins, nArr, color='#660022', marker='o', ls='none')
  return

def showNOfZComp(name1, name2, saveName, verbose=True):
  fig = plt.gcf()
  fig.clf()
  ax = fig.gca()
  
  ## Load
  data1 = loadAscii(name1)
  data2 = loadFitsRedshift(name2)
  
  ## Plot
  plotNOfZ_hist(ax, data1)
  plotNOfZ_cat(ax, data2)
  
  ## Legend
  hand = [plt.Line2D([], [], color='#222288', ls='-'), plt.Line2D([], [], color='#660022', ls='none', marker='o')]
  label = ['Input', 'Reconstructed']
  ax.legend(hand, label, loc=1, fontsize=16)
  
  ## Settings
  ax.set_ylim(-0.01, 0.112)
  ax.tick_params(axis='both', direction='in', labelsize=16)
  ax.set_xlabel('Redshift', size=16)
  ax.set_ylabel('Arbitrary amplitude', size=16)
  
  ## Save
  fig.set_size_inches(6, 6)
  fig.subplots_adjust(left=0.15, right=0.99, bottom=0.09, top=0.99)
  fig.savefig(saveName)
  if verbose:
    print('Saved \"%s\"' % saveName)
  return

################################################################################
## Main

def saveFigure():
  print()
  for i in range(2):
    name1 = '%sinput/mask_type%d.fits' % (DEMO_PATH, i)
    name2 = '%soutput/galCat_ref_type%d.fits' % (DEMO_PATH, i)
    saveName = 'mask_comp_type%d.png' % i
    showFullComp(name1, name2, saveName)
    print()
    
    name1 = '%sinput/nOfZ_hist_type%d.dat' % (DEMO_PATH, i)
    name2 = '%soutput/galCat_ref_type%d.fits' % (DEMO_PATH, i)
    saveName = 'nOfZ_comp_type%d.png' % i
    showNOfZComp(name1, name2, saveName)
    print()
  return

if __name__ == '__main__':
  saveFigure()

###############################################################################


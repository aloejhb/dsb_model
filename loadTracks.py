from __future__ import division
from os import getcwd,path,chdir,walk
from sys import argv
from subprocess import call
from matplotlib import rc
from matplotlib.pyplot import *

from numpy.random import choice
from scipy.stats import spearmanr
from scipy.interpolate import UnivariateSpline
from pylab import *
import pickle
import cPickle

rc('font', family='sans-serif')
rc('font', size=22)
rc('lines', linewidth = 2.)

collist = ['orange','maroon','royalblue','magenta']
on_clist = ['crimson','darkorange','maroon']
off_clist = ['royalblue','magenta','aqua']

    
# sampling frequency
sfreq = 4.

# header = """Frame  XPos  YPos  Foci  Radius Geminin_Int  53BP1_Int mFociArea"""    
# indices ->   0 ,   1,     2,   3,     4,     5,           6,       7


def Plot_FociGemRad(Track, sfreq = sfreq, startFrm = 1):

    frame = Track[:,0]
    if frame[0] < startFrm:
        Track2 = Track[int(startFrm - frame[0]):]
    else:
        Track2 = Track

    # extract observables from array
    tvec = (Track2[:,0]-startFrm)/sfreq #frames to time
    Rs = Track2[:,4] # rad
    Foci = Track2[:,3]
    Geminin = Track2[:,5]
    BP1Int = Track2[:,6]
    fig,axs = subplots(3,1,sharex = True, figsize = (9,9))    

    axs[0].set_ylabel('Number of Foci')

    axs[0].plot(tvec,Foci,'r-',)
    axs[0].set_ylim( (0,max(Foci) + 1) )

    axs[1].set_ylabel('Geminin Intensity a.u.', labelpad = 10)
    axs[1].plot(tvec,Geminin,'g-')
    axs[1].set_ylim(0,16000)


    ymin,ymax = 1.2 * max(Geminin), 0.9 * min(Geminin)
    
    axs[2].set_ylabel('Nuclear Radius (px)', labelpad = 10)
    axs[2].plot(tvec,Rs,'k-')  
    axs[2].set_xlabel('Time (h)' )
    return fig, axs
        

def Load_Tracks(fpath, fn = False):

    tlist = []

    fnames = list(walk(fpath + '/'))[0][2]
    fnames = sorted([f for f in fnames if '.txt' in f])
    print 'Loading..'

    for name in fnames:
        print fpath+name
        A = loadtxt(fpath + name) #numpy function to load txt into numpy.array
        tlist.append(A)

    if fn:
        return tlist, fnames
    else:
        return tlist


# example usage:

# TrackPath = '../data/G1G1/'
# Exp = 'Doxo/'
# track_list = Load_Tracks(TrackPath +'Doxo/')

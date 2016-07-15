from numpy import *
from pylab import *
from model_funcs import plotMnStd
from model_g1 import md_g1
from model_g2 import md_g2

from fit_g1_NCSoff_G1 import model as md_G1 
from fit_s_NCSoff_G1G2 import model as md_G1G2
from fit_g2_NCSoff_G2 import model as md_G2
from fit_g2_NCSoff_G2div import model as md_G2div

import copy

path = '../data/yH2AX/'
d1 = loadtxt(path + 'yH2AX_NCSon.txt')
d2 = loadtxt(path + 'yH2AX_NCSoff.txt')
# time mean(foci) std(foci) median(foci) Q1(foci) Q3(foci)

tvec = [t+5. for t in d1[:,0]]

ion()
figure()

plot(tvec, d1[:,1], 'o', label = 'NCSon', color = 'g')
plot(tvec, d2[:,1], 'o', label = 'NCSoff', color = 'r')
legend(frameon=False)

mds = {}

mds['NCSoff_G1'] = copy.deepcopy(md_G1)
mds['NCSoff_G1G2'] = copy.deepcopy(md_G1G2)
mds['NCSoff_G2'] = copy.deepcopy(md_G2)
mds['NCSoff_G2div'] = copy.deepcopy(md_G2div)

ncell = {'NCSoff_G1':5, 'NCSoff_G1G2':1, 'NCSoff_G2':3, 'NCSoff_G2div': 2}

yhind = {'NCSoff_G1':2, 'NCSoff_G1G2':4, 'NCSoff_G2':4, 'NCSoff_G2div': 4}


yhall = [] 
als = []
for k,m in mds.items():
    al,dum1,dum2 = m.multiSim(tvec, ncell[k])
    als.append(al)
    yhall.append(al[:, :, yhind[k]])

yhall = concatenate(yhall)
yh_mn = yhall.mean(axis = 0)
plot(tvec, yh_mn, color = 'm', label = 'off_sim')
legend(frameon=False)
show()

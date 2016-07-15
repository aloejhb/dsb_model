from numpy import *
from pylab import *
from model_funcs import plotMnStd

from model_g1 import md_g1
from model_s import md_s
from model_g2 import md_g2

from fit_g1_NCS_G1 import pars as p_G1 
from fit_g1_NCS_G1S import pars as p_G1S
from fit_s_NCS_SG2 import pars as p_SG2
from fit_g2_NCS_G2M import pars as p_G2M
from fit_g2_NCS_G2 import pars as p_G2

from fit_g1_NCSoff_G1 import pars as p_oG1 
from fit_g1_NCSoff_G1S import pars as p_oG1S
from fit_s_NCSoff_SG2 import pars as p_oSG2
from fit_g2_NCSoff_G2M import pars as p_oG2M
from fit_s_NCSoff_G2 import pars as p_oG2

import copy

path = '../data/yH2AX/'
d1 = loadtxt(path + 'yH2AX_NCSon.txt')
d2 = loadtxt(path + 'yH2AX_NCSoff.txt')
# time mean(foci) std(foci) median(foci) Q1(foci) Q3(foci)

tvec = [t+5. for t in d1[:,0]]

ion()
figure()

errorbar(tvec, d1[:,1], yerr=d1[:,2], label = 'MYCNon', color = 'g', fmt = 'o')
errorbar(tvec, d2[:,1], yerr=d2[:,2], label = 'MYCNoff', color = 'r', fmt = 'o')
legend(frameon=False)

mds = [copy.deepcopy(md_g1), copy.deepcopy(md_g1), copy.deepcopy(md_s), copy.deepcopy(md_g2), copy.deepcopy(md_g2)] 
mds_off = [copy.deepcopy(md_g1), copy.deepcopy(md_g1), copy.deepcopy(md_s), copy.deepcopy(md_g2), copy.deepcopy(md_s)] 

ps = [p_G1, p_G1S, p_SG2, p_G2M, p_G2]
ps_off = [p_oG1, p_oG1S, p_oSG2, p_oG2M, p_oG2]

for i,md in enumerate(mds):
    md.setAllPars(ps[i])


for i,md in enumerate(mds_off):
    md.setAllPars(ps_off[i])

ncell0 = [2, 2, 5, 3, 8]
ncell_off0 = [16, 2, 4, 3, 5]

'''
ncell = ncell0
ncell_off = ncell_off0
'''

ampf = 5
ncell = []
ncell_off = []
for nc in ncell0:
    ncell.append(nc * 5)

for nc in ncell_off0:
    ncell_off.append(nc * 5)

print ncell
print ncell_off

yhind = [2, 2, 4, 3, 3,]
yhind_off = [2, 2, 4, 3, 4]

etvec = linspace(tvec[0], tvec[-1], 200)

yhall = [] 
als = []
for i,md in enumerate(mds):
    al = md.multiSim(etvec, ncell[i])
    als.append(al)
    yhall.append(al[:, :, yhind[i]])

yhall = concatenate(yhall)
yh_mn = yhall.mean(axis = 0)
yh_st = yhall.std(axis = 0)

plotMnStd(etvec, yh_mn, yh_st, label = 'MYCNon sim', color = 'limegreen')

yhall = [] 
als = []
for i,md in enumerate(mds_off):
    al = md.multiSim(etvec, ncell_off[i])
    als.append(al)
    yhall.append(al[:, :, yhind_off[i]])

yhall = concatenate(yhall)
yh_mn = yhall.mean(axis = 0)
yh_st = yhall.std(axis = 0)

plotMnStd(etvec, yh_mn, yh_st, label = 'MYCNoff sim', color = 'm')

legend(frameon=False, prop={'size':13})
tick_params(labelsize=15)
xlabel('Time (h)', fontsize=16)
ylabel('# yH2AX Foci', fontsize=16)
timestmp = datetime.datetime.now().strftime("-%Y-%m-%d-%H-%M-%s")
savefig('../results/yH2AX_fit/yH2AX_bulk'+timestmp+'.png')
# show()

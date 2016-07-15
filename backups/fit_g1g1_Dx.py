from model_g1g1 import md_g1g1 as model
from data_funcs import *
from fit_funcs import *
from pylab import *
from mpl_toolkits.mplot3d import Axes3D

trackPath = '../data/'
exp = 'G1G1_Doxo'
sfreq = 4.
startFrm = 1
dtxind = 3

track_list, data = loadData(trackPath, exp, dtxind, sfreq, startFrm)


pars = {'D':2.1, 'b':0.45, 'alph':5.,'rNH':0.15, 'rNHup':0.3}
pars.update( {'beta':0.001, 'rHR':0.001} )
# pars.update( {'ttreat':1.5, 'ttreat_dl':0.15, 'twash':5., 'twash_dl':10., 'treg':80., 'treg_dl':2.} )
pars.update( {'ttreat':1.5, 'ttreat_dl':0.15, 'twash':23., 'twash_dl':7.3, 'treg':5.7, 'treg_dl':2.} )

model.setAllPars(pars)

mdxind = 1
trange = (0, 60)

# pltCompr2Data(model, data, mdxind, trange)
# pltCompr2Data(model, data, mdxind, trange, ssim = True, numSim = 50)

par_grid = {'twash':linspace(1., 40., 20), 'treg':linspace(1., 70., 20)}
pnms = par_grid.keys()
lscp_md = copy.deepcopy(model)
# lscp_md.pars['b_to_rNH'] = 3.
# xx, yy, obj = landscape(lscp_md, data, mdxind, par_grid, trange, svfig = True)

'''
figure()
ax = gca(projection = '3d')
ax.plot_surface(xx, yy, obj, cmap = cm.coolwarm)
show()
'''

# without upreg of repair
# fixpars = {'b':0.45, 'beta':0.001, 'rHR':0.001, 'ttreat':2., 'ttreat_dl':0.15, 'treg':98., 'treg_dl':1., 'rNHup': 0.3}
# fitpars = {'D':2., 'alph':4., 'rNH':0.15, 'twash':22., 'twash_dl':4.8}
# bd_dict = {'D':(0.01, 20.), 'alph':(0.1, 10.), 'rNH':(0.05, 1.), 'twash':(3., 60.), 'twash_dl':(1., 30.)}

# with upreg
fixpars = {'b':0.45, 'beta':0.001, 'rHR':0.001, 'ttreat':2., 'ttreat_dl':0.15}
fitpars = {'D':2., 'alph':4., 'rNH':0.15, 'rNHup': 0.3, 'twash':23., 'twash_dl':7.3, 'treg':5.7, 'treg_dl':2.}
bd_dict = {'D':(0.01, 20.), 'alph':(0.1, 10.), 'rNH':(0.05, 1.), 'rNHup':(0.05, 1.), 'twash':(10., 40.), 'twash_dl':(1., 30.), 'treg':(3., 20.), 'treg':(1., 30.)}

fit(model, data, mdxind, fitpars, fixpars, None, trange, 'Nelder-Mead')
# rndFit(model, data, mdxind, fitpars, fixpars, bd_dict, trange, fit_met, numRnd = 10)

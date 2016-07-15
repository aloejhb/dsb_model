from pfuncs_ssa import *
from model_funcs import *
from data_funcs import *
from fit_funcs import *
import copy

ion()

#dtname = 'G2G2_Doxo'
dtname = 'NCS_G2G1s'
print 'Working on ' + dtname + 'data.\n'

trackPath = '../data/'
exps = ('G1G1_Doxo', 'G1G2_Doxo', 'G2G2_Doxo', 'NCS_G1s', 'NCS_G2G1s', 'NCS_G2G2s')
sfreqs = (4., 4., 4, 6., 6., 6.)
startFrms = (1, 1, 1, 200, 200, 200)
dtxind = 3

track_lists = {}
data_list = {}
for i,exp in enumerate(exps):
    trl, dt= loadData(trackPath, exp, dtxind, sfreqs[i], startFrms[i])
    track_lists[exp] = trl
    data_list[exp] = dt

ttime = 100
tvec = linspace(0, ttime, sfreq*ttime + 1)

# pltData((data_list['NCS_G1s'], data_list['NCS_G2G1s'], data_list['NCS_G2G2s']))

mdname = 'g1g2_birepr'

def mk_pfuncs(p):
    if 'twash' not in p:
        p['twash'] = p['ttreat']+p['ttreat_dl']+p['gap_ttreat_twash']
    if 'tcc2' not in p:
        p['tcc2'] = p['tcc1']+p['tcc1_dl']+p['gap_tcc1_tcc2']

    pfuncs = []
    pfuncs.append(spline_pfunc((p['b'], p['D'], p['b']), (p['ttreat'], p['twash']), (p['ttreat_dl'], p['twash_dl']), mk_hfunc(()), ttime))

    pfuncs.append(spline_pfunc((p['alph'], p['alph2'], p['alph']), (p['tcc1'], p['tcc2']), (p['tcc1_dl'], p['tcc2_dl']), mk_hfunc((0,)), ttime))

    pfuncs.append(spline_pfunc((p['beta'], p['beta2'], p['beta']), (p['tcc1'], p['tcc2']), (p['tcc1_dl'], p['tcc2_dl']), mk_hfunc((0,)), ttime))

    pfuncs.append(spline_pfunc((p['rNH'],), (), (), mk_hfunc((1,)), ttime))

    pfuncs.append(spline_pfunc((p['rHR'],), (), (), mk_hfunc((2,)), ttime))

    pfuncs.append(spline_pfunc((p['conv'], p['conv2'], p['conv']), (p['tcc1'], p['tcc2']), (p['tcc1_dl'], p['tcc2_dl']), mk_hfunc((1,)), ttime))

    pfuncs.append(spline_pfunc((p['rHR'],), (), (), mk_hfunc((3,)), ttime))

    return pfuncs

x0 = [0, 1, 0, 0]
stoich = array(([1, 0, 0, 0], [-1, 1, 0, 0], [-1, 0, 1, 0], [0, -1, 0, 0], [0, 0, -1, 0], [0, -1, 0, 1], [0, 0, 0, -1]))

pars = {'D':6., 'b':1., 'alph':0.8, 'alph2':0.4, 'beta':0.01, 'beta2':0.8, 'rNH':0.5, 'rHR':0.4, 'conv':0.01, 'conv2':0.4, 'ttreat':1., 'ttreat_dl':0.5, 'twash':22., 'twash_dl':10., 'tcc1':3., 'tcc1_dl':5., 'tcc2':23., 'tcc2_dl':5.}
pfuncs = mk_pfuncs(pars)

xname = ['gammaH2AX', '53BP1+gammaH2AX', 'RAD51+gammaH2AX', 'RAD51+53BP1+gammaH2AX']
xcolor = ['dodgerblue', 'red', 'green', 'darkorange']
detxcolor = ['cornflowerblue', 'rosybrown', 'lime', 'saddlebrown']

obsinds = [(0,1,2,3), (1,3)]
obsname = ['total gammaH2AX', 'total 53BP1']
obscolor = ['grey', 'm']
detobscolor = ['black', 'blueviolet']

pname = ['DNA_dam', '53BP1_recrt', 'RAD51_recrt', 'NHEJ_repr', 'HR_repr', '53BP1_HRconv', '53BP1_HRrepr']
pcolor = ['blue', 'red', 'green', 'm', 'lime', 'darkorange', 'c']

fuzz = None
# fuzz = {'tcc2':('exponential', 10.)}

model = Model(mdname, x0, stoich, ttime, mk_pfuncs, pars, xname, xcolor, detxcolor, pname, pcolor, obsinds, obsname, obscolor, detobscolor, fuzz)

model_list = {}
md_pars = {}

model_list['Doxo'] = copy.deepcopy(model)

model_list['NCS_G2G2s'] = copy.deepcopy(model)
md_pars['NCS_G2G2s'] = {'D':5.9, 'b':0.45, 'alph':4., 'alph2':0.4, 'beta':0.01, 'beta2':0.2, 'rNH':0.5, 'rHR':0.4, 'conv':0.01, 'conv2':0.4, 'ttreat':5., 'ttreat_dl':0.5, 'twash':7., 'twash_dl':4., 'tcc1':1., 'tcc1_dl':0.5, 'tcc2':90., 'tcc2_dl':5.}
model_list['NCS_G2G2s'].setAllPars(md_pars['NCS_G2G2s'])

# model.pltSim()
# model.pltDetSim()
# model.pltMultiSim()
 
mdxind = 5
trange = (0, 60)

pltCompr2Data(model_list[dtname], data_list[dtname], mdxind, trange)
# pltCompr2Data(model_list[dtname], data_list[dtname], mdxind, trange, ssim = True, numSim = 50)

'''
fixpars = {'b':1.5, 'beta':0.01, 'ttreat':2., 'twash':22.}
# fitpars = {'D':3., 'alph':0.8, 'alph2':0.01, 'beta2':0.8, 'rNH':0.4, 'rHR':0.4, 'conv':0.01, 'conv2':0.2, 'ttreat_dl':5., 'gap_ttreat_twash':15., 'twash_dl':10., 'tcc1':10., 'tcc1_dl':10., 'gap_tcc1_tcc2':10., 'tcc2_dl':10.}
# bd_dict = {'D':(0.01, 20.), 'rNH':(0.5, 1.), 'ttreat_dl':(1., 30.), 'gap_ttreat_twash':(2., 30.), 'twash_dl':(1., 30.), 'tcc1':(1., 30.), 'tcc1_dl':(1., 15.), 'gap_tcc1_tcc2':(1., 30.), 'tcc2_dl':(1., 15.)}
fitpars = {'D':3., 'alph':0.8, 'alph2':0.01, 'beta2':0.8, 'rNH':0.4, 'rHR':0.4, 'conv':0.01, 'conv2':0.2, 'ttreat_dl':5., 'twash_dl':10., 'tcc1':10., 'tcc1_dl':10., 'gap_tcc1_tcc2':10., 'tcc2_dl':10.}
bd_dict = {'D':(fixpars['b'], 20.), 'rNH':(0.5, 1.), 'ttreat_dl':(1., fixpars['twash']-fixpars['ttreat']-0.1), 'twash_dl':(1., 30.), 'tcc1':(1., 30.), 'tcc1_dl':(1., 15.), 'gap_tcc1_tcc2':(1., 30.), 'tcc2_dl':(1., 15.)}

# fit(model, fit_tvec, fit_data_mn, xind, p0, pnames, fixpars, bounds, dtname)
# rndFit(model, fit_tvec, fit_data_mn, xind, pnames, fixpars, bounds, dtname, numRnd = 10)
'''

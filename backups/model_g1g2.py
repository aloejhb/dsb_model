from pfuncs_ssa import *
from model_funcs import *
from data_funcs import *
from fit_funcs import *

# ion()

trackPath = '../data/'
exps = ('G1G1/Doxo/','G1G2/Doxo/', 'G2G2/Doxo/')
dtxind = 3

sfreq = 4

track_lists, data = loadData(trackPath, exps, dtxind)

ttime = 100
tvec = linspace(0, ttime, sfreq*ttime + 1)

# pltData(tvec, data, ('G1G1_Doxo', 'G1G2_Doxo', 'G2G2_Doxo'), ('darkorange', 'red', 'lightgreen'))

mdname = 'g1g2'

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

    return pfuncs

x0 = [0, 0, 0]
stoich = array(([1, 0, 0], [-1, 1, 0], [-1, 0, 1], [0, -1, 0], [0, 0, -1], [0, -1, 1]))

pars = {'D':3., 'b':1.5, 'alph':0.8, 'alph2':0.01, 'beta':0.01, 'beta2':0.8, 'rNH':0.4, 'rHR':0.4, 'conv':0.01, 'conv2':0.1, 'ttreat':2., 'ttreat_dl':5., 'twash':22., 'twash':10., 'tcc1':10., 'tcc1_dl':10., 'tcc2':30., 'tcc2_dl':10.}
pfuncs = mk_pfuncs(pars)

xname = ['gammaH2AX', '53BP1+gammaH2AX', 'RAD51+gammaH2AX']
xcolor = ['dodgerblue', 'red', 'green']
detxcolor = ['cornflowerblue', 'rosybrown', 'lime']

obsinds = [(0,1,2,3)]
obsname = ['total gammaH2AX']
obscolor = ['grey']
detobscolor = ['black']

pname = ['DNA_dam', '53BP1_recrt', 'RAD51_recrt', 'NHEJ_repr', 'HR_repr', '53BP1_conv']
pcolor = ['blue', 'red', 'green', 'm', 'lime', 'darkorange']

fuzz = None
# fuzz = {'tcc2':('exponential', 10.)}

model = Model(mdname, x0, stoich, ttime, mk_pfuncs, pars, xname, xcolor, detxcolor, pname, pcolor, obsinds, obsname, obscolor, detobscolor, fuzz)


# model.pltSim()
# model.pltDetSim()
# model.pltMultiSim(numSim = 5)
 
dtname = 'G1G2_Doxo'
fit_xind = 1
fit_tvec = linspace(0, 60, sfreq*60 + 1)
fit_data_mn = data[dtname][1][:len(fit_tvec)]
fit_data_st = data[dtname][2][:len(fit_tvec)]

# model.pltCompr2Data(fit_tvec, fit_data_mn, fit_data_st, fit_xind, dtname)
# model.pltCompr2Data(fit_tvec, fit_data_mn, fit_data_st, fit_xind, dtname, ssim = True, numSim = 50)

fixpars = {'b':1.5, 'beta':0.01, 'ttreat':2.}
fitpars = {'D':3., 'alph':0.8, 'alph2':0.01, 'beta2':0.8, 'rNH':0.4, 'rHR':0.4, 'conv':0.01, 'conv2':0.2, 'ttreat_dl':5., 'gap_ttreat_twash':15., 'twash_dl':10., 'tcc1':10., 'tcc1_dl':10., 'gap_tcc1_tcc2':10., 'tcc2_dl':10.}
bd_dict = {'D':(0.01, 20.), 'rNH':(0.5, 1.), 'ttreat_dl':(1., 30.), 'gap_ttreat_twash':(2., 30.), 'twash_dl':(1., 30.), 'tcc1':(1., 30.), 'tcc1_dl':(1., 15.), 'gap_tcc1_tcc2':(1., 30.), 'tcc2_dl':(1., 15.)}

pnames = fitpars.keys()
p0 = fitpars.values()

bounds = [(0, 20.)] * len(pnames)
for k,v in bd_dict.items():
    bounds[pnames.index(k)] = v 
 

# fit(model, fit_tvec, fit_data_mn, xind, p0, pnames, fixpars, bounds, dtname)
rndFit(model, fit_tvec, fit_data_mn, xind, pnames, fixpars, bounds, dtname)

from pfuncs_ssa import *
from model_funcs import *

mdname = 's'
def mk_pfuncs(p):
    if 'twash' not in p:
        p['twash'] = p['ttreat']+p['ttreat_dl']+p['gap_ttreat_twash']

    pfuncs = []
    pfuncs.append(spline_pfunc((p['b'], p['D'], p['b']), (p['ttreat'], p['twash']), (p['ttreat_dl'], p['twash_dl']), mk_hfunc(()), ttime))

    pfuncs.append(spline_pfunc((p['alph'], p['alph2']), (p['ts'],), (p['ts_dl'],), mk_hfunc((0,)), ttime))

    pfuncs.append(spline_pfunc((p['beta'], p['beta2']), (p['ts'],), (p['ts_dl'],), mk_hfunc((0,)), ttime))

    pfuncs.append(spline_pfunc((p['rNH'],), (), (), mk_hfunc((1,)), ttime))

    pfuncs.append(spline_pfunc((p['rHR'],), (), (), mk_hfunc((2,)), ttime))

    pfuncs.append(spline_pfunc((p['conv'], p['conv2']), (p['ts'],), (p['ts_dl'],), mk_hfunc((1,)), ttime))

    pfuncs.append(spline_pfunc((p['rHR2'],), (), (), mk_hfunc((3,)), ttime))

    pfuncs.append(spline_pfunc((p['dp'], p['dp2'], p['dp']), (p['ts'], p['tg2']), (p['ts_dl'], p['tg2_dl']), mk_hfunc((1,)), ttime))

    pfuncs.append(spline_pfunc((p['dp_5r'], p['dp2_5r'], p['dp_5r']), (p['ts'], p['tg2']), (p['ts_dl'], p['tg2_dl']), mk_hfunc((3,)), ttime))

    return pfuncs

ttime = 100.
x0 = [0, 1, 0, 0]
stoich = array(([1, 0, 0, 0], [-1, 1, 0, 0], [-1, 0, 1, 0], [0, -1, 0, 0], [0, 0, -1, 0], [0, -1, 0, 1], [0, 0, 0, -1], [1, -1, 0, 0], [1, 0, 0, -1]))
pars = {'D':50., 'b':1., 'alph':0.4, 'alph2':0.4, 'beta':0, 'beta2':0.1, 'rNH':0.5, 'rHR':0.4, 'rHR2':0.4, 'conv':0, 'conv2':0.4}
pars.update({'dp':0, 'dp2':0.8, 'dp_5r':0, 'dp2_5r':0.8})
pars.update({'ttreat':5., 'ttreat_dl':0.05, 'twash':5.6, 'twash_dl':0.05, 'ts':7., 'ts_dl':0.2, 'tg2':9., 'tg2_dl':0.2})

pfuncs = mk_pfuncs(pars)

xname = ['gammaH2AX', '53BP1+gammaH2AX', 'RAD51+gammaH2AX', 'RAD51+53BP1+gammaH2AX']
xcolor = ['dodgerblue', 'red', 'green', 'darkorange']
detxcolor = ['cornflowerblue', 'rosybrown', 'lime', 'saddlebrown']

obsinds = [(0,1,2,3), (1,3)]
obsname = ['total gammaH2AX', 'total 53BP1']
obscolor = ['grey', 'm']
detobscolor = ['black', 'blueviolet']

pname = ['DNA_dam', '53BP1_recrt', 'RAD51_recrt', 'NHEJ_repr', 'HR_repr', '53BP1_HRconv', '53BP1_HRrepr', '53BP1_deplt', '5R_deplt']
pcolor = ['blue', 'red', 'green', 'm', 'lime', 'darkorange', 'c', 'gold', 'saddlebrown']

fuzz = None
# fuzz = {'tcc2':('exponential', 10.)}

md_s = Model(mdname, x0, stoich, ttime, mk_pfuncs, pars, xname, xcolor, detxcolor, pname, pcolor, obsinds, obsname, obscolor, detobscolor, fuzz)


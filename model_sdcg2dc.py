from pfuncs_ssa import *
from model_funcs import *

mdname = 'sdcg2dc'
def mk_pfuncs(p):
    pfuncs = []

    pfuncs.append(spline_pfunc((p['b'], p['D'], p['b']), (p['ttreat'], p['twash']), (p['ttreat_dl'], p['twash_dl']), mk_hfunc(()), ttime))

    pfuncs.append(spline_pfunc((p['alph'], p['alph_s'], p['alph']), (p['ts'], p['tg2']), (p['ts_dl'], p['tg2_dl']), mk_hfunc((0,)), ttime))

    pfuncs.append(spline_pfunc((p['rNH'],), (), (), mk_hfunc((1,)), ttime))

    pfuncs.append(spline_pfunc((0, p['beta']), (p['ts'],), (p['ts_dl'],), mk_hfunc((0,)), ttime))

    pfuncs.append(spline_pfunc((0, p['rHR']), (p['ts'],), (p['ts_dl'],), mk_hfunc((2,)), ttime))

    return pfuncs

ttime = 100.
x0 = [0, 1, 0]
stoich = array(([1, 0, 0], [-1, 1, 0], [0, -1, 0], [-1, 0, 1], [0, 0, -1]))
pars = {'D':20., 'alph':0.4, 'alph_s':0.4, 'rNH':0.5, 'beta':0.1, 'rHR':0.4}
pars.update({'b':0.45})
pars.update({'ttreat':5., 'ttreat_dl':0.05, 'twash':6., 'twash_dl':0.05, 'ts':7., 'ts_dl':0.2, 'tg2':9., 'tg2_dl':0.2})

pfuncs = mk_pfuncs(pars)

xname = ['gammaH2AX', '53BP1+gammaH2AX', 'RAD51+gammaH2AX']
xcolor = ['dodgerblue', 'red', 'green']
detxcolor = ['cornflowerblue', 'rosybrown', 'lime']

obsinds = [(0,1,2)]
obsname = ['total gammaH2AX']
obscolor = ['grey']
detobscolor = ['black']

pname = ['DNA_dam', '53BP1_recrt', 'NHEJ_repr', 'RAD51_recrt', 'HR_repr']
pcolor = ['blue', 'red', 'm', 'green', 'lime']

fuzz = None
# fuzz = {'tcc2':('exponential', 10.)}

md_sdcg2dc = Model(mdname, x0, stoich, ttime, mk_pfuncs, pars, xname, xcolor, detxcolor, pname, pcolor, obsinds, obsname, obscolor, detobscolor, fuzz)


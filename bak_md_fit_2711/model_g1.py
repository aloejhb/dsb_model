from pfuncs_ssa import *
from model_funcs import *

ttime = 100

mdname = 'g1'
def mk_pfuncs(p):
    if 'twash' not in p:
        p['twash'] = p['ttreat']+p['ttreat_dl']+p['gap_ttreat_twash']

    if 'b_to_rNH' in p:
        p['b'] = p['rNH'] * p['b_to_rNH']

    pfuncs = []
    pfuncs.append(spline_pfunc((p['b'], p['D'], p['b']), (p['ttreat'], p['twash']), (p['ttreat_dl'], p['twash_dl']), mk_hfunc(()), ttime))

    pfuncs.append(spline_pfunc((p['alph'],), (), (), mk_hfunc((0,)), ttime))

    pfuncs.append(spline_pfunc((p['rNH'], p['rNHup'],), (p['treg'],), (p['treg_dl'],), mk_hfunc((1,)), ttime))

    return pfuncs

x0 = [0, 1]
stoich = array(([1, 0], [-1, 1], [0, -1]))

pars = {'D':1.16, 'b':0.45, 'alph':0.6, 'rNH':0.0001, 'rNHup':0.2, 'ttreat':2., 'ttreat_dl':0.15, 'twash':22., 'twash_dl':7.3, 'treg':4., 'treg_dl':7.2}

xname = ['gammaH2AX', '53BP1+gammaH2AX']
xcolor = ['dodgerblue', 'red']
detxcolor = ['cornflowerblue', 'rosybrown']

obsinds = [(0,1)]
obsname = ['total gammaH2AX']
obscolor = ['grey']
detobscolor = ['black']

pname = ['DNA_dam', '53BP1_recrt', 'NHEJ_repr']
pcolor = ['blue', 'red', 'm']

fuzz = None
# fuzz = {'twash':('exponential', 10.)}

md_g1 = Model(mdname, x0, stoich, ttime, mk_pfuncs, pars, xname, xcolor, detxcolor, pname, pcolor, obsinds, obsname, obscolor, detobscolor, fuzz)

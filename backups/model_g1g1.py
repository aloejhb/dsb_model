from pfuncs_ssa import *
from model_funcs import *

ttime = 100

mdname = 'g1g1'
def mk_pfuncs(p):
    if 'twash' not in p:
        p['twash'] = p['ttreat']+p['ttreat_dl']+p['gap_ttreat_twash']

    if 'b_to_rNH' in p:
        p['b'] = p['rNH'] * p['b_to_rNH']

    pfuncs = []
    pfuncs.append(spline_pfunc((p['b'], p['D'], p['b']), (p['ttreat'], p['twash']), (p['ttreat_dl'], p['twash_dl']), mk_hfunc(()), ttime))

    pfuncs.append(spline_pfunc((p['alph'],), (), (), mk_hfunc((0,)), ttime))

    pfuncs.append(spline_pfunc((p['beta'],), (), (), mk_hfunc((0,)), ttime))

    pfuncs.append(spline_pfunc((p['rNH'], p['rNHup'],), (p['treg'],), (p['treg_dl'],), mk_hfunc((1,)), ttime))

    pfuncs.append(spline_pfunc((p['rHR'],), (), (), mk_hfunc((2,)), ttime))

    return pfuncs

mdname2 = 'g1g1_step'
def mk_pfuncs2(p):
    if 'twash' not in p:
        p['twash'] = p['ttreat']+p['ttreat_dl']+p['gap_ttreat_twash']

    if 'b_to_rNH' in p:
        p['b'] = p['rNH'] * p['b_to_rNH']

    pfuncs = []
    pfuncs.append(step_pfunc((p['b'], p['D'], p['b']), (p['ttreat'], p['twash']), mk_hfunc(()) ))

    pfuncs.append(step_pfunc((p['alph'],), (), mk_hfunc((0,)) ))

    pfuncs.append(step_pfunc((p['beta'],), (), mk_hfunc((0,)) ))

    pfuncs.append(step_pfunc((p['rNH'], p['rNHup'],), (p['treg'],), mk_hfunc((1,)) ))

    pfuncs.append(step_pfunc((p['rHR'],), (), mk_hfunc((2,)) ))

    return pfuncs


x0 = [0, 1, 0]
stoich = array(([1, 0, 0], [-1, 1, 0], [-1, 0, 1], [0, -1, 0], [0, 0, -1]))

pars = {'D':1.16, 'b':0.45, 'alph':0.6, 'beta':0.001, 'rNH':0.0001, 'rNHup':0.2, 'rHR':0.001, 'ttreat':2., 'ttreat_dl':0.15, 'twash':22., 'twash_dl':7.3, 'treg':4., 'treg_dl':7.2}

xname = ['gammaH2AX', '53BP1+gammaH2AX', 'RAD51+gammaH2AX']
xcolor = ['dodgerblue', 'red', 'green']
detxcolor = ['cornflowerblue', 'rosybrown', 'lime']

obsinds = [(0,1,2)]
obsname = ['total gammaH2AX']
obscolor = ['grey']
detobscolor = ['black']

pname = ['DNA_dam', '53BP1_recrt', 'RAD51_recrt', 'NHEJ_repr', 'HR_repr']
pcolor = ['blue', 'red', 'green', 'm', 'lime']

fuzz = None
# fuzz = {'twash':('exponential', 10.)}

md_g1g1 = Model(mdname, x0, stoich, ttime, mk_pfuncs, pars, xname, xcolor, detxcolor, pname, pcolor, obsinds, obsname, obscolor, detobscolor, fuzz)

md_g1g1_step = Model(mdname2, x0, stoich, ttime, mk_pfuncs2, pars, xname, xcolor, detxcolor, pname, pcolor, obsinds, obsname, obscolor, detobscolor, fuzz, 'step')

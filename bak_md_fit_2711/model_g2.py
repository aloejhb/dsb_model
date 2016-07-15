from pfuncs_ssa import *
from model_funcs import *

mdname = 'g2'
def mk_pfuncs(p):
    if 'twash' not in p:
        p['twash'] = p['ttreat']+p['ttreat_dl']+p['gap_ttreat_twash']

    pfuncs = []
    pfuncs.append(spline_pfunc((p['b'], p['D'], p['b']), (p['ttreat'], p['twash']), (p['ttreat_dl'], p['twash_dl']), mk_hfunc(()), ttime))

    pfuncs.append(spline_pfunc((p['alph2'], ), (), (), mk_hfunc((0,)), ttime))

    pfuncs.append(spline_pfunc((p['beta2'], ), (), (), mk_hfunc((0,)), ttime))

    pfuncs.append(spline_pfunc((p['rNH'],), (), (), mk_hfunc((1,)), ttime))

    pfuncs.append(spline_pfunc((p['rHR'],), (), (), mk_hfunc((2,)), ttime))

    pfuncs.append(spline_pfunc((p['conv2'], ), (), (), mk_hfunc((1,)), ttime))

    pfuncs.append(spline_pfunc((p['rHR2'],), (), (), mk_hfunc((3,)), ttime))

    return pfuncs

mdname2 = 'g2_step'
def mk_pfuncs2(p):
    if 'twash' not in p:
        p['twash'] = p['ttreat']+p['ttreat_dl']+p['gap_ttreat_twash']

    pfuncs = []
    pfuncs.append(step_pfunc((p['b'], p['D'], p['b']), (p['ttreat'], p['twash']), mk_hfunc(()) ))

    pfuncs.append(step_pfunc((p['alph2'],), (), mk_hfunc((0,)) ))

    pfuncs.append(step_pfunc((p['beta2'],), (), mk_hfunc((0,)) ))

    pfuncs.append(step_pfunc((p['rNH'],), (), mk_hfunc((1,)) ))

    pfuncs.append(step_pfunc((p['rHR'],), (), mk_hfunc((2,)) ))

    pfuncs.append(step_pfunc((p['conv2'],), (), mk_hfunc((1,)) ))

    pfuncs.append(step_pfunc((p['rHR2'],), (), mk_hfunc((3,)) ))

    return pfuncs
 
ttime = 100.
x0 = [0, 1, 0, 0]
stoich = array(([1, 0, 0, 0], [-1, 1, 0, 0], [-1, 0, 1, 0], [0, -1, 0, 0], [0, 0, -1, 0], [0, -1, 0, 1], [0, 0, 0, -1]))

pars = {'D':6., 'b':1., 'alph2':0.4, 'beta2':0.8, 'rNH':0.5, 'rHR':0.4, 'rHR2':0.4, 'conv2':0.4, 'ttreat':1., 'ttreat_dl':0.5, 'twash':22., 'twash_dl':10.}
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

md_g2 = Model(mdname, x0, stoich, ttime, mk_pfuncs, pars, xname, xcolor, detxcolor, pname, pcolor, obsinds, obsname, obscolor, detobscolor, fuzz)

md_g2_step = Model(mdname2, x0, stoich, ttime, mk_pfuncs2, pars, xname, xcolor, detxcolor, pname, pcolor, obsinds, obsname, obscolor, detobscolor, fuzz, 'step')

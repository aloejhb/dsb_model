from pfuncs_ssa import *
from model_funcs import *

mdname = 'g2'
def mk_pfuncs(p):
    pfuncs = []

    pfuncs.append(spline_pfunc((p['b'], p['D'], p['b']), (p['ttreat'], p['twash']), (p['ttreat_dl'], p['twash_dl']), mk_hfunc(()), ttime))

    pfuncs.append(spline_pfunc((p['alph'],), (), (), mk_hfunc((0,)), ttime))

    pfuncs.append(spline_pfunc((p['rNH'],), (), (), mk_hfunc((1,)), ttime))

    pfuncs.append(spline_pfunc((p['conv'],), (), (), mk_hfunc((1,)), ttime))

    pfuncs.append(spline_pfunc((p['rHR2'],), (), (), mk_hfunc((2,)), ttime))


    return pfuncs

ttime = 100.
x0 = [0, 1, 0]
stoich = array(([1, 0, 0], [-1, 1, 0], [0, -1, 0], [0, -1, 1], [0, 0, -1]))
pars = {'D':20., 'alph':0.4, 'rNH':0.5, 'conv':0.4, 'rHR2':0.4}
pars.update({'b':0.45})
pars.update({'ttreat':5., 'ttreat_dl':0.05, 'twash':6., 'twash_dl':0.05})

pfuncs = mk_pfuncs(pars)

xname = ['gammaH2AX', '53BP1+gammaH2AX', 'RAD51+53BP1+gammaH2AX']
xcolor = ['dodgerblue', 'red', 'darkorange']
detxcolor = ['cornflowerblue', 'rosybrown', 'saddlebrown']

obsinds = [(0,1,2), (1,2)]
obsname = ['total gammaH2AX', 'total 53BP1']
obscolor = ['grey', 'm']
detobscolor = ['black', 'blueviolet']

pname = ['DNA_dam', '53BP1_recrt', 'NHEJ_repr', '53BP1_HRconv', '53BP1_HRrepr']
pcolor = ['blue', 'red', 'm', 'c', 'gold']

fuzz = None
# fuzz = {'tcc2':('exponential', 10.)}

md_g2 = Model(mdname, x0, stoich, ttime, mk_pfuncs, pars, xname, xcolor, detxcolor, pname, pcolor, obsinds, obsname, obscolor, detobscolor, fuzz)


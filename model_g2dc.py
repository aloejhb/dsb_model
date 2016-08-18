from pfuncs_ssa import *
from model_funcs import *

mdname = 'g2dc'
def mk_pfuncs(p):
    pfuncs = []

    pfuncs.append(spline_pfunc((p['b'], p['D'], p['b']), (p['ttreat'], p['twash']), (p['ttreat_dl'], p['twash_dl']), mk_hfunc(()), ttime))

    pfuncs.append(spline_pfunc((p['alph'],), (), (), mk_hfunc((0,)), ttime))

    pfuncs.append(spline_pfunc((p['rNH'],), (), (), mk_hfunc((1,)), ttime))

    pfuncs.append(spline_pfunc((p['beta'],), (), (), mk_hfunc((0,)), ttime))

    pfuncs.append(spline_pfunc((p['rHR'],), (), (), mk_hfunc((2,)), ttime))


    return pfuncs

ttime = 100.
x0 = [0, 1, 0]
stoich = array(([1, 0, 0], [-1, 1, 0], [0, -1, 0], [-1, 0, 1], [0, 0, -1]))
pars = {'D':20., 'alph':0.4, 'rNH':0.5, 'beta':0.4, 'rHR':0.4}
pars.update({'b':0.45})
pars.update({'ttreat':5., 'ttreat_dl':0.05, 'twash':6., 'twash_dl':0.05})

pfuncs = mk_pfuncs(pars)

xname = ['$\gamma$H2AX', '53BP1+$\gamma$H2AX', 'RAD51+$\gamma$H2AX']
xcolor = ['dodgerblue', 'red', 'green']
detxcolor = ['cornflowerblue', 'rosybrown', 'lime']

obsinds = [(0,1,2)]
obsname = ['Total $\gamma$H2AX']
obscolor = ['grey']
detobscolor = ['black']

pname = ['DNA_dam', '53BP1_recrt', 'NHEJ_repr', 'RAD51_recrt', 'HR_repr']
pcolor = ['blue', 'red', 'm', 'green', 'gold']

fuzz = None
# fuzz = {'tcc2':('exponential', 10.)}

md_g2dc = Model(mdname, x0, stoich, ttime, mk_pfuncs, pars, xname, xcolor, detxcolor, pname, pcolor, obsinds, obsname, obscolor, detobscolor, fuzz)

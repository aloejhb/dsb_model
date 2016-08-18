from pfuncs_ssa import *
from model_funcs import *

mdname = 'sphase'
def mk_pfuncs(p):
    pfuncs = []

    pfuncs.append(spline_pfunc((p['b'],), (), (), mk_hfunc(()), ttime))

    pfuncs.append(spline_pfunc((p['alph'],), (), (), mk_hfunc((0,)), ttime))

    pfuncs.append(spline_pfunc((p['rNH'],), (), (), mk_hfunc((1,)), ttime))

    pfuncs.append(spline_pfunc((p['beta'],), (), (), mk_hfunc((0,)), ttime))

    pfuncs.append(spline_pfunc((p['rHR'],), (), (), mk_hfunc((2,)), ttime))


    return pfuncs

ttime = 100.
x0 = [17, 17, 0]
stoich = array(([1, 0, 0], [-1, 1, 0], [0, -1, 0], [-1, 0, 1], [0, 0, -1]))
pars = {'alph':0.4, 'rNH':0.5, 'beta':0.4, 'rHR':0.4}
pars.update({'b':0.45})

pfuncs = mk_pfuncs(pars)

xname = ['gammaH2AX', '53BP1+gammaH2AX', 'RAD51+gammaH2AX']
xcolor = ['dodgerblue', 'red', 'green']
detxcolor = ['cornflowerblue', 'rosybrown', 'lime']

obsinds = [(0,1,2)]
obsname = ['total gammaH2AX']
obscolor = ['grey']
detobscolor = ['black']

pname = ['DNA_dam', '53BP1_recrt', 'NHEJ_repr', 'RAD51_recrt', 'HR_repr']
pcolor = ['blue', 'red', 'm', 'green', 'gold']

fuzz = None
# fuzz = {'tcc2':('exponential', 10.)}

md_sphase = Model(mdname, x0, stoich, ttime, mk_pfuncs, pars, xname, xcolor, detxcolor, pname, pcolor, obsinds, obsname, obscolor, detobscolor, fuzz)

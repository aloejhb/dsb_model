from pfuncs_ssa import *
from model_funcs import *

mdname = 'sdcg2'
def mk_pfuncs(p):
    pfuncs = []

    pfuncs.append(spline_pfunc((p['b'], p['D'], p['b']), (p['ttreat'], p['twash']), (p['ttreat_dl'], p['twash_dl']), mk_hfunc(()), ttime))

    pfuncs.append(spline_pfunc((p['alph'], p['alph_s'], p['alph']), (p['ts'], p['tg2']), (p['ts_dl'], p['tg2_dl']), mk_hfunc((0,)), ttime))

    pfuncs.append(spline_pfunc((p['rNH'],), (), (), mk_hfunc((1,)), ttime))

    pfuncs.append(spline_pfunc((0, p['beta'], 0), (p['ts'], p['tg2']), (p['ts_dl'], p['tg2_dl']), mk_hfunc((0,)), ttime))

    pfuncs.append(spline_pfunc((0, p['rHR']), (p['ts'],), (p['ts_dl'],), mk_hfunc((2,)), ttime))

    pfuncs.append(spline_pfunc((0, p['conv']), (p['tg2'],), (p['tg2_dl'],), mk_hfunc((1,)), ttime))

    pfuncs.append(spline_pfunc((0, p['rHR2']), (p['tg2'],), (p['tg2_dl'],), mk_hfunc((3,)), ttime))


    return pfuncs

ttime = 100.
x0 = [0, 1, 0, 0]
stoich = array(([1, 0, 0, 0], [-1, 1, 0, 0], [0, -1, 0, 0], [-1, 0, 1, 0], [0, 0, -1, 0], [0, -1, 0, 1], [0, 0, 0, -1]))
pars = {'D':20., 'alph':0.4, 'alph_s':0.4, 'rNH':0.5, 'beta':0.1, 'rHR':0.4, 'conv':0.4, 'rHR2':0.4}
pars.update({'b':0.45})
pars.update({'ttreat':5., 'ttreat_dl':0.05, 'twash':6., 'twash_dl':0.05, 'ts':7., 'ts_dl':0.2, 'tg2':9., 'tg2_dl':0.2})

pfuncs = mk_pfuncs(pars)

xname = ['$\gamma$H2AX', '53BP1+$\gamma$H2AX', 'RAD51+$\gamma$H2AX', 'RAD51+53BP1+$\gamma$H2AX']
xcolor = ['dodgerblue', 'red', 'green', 'darkorange']
detxcolor = ['cornflowerblue', 'rosybrown', 'lime', 'saddlebrown']

obsinds = [(0,1,2,3), (1,3)]
obsname = ['total $\gamma$H2AX', 'total 53BP1']
obscolor = ['grey', 'm']
detobscolor = ['black', 'blueviolet']

pname = ['DNA_dam', '53BP1_recrt', 'NHEJ_repr', 'RAD51_recrt', 'HR_repr', '53BP1_HRconv', '53BP1_HRrepr']
pcolor = ['blue', 'red', 'm', 'green', 'lime', 'c', 'gold']

fuzz = None
# fuzz = {'tcc2':('exponential', 10.)}

md_sdcg2 = Model(mdname, x0, stoich, ttime, mk_pfuncs, pars, xname, xcolor, detxcolor, pname, pcolor, obsinds, obsname, obscolor, detobscolor, fuzz)


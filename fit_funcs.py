from model_funcs import *
from matplotlib.pyplot import savefig
from matplotlib import cm
from scipy.optimize import minimize
from numpy.random import uniform
from scipy.stats import scoreatpercentile as perc
import os
import copy
import re

def sumSqrRes(data_mn, md_mn):
    sqrvec = []
    for i in range(len(data_mn)):
        sqr = pow(data_mn[i] - md_mn[i], 2)
        sqrvec.append(sqr)

    return sum(sqrvec)

def objFunc(p0, model, tvec, data_mn, mdxind, pnames, fixpars):
    fitpars = dict(zip(pnames, p0))
    fitpars.update(fixpars)
    md = copy.deepcopy(model)
    try:
        md.setAllPars(fitpars)
    except Exception as e:
        print 'Warning: parameter values are invalid!'
        print str(e) + '\n'
        return 1000000
    res = md.detSim(tvec)
    sumsqr = sumSqrRes(data_mn, res[:,mdxind])

    print fitpars
    print sumsqr
    print ''

    return sumsqr

def fit(model, data, mdxind, fitpars, fixpars, bd_dict = None, trange = None, fit_met = 'L-BFGS-B', path = '../results/fit/'):

    pnames = fitpars.keys()
    p0 = fitpars.values()

    if bd_dict:
        bounds = [(0, 10.)] * len(pnames)
        for k,v in bd_dict.items():
            bounds[pnames.index(k)] = v 
    else:
        bounds = None
    

    start_time = datetime.datetime.now()

    if  trange:
        sel = logical_and(data.tvec>=trange[0], data.tvec<=trange[1])
        tvec = data.tvec[sel]
        data_mn = data.mn[sel]
    else:
        tvec = data.tvec
        data_mn = data.mn

    if fit_met == 'Nelder-Mead':
        opt = {'xtol': 1e-4, 'ftol': 1e-2, 'maxfev': 10000, 'disp':True}
        bounds = None
    else:
        opt = {'disp':True, 'factr':0.1}

    res = minimize(objFunc, p0, args = (model, tvec, data_mn, mdxind, pnames, fixpars), bounds = bounds, method = fit_met, options = opt)

    end_time = datetime.datetime.now()

    initpars = dict(zip(pnames, p0))
    respars = dict(zip(pnames, res.x))
    respars.update(fixpars)

    # write to file
    filename = path + '/' + model.name + '_vs_' + data.name + start_time.strftime("-%Y-%m-%d-%H-%M-%s")
    f = open(filename + '.txt', 'w')
    f.write('Fit result:\n')
    f.write(str(res) + '\n')
    f.write('Fixed pars:\n')
    f.write(str(fixpars) + '\n')
    f.write('Initial pars:\n')
    f.write(str(initpars) + '\n')
    f.write('Result pars:\n')
    f.write(str(respars) + '\n')
    for k in sorted(respars.keys()):
        f.write('{0}\t{1:.4f}\n'.format(k,respars[k]))
    f.write('Time consumed: ' + str(end_time - start_time))
    f.close()

    model.setAllPars(respars)
    pltCompr2Data(model, data, mdxind, trange = trange)
    savefig(filename + '.pdf')

    return res

def rndFit(model, data, mdxind, fitpars, fixpars, bd_dict, trange = None, fit_met = 'L-BFGS-B', numRnd = 2):
    start_time = datetime.datetime.now()
    path = '../results/rndfit/' + model.name + '_vs_' + data.name + start_time.strftime("-%Y-%m-%d-%H-%M-%s/")
    if not os.path.exists(path):
        os.makedirs(path)

    bds = dict.fromkeys(fitpars.keys(), (0, 10.))
    bds.update(bd_dict)

    for nr in range(numRnd):
        rnd_fitpars = fitpars.copy()
        for pn in rnd_fitpars.keys():
            bd = bds[pn]
            rnd_fitpars[pn] = uniform(bd[0], bd[1])

        fit(model, data, mdxind, rnd_fitpars, fixpars, bds, trange, fit_met, path)

def pltCompr2Data(model, data, mdxind, trange = None, ssim = False, numSim = 200, svfig = False, opath = '../results/manual_fit/'):
    
    if trange:
        sel = logical_and(data.tvec>=trange[0], data.tvec<=trange[1])
        tvec = data.tvec[sel]
        data_mn = data.mn[sel]
        # data_st = data.st[sel]
        data_q1 = data.q1[sel]
        data_q3 = data.q3[sel]
    else:
        # tvec, data_mn, data_st = data.tvec, data.mn, data.st
        tvec, data_mn, data_q1, data_q3 = data.tvec, data.q1, data.q3

    # figname = model.name + '_vs_' + data.name
    figname = data.name
    fig, ax1, ax2 = setNewfig(figname)
    sca(ax1)
    ax1.set_ylabel( '# 53BP1 Foci' , fontsize = 16)
    # plotMnStd(tvec, data_mn, data_st, label = 'data', color = 'royalblue')
    plotPerc(tvec, data_mn, data_q1, data_q3, label = 'data', color = 'royalblue')

    if ssim:
        xall = model.multiSim(tvec, numSim = numSim)
        xall2 = xall[:,:,mdxind]
        x_mn = mean(xall2, axis = 0) 

        q1 = []
        q3 = []
        for cl in xall2.T:
            cl2 = [c for c in cl if not isnan(c)]
            q1.append(perc(cl2, 25))
            q3.append(perc(cl2, 75))
        q1 = array(q1)
        q3 = array(q3)

        # plotMnStd(tvec, x_mn[:,mdxind], x_st[:,mdxind], label = 'stoich. sim.', color = 'darkorange')
        plotPerc(tvec, x_mn, q1, q3, label = 'stoich. sim.', color = 'darkorange')

    detsol = model.detSim(tvec)
    plot(tvec, detsol[:,mdxind],'--', color = 'saddlebrown',lw = 1.5, label = 'det. sim.')

    legend(frameon = False, prop = {'size':12})
    if re.match('.*off.*', data.name):
        ylim(0,20)
    else:
        ylim(0,30)

    sca(ax2)
    model.pltCfuncs(tvec)

    if ssim:
        fitness = sumSqrRes(data_mn, x_mn)
    else:
        fitness =  sumSqrRes(data_mn, detsol[:,mdxind])
    
    if svfig:
        # timestmp = datetime.datetime.now().strftime("-%Y-%m-%d-%H-%M-%s")
        timestmp = ''
        savefig(opath + figname + timestmp + '.png')
        f = open(opath + figname + timestmp + '.txt', 'w')
        f.write('Parameters:\n')
        f.write(model.printPars(tofile = True))
        f.write('\nFitness:\n')
        f.write(str(fitness))
        

    return fitness

def landscape(model, data, mdxind, par_grid, trange = None, svfig = False, svpath = '../results/landscape/'):


    if  trange:
        sel = logical_and(data.tvec>=trange[0], data.tvec<=trange[1])
        tvec = data.tvec[sel]
        data_mn = data.mn[sel]
    else:
        tvec = data.tvec
        data_mn = data.mn

    pnms = [pg[0] for pg in par_grid]
    grids = [pg[1] for pg in par_grid]

    fixpars = model.pars.copy()
    for pn in pnms:
        del fixpars[pn]

    figure()

    if len(par_grid) == 1:

        xx = grids[0]
        yy = None
        obj = zeros(xx.shape)
        for i in range(len(xx)):
            obj[i] = objFunc((xx[i],), model, tvec, data_mn, mdxind, pnms, fixpars)

        figttl = '{}_vs_{}_{}'.format(model.name, data.name, pnms[0])
        title(figttl, fontsize = 18)
        xlabel(pnms[0])
        ylabel('obj. func.')

        plot(xx, obj)

    elif len(par_grid) == 2:

        xx, yy = meshgrid(grids[0], grids[1])
        obj = zeros(xx.shape)

        for i in ndindex(xx.shape):
            obj[i] = objFunc((xx[i], yy[i]), model, tvec, data_mn, mdxind, pnms, fixpars)

        figttl = '{}_vs_{}_{}_{}'.format(model.name, data.name, pnms[0], pnms[1])
        title(figttl, fontsize = 18)
        xlabel(pnms[0])
        ylabel(pnms[1])

        imshow(obj, vmin=obj.min(), vmax=5000., cmap=cm.coolwarm, extent=(xx.min(), xx.max(), yy.min(), yy.max()), aspect = 'auto', origin = 'lower')
        colorbar()

    else:
        print 'landscape can accept either one or two varied parameters!'
        return
    

    if svfig:
        timestmp = datetime.datetime.now().strftime("-%Y-%m-%d-%H-%M-%s")
        savefig(svpath + '/{}{}.pdf'.format(figttl, timestmp))
        filename = svpath + '/{}{}.txt'.format(figttl, timestmp) 
        f = open(filename, 'w')
        f.write('Fixed pars:\n')
        f.write(str(fixpars) + '\n\n')
        f.write('Varied pars:\n')
        for pg in par_grid:
            f.write('{0}\t{1:.4f}\t{2:.4f}\n'.format(pg[0], pg[1][0], pg[1][-1]))
        f.close()

    return xx, yy, obj

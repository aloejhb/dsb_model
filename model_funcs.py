from pfuncs_ssa import *
from scipy.integrate import odeint
from pylab import *
from matplotlib import gridspec
from numpy import random
from scipy.stats import scoreatpercentile as perc
import sys

def mk_hfunc(inds):
    def hfunc(x):
        h = 1
        for i in inds:
            h *= x[i]
        return h
    return hfunc

def convertJumps(sol, t, tvec):
    j = 0
    xvec = []
    for i,tv in enumerate(tvec):
        while tv >= t[j]:
            j = j + 1
        x_now = sol[j - 1,:] 
        xvec.append(x_now)
    return array(xvec)

def setNewfig(title):
    fig = figure()
    fig.suptitle(title, fontsize=18)
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 

    ax1 = fig.add_subplot(gs[0])
    setp(ax1.get_xticklabels(), visible=False)
    setp(ax1.get_yticklabels(), fontsize = 15)
    ax1.set_ylabel( '# Foci' , fontsize = 16)

    ax2 = fig.add_subplot(gs[1], sharex = ax1)
    setp(ax2.get_xticklabels(), fontsize = 15)
    setp(ax2.get_yticklabels(), fontsize = 10)
    ax2.set_xlabel( 'Time (h)' , fontsize = 16)
    ax2.set_ylabel('cfunc(t)', fontsize = 16)
    return fig, ax1, ax2

def plotMnStd(tvec, mn, st, label, color = None):
    if len(mn) != len(tvec):
        print "The length of mn and tvec should be the same"
        return

    plot(tvec, mn, label = label, color = color, alpha = 0.7, lw = 2.5)
    if st is not None:
        if len(st) != len(mn):
            print "The length of st and mn should be the same"
            return
        fill_between(tvec, mn - st, mn + st, color = color, alpha = 0.3)

    # legend(frameon = False, prop={'size':12})

def plotPerc(tvec, mn, q1, q3, label, color = None):
    if len(mn) != len(tvec):
        print "The length of mn and tvec should be the same"
        return

    plot(tvec, mn, label = label, color = color, alpha = 0.7, lw = 2.5)
    if q1 is not None:
        if len(q1) != len(mn):
            print "The length of st and mn should be the same"
            return
        fill_between(tvec, q1, q3, color = color, alpha = 0.3)

    # legend(frameon = False, prop={'size':12})

class Model:
    # I can:
    #   accept parameters
    #   simulate
    #   plot
    def __init__(self, name, x0, stoich, ttime, mk_pfuncs, pars, xname, xcolor, detxcolor, pname, pcolor, obsinds = [], obsname = [], obscolor = [], detobscolor = [], fuzz = None, pfunc_typ = 'spline'):
        self.name = name

        self.x0 = x0
        self.stoich = stoich
        self.ttime = ttime
        self.pars = pars
        self.mk_pfuncs = mk_pfuncs
        self.pfuncs = mk_pfuncs(pars)

        self.xname = xname
        self.xcolor = xcolor
        self.detxcolor = detxcolor

        self.pname = pname
        self.pcolor = pcolor

        self.obsinds = obsinds
        self.xname.extend(obsname)
        self.xcolor.extend(obscolor)
        self.detxcolor.extend(detobscolor)

        self.fuzz = fuzz
        
        self.pfunc_typ = pfunc_typ


    def sim(self, ttime, pfuncs = None):
        if pfuncs is None:
            pfuncs = self.pfuncs

        if self.pfunc_typ == 'spline':
            sol, t = spline_switch_ssa(self.x0, ttime, pfuncs, self.stoich)
        else:
            sol, t = multi_switch_ssa(self.x0, ttime, pfuncs, self.stoich)

        for inds in self.obsinds:
            obs = array([sol[:,ind] for ind in inds]).sum(axis = 0)
            sol = c_[sol, obs]
        return sol, t

    def detSim(self, tvec):
        ode = self.mk_ode()
        detsol = odeint(ode, self.x0, tvec)

        for inds in self.obsinds:
            obs = array([detsol[:,ind] for ind in inds]).sum(axis = 0)
            detsol = c_[detsol, obs]

        return detsol

    def multiSim(self, tvec, numSim = 200):
        xall = []
        ## process bar
        pbstr = ''
        pbstep = floor(float(numSim)/float(50))
        if pbstep == 0:
            pbstep = 1
        print "Start the simulations ..."

        for i in range(numSim):
            ## process bar
            if (i+1) % pbstep == 0:
                pbstr += '='
                sys.stdout.write('\r[' + str(int(i*100/numSim)) + ' %]' + pbstr)
                sys.stdout.flush() 

            if self.fuzz:
                pfuncs = self.mk_rnd_pfuncs()
            else:
                pfuncs = self.pfuncs

            sol, t = self.sim(tvec[-1], pfuncs) 
            xvec = convertJumps(sol, t, tvec)
            xall.append(xvec)
        xall = array(xall)
        print ""

        return xall


    def setPars(self, newpars):
        for k,v in newpars.items():
            if v < 0:
                raise Exception('Value of parameters should not be negative', k)
        self.pars.update(newpars)
        self.pfuncs = self.mk_pfuncs(self.pars)

    def setAllPars(self, newpars):
        for k,v in newpars.items():
            if v < 0:
                raise Exception('Value of parameters should not be negative', k)

        self.pars = newpars
        self.pfuncs = self.mk_pfuncs(self.pars)

    def printPars(self, tofile = False):
        pstr = ''
        for k in sorted(self.pars.keys()):
            pstr += '{0}\t{1:.4f}\n'.format(k,self.pars[k])
        if tofile:
            return pstr
        print pstr
        


    def pltSim(self, ttime = None):
        if ttime is None:
            ttime = self.ttime

        if self.fuzz:
            pfuncs = self.mk_rnd_pfuncs()
        else:
            pfuncs = self.pfuncs

        sol, t = self.sim(ttime, pfuncs)

        fig, ax1, ax2 = setNewfig(self.name + ' single sim')
        sca(ax1)
        for i in range(len(sol[0,:])):
            step(t,sol[:,i],label = self.xname[i], color = self.xcolor[i], where = 'post', alpha = 0.8, lw = 2.5)
        legend(frameon = False, prop={'size':12})
        sca(ax2)
        self.pltCfuncs(t, pfuncs)

    def pltDetSim(self, tvec = None):
        if tvec is None:
            tvec = linspace(0, self.ttime, 4*self.ttime + 1)
        detsol = self.detSim(tvec)
        fig, ax1, ax2 = setNewfig(self.name + ' det sim')
        sca(ax1)
        for i in range(len(detsol[0,:])):
            plot(tvec, detsol[:,i],'--', color = self.detxcolor[i],lw = 1.5, label = 'det ' + self.xname[i])
        legend(frameon = False, prop={'size':12})
        sca(ax2)
        self.pltCfuncs(tvec)

    def pltMultiSim(self, tvec = None, numSim = 200):
        if tvec is None:
            tvec = linspace(0, self.ttime, 4*self.ttime + 1)

        xall = self.multiSim(tvec, numSim = numSim)

        detsol = self.detSim(tvec)

        fig, ax1, ax2 = setNewfig(self.name + ' multi sim')
        sca(ax1)
        for i in range(len(detsol[0,:])):
            xall2 = xall[:,:,i]
            mn = mean(xall2, axis = 0) 

            q1 = []
            q3 = []
            for cl in xall2.T:
                cl2 = [c for c in cl if not isnan(c)]
                q1.append(perc(cl2, 25))
                q3.append(perc(cl2, 75))
            q1 = array(q1)
            q3 = array(q3)

            plotPerc(tvec, mn, q1, q3, label = self.xname[i], color = self.xcolor[i])
            plot(tvec, detsol[:,i],'--', color = self.detxcolor[i],lw = 1.5, label = 'det ' + self.xname[i])
        legend(frameon = False, ncol=2,  prop={'size':10})
        sca(ax2)
        self.pltCfuncs(tvec)
        
        return fig, ax1, ax2



    def mk_rnd_pfuncs(self):
        rndpars = self.pars.copy()
        for k,v in self.fuzz.items():
            if v[0] == 'exponential':
                rndpars[k] = rndpars[k] + random.exponential(v[1])
            elif v[0] == 'normal':
                rndpars[k] = random.normal(rndpars[k], v[1])
            elif v[0] == 'lognormal':
                rndpars[k] = rndpars[k] + random.lognormal(0, v[1])
        return self.mk_pfuncs(rndpars)

    def mk_ode(self):
        def ode(x,t):
            rates = self.stoich.transpose() * array([pf(t, x) for pf in self.pfuncs])
            rates = sum(rates, axis = 1)
            return rates
        return ode

    def pltCfuncs(self, tvec, pfuncs = None):
        if pfuncs is None:
            pfuncs = self.pfuncs

        etvec = linspace(tvec[0], tvec[-1], 1000)
        for i,pf in enumerate(pfuncs):
            cf = vectorize(pf.cfunc)
            plot(etvec, cf(etvec), color = self.pcolor[i], label = self.pname[i], alpha = 0.7)
        legend(mode = 'expand', ncol = 4, bbox_to_anchor = (0,1,1,0.35),borderaxespad=0., prop={'size':8})
        ylim(0, 1.)


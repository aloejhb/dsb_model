from numpy import *
from numpy.random import random as U
from scipy.interpolate import UnivariateSpline as spline
from switch_splines import mk_SwitchSpline
tauMax = 10.
Ntau = 20.

def _draw(ps):

        s = 0
        Vflat = []
        for p in ps:
            s = s + p
            Vflat.append(s)
        
        r = U()
        for k,Xi in enumerate(Vflat):
            if Xi > r:
                return k

class step_pfunc:
    def __init__(self, level, switch, hfunc):
        self.level = level
        self.switch = switch
        self.hfunc = hfunc
        self.step_now = 0 
    def cfunc(self,t):
        for i, tsw in enumerate(self.switch):
            if t < tsw:
                return self.level[i]
        return self.level[-1]
    def __call__(self,t,X):
        return self.cfunc(t)*self.hfunc(X)
    def p_now(self, X):
        return self.level[self.step_now] * self.hfunc(X)

class spline_pfunc:
    def __init__(self, level, switch, delay, hfunc, ttime, tauMax = tauMax):
        self.level = level
        self.switch = switch
        self.delay = delay
        self.hfunc = hfunc
        if len(level) == 1:
            self.cfunc = lambda x: level[0]
            self.intCfunc = lambda x: level[0] * x
        else:
            if len(level) == 2:
                cfunc = mk_SwitchSpline(level1 = level[0], level2 = level[1], switch1 = switch[0], delay1 = delay[0], maxT = ttime + tauMax)
            else:
                cfunc = mk_SwitchSpline(level1 = level[0], level2 = level[1], level3 = level[2], switch1 = switch[0], switch2 = switch[1], delay1 = delay[0], delay2 = delay[1], maxT = ttime + tauMax)

            self.cfunc = cfunc
            self.intCfunc = cfunc.antiderivative()
    def __call__(self, t, x):
        return self.cfunc(t) * self.hfunc(x)
    def intFunc(self, t, x):
        def ifunc(tau):
           return (self.intCfunc(t+tau) - self.intCfunc(t)) * self.hfunc(x)
        return vectorize(ifunc)

def multi_switch_ssa(x0, ttime, pfuncs, stoich):
    
    # initiation
    t = 0
    x = array(x0)
    time = [0]
    sol = [x0]
    nextSwitch = 0
  
    # take out switches time from each reaction, and record reaction index 
    tswitch = []
    rswitch = []
    for i, pf in enumerate(pfuncs):
        tswitch.extend(pf.switch)
        rswitch.extend([i] * len(pf.switch))

    # sort tswitch and according reaction index 
    enum_tsw_sort = sorted(enumerate(tswitch), key = lambda x: x[1])
    tsw_sort = [et[1] for et in enum_tsw_sort]
    rsw_sort = [rswitch[et[0]] for et in enum_tsw_sort]

    # Append total time to tswitch to include the situation when tau > t + ttime
    tsw_sort.append(ttime+1)

    ak = [pf.p_now(x) for pf in pfuncs]
    # Start simulation
    while t < ttime:
        # Draw random number for tau calculation
        r1 = U()

        # Loop over sorted switches to find which switches can be jumped over
        # within this step of simulation
        Asw = 0
        S = 0
        dS = 0
        # for i in range(nextSwitch, len(tsw_sort)):
        i = nextSwitch
        while i < len(tsw_sort):

            delta = tsw_sort[i] - t
            ak = [pf.p_now(x) for pf in pfuncs]
            a0 = sum(ak) 
            if i == nextSwitch:
                Asw = Asw + a0 * delta
            else:
                Asw = Asw + a0 * (tsw_sort[i] - tsw_sort[i-1])
            S = S + dS
            
            if log(1/r1) < Asw:
                break

            if i < len(tsw_sort) - 1:
                # Jump over this switch, adding a shift term to S and update the propensity
                pfjump = pfuncs[rsw_sort[i]]
                dS = (pfjump.level[pfjump.step_now] - pfjump.level[pfjump.step_now + 1]) * pfjump.hfunc(x) * delta
                pfjump.step_now += 1

                # Jump over all switches if they happen at the same time
                while i < len(tsw_sort) - 1 and tsw_sort[i + 1] == tsw_sort[i]:
                    i += 1
                    if i < len(rsw_sort):
                        pfjump = pfuncs[rsw_sort[i]]
                        dS = (pfjump.level[pfjump.step_now] - pfjump.level[pfjump.step_now + 1]) * pfjump.hfunc(x) * delta
                        pfjump.step_now += 1
            i += 1
       
        # Record the switch index that is jumped.
        # nextSwitch is the index of the first switch for the next step of simulation
        nextSwitch = i

        # If total propensity is not 0, 
        # calculate tau and select which reation happens and update x
        if a0 != 0:
            tau = 1 / a0 * (log(1/r1) - S)
            t = t + tau
            mu = _draw(ak/a0)
            x = x + stoich[mu]
        # If total propensity is 0,
        # nothing happens until the next switch or end of simulation
        else:
            t = tsw_sort[nextSwitch]
            if nextSwitch < len(tsw_sort) - 1:
                pfjump = pfuncs[rsw_sort[i]]
                pfjump.step_now += 1
            x = x
        time.append(t)
        sol.append(x)
            
    for pf in pfuncs:    
        pf.step_now = 0 

    return array(sol), array(time)

def spline_switch_ssa(x0, ttime, pfuncs, stoich, tauMax = tauMax, Ntau = Ntau):

    # Initiation
    t = 0
    x = array(x0)
    time = [0]
    sol = [x0]
    tauMax = float(tauMax)

    while t < ttime:

        r1 = U()
        R = log(1/r1)

        Amax = A0(tauMax)
        if R > Amax:
        tauvec = linspace(0, tauMax, Ntau)

        Ak = array([pf.intFunc(t, x)(tauvec) for pf in pfuncs])
        A0vec = sum(Ak, axis = 0)
        revA0 = spline(A0vec, tauvec, s = 0)
        tau = revA0(-log(r1))
        if tau > tauMax2:
            print 'Warning, tauMax too small!!'
    
        
        ak = [pf(t+tau, x) for pf in pfuncs]
        a0 = sum(ak)
        mu = _draw(ak/a0)

        t = t + tau
        x = x + stoich[mu]

        time.append(t)
        sol.append(x)

    return array(sol), array(time)





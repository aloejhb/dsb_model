from numpy import *
from numpy.random import random as U
from scipy.interpolate import UnivariateSpline as spline
from switch_splines import mk_SwitchSpline
Ntau = 50.
tauMax = 10.

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
    def intFunc_bak(self, t, x):
        def ifunc(tau):
           return (self.intCfunc(t+tau) - self.intCfunc(t)) * self.hfunc(x)
        return vectorize(ifunc)
    def intFunc(self, t, x, tau):
       return (self.intCfunc(t+tau) - self.intCfunc(t)) * self.hfunc(x)

def multi_switch_ssa(x0, ttime, pfuncs, stoich, react_count=False):
    
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

def a0_func(pfuncs, t, x):
    return sum([pf(t, x) for pf in pfuncs])

def A0_func(pfuncs, t, x, tau):
    return sum([pf.intFunc(t, x, tau) for pf in pfuncs])

# Function for simulating reactions by modified Gillespie algorithm using spline
# Input:
#   x0: initial particle num of species
#   ttime: total time for simulation
#   pfuncs: array of spline_pfunc objects, containing propensity functions for each reaction
#   stoich: stoichiometric matrix of the system of reactions
#   tauMax: ??
#   Ntau: ??
# Output:
#   sol: particle num of species at each time point of reaction
#   time: array of time points for reactions to take place
def spline_switch_ssa(x0, ttime, pfuncs, stoich, tauMax=tauMax, Ntau=Ntau, react_count=False):

    # Initiation
    t = 0 # curren time
    x = array(x0) # current particle num of species
    time = [0] # time: array of time points for reactions to take place
    sol = [x0] # particle num of species at each time point of reaction
    rinds = [None]
    tauMax = float(tauMax)

    tsw = []
    for pf in pfuncs:
        tsw.extend(pf.switch)
    tsw = sorted(list(set(tsw)))
    tsw = [sw for sw in tsw if sw < ttime]
    tsw.append(ttime)

    while t < ttime:

        r1 = U()
        R = -log(r1)

        tausw = [sw - t for sw in tsw if sw - t > 0]
        A0sw = [A0_func(pfuncs, t, x, tas) for tas in tausw]

        for i,A0s in enumerate(A0sw):
            if R < A0s:
                swi = i
                break

        if swi == len(tausw):
            time.append(ttime)
            sol.append(x)
            break

        if swi == 0:
            x1, y1 = 0, 0
        else:
            x1, y1 = tausw[swi-1], A0sw[swi-1]

        x2, y2 = tausw[swi], A0sw[swi]

        ### still have problem of determining taumax when propensity is really large
        ymax = 10.
        ntau = Ntau 
        taumax = ymax / (y2 - y1) * (x2 - x1)

        ext = 1.
        tau2 = min(x1 + taumax, x2 + ext)
        tauvec = linspace(x1, tau2, ntau)
        A0vec = [A0_func(pfuncs, t, x, tau) for tau in tauvec]
        spl = spline(A0vec, tauvec, s=0)
        tau = spl(R)
        
        if x1 + taumax == tau2 and tau > tau2 :
            print "Warning: taumax too small!"
            
     

        ak = [pf(t+tau, x) for pf in pfuncs]
        a0 = sum(ak)
        mu = _draw(ak/a0)

        t = t + tau
        x = x + stoich[mu]

        time.append(t)
        sol.append(x)

        if react_count:
                rinds.append(mu)


    return array(sol), array(time), array(rinds)





def spline_switch_ssa_bak(x0, ttime, pfuncs, stoich, tauMax = tauMax, Ntau = Ntau):

    # Initiation
    t = 0
    x = array(x0)
    time = [0]
    sol = [x0]
    tauMax = float(tauMax)

    while t < ttime:

        r1 = U()

        tauvec = linspace(0, tauMax, Ntau)
        Ak = array([pf.intFunc(t, x)(tauvec) for pf in pfuncs])
        A0vec = sum(Ak, axis = 0)
        revA0 = spline(A0vec, tauvec, s = 0)
        tau = revA0(-log(r1))
        if tau > tauMax:
            print 'Warning, tauMax too small!!'
    
        
        ak = [pf(t+tau, x) for pf in pfuncs]
        a0 = sum(ak)
        mu = _draw(ak/a0)

        t = t + tau
        x = x + stoich[mu]

        time.append(t)
        sol.append(x)

    return array(sol), array(time)





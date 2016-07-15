from model_funcs import *
from pfuncs_ssa import *

mdname = 'toy'
ttime = 100.
def mk_pfuncs(pars):
    p1 = spline_pfunc((0.5, 150., 1.), (5., 5.3), (0.1, 0.1), mk_hfunc(()), ttime)
    p2 = spline_pfunc((0.4, 1.), (7.,), (10.,), mk_hfunc((0,)), ttime)
    return [p1, p2]

pfuncs = mk_pfuncs({})

x0 = [10]
stoich = array(([1], [-1]))
pars = {}
# spline_switch_ssa(x0, ttime, pfuncs, stoich)

xname = ['foci']
xcolor = ['dodgerblue']
detxcolor = ['cornflowerblue']

pname = ['dam', 'repair']
pcolor = ['blue', 'red']

model = Model(mdname, x0, stoich, ttime, mk_pfuncs, pars, xname, xcolor, detxcolor, pname, pcolor)

model.pltSim()
model.pltMultiSim(numSim = 50)
show()


tsw = []
for pf in pfuncs:
    tsw.extend(pf.switch)
    nd = [pf.delay[i] + sw for i,sw in enumerate(pf.switch)]
    tsw.extend(nd)
tsw = sorted(list(set(tsw)))
tsw = [sw for sw in tsw if sw < ttime]
tsw.append(ttime)

t = 0
x = array(x0)
tausw = [sw - t for sw in tsw if sw - t > 0]
A0sw = [A0_func(pfuncs, t, x, tas) for tas in tausw]

etau = linspace(0, ttime-t, 10000)
eA0 = [A0_func(pfuncs, t, x, ta) for ta in etau]

ion()
figure()
plot(etau, eA0, label = 'A0_eval')
plot(tausw, A0sw, '*', label = 'A0sw')

taumax = 10.

# r1 = U()
# R = -log(r1)
R = 950.

swi = 0
for i,A0s in enumerate(A0sw):
    if R < A0s:
        swi = i
        break

if swi == 0:
    x1, y1 = 0, 0
else:
    x1, y1 = tausw[swi-1], A0sw[swi-1]
x2, y2 = tausw[swi], A0sw[swi]

plot((x1,x2), (y1,y2), 'h')
axhline(R)

def met1():
    yintvl = 5.
    npoint = int((y2 - y1)/yintvl)

    # add to switch end
    tau2 = min(x1 + taumax, x2)
    tauvec = linspace(x1, tau2, npoint)
    A0vec = [A0_func(pfuncs, t, x, tau) for tau in tauvec]
    spl = spline(A0vec, tauvec, s=0)
    eA0vec = linspace(A0vec[0], A0vec[-1], npoint*100)
    tau = spl(R)
    plot(tau, R, 'H')
    # plot(tauvec, A0vec, 'd')
    plot(spl(eA0vec), eA0vec)
    # plot(spl(A0vec), A0vec, 'D') 

def met2():
    taut = (R - y1) * (x1 - x2)/(y1 - y2) + x1
    A0t = A0_func(pfuncs, t, x, taut)

    plot(taut, A0t, 'o', label = 'A0t')


    ds = 0.01
    eps = 0.5
    stn = 0
    while abs(R - A0t) > eps:
        a0t = a0_func(pfuncs, t+taut, x)
        dtau = ds / a0t
        taut += - sign(A0t - R)*dtau
        A0t = A0_func(pfuncs, t, x, taut)
        stn += 1

    print stn
    tau = taut
    plot(tau, R, 'H')



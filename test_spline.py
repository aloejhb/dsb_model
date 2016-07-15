from switch_splines import mk_SwitchSpline
from scipy.interpolate import UnivariateSpline as spline
from numpy import *
from pylab import *
s1 = mk_SwitchSpline(0.45, 5., 0.01, 100., 5.5, 0.01, 0.45, smord=1)
s2 = mk_SwitchSpline(0.2, 5., 2., 0.4)
def s3(x):
    return s1(x) + s2(x)

etvec = linspace(0, 100, 1001)
is1 = s1.antiderivative()
is2 = s2.antiderivative()
def is3(x):
    return is1(x) + is2(x)

tauMax = 10.
Ntau = 20.
tauvec = linspace(0, tauMax, Ntau)
etauvec = linspace(0, tauMax, Ntau*10)

A0vec = is3(tauvec)
eA0vec = linspace(A0vec.min(), A0vec.max(), 100000)
revA0 = spline(A0vec, tauvec, s=0)
revA02 = spline(A0vec[:9], tauvec[:9], s=0)
zA0 = spline(tauvec, A0vec, s=0)

ion()

figure()
plot(tauvec, A0vec, 'o', label = 'A0vec')
plot(etauvec, is3(etauvec), label = 'is3')

plot(revA0(A0vec), A0vec, '*', label = 'revA0')
plot(revA0(eA0vec), eA0vec, label = 'revA0')

# plot(revA0(A0vec[:9]), A0vec[:9], '*', label = 'revA01')
# plot(revA0(eA0vec[:100000]), eA0vec[:100000], label = 'revA0')

# plot(etauvec, zA0(etauvec), label = 'zA0')

legend(frameon = False)

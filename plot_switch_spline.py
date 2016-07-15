from pylab import *
from switch_splines import mk_SwitchSpline

maxT = 100.

level1 = 0.5
level2 = 2.
level3 = 1.

switch1 = 30.
switch2 = 60.

delay1 = 10.
delay2 = 15.

s1,x1,y1 = mk_SwitchSpline(level1,switch1,delay1,level2,maxT=maxT, debug = True)
s2,x2,y2 = mk_SwitchSpline(level1,switch1,delay1,level2,switch2,delay2,level3,maxT, debug = True)

xgrid = linspace(0, maxT, 1000)

close('all')
plot(xgrid, s2(xgrid), linewidth=2.0, color='lightgrey')

plot(x2[::5], y2[::5], 'o', ms=4, color='k')


# vlclr = 'lightgreen'
# plot([switch1,switch1], [0,level1], color=vlclr)
# plot([switch1+delay1,switch1+delay1], [0,level2], color=vlclr)
# plot([switch2,switch2], [0,level2], color=vlclr)
# plot([switch2+delay2,switch2+delay2], [0,level3], color=vlclr)

text(switch1, level1, r'$(\tau_k^1,v_k^0)$', ha='right', va = 'bottom')
text(switch1+delay1, level2, r'$(\tau_k^1+\delta_k^1,v_k^1)$', ha='right')
text(switch2, level2, r'$(\tau_k^2,v_k^1)$')
text(switch2+delay2, level3, r'$(\tau_k^2+\delta_k^2,v_k^2)$', va='bottom')

tick_params(axis='both', which='both', top='off', bottom='off', right='off', left='off',  labelleft='off', labelbottom='off') 

ax = gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)

ax.set_xlabel(r'$t$', fontsize=15)
ax.set_ylabel(r'$c_k(t)$', fontsize=15)

ax.set_ylim(0,2.5)
savefig('../protocol/figures/spline.pdf')

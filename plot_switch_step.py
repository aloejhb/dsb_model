from pylab import *
from pfuncs_ssa import step_pfunc
maxT = 100.

level = [0.5, 2., 1.]

switch = [30., 60.]

sp = step_pfunc(level, switch, lambda x: 1)
sf = vectorize(sp.cfunc)

xgrid = linspace(0, maxT, 1000)

close('all')
plot(xgrid, sf(xgrid), linewidth=2.0, color='lightgrey')


text(switch[0], level[0], r'$(\tau_k^1,v_k^0)$', ha='right', va = 'bottom')
text(switch[0], level[1], r'$(\tau_k^1,v_k^1)$', ha='right')
text(switch[1], level[1], r'$(\tau_k^2,v_k^1)$')
text(switch[1], level[2], r'$(\tau_k^2,v_k^2)$', va='bottom')

tick_params(axis='both', which='both', top='off', bottom='off', right='off', left='off',  labelleft='off', labelbottom='off') 

ax = gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)

ax.set_xlabel(r'$t$', fontsize=15)
ax.set_ylabel(r'$c_k(t)$', fontsize=15)

ax.set_ylim(0,2.5)
savefig('../figures/step.pdf')

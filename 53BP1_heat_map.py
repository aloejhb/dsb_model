
from data_funcs import *
from pylab import *
from matplotlib.pyplot import savefig
from cc_strat import tag_tr
from collections import Counter

# ipath = '../results/replot_data/'
# exp = 'NCS_G1G2'
# exp_path = exp.replace('_', '/')
# if exp_path[-1] != '/':
#    exp_path +='/'
ipath = '../data/'
exp_path = 'NCSoffpool/'

opath = '../results/heat_map/'
figname = 'NCSoff_heat_map'
figtitle = '(MCYN-off 100nM NCS at 5h)'
sfreq = 6 
startFrm = 200
treatFrm = 230
trl = Load_Tracks(ipath + exp_path)

ttreat = (treatFrm - startFrm) / sfreq
# sort by tag, gm_treat, tup, tdiv
trinfo = [list(tag_tr(tr, ttreat)) for tr in trl]
ztr = zip(trinfo, trl)
tag_order = ['G1', 'G1S', 'SG2', 'G2', 'G2M']
def sort_info(z):
    zi = z[0]
    to = tag_order.index(zi[0])
    return (to, zi[1:])
def sort_tag(z):
    to = tag_order.index(z[0])
    return to
trinfo_srt = [x for (x,y) in sorted(ztr, key=sort_info)]

trl_srt = [y for (x,y) in sorted(ztr, key=sort_info)]
fc = getDataMatrix(exp, trl_srt, 3, sfreq, startFrm)
gm = getDataMatrix(exp, trl_srt, 5, sfreq, startFrm)
tvec = fc.tvec
fca = fc.all
gma = gm.all
by = 1
bx = 3
ny = len(fca)
nx = len(fca[0,:])
my = ny % by
mx = nx % bx

if mx:
    fca2 = fca[:, :-mx]
    gma2 = gma[:, :-mx]
if my > 1:
    trinfo_srt2 = trinfo_srt[:-my]
    fca2 = fca2[:-my, :]
    gma2 = gma2[:-my, :]
else:
    trinfo_srt2 = trinfo_srt
    
ny2 = len(fca2)/by
nx2 = len(fca2[0,:])/bx

tag_ls = [tri[0] for tri in trinfo_srt2]
cttag = Counter(tag_ls)
cttag_srt = sorted(list(cttag.items()), key=sort_tag)

fc2 = nanmean(nanmean(fca2.reshape(ny2,by,nx2,bx), axis = 1), axis = 2)
gm2 = nanmean(nanmean(gma2.reshape(ny2,by,nx2,bx), axis = 1), axis = 2)

# ion()
fig1 = figure(1)
figure(1)
imshow(fc2, vmin=0, vmax=40, extent=(tvec.min(), tvec.max(), 1, len(fca2)), aspect = 'auto', origin = 'lower', interpolation = 'none', cmap=cm.coolwarm)
cb = colorbar()
cb.ax.tick_params(labelsize=14)
tick_params(labelsize=14)
xlabel('Time (h)', size=14)
title('# 53BP1 Foci ' + figtitle, size=14)

fig2 = figure(2)
figure(2)
imshow(gm2, vmin=nanmin(gm2), vmax=nanmax(gm2), extent=(tvec.min(), tvec.max(), 1, len(gma2)), aspect = 'auto', origin = 'lower', interpolation = 'none', cmap=cm.YlGn)
cb = colorbar()
cb.ax.tick_params(labelsize=14)
tick_params(labelsize=14)
xlabel('Time (h)', size=14)
title('Geminin Intensity [AU] ' + figtitle, size=14)

cpos = 0
for i,tg in enumerate(cttag_srt):
    ytxt = 0.5 * tg[1] + cpos
    yl = tg[1] + cpos
    figure(1)
    # axhline(yl, ls = '--', linewidth = 2., color = 'darkgreen')
    axhline(yl, ls = '--', linewidth = 3., color = 'firebrick')
    text(-10, ytxt, tg[0], color = 'firebrick', size = 14)

    figure(2)
    #mediumspringgreen
    axhline(yl, ls = '--', linewidth = 3., color = 'firebrick')
    text(-10, ytxt, tg[0], color = 'firebrick', size = 14)

    cpos = yl
fig1.savefig(opath + figname + '_Foci' + '.pdf')
fig2.savefig(opath + figname + '_GMNN' + '.pdf')
#show()
close('all')

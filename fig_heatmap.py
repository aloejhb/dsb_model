from data_funcs import *
from pylab import *
from matplotlib.pyplot import savefig
from cc_strat import tag_tr
from collections import Counter

ipath = '../data/'
opath = '../figures/'
figname = 'heatmap_zusammen'

sfreq = 6
startFrm = 200
treatFrm = 230
ttreat = (treatFrm - startFrm) / sfreq

exp_paths = ['NCSpool/', 'NCSoffpool/']
tag_order = ['G1', 'G1S', 'SG2', 'G2', 'G2M']

def sort_info(z):
    zi = z[0]
    to = tag_order.index(zi[0])
    return (to, zi[1:])

def sort_by_tag(trl):
    # sort by tag, gm_treat, tup, tdiv
    trinfo = [list(tag_tr(tr, ttreat)) for tr in trl]
    ztr = zip(trinfo, trl)
    ztr_srt = sorted(ztr, key=sort_info)

    trinfo_srt = [x for (x,y) in ztr_srt]
    trl_srt = [y for (x,y) in ztr_srt]
    return trl_srt, trinfo_srt

def bin_dtmat(mat, bx=3, by=1):
    bx = int(bx)
    by = int(by)
    nx = len(mat[0,:])
    ny = len(mat)
    mx = nx % bx
    my = ny % by

    if mx:
        binmat = mat[:, :-mx]
    if my:
        binmat = mat[:-my, :]

    nx2 = nx/bx
    ny2 = ny/by

    binmat = nanmean(nanmean(binmat.reshape(ny2,by,nx2,bx), axis = 1), axis = 2)

    return binmat

def tag_count(tags):
    tag_dict = {}
    for tag in tags:
        if tag in tag_dict:
            tag_dict[tag] += 1
        else:
            tag_dict[tag] = 1
            
    tag_cnt = []
    for tag in tag_order:
        tag_cnt.append(tag_dict[tag])

    return tag_cnt
    

def get_fc_gm(exp_path):
    trl = Load_Tracks(ipath + exp_path)
    # sort by tag
    trl_srt, trinfo_srt = sort_by_tag(trl)
    tag_cnt = tag_count([tri[0] for tri in trinfo_srt])
    # get data matrix
    fc = getDataMatrix(exp, trl_srt, 3, sfreq, startFrm)
    gm = getDataMatrix(exp, trl_srt, 5, sfreq, startFrm)
    tvec = fc.tvec
    fca = fc.all
    gma = gm.all
    # bin
    fca = bin_dtmat(fca, bx=3, by=1)
    gma = bin_dtmat(gma, bx=3, by=1)
    
    return tvec, fca, gma, tag_cnt


def plot_tag(tag_cnt):
    cpos = 0
    for i,tag in enumerate(tag_order):
        ytxt = 0.5 * tag_cnt[i] + cpos
        yl = tag_cnt[i] + cpos
        if i < len(tag_order)-1:
            axhline(yl, ls = '--', inewidth = 2., color = 'firebrick')
        text(-8, ytxt, tag, color = 'firebrick', size = 9.)
        cpos = yl

rcParams['xtick.labelsize'] = 9.
rcParams['ytick.labelsize'] = 9.

fig, axs = subplots(2,2)

for i,exp_path in enumerate(exp_paths):
    tvec, fca, gma, tag_cnt = get_fc_gm(exp_path)
    #tvec = arange(1,60)
    #fca = zeros((10,10))
    #gma = zeros((10,10))
    sca(axs[0,i])
    imshow(fca, vmin=0, vmax=40, extent=(tvec.min(), tvec.max(), 1, len(fca)), aspect = 'auto', origin = 'lower', interpolation = 'none', cmap=cm.coolwarm)
    colorbar()
    plot_tag(tag_cnt)
    
    sca(axs[1,i])
    imshow(gma, vmin=0, vmax=14000, extent=(tvec.min(), tvec.max(), 1, len(gma)), aspect = 'auto', origin = 'lower', interpolation = 'none', cmap=cm.YlGn)
    colorbar()
    plot_tag(tag_cnt)

for ax in axs.flat:
    ax.get_yaxis().set_ticks([])

axs[0,0].text(0, 103, 'MYCN-on', fontsize=14)
axs[0,1].text(0, 160, 'MYCN-off', fontsize=14)

    
fig.savefig(opath + figname + '.pdf')
close()



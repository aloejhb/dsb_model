from data_funcs import *
from pylab import *
from matplotlib.pyplot import savefig

figname = 'foci_dynm_nuc_rad'
opath = '../figures/'

trackPath = '../results/replot_data/'
exps = ('NCS', 'NCSoff')

labs = ('G1', 'G1S', 'SG2', 'G2', 'G2M')
labs_show = ('g1', 'g1s', 'sg2', 'g2', 'g2m')
colors = ('m', 'royalblue', 'gold', 'red', 'mediumseagreen')

sfreq = 6 
startFrm = 200
treatFrm = 230
ttreat = (treatFrm-startFrm)/sfreq

rc('axes', labelsize = 15)
rc('xtick', labelsize = 15)
rc('ytick', labelsize = 15)

fig,axs = subplots(2,2, sharex='col', figsize = (10,7))


for i,exp in enumerate(exps):
    track_lists = []
    foci_list = []
    rad_list = []
    area_list = []

    for lab in labs:
#    for lab in [labs[0]]:
        exp_path = exp + '/' + lab + '/'

        trl = Load_Tracks(trackPath + exp_path)
        track_lists.append(trl)

        fc = getDataMatrix(exp, trl, 3, sfreq, startFrm)
        foci_list.append(fc)

        # ar = getDataMatrix(exp, trl, 7, sfreq, startFrm)
        # area_list.append(ar)

        rd = getDataMatrix(exp, trl, 4, sfreq, startFrm)
        rad_list.append(rd)


    sca(axs[0, i])
    pltData(foci_list, labs_show, colors, newfig = False)
    axvline(ttreat, linestyle='--', color='grey')
    axs[0, i].legend_.remove()
    ylim(0,30)
    
    sca(axs[1, i])
    pltData(rad_list, labs_show, colors, newfig = False)
    axvline(ttreat, linestyle='--', color='grey')
    xlabel('Time (h)')
    axs[1, i].legend_.remove()
    ylim(15,30)
    xlim(0, 55)

    axs[0, i].legend(mode = 'expand', ncol = 5, loc = 3, bbox_to_anchor = (0., 1.02, 1., 1.02), borderaxespad=0., prop={'size':10})
    axs[0, i].set_ylabel('# 53BP1 foci')
    axs[1, i].set_ylabel('Mean nuclear radius [AU]')


axs[0,0].text(0, 35, 'MYCN-on', fontsize=15)
axs[0,1].text(0, 35, 'MYCN-off', fontsize=15)
#tight_layout()
savefig(opath + figname + '.pdf')
show()


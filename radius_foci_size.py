from data_funcs import *
from pylab import *
from matplotlib.pyplot import savefig

# trackPath = '../data/'
# expname = 'NCSpool'
trackPath = '../results/replot_data/'
expname = 'rrNCSoff'
figname = expname + '_cell_rad_foci_size'

# colors = ('b',)
# exps = ('NCSpool',)
colors = ('m', 'royalblue', 'gold', 'red', 'mediumseagreen')

exps = ('NCSoff_G1', 'NCSoff_G1S', 'NCSoff_SG2', 'NCSoff_G2', 'NCSoff_G2M')
# exps = ('NCS_G1', 'NCS_G1S', 'NCS_SG2', 'NCS_G2', 'NCS_G2M')
sfreq = 6 
startFrm = 200
# exps = ('G1G1_Doxo', 'G1G2_Doxo', 'G2G2_Doxo')
# sfreq = 4
# startFrm = 1

track_lists = []
foci_list = []
rad_list = []
area_list = []
for i,exp in enumerate(exps):

    exp_path = exp.replace('_', '/')
    if exp_path[-1] != '/':
        exp_path +='/'

    trl = Load_Tracks(trackPath + exp_path)
    track_lists.append(trl)

    fc = getDataMatrix(exp, trl, 3, sfreq, startFrm)
    foci_list.append(fc)

    ar = getDataMatrix(exp, trl, 7, sfreq, startFrm)
    area_list.append(ar)

    rd = getDataMatrix(exp, trl, 4, sfreq, startFrm)
    rad_list.append(rd)

fig,axs = subplots(3,1,sharex = True, figsize = (8,9))    
for ax in axs:
    setp(ax.get_yticklabels(), fontsize = 15)

sca(axs[0])
pltData(foci_list, colors, newfig = False)
ylabel('# foci', fontsize = 16)
axs[0].legend_.remove()
legend(mode = 'expand', ncol = 3, loc = 3, bbox_to_anchor = (0., 1.02, 1., 1.02), borderaxespad=0., prop={'size':12})

sca(axs[1])
pltData(area_list, colors, newfig = False)
ylabel('mean foci area a.u.', fontsize = 16)
axs[1].legend_.remove()
ylim(0,20)

sca(axs[2])
pltData(rad_list, colors, newfig = False)
ylabel('mean nuclear radius a.u.', fontsize = 16)
xlabel('time a.u.', fontsize = 16)
axs[2].legend_.remove()
setp(axs[2].get_xticklabels(), fontsize = 15)
ylim(15,28)


xlim(0, 60)
savefig('../results/radius_foci_size/' + figname + '.pdf')
show()


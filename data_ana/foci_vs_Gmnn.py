from stratif_by_Gmnn import *
from loadTracks import *
import os

trackPath = '../data/'
# exps = ('G1G1/Doxo/', 'G1G1/Doxo+ATRi/', 'G1G1/Doxo+ATMi/', 'G1G2/Doxo/', 'G1G2/Doxo+ATRi/', 'G1G2/Doxo+ATMi/')
exp = 'G1G1/Doxo/'
exp2 = exp.replace('/', '_')
track_list = Load_Tracks(trackPath + exp)

fnames = list(walk(trackPath + exp + '/'))[0][2]
fnames = sorted([f for f in fnames if '.txt' in f])
tnames = [f.replace('.txt', '') for f in fnames]

peakPath = '../results/peak_detect/'
peaks_list = Load_Tracks(peakPath + exp)

figPath = '../results/foci_vs_Gmnn/'
if not os.path.exists(figPath + exp):
    os.makedirs(figPath + exp)

gmnnTh = 1000.
firUps, tcross = getTcross(track_list, gmnnTh)

# for i in range(5):
for i in range(len(track_list)):
    fig, axs = Plot_FociGemRad(track_list[i])
    peaks = peaks_list[i]
    if peaks.ndim < 2:
        peakX = peaks[1]
        peakY = peaks[2]
    else:
        peakX = peaks[:,1]
        peakY = peaks[:,2]
    axs[0].plot(peakX, peakY, '*', markersize = 12)
    axs[1].axhline(y = gmnnTh)
    tc = tcross[i][0]
    axs[1].plot(tc, [gmnnTh] * len(tc), 'o')
    axs[2].set_xlim(0,60)
    savefig(figPath + exp + tnames[i] + '.pdf')
    close()

firPeaks = []
for peaks in peaks_list:
    if peaks.ndim < 2:
        firPk = peaks[1]
    else:
        firPk = peaks[0][1]
    firPeaks.append(firPk)    

# ion()
figure()
plot(firUps, firPeaks, 'o')
xlabel('time of Gmnn cross ' + str(gmnnTh))
ylabel('1st peak of #foci')
show()

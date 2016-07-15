from data_funcs import *
trackPath = '../data/'
exps = ('G1G1/Doxo/', 'G1G2/Doxo/', 'G2G2/Doxo/')
dtxind = 3

sfreq = 4

track_lists, data = loadData(trackPath, exps, dtxind)

ttime = 60
tvec = linspace(0, ttime, sfreq*ttime + 1)

ion()
pltData(tvec, data, ('G1G1_Doxo', 'G1G2_Doxo', 'G2G2_Doxo'), ('darkorange', 'red', 'lightgreen'))



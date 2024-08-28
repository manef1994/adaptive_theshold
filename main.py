import bisect
import math
import statistics
import neurokit2 as nk
from shapely.geometry import LineString
import numpy as np
from matplotlib import *
from pathlib import Path
import matplotlib.pyplot as plt
# import entropy as ent
import wfdb
import antropy as ant
from shapely import Polygon

import statistics
from mne.io import read_raw_edf

# ===================================
" Reading files functions "
# =====================================================
def read_wfdb(path, start, end):

    record = wfdb.rdsamp(path, sampfrom= start, sampto= end)
    # record = wfdb.rdsamp(path)
    # annotation = wfdb.rdann(path, 'ecg')

    # annotation = wfdb.rdann(path, 'ari')
    annotation = wfdb.rdann(path, 'ari', sampfrom= start, sampto= end)
    print('annotion')
    # print(annotation.sample)
    # print(annotation.symbol)

    # Read an annotation as an Annotation object
    sig = record[0]
    rr = record[1]
    print('rr\t', rr)

    fs = rr['fs']
    print('fs equals to: ', fs)

    i = 0
    # l_n = end - start
    start = 0
    end = rr['sig_len']
    l_n = rr['sig_len']
    signal = []

    while i < l_n:
        an = sig[i]
        signal.append(an)
        i = i + 1

    annotation.sample = annotation.sample - start
    anno = ["|", "N", "L", "R", "B", "A", "a", "J", "S", "V", "r", "F", "e", "j", "n", "E", "/", "f", "Q", "?"]
    peak = []

    for i in range(len(annotation.symbol)):

        if (annotation.symbol[i] in anno):
            peak.append(annotation.sample[i])

    print('peaks length\t', len(peak))
    return signal, peak, fs

def reading_EDF(path, ID):

    file_path = path + ID + ".edf"
    edf = read_raw_edf(file_path, preload=False, stim_channel=None, verbose=False)
    xx = edf.ch_names
    index = xx.index("EEG T5")
    fs = edf.info['sfreq']
    fs = int(fs)
    signal_input = edf[index]
    signal = signal_input[0]

    signal_input = signal[0]

    return signal_input


# =====================================================
" Data processing function including:   "
"     Features extraction               "
"     computing the Standard deviation  "
"     computing the thresholdi value    "
# =====================================================

def features(peaks, fs):

    NNRi, approximate, = ([] for i in range(2))
    start = 0
    end = 120 * fs

    while True:

        go = bisect.bisect_left(peaks, start)
        out = bisect.bisect_left(peaks, end)

        peaks_in = peaks[go:out]
        RRi = np.diff(peaks_in)

        NNRi.append(len(RRi))
        approximate.append(ant.app_entropy(RRi))

        start = start + (10 * fs)
        end = end + (10 * fs)

        if (end > peaks[-1]):
            break

    return approximate, NNRi

def std_compute(xx):
    i = 0
    j = 6
    STD = []

    while (True):
        go_in = []

        go_in = xx[i:j]
        std = np.std(go_in)
        STD.append(std)

        j = j + 6
        i = i + 6

        if j > len(xx):
            break
    return STD

def thresholding(array, thresh):

    i = 0
    j = 5
    Thresh = []

    while (True):

        if i == 0:

            table = [array[0], array[1], array[2], array[3], array[4]]
            mean = np.mean(table)
            percent = (mean / 100) * thresh

            Thresh.append(mean + percent)
            Thresh.append(mean + percent)
            Thresh.append(mean + percent)
            Thresh.append(mean + percent)
            Thresh.append(mean + percent)

        else:
            go_in = array[i:j]
            mean = np.mean(go_in)

            percent = (mean / 100) * thresh
            Thresh.append(mean + percent)

        i = i + 1
        j = j + 1


        if j > len(array):
            break


    return Thresh

# =============================================================================================================
# == reading the Siena scalp EEG Database

# =====================================================
"Reading the data"
# =====================================================

path = "C:\\Users\\Manef\\Desktop\\tests_22_02_2024\\tests\\Peaks_RR\\PN06\\"
ID = "PN06-2"
fs = 512
with open(path + "peaks-" + ID + ".txt", 'r') as file1:
    peaks = [float(i) for line in file1 for i in line.split('\n') if i.strip()]

# =====================================================
"define the percentage to use to compute the theshold values"
# =====================================================

percentage = 60

# =====================================================
"Here ewe define how much ictal period we have in the signal (the period after the epileptic seizure )"
# =====================================================

pre_ict = 10

# =====================================================
" Compute the ApEN and the NRRi Features for the input ECG signal"
# =====================================================

app, NNRi = features(peaks,fs)


# =====================================================
"Computing the standard deviation fro the computed features"
# =====================================================

STD_app = std_compute(app)
STD_NNRi = std_compute(NNRi)

thresh_AP = thresholding(STD_app , percentage)
thresh_nn = thresholding(STD_NNRi, percentage)

arr=[]
arr_n=[]
arr1 = [0] * len(thresh_AP)

# =====================================================
" Extraction of the intersections of the STD curve of the ApEn and the NNRi features computed with the threshold values"
# =====================================================

tt = np.arange(0, len(STD_app))

first_line = LineString(np.column_stack((tt,STD_app)))
second_line = LineString(np.column_stack((tt,thresh_AP)))
intersection = first_line.intersection(second_line)
x, y = LineString(intersection).xy
print(" x \t", x)



first_line = LineString(np.column_stack((tt, STD_NNRi)))
second_line = LineString(np.column_stack((tt, thresh_nn)))
intersection_n = first_line.intersection(second_line)
x_n, y_n = LineString(intersection_n).xy
print((x_n))

# =====================================================
"Plotting the results"
# =====================================================

fig, axs = plt.subplots(2, 1)
axs[0].plot(app, label="approximate entropy of the input signal", marker='o')
axs[1].plot(STD_app, label="STD of approximate entropy", marker='o')

axs[0].axvline(x=len(app) - 60, color='red', linestyle='--')
axs[1].axvline(x=len(STD_app) - 10, color='red', linestyle='--')

axs[0].set_title('ApEn based curves ',fontsize=24, y=1)

axs[0].set_xlabel('sample per min')
axs[0].set_ylabel('entropy value')
axs[1].set_xlabel('sample per min')
axs[1].set_ylabel('entropy value')
axs[0].grid(True)
axs[1].grid(True)
axs[0].legend()
axs[1].legend()

ymin = 0.1
ymax = 1.3
axs[0].set_ylim([ymin,ymax])

ymin_1 = 0
ymax_1 = 0.25
axs[1].set_ylim([ymin_1,ymax_1])

if intersection.geom_type == 'MultiPoint':
    axs[1].plot(*LineString(intersection).xy, 'o')
elif intersection.geom_type == 'Point':
    axs[1].plot(*intersection.xy, 'o')

# =====================================================
"plotting NNRi feature curves"
# =====================================================

fig, axs = plt.subplots(2, 1)
axs[0].plot(NNRi, label=" NNRi curve", marker='o')
axs[1].plot(STD_NNRi, label="NRRi STD curve", marker='o')

axs[0].axvline(x=len(NNRi) - 60, color='red', linestyle='--')
axs[1].axvline(x= (len(STD_NNRi) - pre_ict), color='red', linestyle='--')

axs[0].set_title('NNRi based curves ',fontsize=24, y=1)

axs[0].set_xlabel('sample per min')
axs[0].set_ylabel('entropy value')
axs[1].set_xlabel('sample per min')
axs[1].set_ylabel('entropy value')
axs[0].grid(True)
axs[1].grid(True)
axs[0].legend()
axs[1].legend()

ymin = 70 ; ymax = 275
axs[0].set_ylim([ymin,ymax])

ymin_1 = 0; ymax_1 = 14
axs[1].set_ylim([ymin_1,ymax_1])

if intersection_n.geom_type == 'MultiPoint':
    axs[1].plot(*LineString(intersection_n).xy, 'o')
elif intersection_n.geom_type == 'Point':
    axs[1].plot(*intersection_n.xy, 'o')

# =====================================================
"plotting both features STD curves"
# =====================================================

fig, axs = plt.subplots(2, 1)
axs[0].plot(STD_app, label="STD de ApEn", marker='o')
axs[1].plot(STD_NNRi, label="STD de NNRi", marker='o')

axs[0].set_title('STD curves ',fontsize=24, y=1)

axs[0].axvline(x= (len(STD_app) - pre_ict), color='red', linestyle='--')
axs[1].axvline(x= (len(STD_NNRi) - pre_ict), color='red', linestyle='--')

axs[0].set_xlabel('sample per min')
axs[0].set_ylabel('entropy value')
axs[1].set_xlabel('sample per min')
axs[1].set_ylabel('entropy value')
axs[0].grid(True)
axs[1].grid(True)
axs[0].legend()
axs[1].legend()
axs[0].plot(thresh_AP, color='g', linestyle='--')
axs[1].plot(thresh_nn, color='g', linestyle='--')

if intersection.geom_type == 'MultiPoint':
    axs[0].plot(*LineString(intersection).xy, 'o')

elif intersection.geom_type == 'Point':
    axs[0].plot(*intersection.xy, 'o')


if intersection_n.geom_type == 'MultiPoint':
    axs[1].plot(*LineString(intersection_n).xy, 'o')

elif intersection_n.geom_type == 'Point':
    axs[1].plot(*intersection_n.xy, 'o')

plt.show()
# =====================================================
"saving the results"
# =====================================================

# path = "/home/ftay/Fabrice/Features/PN09/"
# # path = "C:\\Users\\Ftay\\Desktop\\PhD\\tests\\siena_vs_fantasia\\Peaks_RR\\PN013"
# Path(path).mkdir(parents=True, exist_ok=True)
#
# np.savetxt(path + 'app-'+ ID + '.txt', np.array(app))
# np.savetxt(path + 'NRRi-'+ ID + '.txt', np.array(NNRi))


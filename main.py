import bisect
import math
import statistics

import pyentrp as pyentrp
# import neurokit2 as nk
from shapely.geometry import LineString
import numpy as np
from matplotlib import *
from pathlib import Path
import matplotlib.pyplot as plt
import entropy as ent
import wfdb
import antropy as ant
import statistics
from mne.io import read_raw_edf
import EntropyHub
import entropy as ent

import warnings
warnings.filterwarnings("ignore")

def read_wfdb(path):

    # record = wfdb.rdsamp(path, sampfrom= start, sampto= end)
    record = wfdb.rdsamp(path)
    # annotation = wfdb.rdann(path, 'ecg')

    # annotation = wfdb.rdann(path, 'ari', sampfrom= start, sampto= end)
    annotation = wfdb.rdann(path, 'ari')
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
def read_wfdb_interictal(path, start, end):

    record = wfdb.rdsamp(path, sampfrom= start, sampto= end)
    annotation = wfdb.rdann(path, 'ari', sampfrom= start, sampto= end)
    print('annotion')

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
def features(peaks, fs):

    NN50, approximate, = ([] for i in range(2))
    start = 0
    end = 120 * fs

    # print("length of peaks:\t", len(peaks))

    while True:

        go = bisect.bisect_left(peaks, start)
        out = bisect.bisect_left(peaks, end)

        RRi = []
        peaks_in = []
        ff = 0
        peaks_in = peaks[go:out]

        for i in range(len(peaks_in) - 1):
            new = peaks_in[i + 1] - peaks_in[i]
            new = (new / fs)

            RRi.append(new)
            if (new > 0.05):
                ff = ff+1

        NN50.append(ff)

        approximate.append(ant.app_entropy(RRi))

        start = start + (10 * fs)
        end = end + (10 * fs)
        i += 1

        if (end > peaks[-1]):
            break

    return approximate, NN50
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
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]
def thresholding(array, thresh):
    i = 0
    j = 5
    Thresh = []

    while (True):
        go_in = []

        go_in = array[i:j]
        mean = np.mean(go_in)

        percent = (mean / 100) * thresh
        Thresh.append(mean + percent)

        j = j + 1
        i = i + 1

        if j > len(array):
            break


    return Thresh
def thresh_values(array):

    Thresh = []
    per = []
    x = 0

    for i in range (60, 100):
        per.append(i)
        it_n = i

    while (True):

        mean = np.mean(array)
        percent = (mean/100) * per[x]
        Thresh.append(mean + percent)

        x = x + 1

        if x > len(per):
            break

    return Thresh
def thresh(array, thresh):

    mean = np.mean(array)
    percent = (mean / 100) * thresh

    return (mean + percent)
def code1(array1, array2):

    print("Workin on the inter-ictal period")
    per = []
    T = []

    for i in range (40, 101):
        per.append(i)

    # print("per values are \t", per)
    # print('length of the input array\t', )

    for i in range(len(per)):
        ff = 0
        seg = 0

        ii = 0
        j = 5

        print("testing precentage\t", per[i])

        while (True):

            # print("length\t", len(array1) - j)

            go_in_A = array1[ii:j]
            go_in_N = array2[ii:j]

            print("segment", seg)
            seg = seg + 1

            thresh_A = thresh(go_in_A, per[i])
            thresh_N = thresh(go_in_N, per[i])

            if (ff >= 2):
                ff = 0
                print("break")
                break

            if ((thresh_A <= go_in_A[-1]) and (thresh_N <= go_in_N[-1])):
                ff = ff + 1
                print("above the threshold")

            if j > len(array1):
                # print("condition on")
                seg = 0

                if (ff < 2):
                    T.append(per[i])

                break

            ii = ii + 1
            j = j + 1



    print("the percentage extracted are :\t", T)
    return T
def code2(array1, array2, T):

    print("Working on the pre-ictal period")

    index = []
    final = []
    # ii = 0
    # j = 5
    seg = 0
    test = 0
    # == testing on the pre-ictal period
    for i in range(len(T)):

        seg = 0
        ii = 0
        j = 5

        cond = False
        print("###################################")

        while (True):

            ff = 0
            go_in_A = array1[ii:j]
            go_in_N = array2[ii:j]
            seg = seg + 1

            print("pourcentage working on\t", T[i])
            print("segments number\t", seg)

            rrr = T[i]

            thresh_A = thresh(go_in_A, rrr)

            thresh_N = thresh(go_in_N, rrr)

            if ((thresh_A <= go_in_A[-1]) and (thresh_N <= go_in_N[-1])):
                ff = ff + 1

                print("thresh A\t", thresh_A)
                print("go_in_a \t", go_in_A[-1])

                print("############################")

                print("thresh N\t", thresh_N)
                print("go_in_N \t", go_in_N[-1])

                cond = True
                print("condition true")
                # print("###################################")
                final.append(T[i])
                index.append(seg)
                # break

            if (cond == True):
                seg = 0
                cond = False
                break

            if j > len(array1) or (ff > 0):
                print("break")

                if (ff > 0 ):
                    final.append(T[i])
                    # index.append(test)
                    seg = 0
                    # print("segment equals to\t", seg)

                ff = 0
                seg = 0

                break

            ii = ii + 1
            j = j + 1

    print("finals are\t", final)
    print("index are\t", index)
    if not final:
        final = T

    return final, index
def threshold_finding(array1, array2, array3, array4):

    T = code1(array1, array2)

    Y = code2(array3, array4, T)

    print("Y equals to \t", Y)
    teste = Y[0]
    index = Y[1]
    # index = 0
    mini = index[0]

    for i in range(0, len(index)):
        # Compare elements of array with min
        if (index[i] < mini):
            mini = index[i];


    iidex = []

    pourc = Y[0]
    xtt = []
    for i in range (len(teste)):
        if(index[i] == mini):
            xtt.append(pourc[i])
            iidex.append(index[i])

    return xtt, iidex
def threshold_finding_pre(array1, array2, T):

    print("thresh pre")
    print("T\t", T)

    XR = T

    Y = code2(array1, array2, XR)

    if (len(Y[0]) > 0):

        print("Y equals to \t", Y)
        teste = Y[1]

        mini = teste[0]

        for i in range(0, len(teste)):
            # Compare elements of array with min
            if (teste[i] <= mini):
                mini = teste[i];

        pourc = Y[0]
        xtt = []
        for i in range(len(teste)):
            if (teste[i] == mini):
                xtt.append(pourc[i])
    else:
        xtt = T
        print("no detection was made for this acquisition !!")

    return xtt
def reading_EDF(path, ID):

    with open(path + "peaks-" + ID + ".txt", 'r') as file1:
        peaks = [float(i) for line in file1 for i in line.split('\n') if i.strip()]

    return peaks
def test_pre(path, ID, fs, T):

    peak = reading_EDF(path, ID)

    app_pre, NN50_pre = features(peak, fs)

    STD_app = std_compute(app_pre)
    STD_NN50 = std_compute(NN50_pre)


    print("lenghgtt equals to\t", len(STD_app))

    print("T equals to \t", T)

    XX = threshold_finding_pre(STD_app, STD_NN50, T)

    return XX

# =============================================================================================================

# path = "C:\\Users\\Ftay\\Desktop\\PhD\\tests\\siena_vs_fantasia\\unhealthy\\"
# path = "E:\\data\\tests\\siena_vs_fantasia\\Peaks_RR\\PN06\\"
# ID = "PN06-2"

# == inter-ictal period
# path = "E:\\data\\tests\\Peaks_RR\\Interictal\\PN17\\"
# path = "E:\\data\\tests\\Peaks_RR\\siena\\Interictal\\PN10\\"
# path = "E:\\data\\tests\\Peaks_RR\\fabrice\\inter-ictal\\PN07\\"
# path = "E:\\data\\tests\\siena_vs_fantasia\\Peaks_RR\\Interictal\\PN00\\"
# fs = 512
#
# ID = "PN10-2-0_60"
# ID = "PN00-3-30"

# peaks = reading_EDF(path, ID)
# # = =
# ID = "PN03-06-inter-ictal"
# peaks1 = reading_EDF(path, ID)
# # ==
# ID = "PN17-2-0_60"
# peaks2 = reading_EDF(path, ID)
# = =
# ID = "PN13-2-0_60"
# peaks3 = reading_EDF(path, ID)
# # = =
# ID = "PN13-2-30"
# peaks4 = reading_EDF(path, ID)
# # = =
# ID = "PN13-3-30"
# peaks5 = reading_EDF(path, ID)

# NN50, app, NN501, app1, NN502, app2 = ([] for i in range(6))
# app, NN50 = features(peaks,fs)
# app1, NN501 = features(peaks1,fs)
# app2, NN502 = features(peaks2,fs)
# app3, NN503 = features(peaks3,fs)
# app4, NN504 = features(peaks4,fs)
# app5, NN505 = features(peaks5,fs)

# app = np.append(app, app1)
# app = np.append(app, app2)
# app = np.append(app, app3)
# app = np.append(app, app4)
# app = np.append(app, app5)

# NN50 = np.append(NN50, NN501)
# NN50 = np.append(NN50, NN502)
# NN50 = np.append(NN50, NN503)
# NN50 = np.append(NN50, NN504)
# NN50 = np.append(NN50, NN505)
########################################################################################################################
# == pre-ictal periods
path = "E:\\data\\tests\\Peaks_RR\\fabrice\\pre-ictal\\PN06\\"
# path = "E:\\data\\tests\\Peaks_RR\\fabrice\\inter-ictal\\PN06\\"
# path = "E:\\data\\tests\\siena_vs_fantasia\\Peaks_RR\\PN00\\"
# fs = 512
fs = 256

ID = "211109C-AEX_0004"
# ID = "PN06-03-inter-ictal"
# ID = "signal-PN08-180314D-AEX_0001"

percentage = 53

seizure = 13 + (35/60)

peaks_pre = reading_EDF(path, ID)

# path = "E:\\data\\tests\\Peaks_RR\\fabrice\\pre-ictal\\PN06\\"
# ID = "211109C-AEX_0004"
#
# pre = reading_EDF(path, ID)
#
# peaks_pre = np.append(peaks_pre, pre)
#
# print("lenght of the pre-ictal period is:\t", )

app_pre, NN50_pre = features(peaks_pre,fs)

# =============================================================================================================

# path = 'E://data//post-ictal-heart-rate-oscillations-in-partial-epilepsy-1.0.0//sz07'
# ID = "sz03"
# fs = 200
# # # =====================
# #  01:02:43
#
# signal_input, peaks, fs = read_wfdb(path)
# sig_len = len(signal_input)
# print("sig_len\t", sig_len)
# print("length of the input signal is \t", ((len(signal_input)/fs)/60))
#
# start = (fs * 60 * 60 * 0) + (fs * 60 * 0) + (34 * fs)
# end = (fs * 60 * 60 * 0) + (fs * 60 * 24) + (34 * fs)
#
# start = (fs * 60 * 60 * 0) + (fs * 60 * 24) + (36 * fs)
# # # end = (fs * 60 * 60 * 1) + (fs * 60 * 55) + (43 * fs)
# print("start\t", start)
#
# # signal_input, peaks, fs = read_wfdb_interictal(path,start, end)
# signal_input, peaks, fs = read_wfdb_interictal(path,start, sig_len)
#
# peaks = [x - start for x in peaks]
# app, NN50 = features(peaks,fs)
#
# # === The pre-ictal period
#
# ## 02:55:51
# start = (fs * 60 * 60 * 0) + (fs * 60 * 0) + (36 * fs)
# end = (fs * 60 * 60 * 0) + (fs * 60 * 24) + (36 * fs)
#
# # start= (fs * 60 * 60 * 1) + (fs * 60 * 4) + (45 * fs)
# # end = (fs * 60 * 60 * 2) + (fs * 60 * 14) + (45 * fs)
#
# signal_input, peaks_re, fs = read_wfdb_interictal(path, start, end)
#
# peaks_re = [x - start for x in peaks_re]
#
# print("pre peak equals to\t", peaks_re)
#
# NN50_pre, app_pre = ([] for i in range(2))
# app_pre, NN50_pre = features(peaks_re,fs)

# =============================================================================================================
# # == loading epileptic patients acquisitions taken from INS " DAvid Olivier"
# path = "E:\\data\\tests\\INS - David\\"
# fs = 512
# ID = "0252GRE-EEG_13"
# ID = "0252GRE-EEG_14"
# with open(path + "peaks-" + ID + ".txt", 'r') as file1:
#     peaks_pre = [float(i) for line in file1 for i in line.split('\n') if i.strip()]
#
# peaks_pre = reading_EDF(path, ID)
# app_pre, NN50_pre = features(peaks_pre,fs)
#
# # ============================================================================================================
#
# path = "E:\\data\\tests\\INS - David\\inter-ictal\\"
# ID = "0252GRE-end-EEG-13"
# ID = "0252GRE-end-EEG-14"
# # ID = "0047GRE-60-125"
# # ID = "0047GRE-125-194"
# fs = 512
# # == load the epileptic patient
# # with open(path + "peaks-" + ID + "-end-EEG-14.txt", 'r') as file1:
#
# path = "E:\\data\\tests\\Peaks_RR\\siena\\pre-ictal\\PN06\\"
# ID = "PN06-2"
#
# with open(path + "peaks-" + ID + ".txt", 'r') as file1:
#     peaks = [float(i) for line in file1 for i in line.split('\n') if i.strip()]

# app, NN50 = features(peaks,fs)
# seizure = 6.21
# seizure = 5 + (49/60)

# =============================================================================================================
pre_ictal = []
STD = []
STD_healthy = []

# STD_app = std_compute(app)
# STD_NN50 = std_compute(NN50)

STD_app_pre = std_compute(app_pre)
STD_NN50_pre = std_compute(NN50_pre)

# percentage = 62

# STD_app_test = STD_app_pre[0: (len(STD_app_pre)-10)]
# STD_NN50_test = STD_NN50_pre[0: (len(STD_NN50_pre)-10)]

# XX = threshold_finding(STD_app, STD_NN50, STD_app_pre, STD_NN50_pre)
#
# IDs = ["PN10-2", "PN10-3", "PN10-4", "PN10-5", "PN10-6", "PN10-7", "PN10-8", "PN10-9", "PN10-10"]
# #
# T = XX[0]
#
# for i in range(len(IDs)):
#
#     # path = "E:\\data\\tests\\Peaks_RR\\fabrice\\pre-ictal\\PN10\\"
#     path = "E:\\data\\tests\\siena_vs_fantasia\\Peaks_RR\\PN10"
#     print("######################################################################")
#     print("working on the acquisition number\t", IDs[i])
#
#     XV = test_pre(path, IDs[i], fs, T)
#
#     T = XV
#
# print("final T equals to", T)
# print("final pourcentage to use\t", T[0])
# == working on the pre-ictal or inter-ictal
#

app = app_pre
NN50 = NN50_pre

STD_app = STD_app_pre
STD_NN50 = STD_NN50_pre

# print("seizure\t", len(STD_app) - seizure)

thresh_AP = thresholding(STD_app , percentage)
thresh_nn = thresholding(STD_NN50, percentage)

first = 0
tt = [first]
i = 1
while i <= len(app) - 1:
    first = first + 10
    tt.append(first)
    i = i + 1

first = 0
tt_healthy = [first]
i = 1
while i <= len(STD_app) - 1:
    first = first + 1
    tt_healthy.append(first)
    i = i + 1

# == computing the adaptive threshold values

# thresh_AP = thresh(STD_app, 82)
# thresh_nn = thresh(STD_NN50, 82)

first = 0
tt = [first]
i = 1
while i <= len(app) - 1:
    first = first + 10
    tt.append(first)
    i = i + 1

first = 0
tt_healthy = [first]
i = 1
while i <= len(STD_app) - 1:
    first = first + 1
    tt_healthy.append(first)
    i = i + 1

# ____ woerking on the approximate thresholding
x = numpy.zeros(len(STD_app), dtype=float, order='C')
xy = numpy.zeros(len(app), dtype=float, order='C')

arr=[]
arr1=[]

ln = int(len(STD_app))
xx = 0

i = 0

while(True):

    if (i == 0):
        arr.append(thresh_AP[i])
        arr.append(thresh_AP[i])
        arr.append(thresh_AP[i])
        arr.append(thresh_AP[i])
        arr.append(thresh_AP[i])

        i = i + 1

        arr1.append(0)
        arr1.append(0)
        arr1.append(0)
        arr1.append(0)
        arr1.append(0)
    else:
        arr.append(thresh_AP[i])
        arr1.append(0)

        i = i + 1
        # print("i value\t", i)
        # print("left\t", len(STD_app) - i)

    if (i == len(STD_app) - 4):
        break

first_line = LineString(np.column_stack((tt_healthy,STD_app)))
second_line = LineString(np.column_stack((tt_healthy,arr)))
intersection = first_line.intersection(second_line)
x, y = LineString(intersection).xy
print(" x \t", x)

arr_n=[]
arr1_n=[]

ln = int (len(STD_NN50))

# for i in range (ln):
#     if (i == 0):
#         arr_n.append(thresh_nn[i])
#         arr_n.append(thresh_nn[i])
#         arr_n.append(thresh_nn[i])
#         arr_n.append(thresh_nn[i])
#         # arr_n.append(thresh_nn[i])
#
#         arr1_n.append(0)
#         arr1_n.append(0)
#         arr1_n.append(0)
#         arr1_n.append(0)
#         # arr1_n.append(0)
#
#     else:
#         arr.append(thresh_AP[i])
#         arr1.append(0)

i = 0
while(True):

    if (i == 0):
        arr_n.append(thresh_nn[i])
        arr_n.append(thresh_nn[i])
        arr_n.append(thresh_nn[i])
        arr_n.append(thresh_nn[i])
        arr_n.append(thresh_nn[i])

        i = i + 1
        # print("i value\t", i)
        arr1_n.append(0)
        arr1_n.append(0)
        arr1_n.append(0)
        arr1_n.append(0)
        arr1_n.append(0)
    else:
        arr_n.append(thresh_nn[i])
        arr1_n.append(0)

        i = i + 1
        # print("i value\t", i)
        # print("left\t", len(STD_NN50) - i)

    if (i == len(STD_NN50) - 4):
        break

first_line_n = LineString(np.column_stack((tt_healthy, STD_NN50)))
second_line_n = LineString(np.column_stack((tt_healthy, arr_n)))
intersection_n = first_line_n.intersection(second_line_n)
x_n, y_n = LineString(intersection_n).xy
print((x_n))

round_app = []
round_NN50 = []

# =============================================================================================================
# == plotting the results

fig, axs = plt.subplots(2, 1)
axs[0].plot(tt, app, label="approximate entropy feature", marker='o')
axs[1].plot(tt_healthy,STD_app, label="STD of approximate entropy feature", marker='o')
axs[0].axvline(x=tt[-1] - 600, color='red', linestyle='--')
axs[1].axvline(x=len(STD_app) - 5, color='red', linestyle='--')

axs[0].set_title('Acquisition ' + ID, fontsize=24, y=1)

axs[1].plot(tt_healthy, arr, color='blue', linestyle='--')

axs[0].set_xlabel('Time in min')
axs[0].set_ylabel('entropy value')
axs[1].set_xlabel('Time in min')
axs[1].set_ylabel('STD value')
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


# == plotting NN50
fig, axs = plt.subplots(2, 1)
axs[0].plot(tt_healthy,STD_app, label="STD of approximate entropy feature", marker='o')
axs[1].plot(tt_healthy, STD_NN50, label="STD of NRRi feature", marker='o')

axs[0].set_title('Acquisition ' + ID,fontsize=24, y=1)

# axs[0].axvline(x= len(STD_app) - 10, color='red', linestyle='--')
# axs[1].axvline(x= len(STD_NN50) - 10, color='red', linestyle='--')

axs[0].axvline(x= (len(STD_app) - seizure), color='red', linestyle='--')
axs[1].axvline(x= (len(STD_app) - seizure), color='red', linestyle='--')

print("szireu time\t", len(STD_app) - seizure)
# axs[1].axvline(x= (len(STD_NN50) - seizure), color='red', linestyle='--')

axs[0].set_xlabel('Time in min')
axs[0].set_ylabel('STD value')
axs[1].set_xlabel('Time in min')
axs[1].set_ylabel('STD value')
axs[0].grid(True)
axs[1].grid(True)
axs[0].legend()
axs[1].legend()

axs[0].plot(tt_healthy, arr, color='g', linestyle='--')
axs[1].plot(tt_healthy, arr_n, color='g', linestyle='--')

# ymin_1 = 0; ymax_1 = 0.7
# axs[0].set_ylim([ymin_1,ymax_1])
#
# ymin_1 = 0; ymax_1 = 0.70
# axs[1].set_ylim([ymin_1,ymax_1])

if intersection.geom_type == 'MultiPoint':
    axs[0].plot(*LineString(intersection).xy, 'o')

elif intersection.geom_type == 'Point':
    axs[0].plot(*intersection.xy, 'o')


if intersection_n.geom_type == 'MultiPoint':
    axs[1].plot(*LineString(intersection_n).xy, 'o')

elif intersection_n.geom_type == 'Point':
    axs[1].plot(*intersection_n.xy, 'o')


# == plotting NN50
fig, axs = plt.subplots(2, 1)
axs[0].plot(tt, NN50, label=" NRRi feature", marker='o')
axs[1].plot(tt_healthy, STD_NN50, label="STD of NRRi feature", marker='o')
axs[0].set_title('Subject ' + ID,fontsize=24, y=1)

axs[0].axvline(x=tt[-1] - 600, color='red', linestyle='--')
axs[1].axvline(x= len(STD_NN50) - 10, color='red', linestyle='--')

axs[0].set_xlabel('Time in min')
axs[0].set_ylabel('NRRi value')
axs[1].set_xlabel('Time in min')
axs[1].set_ylabel('STD value')
axs[0].grid(True)
axs[1].grid(True)
axs[0].legend()
axs[1].legend()

axs[1].plot(tt_healthy, arr_n, color='blue', linestyle='--')

ymin = 70 ; ymax = 275
axs[0].set_ylim([ymin,ymax])

ymin_1 = 0; ymax_1 = 14
axs[1].set_ylim([ymin_1,ymax_1])

if intersection_n.geom_type == 'MultiPoint':
    axs[1].plot(*LineString(intersection_n).xy, 'o')
elif intersection_n.geom_type == 'Point':
    axs[1].plot(*intersection_n.xy, 'o')

plt.show()



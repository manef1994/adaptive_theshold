import bisect
import math
import statistics
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
            # print(new)
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
    ff = 0
    index = []
    per = []
    T = []
    it_n = 90 - 70

    for i in range (70, 100):
        per.append(i)

    # print("per values are \t", per)
    # print('length of the input array\t', )

    for i in range(len(per)):
        ff = 0
        seg = 0

        ii = 0
        j = 5

        # print("testing precentage\t", per[i])

        while (True):

            # print("length\t", len(array1) - j)

            go_in_A = array1[ii:j]
            go_in_N = array2[ii:j]

            # print("segment", seg)
            seg = seg + 1

            thresh_A = thresh(go_in_A, per[i])
            thresh_N = thresh(go_in_N, per[i])

            for y in range(len(go_in_N)):

                if (ff > 1):
                    ff = 0
                    break

                if ((thresh_A <= go_in_A[y]) and (thresh_N <= go_in_N[y])):
                    ff = ff + 1
                    # print("above the threshold")

            if j > len(array1):
                # print("condition on")
                seg = 0

                if (ff == 0):
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



            for y in range(len(go_in_N)):


                if ((thresh_A <= go_in_A[y]) and (thresh_N <= go_in_N[y])):


                    ff = ff + 1

                    print("thresh A\t", thresh_A)
                    print("go_in_a \t", go_in_A[y])

                    print("############################")

                    print("thresh N\t", thresh_N)
                    print("go_in_N \t", go_in_N[y])



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

    return final

def threshold_finding(array1, array2, array3, array4):

    T = code1(array1, array2)

    Y = code2(array3, array4, T)

    print("Y equals to \t", Y)
    teste = Y

    index =0

    mini = teste[0]

    for i in range(0, len(teste)):
        # Compare elements of array with min
        if (teste[i] <= mini):
            mini = teste[i];

    pourc = Y
    xtt = []
    for i in range (len(teste)):
        if(teste[i] == mini):
            xtt.append(pourc[i])

    return xtt, index

def threshold_finding_pre(array1, array2, T):

    print("thresh pre")
    print("T\t", T)

    Y = code2(array1, array2, T)

    print("Y equals to \t", Y)
    teste = Y

    mini = teste[0]

    for i in range(0, len(teste)):
        # Compare elements of array with min
        if (teste[i] <= mini):
            mini = teste[i];

    pourc = Y
    xtt = []
    for i in range (len(teste)):
        if(teste[i] == mini):
            xtt.append(pourc[i])

    return xtt

def reading_EDF(path, ID):

    with open(path + "peaks-" + ID + ".txt", 'r') as file1:
        peaks = [float(i) for line in file1 for i in line.split('\n') if i.strip()]

    return peaks

def test_pre(path, ID, fs, T):

    peak = reading_EDF(path, ID)

    app_pre, NN50_pre = features(peak, fs)

    print("lenghgtt equals to\t", len(app_pre))

    print("T equals to \t", T)

    XX = threshold_finding_pre(app_pre, NN50_pre, T)

    return XX

# =============================================================================================================
# # == loading epileptic patients acquisitions taken from INS " DAvid Olivier"
# path = "C:\\Users\\Ftay\Desktop\\PhD\\tests\\INS - David\\"
# fs = 512
# ID = "0252GRE-EEG_13"
# with open(path + "peaks-" + ID + ".txt", 'r') as file1:
#     peaks = [float(i) for line in file1 for i in line.split('\n') if i.strip()]

# path = "C:\\Users\\Ftay\Desktop\\PhD\\tests\\INS - David\\inter-ictal\\"
# ID = "0252GRE"
# fs = 512
# # == load the epileptic patient
# with open(path + "peaks-" + ID + "-end-EEG-14.txt", 'r') as file1:
#     peaks = [float(i) for line in file1 for i in line.split('\n') if i.strip()]

# thresh_AP = 0.06
# thresh_nn = 6

# =============================================================================================================
# path = "C:\\Users\\Ftay\\Desktop\\PhD\\tests\\siena_vs_fantasia\\Peaks_RR\\PN00\\"
# path = "C:\\Users\\Ftay\\Desktop\\PhD\\tests\\siena_vs_fantasia\\unhealthy\\"
# path = "E:\\data\\tests\\siena_vs_fantasia\\Peaks_RR\\PN06\\"
# ID = "PN06-2"

# == inter-ictal period
path = "E:\\data\\tests\\Peaks_RR\\Interictal\\PN06\\"
fs = 512
# = =
ID = "PN06-3"
peaks = reading_EDF(path, ID)
# = =
ID = "PN06-4"
peaks1 = reading_EDF(path, ID)
# = =
ID = "PN06-5"
peaks2 = reading_EDF(path, ID)
# = =

NN50, app, NN501, app1, NN502, app2 = ([] for i in range(6))
app, NN50 = features(peaks,fs)
app1, NN501 = features(peaks1,fs)
app2, NN502 = features(peaks2,fs)

app = np.append(app, app1)
app = np.append(app, app2)

NN50 = np.append(NN50, NN501)
NN50 = np.append(NN50, NN502)
########################################################################################################################
# == pre-ictal periods
path = "E:\\data\\tests\\siena_vs_fantasia\\Peaks_RR\\PN06\\"
fs = 512

ID = "PN06-3"
peaks_pre = reading_EDF(path, ID)

NN50_pre, app_pre = ([] for i in range(2))
app_pre, NN50_pre = features(peaks_pre,fs)

# thresh_AP = 0.14
# thresh_nn = 0.04

# =============================================================================================================
# NN50, app = ([] for i in range(2))
# NN50_pre, app_pre = ([] for i in range(2))

pre_ictal = []
STD = []
STD_healthy = []

STD_app = std_compute(app)
STD_NN50 = std_compute(NN50)

STD_app_pre = std_compute(app_pre)
STD_NN50_pre = std_compute(NN50_pre)

XX = threshold_finding(STD_app, STD_NN50, STD_app_pre,STD_NN50_pre)

IDs = ["PN06-3",  "PN06-4", "PN06-5"]

T = XX[0]
path = "E:\\data\\tests\\siena_vs_fantasia\\Peaks_RR\\PN06\\"
fs = 512

# for i in range(len(IDs)):
#
#     print("######################################################################")
#     print("working on the acquisition number\t", IDs[i])
#
#     XV = test_pre(path, IDs[i], fs, T)
#
#     T = XV

print("final T equals to", T)
# == working on the pre-ictal or inter-ictal


app = app_pre
NN50 = NN50_pre

STD_app = STD_app_pre
STD_NN50 = STD_NN50_pre


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

percentage = 73
thresh_AP = thresholding(STD_app , percentage)
thresh_nn = thresholding(STD_NN50, percentage)

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
# x, y = LineString(intersection).xy
# print(" x \t", x)

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
# x_n, y_n = LineString(intersection_n).xy
# print((x_n))

round_app = []
round_NN50 = []

# =============================================================================================================
# == plotting the results

fig, axs = plt.subplots(2, 1)
axs[0].plot(tt, app, label="approximate entropy of healthy signal", marker='o')
axs[1].plot(tt_healthy,STD_app, label="STD of approximate entropy", marker='o')
axs[0].axvline(x=tt[-1] - 600, color='red', linestyle='--')
axs[1].axvline(x=len(STD_app) - 10, color='red', linestyle='--')

axs[0].set_title('Subject ' + ID,fontsize=24, y=1)

axs[1].plot(tt_healthy, arr, color='blue', linestyle='--')

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


# == plotting NN50
fig, axs = plt.subplots(2, 1)
axs[0].plot(tt_healthy, STD_app, label="STD of ApEn", marker='o')
axs[1].plot(tt_healthy, STD_NN50, label="STD of NN50", marker='o')

axs[0].set_title('Subject ' + ID,fontsize=24, y=1)

# axs[0].axvline(x= len(STD_app) - 10, color='red', linestyle='--')
# axs[1].axvline(x= len(STD_NN50) - 10, color='red', linestyle='--')

axs[0].axvline(x= (len(STD_app) - 10), color='red', linestyle='--')
axs[1].axvline(x= (len(STD_NN50) - 10), color='red', linestyle='--')

axs[0].set_xlabel('sample per min')
axs[0].set_ylabel('entropy value')
axs[1].set_xlabel('sample per min')
axs[1].set_ylabel('entropy value')
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
axs[0].plot(tt, NN50, label=" NN50 curve", marker='o')
axs[1].plot(tt_healthy, STD_NN50, label="STD of NN50", marker='o')
axs[0].set_title('Subject ' + ID,fontsize=24, y=1)

axs[0].axvline(x=tt[-1] - 600, color='red', linestyle='--')
axs[1].axvline(x= len(STD_NN50) - 10, color='red', linestyle='--')

axs[0].set_xlabel('sample en seconde')
axs[0].set_ylabel('valeur de entropie')
axs[1].set_xlabel('sample en seconde')
axs[1].set_ylabel('valeur de entropie')
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

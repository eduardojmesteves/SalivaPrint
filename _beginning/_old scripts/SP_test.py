import math, collections, os, re, sys, glob, codecs
from pprint import pprint
import pandas as pd
import numpy as np
from scipy import signal
from matplotlib import pylab as plt
import matplotlib.pyplot as plot
from scipy.signal import find_peaks
from scipy import stats

# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------
# Script 1
# Capillary electrophoresis
# Data extraction from XML file (with .txt file associated)
# Data Export

profiles = {}
samples = 0

#Get the Raw data points
for file in glob.glob("*.xml"):
    profile_file = file

file = open(profile_file)
content = file.read()
file.close()
content = content[content.find('<RawDataPoints>'):]

# Gets the fluorescence points.
# aux = re.search(r"<Signal>(.*?)</Signal>", content)
aux = {}
ix = 0
key = 0
matches = re.finditer(r"<Signal>(.*?)</Signal>", content)

for matchNum, match in enumerate(matches):
    for groupNum in range(0, len(match.groups())):
        groupNum = groupNum + 1
        group = match.group(groupNum)
        if (ix % 1380 == 0):
            key += 1
        if key in aux.keys():
            aux[key].append(float(group))
        else:
            aux[key] = [float(group)]
        ix += 1
    matchNum = matchNum + 1

meta_file = open(profile_file.split('.')[0] + ".csv", 'w')
inds = aux.keys()

for i in range(0, 1380):
    vals = []
    for ind in inds:
        vals.append(aux[ind][i])

#Get the Sample IDs
ind_names = []
f1 = codecs.open(profile_file.split('.')[0] + "_RunSummary.txt", encoding='utf-16-le')
for line in f1:  # note, the readlines is not really needed
    if samples != 0:
        ind_names.append(line.split('\t')[1])  # the comma strips the trailing newline in case that's bothering you
    samples += 1


# Calculate the Molecular Weight
meta_file.write(',')
for i in range(466,870): # 471 nao é erro, tem de ser mais 1 do que 470
    if i < 869:
        meta_file.write(str(0.00042006*(i*i)-0.28216*i+49.225)+',')
    else:
        meta_file.write(str(0.00042006*(i*i)-0.28216*i+49.225)+'\n')

ix = 0

# write data
for ind in aux.keys():
    meta_file.write(ind_names[ix]+',')
    for i in range(465, 869):
        if i < 868:
            meta_file.write(str(aux[ind][i]) + ',')
        else:
            meta_file.write(str(aux[ind][i]) + '\n')
    ix +=1

meta_file.flush()
meta_file.close()

# Read in the file
with open(profile_file.split('.')[0]+".csv", 'r') as file :
  filedata = file.read()

# Write the file out again
with open(profile_file.split('.')[0]+".csv", 'w') as file:
  file.write(filedata)

# Read the file into a data frame
df = pd.read_csv(profile_file.split('.')[0] + ".csv",index_col=0)
print(df)

# # create the transpose
df_transposed = df.T
# save data to csv
print(df_transposed)
# df_transposed.to_csv('/home/eduardo/Desktop/t1.csv')

# # # -----------------------------------------------------------------------------------------------------------------------
# # # -----------------------------------------------------------------------------------------------------------------------
# # # -----------------------------------------------------------------------------------------------------------------------
# # # Script2
# # # Profile alignement

# # ### Peal Alignment ###

# # ###PROFILE IMPORT

# data = pd.read_csv('/home/eduardo/Desktop/t1.csv', index_col=0, skiprows=1)
# print(data)

# ID = data.ID
# print(ID)

# data = data.drop(columns=['ID'])
# # data.to_csv('/home/eduardo/Desktop/align.csv')
# # print(data)

# all_profiles = list(data.columns.values)
# master = all_profiles[0]

# master_profile = data[master]

# all_profiles_minus_actual = all_profiles.copy()
# all_profiles_minus_actual.remove(master)

# targets = data[all_profiles_minus_actual]

# ###Peak alignment

# f1 = plt.figure(0)
# for target in targets:
#     plt.plot(mw, targets[target])

# f2 = plt.figure(1)
# for target in targets:
#     dx = np.mean(np.diff(mw.values))
#     shift = (np.argmax(signal.correlate(master_profile, targets[target])) - len(targets[target])) * dx
#     plt.plot(mw + shift, targets[target])
# plt.show()


# # -----------------------------------------------------------------------------------------------------------------------
# # -----------------------------------------------------------------------------------------------------------------------
# # -----------------------------------------------------------------------------------------------------------------------
# # Script3
# # Area Reckon

# # def sum_fluorescence_per_ind(profile, peak_intervals):
# #     sum1 = 0
# #     sum2 = 0
# #     sum3 = 0
# #     sum4 = 0
# #     sum5 = 0
# #     sum6 = 0

# #     to_consider1 = list(range(int(peak_intervals[0][0]), int(peak_intervals[0][1])))
# #     to_consider2 = list(range(int(peak_intervals[1][0]), int(peak_intervals[1][1])))
# #     to_consider3 = list(range(int(peak_intervals[2][0]), int(peak_intervals[2][1])))
# #     to_consider4 = list(range(int(peak_intervals[3][0]), int(peak_intervals[3][1])))
# #     to_consider5 = list(range(int(peak_intervals[4][0]), int(peak_intervals[4][1])))
# #     to_consider6 = list(range(int(peak_intervals[5][0]), int(peak_intervals[5][1])))

# #     for i in range(0, len(profile)):
# #         if i in to_consider1:
# #             sum1 += profile[i]

# #     for i in range(0, len(profile)):
# #         if i in to_consider2:
# #             sum2 += profile[i]

# #     for i in range(0, len(profile)):
# #         if i in to_consider3:
# #             sum3 += profile[i]

# #     for i in range(0, len(profile)):
# #         if i in to_consider4:
# #             sum4 += profile[i]

# #     for i in range(0, len(profile)):
# #         if i in to_consider5:
# #             sum5 += profile[i]

# #     for i in range(0, len(profile)):
# #         if i in to_consider6:
# #             sum6 += profile[i]

# #     return [sum1, sum2, sum3, sum4, sum5, sum6]


# # def find_nearest(array, value):
# #     idx = np.searchsorted(array, value, side="left")
# #     if idx > 0 and (idx == len(array) or math.fabs(value - array[idx - 1]) < math.fabs(value - array[idx])):
# #         return array[idx - 1]
# #     else:
# #         return array[idx]


# # def find_half_peak_from_top(peak, profile, index_top):
# #     margin_bottom = 0
# #     margin_top = len(profile)

# #     if (peak == '1'):
# #         margin_bottom = 2
# #         margin_top = peaks[str(int(peak) + 1)][0]
# #     if (peak == '2'):
# #         margin_bottom = peaks[str(int(peak) - 1)][1]
# #         margin_top = peaks[str(+int(peak) + 1)][0]
# #     if (peak == '3'):
# #         margin_bottom = peaks[str(int(peak) - 1)][1]
# #         margin_top = peaks[str(+int(peak) + 1)][0]
# #     if (peak == '4'):
# #         margin_bottom = peaks[str(int(peak) - 1)][1]
# #         margin_top = peaks[str(+int(peak) + 1)][0]
# #     if (peak == '5'):
# #         margin_bottom = peaks[str(int(peak) - 1)][1]
# #         margin_top = peaks[str(+int(peak) + 1)][0]
# #     if peak == '6':
# #         margin_bottom = peaks[str(int(peak) - 1)][1]
# #         # margin_top = peaks[str(+int(peak)+1)][0]
# #         margin_top = 340

# #     half_value = profile[index_top] / 2

# #     left_half = profile[:index_top]
# #     right_half = profile[index_top:]

# #     closest_right = 0
# #     distance_right = 9999999999

# #     # go_right
# #     for v in range(index_top, len(profile)):
# #         if profile[v] - half_value < distance_right:
# #             distance_right = profile[v] - half_value
# #             closest_right = v
# #         if (profile[v] < half_value) or (v > margin_top):
# #             break

# #     closest_left = 0
# #     distance_left = 9999999999

# #     # go_left
# #     for v in range(index_top, 0, -1):
# #         if profile[v] - half_value < distance_left:
# #             distance_left = profile[v] - half_value
# #             closest_left = v
# #         if (profile[v] < half_value) or (v < margin_bottom):
# #             break

# #     return [closest_left, closest_right]


# # f = open('/home/eduardo/Desktop/area1.csv')
# # d = {}

# # lines = f.readlines()

# # for line in lines:
# #     line = line.strip()
# #     line = line.replace(',', '.')
# #     aux = line.split(';')
# #     k = aux[0]
# #     val = aux[1:]
# #     values = list(map(float, val))
# #     d[k] = values

# # peaks = {'1': [48, 56], '2': [112, 127], '3': [136, 161], '4': [178, 199], '5': [241, 249], '6': [262, 276]}

# # for ind in d.keys():

# #     # print(ind)
# #     plot.figure()

# #     plot.plot(d[ind])

# #     peak_intervals = []

# #     for peak in peaks.keys():
# #         # Janela inicial
# #         [init, fim] = peaks[peak]
# #         window = d[ind][init - 1:fim]
# #         xs = range(init - 1, fim)
# #         index_max_window = window.index(max(window))
# #         plot.plot(xs, window, 'r')

# #         # Maximums
# #         index_max = init - 1 + index_max_window
# #         plot.plot(index_max, window[index_max_window], 'go')

# #         boundaries = find_half_peak_from_top(peak, d[ind], index_max)

# #         plot.plot(boundaries[0], d[ind][boundaries[0]], 'b*')
# #         plot.plot(boundaries[1], d[ind][boundaries[1]], 'y*')

# #         peak_intervals.append(boundaries)

# #     # print(peak_intervals)
# #     print(str(sum_fluorescence_per_ind(d[ind], peak_intervals)))

# #     plot.show()
# #     # plot.savefig ("area"+ind+".png")
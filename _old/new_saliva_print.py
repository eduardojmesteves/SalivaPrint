#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os, csv, re, glob, math, collections, sys, codecs
import xml.etree.ElementTree as ET
from matplotlib import pylab as plt
import matplotlib.pyplot as plot
import pandas as pd
from io import BytesIO
import numpy as np
import pandas as pd
from scipy import signal
from scipy.signal import find_peaks
from scipy import stats

# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------
# Script1
# Createa csv data from xml source
full_xml_path = os.path.abspath(os.path.join(".", "test.xml"))
full_csv_path = os.path.abspath(os.path.join(".", "Samples.csv"))


tree = ET.parse(full_xml_path)
root = tree.getroot()

file_header = ['Sample']
file_header_vals = []
sample_well_info = []

# Set File Header
for i in range(476, 869):
    file_header.append(str(round((0.00042006 * (i * i) - 0.28216 * i + 49.225), 1)))
    file_header_vals.append(str(round((0.00042006 * (i * i) - 0.28216 * i + 49.225), 1)))

# get all well id's and sample names from file
for sample in root.findall('Results'):
    # write sample name
    well_id = sample.find('WellId').text

    if well_id != '0':
        sample_name = sample.find('SampleName').text
        sample_well_info.append({'Sample': sample_name, 'well_id': well_id})

# create file with headers and data
with open(full_csv_path, 'w', newline='') as outcsv:
    writer = csv.DictWriter(outcsv, fieldnames = file_header)
    writer.writeheader()
    xmlstr = ET.tostring(root, encoding='utf8', method='xml').decode("utf-8")

    for elem in sample_well_info:
        
        idx = 0
        row = {'Sample': elem['Sample']}
        regex= '<RawDataPoints><WellId>%s<\/WellId><Signal>(.*?)<\/Signal><\/RawDataPoints>'% elem['well_id']
        raw_data_points = re.findall(regex, xmlstr)
        
        # write data for each column
        for count, raw_data_point in enumerate(raw_data_points, start=0):
            if 869 >= count >= 477:
                row[file_header_vals[idx]] = raw_data_point
                idx+=1
            
            if count > 869:
                break
        
        writer.writerow(row)
'''
# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------
# Script2
# Profile alignement

### Peal Alignment ###

###PROFILE IMPORT
"""
df = pd.read_csv('Samples.csv')
data = df.drop(columns=['Sample'])
data = data.T

all_profiles = list(data.columns.values)
master = all_profiles[0]
master_profile = data[master]

all_profiles_minus_actual = all_profiles.copy()
all_profiles_minus_actual.remove(master)
targets = data[all_profiles_minus_actual]

###Peak alignment

f1 = plt.figure(0)
for target in targets:
    plt.plot(samples, targets[target])

f2 = plt.figure(1)
for target in targets:
    dx = np.mean(np.diff(samples.values))
    shift = (np.argmax(signal.correlate(master_profile, targets[target])) - len(targets[target])) * dx
    plt.plot(samples + shift, targets[target])
plt.show()
"""
# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------
# Script3
# Area Reckon (manual interval decision)

def sum_fluorescence_per_ind(profile, peak_intervals):
    sum1 = 0
    sum2 = 0
    sum3 = 0
    sum4 = 0
    sum5 = 0
    sum6 = 0

    to_consider1 = list(range(int(peak_intervals[0][0]), int(peak_intervals[0][1])))
    to_consider2 = list(range(int(peak_intervals[1][0]), int(peak_intervals[1][1])))
    to_consider3 = list(range(int(peak_intervals[2][0]), int(peak_intervals[2][1])))
    to_consider4 = list(range(int(peak_intervals[3][0]), int(peak_intervals[3][1])))
    to_consider5 = list(range(int(peak_intervals[4][0]), int(peak_intervals[4][1])))
    to_consider6 = list(range(int(peak_intervals[5][0]), int(peak_intervals[5][1])))

    for i in range(0, len(profile)):
        if i in to_consider1:
            sum1 += profile[i]

    for i in range(0, len(profile)):
        if i in to_consider2:
            sum2 += profile[i]

    for i in range(0, len(profile)):
        if i in to_consider3:
            sum3 += profile[i]

    for i in range(0, len(profile)):
        if i in to_consider4:
            sum4 += profile[i]

    for i in range(0, len(profile)):
        if i in to_consider5:
            sum5 += profile[i]

    for i in range(0, len(profile)):
        if i in to_consider6:
            sum6 += profile[i]

    return [sum1, sum2, sum3, sum4, sum5, sum6]


def find_nearest(array, value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx - 1]) < math.fabs(value - array[idx])):
        return array[idx - 1]
    else:
        return array[idx]


def find_half_peak_from_top(peak, profile, index_top):
    margin_bottom = 0
    margin_top = len(profile)

    if (peak == '1'):
        margin_bottom = 2
        margin_top = peaks[str(int(peak) + 1)][0]
    if (peak == '2'):
        margin_bottom = peaks[str(int(peak) - 1)][1]
        margin_top = peaks[str(+int(peak) + 1)][0]
    if (peak == '3'):
        margin_bottom = peaks[str(int(peak) - 1)][1]
        margin_top = peaks[str(+int(peak) + 1)][0]
    if (peak == '4'):
        margin_bottom = peaks[str(int(peak) - 1)][1]
        margin_top = peaks[str(+int(peak) + 1)][0]
    if (peak == '5'):
        margin_bottom = peaks[str(int(peak) - 1)][1]
        margin_top = peaks[str(+int(peak) + 1)][0]
    if peak == '6':
        margin_bottom = peaks[str(int(peak) - 1)][1]
        margin_top = peaks[str(+int(peak)+1)][0]
        margin_top = 340

    half_value = profile[index_top] / 2

    left_half = profile[:index_top]
    right_half = profile[index_top:]

    closest_right = 0
    distance_right = 9999999999

    # go_right
    for v in range(index_top, len(profile)):
        if profile[v] - half_value < distance_right:
            distance_right = profile[v] - half_value
            closest_right = v
        if (profile[v] < half_value) or (v > margin_top):
            break

    closest_left = 0
    distance_left = 9999999999

    # go_left
    for v in range(index_top, 0, -1):
        if profile[v] - half_value < distance_left:
            distance_left = profile[v] - half_value
            closest_left = v
        if (profile[v] < half_value) or (v < margin_bottom):
            break

    return [closest_left, closest_right]


f = open('Samples.csv')
d = {}

lines = f.readlines()

for line in lines:
    line = line.strip()
    aux = line.split(',')
    k = aux[0]
    val = aux[1:]
    values = list(map(float, val))
    d[k] = values

peaks = {'1': [48, 57], '2': [112, 129], '3': [136, 166], '4': [178, 201], '5': [241, 251], '6': [263, 278]}

for ind in d.keys():

    # print(ind)
    plot.figure()

    plot.plot(d[ind])

    peak_intervals = []

    for peak in peaks.keys():
        # Janela inicial
        [init, fim] = peaks[peak]
        print([init, fim])
        print(d[ind])
        window = d[ind][init - 1:fim]
        print(window)
        xs = range(init - 1, fim)
        index_max_window = window.index(max(window))
        print(index_max_window)
        plot.plot(xs, window, 'r')

        # Maximums
        index_max = init - 1 + index_max_window
        plot.plot(index_max, window[index_max_window], 'go')

        boundaries = find_half_peak_from_top(peak, d[ind], index_max)

        plot.plot(boundaries[0], d[ind][boundaries[0]], 'b*')
        plot.plot(boundaries[1], d[ind][boundaries[1]], 'y*')

        peak_intervals.append(boundaries)

    # print(peak_intervals)
    print(str(sum_fluorescence_per_ind(d[ind], peak_intervals)))

    plot.show()
    # plot.savefig ("area"+ind+".png")

# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------
# Script4
# Signature finding (application of Kruskal-Wallis statistical test to find the Molecular Weight where both groups are more distinct)

f = open("merged.csv", "r")
ix = 0
data = {}
for line in f.readlines():
    if (ix != 0):
        # print(ix)
        info = line.split(',')
        ind = info[0]
        gender = info[1]
        values = info[2:401]
        # values = [w.replace(',', '.') for w in mw]
        values = list(map(float, values))
        print(values)
        soma_total = np.sum(values)
        values = np.true_divide(values, 1000 / soma_total)
        data[ind] = [gender, [values]]
    ix += 1

group1 = {}
group2 = {}
ind = 0

for key, value in data.items():
    gender = value[0]
    if (gender == 'Female'):
        group1[key] = value[1][0]
    elif (gender == 'Male'):
        if ind > 470:
            continue
        else:
            group2[key] = value[1][0]
        ind += 1

plot.figure(1)
plot.subplot(2, 1, 1)
for key, value in group1.items():
    plot.plot(value, 'b')
for key, value in group2.items():
    plot.plot(value, 'r')


# kruskall wallis
def getGroupMWs(group, mw_ix):
    values = []
    for key, value in group.items():
        values.append(value[mw_ix])
    return values


kruskal_values = []
p_value = []

for v in range(0, 399):
    g1 = getGroupMWs(group1, v)
    g2 = getGroupMWs(group2, v)
    stats.kruskal(g1, g2)

    kruskal, pvalue = stats.kruskal(g1, g2)
    kruskal_values.append(kruskal)
    p_value.append(pvalue)

kruskal_values = np.array(kruskal_values)
plot.plot(kruskal_values, 'g')
plot.plot(p_value, 'y')

plot.subplot(2, 1, 2)
plot.plot(kruskal_values, 'g')
plot.plot(p_value, 'y')

peaks, _ = find_peaks(kruskal_values, height=0, prominence=10)
plot.plot(peaks, kruskal_values[peaks], "x")
plot.show()

plot.figure(2)
# g1 = getGroupMWs(group1, peaks[0])
# g2 = getGroupMWs(group2, peaks[0])

axisg1 = np.zeros(len(g1))
axisg2 = np.ones(len(g2))

plot.plot(axisg1, g1, 'ob')
plot.plot(axisg2, g2, 'or')
plot.show()
'''
import matplotlib.pyplot as plt
import numpy as np
import math

def sum_fluorescence_per_ind(profile, peak_intervals):

    sum1 = 0
    sum2 = 0
    sum3 = 0
    sum4 = 0
    sum5 = 0
    sum6 = 0

    to_consider1 = list(range(int(peak_intervals[0][0]),int(peak_intervals[0][1])))
    to_consider2 = list(range(int(peak_intervals[1][0]),int(peak_intervals[1][1])))
    to_consider3 = list(range(int(peak_intervals[2][0]),int(peak_intervals[2][1])))
    to_consider4 = list(range(int(peak_intervals[3][0]),int(peak_intervals[3][1])))
    to_consider5 = list(range(int(peak_intervals[4][0]),int(peak_intervals[4][1])))
    to_consider6 = list(range(int(peak_intervals[5][0]),int(peak_intervals[5][1])))

    for i in range(0,len(profile)):
        if i in to_consider1:
            sum1 += profile[i]

    for i in range(0,len(profile)):
        if i in to_consider2:
            sum2 += profile[i]

    for i in range(0,len(profile)):
        if i in to_consider3:
            sum3 += profile[i]

    for i in range(0,len(profile)):
        if i in to_consider4:
            sum4 += profile[i]

    for i in range(0,len(profile)):
        if i in to_consider5:
            sum5 += profile[i]

    for i in range(0,len(profile)):
        if i in to_consider6:
            sum6 += profile[i]

    return [sum1, sum2, sum3, sum4, sum5, sum6]

def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]


def find_half_peak_from_top(peak,profile,index_top):
    margin_bottom = 0;
    margin_top = len(profile)

    if (peak == '1'):
        margin_bottom = 2
        margin_top = peaks[str(int(peak)+1)][0]
    if (peak == '2'):
        margin_bottom = peaks[str(int(peak)-1)][1]
        margin_top = peaks[str(+int(peak)+1)][0]
    if (peak == '3'):
        margin_bottom = peaks[str(int(peak)-1)][1]
        margin_top = peaks[str(+int(peak)+1)][0]
    if (peak == '4'):
        margin_bottom = peaks[str(int(peak)-1)][1]
        margin_top = peaks[str(+int(peak)+1)][0]
    if (peak == '5'):
        margin_bottom = peaks[str(int(peak)-1)][1]
        margin_top = peaks[str(+int(peak)+1)][0]
    if (peak == '6'):
        margin_bottom = peaks[str(int(peak)-1)][1]
        #margin_top = peaks[str(+int(peak)+1)][0]
        margin_top = 340

    half_value = profile[index_top]/2

    left_half = profile[:index_top]
    right_half = profile[index_top:]

    closest_right = 0
    distance_right = 9999999999

    #go_right
    for v in range(index_top,len(profile)):
        if profile[v] - half_value < distance_right:
            distance_right = profile[v] - half_value
            closest_right = v
        if (profile[v] < half_value) or (v > margin_top):
            break;

    closest_left = 0
    distance_left = 9999999999

    # go_left
    for v in range(index_top, 0, -1):
        if profile[v] - half_value < distance_left:
            distance_left = profile[v] - half_value
            closest_left = v
        if (profile[v] < half_value) or (v < margin_bottom):
            break;

    return [closest_left,closest_right]

f = open('data.csv')
d = {}

lines = f.readlines()

for line in lines:
    line = line.strip()
    # line = line.replace(',','.')
    aux = line.split(',')
    k = aux[0]
    val = aux[1:]
    values = list(map(float,val))
    d[k] = values


peaks = {'1': [48, 57], '2': [112, 129], '3': [136, 166], '4': [178, 201], '5': [241, 251], '6': [263, 278]}


for ind in d.keys():
    
    #print(ind)
    plt.figure()

    plt.plot(d[ind])

    peak_intervals = []

    for peak in peaks.keys():
        # Janela inicial
        [init,fim] = peaks[peak]
        window = d[ind][init-1:fim]
        xs = range(init-1,fim)
        index_max_window = window.index(max(window))
        plt.plot(xs , window, 'r')

        # Maximums
        index_max = init-1 + index_max_window
        plt.plot(index_max, window[index_max_window],'go')

        boundaries = find_half_peak_from_top(peak,d[ind],index_max)

        plt.plot(boundaries[0], d[ind][boundaries[0]],'b*')
        plt.plot(boundaries[1], d[ind][boundaries[1]],'y*')

        peak_intervals.append(boundaries)

    #print(peak_intervals)
    print(str(sum_fluorescence_per_ind(d[ind],peak_intervals)))

    plt.show()
    #plt.savefig ("area"+ind+".png")



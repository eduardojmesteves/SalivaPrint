import matplotlib.pyplot as plt
import math
from scipy.fftpack import fft
from scipy.fftpack import ifft
from scipy.signal import find_peaks
import numpy as np
from scipy.optimize import curve_fit
from scipy import stats

f = open("genero.csv", "r")
ix = 0
data = {}
for line in f.readlines():
    if (ix != 0):
        #print(ix)
        info = line.split(';')
        ind = info[0]
        gender = info[1]
        mw = info[2:400]
        values = [w.replace(',', '.') for w in mw]
        values = list(map(float, values))
        soma_total = np.sum(values)
        values = np.true_divide(values, 1000/soma_total)
        data[ind] = [gender, [values]]
    ix+=1

group1 = {}
group2 = {}
ind = 0;

for key, value in data.items():
    gender = value[0]
    if (gender == 'Feminino'):
        group1[key] = value[1][0]
    elif (gender == 'Masculino'):
        if ind>470:
            continue
        else:
            group2[key] = value[1][0]
        ind+=1
print(len(group1.items()))
print(len(group2.items()))

plt.figure(1)
plt.subplot(2, 1, 1)
for key, value in group1.items():
    plt.plot(value,'b')
for key, value in group2.items():
    plt.plot(value,'r')

#kruskall wallis
def getGroupMWs(group, mw_ix):
    values = []
    for key, value in group.items():
        values.append(value[mw_ix])
    return values;

kruskal_values = []
p_value = []

for v in range(0,398):
    
    g1 = getGroupMWs(group1,v)
    g2 = getGroupMWs(group2,v)
    stats.kruskal(g1,g2)
    
    kruskal, pvalue = stats.kruskal(g1,g2)
    kruskal_values.append(kruskal)
    p_value.append(pvalue)

kruskal_values = np.array(kruskal_values)
plt.plot(kruskal_values,'g')
plt.plot(p_value,'y')

plt.subplot(2, 1, 2)
plt.plot(kruskal_values,'g')
plt.plot(p_value,'y')

peaks, _ = find_peaks(kruskal_values, height=0, prominence=10)
print(kruskal_values[peaks])
plt.plot(peaks, kruskal_values[peaks], "x")
plt.show()


plt.figure(2)
g1 = getGroupMWs(group1,peaks[0])
g2 = getGroupMWs(group2,peaks[0])

print(g1)
print(g2)
axisg1 = np.zeros(len(g1))
axisg2 = np.ones(len(g2))

plt.plot(axisg1,g1,'ob')
plt.plot(axisg2,g2,'or')
plt.show()



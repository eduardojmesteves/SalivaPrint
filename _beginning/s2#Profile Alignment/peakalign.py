### Peal Alignment ###
import pandas as pd
import numpy as np
from scipy import signal
from matplotlib import pylab as plt



###PROFILE IMPORT
data = pd.read_csv('/home/eduardo/Documents/salivaprint_py/s2 #Profile alignement/data1.csv')

mw = data.mw

data = data.drop(columns=['ix', 'mw'])

all_profiles = list(data.columns.values)

master = all_profiles[0]

master_profile = data[master]

all_profiles_minus_actual = all_profiles.copy()
all_profiles_minus_actual.remove(master)

targets = data[all_profiles_minus_actual]

###Peak alignment

f1 = plt.figure(0)
for target in targets:
    plt.plot(mw, targets[target])

f2 = plt.figure(1)
for target in targets:
    dx = np.mean(np.diff(mw.values))
    shift = (np.argmax(signal.correlate(master_profile, targets[target])) - len(targets[target])) * dx
    plt.plot(mw + shift, targets[target])
plt.show()



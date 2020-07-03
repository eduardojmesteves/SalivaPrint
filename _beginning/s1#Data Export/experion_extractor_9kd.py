import re
import math
import collections
import os
import pandas as pd
import csv
# def mean(data):
#     """Return the sample arithmetic mean of data."""
#     n = len(data)
#     if n < 1:
#         raise ValueError('mean requires at least one data point')
#     return sum(data)/n # in Python 2 use sum(data)/float(n)
# def _ss(data):
#     """Return sum of square deviations of sequence data."""
#     c = mean(data)
#     ss = sum((x-c)**2 for x in data)
#     return ss
# def stddev(data, ddof=0):
#     """Calculates the population standard deviation
#     by default; specify ddof=1 to compute the sample
#     standard deviation."""
#     n = len(data)
#     if n < 2:
#         raise ValueError('variance requires at least two data points')
#     ss = _ss(data)
#     pvar = ss/(n-ddof)
#     return pvar**0.5

import sys

import glob, os
for file in glob.glob("*.xml"):
    filename = file
    
profiles = {}
ind_names = []

i = 0;
import codecs
f1 = codecs.open(filename.split('.')[0]+"_RunSummary.txt", encoding='utf-16-le')
for line in f1:  # note, the readlines is not really needed
    if i != 0:
        ind_names.append(line.split('\t')[1]);  # the comma strips the trailing newline in case that's bothering you
    i+=1

file = open(filename)
content = file.read()
file.close()

content = content[content.find('<RawDataPoints>'):]

# Gets the points.
# aux = re.search(r"<Signal>(.*?)</Signal>", content)
aux = {}
ix = 0
key = 0
matches = re.finditer(r"<Signal>(.*?)</Signal>", content)

for matchNum, match in enumerate(matches):
    matchNum = matchNum + 1
    for groupNum in range(0, len(match.groups())):
        groupNum = groupNum + 1
        group = match.group(groupNum)
        if (ix % 1380==0):
            key+=1;
        if key in aux.keys():
            aux[key].append(float(group))
        else:
            aux[key] = [float(group)]
        ix+=1;

meta_file = open(filename.split('.')[0]+".csv",'w')
medias = []
inds = aux.keys()

for i in range(0,1380):
    vals = []
    for ind in inds:
        vals.append(aux[ind][i])
    medias.append(sum(vals)/float(len(vals)))

# # -- - - stdevs
# dvs = []
# for i in range(0,1380):
#     vals = []
#     for ind in inds:
#         vals.append(aux[ind][i])
#     dvs.append(stddev(vals))

ix = 0;
#------ write data
for ind in aux.keys():
    meta_file.write(ind_names[ix]+',')
    for i in range(465, 869):
        if i < 868:
            meta_file.write(str(aux[ind][i]) + ',')
        else:
            meta_file.write(str(aux[ind][i]) + '\n')
    ix+=1
def truncate(n, decimals=3)

# meta_file.write('Average,')
# for i in range(465,869):
#     if i < 868:
#         meta_file.write(str(medias[i])+',')
#     else:
#         meta_file.write(str(medias[i])+'\n')

# meta_file.write('Std_Dev,')
# for i in range(465,869):
#     if i < 868:
#         meta_file.write(str(dvs[i])+',')
#     else:
#         meta_file.write(str(dvs[i])+'\n')

# Indexes
# meta_file.write('Index,')
# for i in range(466,869): # 471 nao é erro, tem de ser mais 1 pk em python começa em 0 e é preciso para os pesos ficarem certos
#     if i < 868:
#         meta_file.write(str(i)+',')
#     else:
#         meta_file.write(str(i)+'\n')

# MOL WEIGHT
meta_file.write('mw,')
for i in range(466,869): # 471 nao é erro, tem de ser mais 1 do que 470
    if i < 868:
        meta_file.write(str(0.00042006*(i*i)-0.28216*i+49.225)+',')
    else:
        meta_file.write(str(0.00042006*(i*i)-0.28216*i+49.225)+'\n')
    def truncate(n, decimals=3)

meta_file.flush()
meta_file.close()

# Read in the file
with open(filename.split('.')[0]+".csv", 'r') as file :
  filedata = file.read()

# Replace the target string
# filedata = filedata.replace('.', ',')

# Write the file out again
with open(filename.split('.')[0]+".csv", 'w') as file:
  file.write(filedata)




#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import re

colourDict = {}
df = pd.read_table("/t1-data1/WTSA_Dev/jkerry/CaptureC/CRISPR-OT/AllOligos_Info.txt",sep="\t")

## Define colour codes
Counter = 0.3
for i in df['Site'].unique():
    colourDict[i] = Counter
    Counter=Counter+0.3

## Add colour list codes to dataframe    
CodeList = []
for i in df['Site']:
    CodeList.append(colourDict[i])
df['Colour code'] = CodeList

## Add sizes/(number from cut-site) to dataframe
NumberList = []
for i in df['Location']:
    NumberList.append(16-int(re.search(r'\d+',i).group(0)))
    #print(numbers)
df['Point Size'] = NumberList

fig, ax = plt.subplots()
plt.scatter(df['Density score'],df['Repeat length'],c=cm.rainbow(df['Colour code']),s=[x*3 for x in df['Point Size']])
ax.set_xlim(0,max(df['Density score'])+10)
ax.set_xlabel('STAR Density Score')
ax.set_ylabel('Repeat Length')
plt.savefig('Distribution.png')
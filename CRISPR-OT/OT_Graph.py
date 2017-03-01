#!/usr/bin/env python

from __future__ import division
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import re,getopt,sys

def usage():
    print("usage: OT_Graph.py -o <oligo info file> -s <site name>")


oligo_file = ""
key = ""

try:
    opts, args = getopt.getopt(sys.argv[1:], 'o:s:h',)
except getopt.GetoptError:
    usage()
    sys.exit(2)
    
if not opts:
    usage()
    sys.exit(2)
else:
    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit(2)
        elif opt == '-o':
            oligo_file = arg
        elif opt == '-s':
            key = arg
        else:
            usage()
            sys.exit(2)

df = pd.read_table(oligo_file,sep="\t")
group_df = df.groupby('Site')
    
SiteList = []
for i in group_df.get_group(key).index.values:
    if df.iloc[i]['Location'][:2]=="Up":
        SiteList.append(-1*int(df.iloc[i]['Distance from site (bp)']))
    else:
        SiteList.append(int(df.iloc[i]['Distance from site (bp)']))


fig, ax = plt.subplots()
ColourList = np.where(group_df.get_group(key)['Repeat length']<(len(group_df.get_group(key).reset_index()['Sequence'][0])/4),'b','r')
plt.scatter(SiteList,group_df.get_group(key)['Density score'],c=ColourList,s=[x+10 for x in group_df.get_group(key)['Repeat length']])

ax.set_xlabel('Distance from off-target site (bp)')
ax.set_ylabel('STAR Density Score',rotation=0)

ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.tick_params(axis='both', labelsize=8)
plt.title(key,y=1.07, fontsize=14, fontweight='bold')
ax.yaxis.set_label_coords(0.5,1.02)
ax.set_xlim(min(SiteList)-10,max(SiteList)+10)
ax.set_ylim(0,max(group_df.get_group(key)['Density score'])+10)
plt.plot([100,100],[0,max(group_df.get_group(key)['Density score'])+10], 'k--', lw=1)
plt.plot([-100,-100],[0,max(group_df.get_group(key)['Density score'])+10], 'k--', lw=1)
plt.plot([min(SiteList)-10,max(SiteList)+10],[50,50], 'r--', lw=1)
#red_patch = mpatches.Patch(color='red',label='Repeat length>=cutoff')
#blue_patch = mpatches.Patch(color='blue',label='Repeat length<cutoff')
#first_legend = plt.legend(handles=[red_patch,blue_patch],loc=2,bbox_to_anchor=(1.01,1),borderaxespad=0.)
#plt.gca().add_artist(first_legend)
for i, txt in enumerate(group_df.get_group(key)['Location']):
    ax.annotate(txt,(SiteList[i],group_df.get_group(key).reset_index()['Density score'][i]), fontsize=8)
plt.subplots_adjust(bottom=0.1) 
#plt.show()
plt.savefig("ExampleGraph.png")

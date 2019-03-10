from __future__ import division
import numpy as np
from matplotlib import pyplot as plt
import seaborn.apionly as sns
import matplotlib as mpl
import os

os.environ["PROJ_LIB"] = 'C:\\ProgramData\\Anaconda3\\envs\\test\\Library\\share'

from mpl_toolkits.basemap import Basemap

# load data
data = np.loadtxt('Profitability_AfricaSites_fgls_lev_bestmod.csv',delimiter=',',skiprows=1,usecols=[0,61,62,63,64,73,74,77,80]) 
#columns correspond to [cellID, class_nono_70, class_yesyes_70, class_type1_70, class_type2_70, lat, lon, prneeded_robust, prneeded_naive]
#prneeded code:
# 1 = Profitable if F/M <= 0 (i.e. never)
# 2 = Profitable if 0 < F/M < 5
# 3 = Profitable if 5 < F/M < 10
# 4 = Profitable if F/M >= 10

# simulated VCR distributions at different sites under naive and robust methods
naiveCDFs = np.loadtxt('AfricaSites_Sim_Shortlist_naive_fgls_lev_bestmod.csv',delimiter=',',skiprows=1)
robustCDFs = np.loadtxt('AfricaSites_Sim_Shortlist_robust_fgls_lev_bestmod.csv',delimiter=',',skiprows=1)

# locations of sites with simulated VCRs
lat = np.arange(-39.875,40.125,0.25)
lon = np.arange(-19.875,55.125,0.25)
selectSites = np.array([598.0,3200.0,3593.0,8271.0,9166.0,9331.0,10168.0,10760.0,10976.0,13964.0])
selectData = np.zeros([len(selectSites),np.shape(data)[1]])
for i in range(len(selectSites)):
    row = [k for k in range(len(data[:,0])) if data[k,0] == selectSites[i]]
    selectData[i,:] = data[row,:]
    
# color of line for each site and where to label it on the map
siteColors = ['#e31a1c','#fb9a99','#ff7f00','#fdbf6f','#6a3d9a','#cab2d6','#1f78b4','#a6cee3','#33a02c','#b2df8a']
label_ys = np.array(selectData[:,5] + 1.5)
label_ys[6] = selectData[6,5] - 1.75
label_xs = np.array(selectData[:,6])
label_xs[5] = selectData[5,6] - 1.0
label_xs[4] = selectData[4,6] + 1.0

def makeFigure2(simVCR, selectSites, selectData, label_xs, label_ys, colors, estimator):
    
    sns.set()
    fig = plt.figure()
    ax = fig.add_subplot(1,2,1)
    
    # plot basemap and countries
    m = Basemap(llcrnrlat=-40,urcrnrlat=40,llcrnrlon=-20,urcrnrlon=55,resolution='l')
    m.arcgisimage(service='World_Shaded_Relief')
    m.drawcountries(color='0.6', linewidth=0.5)
    m.drawparallels(np.arange(-40,41,20), labels=[1,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5', fontsize=18)
    m.drawmeridians(np.arange(-20,56,15), labels=[1,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5', fontsize=18)
    
    # plot select sites
    x, y = m(selectData[:,6],selectData[:,5])
    m.scatter(x, y, facecolor=colors, edgecolor='k', s=50, marker='s', alpha=1.0)
    
    # label select sites
    for i in range(len(colors)):
        plt.text(label_xs[i], label_ys[i], str(int(selectSites[i])), fontsize=22, \
            ha='center', va='center', color='k')
    
    # plot CDFs from each site
    ax = fig.add_subplot(1,2,2)
    n = np.shape(simVCR)[0]
    M = np.array(range(1,n+1))
    P = M/(n+1)
    for i in range(len(colors)):
        ax.step(np.sort(simVCR[:,i]), P, c=colors[i], label=str(int(selectSites[i])), linewidth=2)
        
    ax.set_xlim([-0.5,3.25])
    ax.set_ylim([0,1])
    ax.plot([0,0],[0,1],linewidth=2, linestyle='--', c='k')
    
    ax.tick_params(axis='both',labelsize=22)
    ax.set_xlabel('Yield Difference (t/ha)',fontsize=26)
    ax.set_ylabel('Cumulative Distribution',fontsize=26)
    handles, labels = plt.gca().get_legend_handles_labels()
    ax.legend(handles, labels, loc='lower right', ncol=1, fontsize=22)
    fig.set_size_inches([24, 12.05])
    fig.savefig('Figures/Figure2_' + estimator + '.pdf')
    fig.savefig('Figures/Figure2_' + estimator + '.png')
    fig.clf()
    
    return None
    
def makeFigure3(data):
    
    counts = [np.sum(data[:,1]), np.sum(data[:,2]), np.sum(data[:,3]), np.sum(data[:,4])]
    labels = ['Never Profitable','Always Profitable','Underutilize\nw/ naive','Overutilize\nw/ naive']
    explode = (0, 0, 0.1, 0.1)
    colors = ['#e41a1c','#4daf4a','#377eb8','#ff7f00']
    mpl.rcParams['font.size'] = 16
    
    sns.set_style("dark")
    fig, ax = plt.subplots()
    patches, texts, autotexts = ax.pie(counts, explode=explode, labels=labels, \
        autopct='%1.1f%%', shadow=False, startangle=90, colors=colors)
    for i in range(len(colors)):
        texts[i].set_fontsize(16)

    ax.axis('equal')
    fig.set_size_inches([11.7375, 8.4375])
    fig.savefig('Figures/Figure3.pdf')
    fig.savefig('Figures/Figure3.png')
    fig.clf()
    
    return None
    
def makeFigure4(data, lon, lat):

    C = np.array([[228, 26, 28],[77, 175, 74],[55, 126, 184],[255, 127, 0]])/255.0
    labels = ['Never\nProfitable','Always\nProfitable','Underutilize\nw/ naive','Overutilize\nw/ naive']
    
    sns.set()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    # plot basemap and countries
    m = Basemap(llcrnrlat=-40,urcrnrlat=40,llcrnrlon=-20,urcrnrlon=55,resolution='l')
    m.arcgisimage(service='World_Shaded_Relief')
    m.drawcountries(color='0.6', linewidth=0.5)
    m.drawparallels(np.arange(-40,41,20), labels=[1,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5', fontsize=18)
    m.drawmeridians(np.arange(-20,56,15), labels=[1,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5', fontsize=18)
    
    x, y = np.meshgrid(lon, lat)
    array = np.empty((80*4,75*4))
    array[:] = np.NAN;
    
    for i in range(np.shape(data)[0]):
        row = int((data[i,5]-np.min(lat))/0.25)
        col = int((data[i,6]-np.min(lon))/0.25)
        if data[i,1] == 1:
            array[row,col] = 0
        elif data[i,2] == 1:
            array[row,col] = 1
        elif data[i,3] == 1:
            array[row,col] = 2
        else:
            array[row,col] = 3
    
    array_mask = np.ma.masked_invalid(array)
    x, y = np.meshgrid(lon, lat)
    x, y = m(x,y)
    cmap = mpl.colors.ListedColormap(C)
    m.pcolormesh(x, y, array_mask, cmap=cmap, rasterized=False, edgecolor='0.6', linewidth=0)
        
    formatter = plt.FuncFormatter(lambda val, loc: labels[val])
        
    cbar = plt.colorbar(ticks=[0,1,2,3], format=formatter)
    plt.clim(-0.5, 3.5)
    cbar.solids.set_edgecolor("face")
    cbar.ax.tick_params(labelsize=22)
    fig.set_size_inches([24, 12.05])
    fig.savefig('Figures/Figure4.pdf')
    fig.savefig('Figures/Figure4.png')
    fig.clf()
    
    return None
    
def makeFigure5(data, lon, lat, estimator):
    
    C = np.array([[178, 24, 43],[214, 96, 77],[135, 135, 135],[26, 26, 26]])/255.0
    labels = ['Never\nProfitable','Profitable if\nF:M $\leq$ 5','Profitable if\nF:M $\leq$ 10','Profitable even if\nF:M > 10']
    
    sns.set()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    # plot basemap and countries
    m = Basemap(llcrnrlat=-40,urcrnrlat=40,llcrnrlon=-20,urcrnrlon=55,resolution='l')
    m.arcgisimage(service='World_Shaded_Relief')
    m.drawcountries(color='0.6', linewidth=0.5)
    m.drawparallels(np.arange(-40,41,20), labels=[1,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5', fontsize=18)
    m.drawmeridians(np.arange(-20,56,15), labels=[1,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5', fontsize=18)
    
    array = np.empty((80*4,75*4))
    array[:] = np.NAN;
    
    for i in range(np.shape(data)[0]):
        row = int((data[i,5]-np.min(lat))/0.25)
        col = int((data[i,6]-np.min(lon))/0.25)
        if estimator == 'naive':
            array[row,col] = int(data[i,8]-1)
        else:
            array[row,col] = int(data[i,7]-1)
        
    array_mask = np.ma.masked_invalid(array)
    x, y = np.meshgrid(lon, lat)
    x, y = m(x,y)
    cmap = mpl.colors.ListedColormap(C)
    m.pcolormesh(x, y, array_mask, cmap=cmap, rasterized=False, edgecolor='0.6', linewidth=0)
        
    formatter = plt.FuncFormatter(lambda val, loc: labels[val])
        
    cbar = plt.colorbar(ticks=[0,1,2,3], format=formatter)
    plt.clim(-0.5, 3.5)
    cbar.solids.set_edgecolor("face")
    cbar.ax.tick_params(labelsize=22)
    fig.set_size_inches([24, 12.05])
    fig.savefig('Figures/Figure5_' + estimator + '.pdf')
    fig.savefig('Figures/Figure5_' + estimator + '.png')
    fig.clf()
    
    return None
    
makeFigure2(naiveCDFs, selectSites, selectData, label_xs, label_ys, siteColors, 'naive')
makeFigure2(robustCDFs, selectSites, selectData, label_xs, label_ys, siteColors, 'robust')
makeFigure3(data)
makeFigure4(data, lon, lat)
makeFigure5(data, lon, lat, 'naive')
makeFigure5(data, lon, lat, 'robust')

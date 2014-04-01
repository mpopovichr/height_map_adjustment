__author__ = 'mpopovic'


import sqlite3 as lite
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pandas.io.sql as psql
import pandas.io.parsers as pp
import matplotlib.image as mpimg
import tifffile as tiff
import scipy.signal as signal
from matplotlib.colors import LogNorm


inPath = '/Users/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/height_maps/'
inFile = 'HeightMa.png'
inFile = 'HM_Stitch_Time_196.tif'

a = tiff.imread(inPath+inFile)


tiff.imshow(a)
plt.show()

gx, gy = np.gradient(a)
gg = np.sqrt(gx**2+gy**2)
gx_cut = gx
gx.max()
bin = np.arange(18)
hist, bins = np.histogram(gx, bin)

tiff.imshow(smooth_30_gx[500:2000,500:1500], cmap='hot')
tiff.imshow(np.exp(smooth_20_gx), cmap='gist_rainbow')
tiff.imshow(np.exp(smooth_20_gy), cmap='gist_rainbow')
tiff.imshow(a, cmap='gist_rainbow')
plt.show()

c = np.exp(smooth_20_gx[500:1200,500:1200])
tiff.imshow(c, cmap='gist_rainbow')
plt.show()

Nsmooth = 30
smooth = np.ones((Nsmooth,Nsmooth))/Nsmooth**2


smooth_30_gx = signal.convolve(gx, smooth)
smooth_30_gy = signal.convolve(gy, smooth)
smooth_20_gx = signal.convolve(gx, smooth)
smooth_20_gy = signal.convolve(gy, smooth)

dbPath = '/Users/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/DB/'
dbName = 'WT_25deg_111102'


con = lite.connect(dbPath+dbName+'/'+dbName+'.sqlite')
df = psql.frame_query('SELECT * FROM cells WHERE (frame==196 AND cell_id>10000)', con)

cell_gx, cell_gy = [],[]
for i in range(len(df['cell_id'])):
    side = np.int(np.sqrt(df['area'][i]))
    pos_x, pos_y = np.int(df['center_x'][i]), np.int(df['center_y'][i])
    cell_gx.append(np.mean(gx[pos_y-side/2:pos_y+side/2,pos_x-side/2:pos_x+side/2].flatten()))
    cell_gy.append(np.mean(gy[pos_y-side/2:pos_y+side/2,pos_x-side/2:pos_x+side/2].flatten()))

cell_gx, cell_gy = np.array(cell_gx)/0.208, np.array(cell_gy)/0.208
df['cell_gx'] = cell_gx
imPath = '/Users/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/raphael_data/111102_SEG/segmentation_data/111102_segmentation_v2/'
imName = 'Optimized_projection'
wing_im = plt.imread(imPath+imName+'_196.png')
area = df['area']
corr_factor = np.sqrt((1.+cell_gy**2.)/(1+cell_gx**2.))
pt_size = np.abs([x*40./2000. for x in area])
cm = plt.cm.get_cmap('gist_rainbow')
fig = plt.figure()
p = plt.scatter(df['center_x'],df['center_y'],c=corr_factor,s=pt_size, cmap=cm,vmin=1, vmax=1.05, lw=0)
plt.imshow(wing_im)
ax = plt.gca()
cbar = fig.colorbar(p, ax = ax, shrink = .905, orientation = 'horizontal', pad = -0.03)
plt.show()
np.ones((10,10)).flatten()
len(cell_gx)
gx.max()
df['gx'] = np.mean(a[np.floor(df['center_x']),np.int(df['center_y'])])
np.round(df['center_y'])
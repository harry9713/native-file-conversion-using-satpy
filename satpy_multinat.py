#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  1 00:04:57 2021

@author: edith005

This python script loads the Spinning Enhanced Visible Infra-Red Imager(SEVIRI) data from 
Meteosat Second Generation(MSG-1) satellite in native(.nat) format using satpy and 
pyresample. Satpy loads the data as scene object using its pre-built caliberation functions
for the specific reader of SEVIRI level 1.5 image data.

"""


#%%
import glob
import datetime
import cartopy.crs as ccrs
import cftime
import numpy as np
import xarray as xr
from netCDF4 import Dataset
from pyresample.geometry import AreaDefinition, SwathDefinition
from pyresample.kd_tree import resample_gauss, resample_nearest
from satpy import MultiScene, Scene
from satpy.multiscene import timeseries
from satpy.resample import get_area_def

#%%Data Loading
pdir ='/DataShared/DataFiles/TB data/Met7_SEVIRI_20180812_20180819/native files/MSG1-SEVI-MSG15-0100-NA-20180812*.nat'
natfiles = glob.glob(pdir)


#create Multiscene Object to load multiple data files simultaneously.
#here, the reader is specified as "seviri_l1b_native" to handle native files
#and loading the data for a specific channel, here 10.8um IR channel is used.
#then saves the data to netCDF file for re-sampling to a uniform grid.


msc = MultiScene.from_files(natfiles,reader='seviri_l1b_native')
msc.load(['IR_108'])
msc.save_datasets(base_dir='/DataShared/DataFiles/TB data/Met7_SEVIRI_20180812_20180819/converted nc satpy',writer='cf')

#%%Data processing
pdir_1 = '/DataShared/DataFiles/TB data/Met7_SEVIRI_20180812_20180819/converted nc satpy/IR_108_20180812*.nc'
ncfiles = glob.glob(pdir_1)
ncfiles.sort()
te=[]
for i in range(len(ncfiles)):
#---opens the data files one by one using xarray---------------------------
    ds1 = xr.open_dataset(ncfiles[i])
#---Finding the lat-lon points for indian region---------------------------
#---here latitude and longitude are matrices of----------------------------
#---size 3712x3712 so selecting a central row/column ----------------------
#---to ge the lat-lon target-grid-------  --------------  -----------------
    ln = ds1.longitude.data[2500,:]
    lo = ln[(ln<=100)&(ln>=60)]
    lt = ds1.latitude.data[:,2500]
    la = lt[(lt>=0)&(lt<=30)]
#---creating a target uniform meshgrid------------------------------------
    X,Y = np.meshgrid(lo[::-1],la[::-1])
    
    
    
#---creating swath defenition from the original grid----------------------
    origin_grid = SwathDefinition(lons=ds1.longitude, lats=ds1.latitude)
#---creating swath defenition for the uniform target grid-----------------    
    target_grid = SwathDefinition(lons=X, lats=Y)
#---Extracting the entire data and resampling using nearest method in pyresample    
    data = ds1.IR_108.to_masked_array()
    data_near = resample_nearest(origin_grid, data, target_grid, radius_of_influence=20000)
    # data_gaus = resample_gauss(origin_grid, data, target_grid, radius_of_influence=20000, sigmas=1000)


#---This step extracts the lat-lon from the new data created--------------    
    lons,lats = target_grid.get_lonlats()
    latt = lats[:,450]
    lonn = lons[450,:]
    # X1,Y1 = np.meshgrid(lonn,latt)
#---creating the data array for store the data at different time steps----    
    if i==0:
        data_new = np.empty(shape=(len(ncfiles),len(latt),len(lonn)))        
    data_new[i,:,:] = data_near


#---Extracting end_time from the data    
    ti = ds1.IR_108.end_time
    tid = datetime.datetime.strptime(ti[:19],'%Y-%m-%d %H:%M:%S')
    print('end_time: ',tid)
#----Converting it to cftime object---------------------------------------    
    t = cftime.date2num(dates=tid, units='seconds since 1970-1-1', calendar='standard')
    print('cftime: ',t)
    te.append(int(t))

#%%Dataset creation.
nwti = datetime.datetime.now()
initime = datetime.datetime.utcfromtimestamp(te[0])
fintime = datetime.datetime.utcfromtimestamp(te[-1])
Fill_Value =-1252.685
root = Dataset("/DataShared/DataFiles/TB data/Met7_SEVIRI_20180812_20180819/nc_processed/"+\
    f"IR_108_NAT_combined_{initime}_{fintime}.nc4","w",format="NETCDF4")
root.title = 'tb_test'
root.history = 'Generated from the MET-8 native files \n' +\
                nwti.strftime('created on %d-%m-%Y %H:%M\n')+\
                    'generating code: read_met'

                                     
# dimensions
lat = root.createDimension("lat",len(latt))
lon = root.createDimension("lon",len(lonn))
#lev = root.createDimension(dd['levt'][i],len(dat.attr['levs']))
time = root.createDimension("time",len(te))
lats = root.createVariable("lat","f4",("lat"),fill_value=False)
lons = root.createVariable("lon","f4",("lon"),fill_value=False)
times = root.createVariable("time","int32",("time"),fill_value=False)
times.units = "seconds since 1970-1-1"
lats.long_name = 'latitude'
lats.units = "degrees_north"
lons.long_name = 'latitude'
lons.units = "degrees_east"
lats[:]=latt
lons[:]=lonn
times[:]=te


Tb = root.createVariable('Tb',"f4",("time","lat","lon"),fill_value=Fill_Value,zlib=True)
Tb.long_name = 'Brigtness temperature'
Tb.units = 'K'
Tb[:,:,:] = data_new
root.close()
print('Dataset saved.')
#%%test plot


import plotter_class
import xarray as xr
ds_proc = xr.open_dataset('/DataShared/DataFiles/TB data/Met7_SEVIRI_20180812_20180819/nc_processed/'+\
    'IR_108_NAT_combined_2018-08-12 00:30:11_2018-08-13 00:00:10.nc4')
dat = plotter_class.sat()
data_ind_IR = 90
data_ind_ncep = 46
dat._plotter_(ds_proc.Tb[data_ind_IR],t=str(ds_proc.time[data_ind_IR].data)[:19],bbox=[60,98,0,29])
ds_ncep = xr.open_dataset('/DataShared/DataFiles/TB data/NCEP_merg_IR/20180812-19/oneweek.nc')
dat._plotter_(ds_ncep.Tb[data_ind_ncep],t=str(ds_ncep.time[data_ind_ncep].data)[:19],bbox=[60,98,0,29])

# %%

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 10:38:03 2021

@author: ajil
"""


import numpy as np

import matplotlib.pyplot as plt
#path='/home/ajil/Desktop/GRIDSAT-B1.1996.01.01.03.v02r01.nc'
import matplotlib.colors as colors
import matplotlib.colorbar as clb
import cartopy.crs as ccrs
from cartopy import feature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

class sat(object):
    ''' Class for plotting different satellite datas '''
    
    listcolors=['#161d58','#253494','#2850a6','#2c7fb8','#379abe','#41b6c4',
            '#71c8bc','#a1dab4','#d0ecc0','#ffffcc','#fef0d9','#fedeb1',
            '#fdcc8a','#fdac72','#fc8d59','#ef6b41','#e34a33','#cb251a',
            '#b30000','#7f0000']
    mymap=colors.ListedColormap(listcolors)
    bbox=[60,100,0,30]

    def __init__(self):
        pass
        #self.file=path
        #print(self.convert_GridSat['IR0'])
   
#%%
   
    def _plotter_(self,Tb_arr,t='',sat_name=None,bbox=bbox,cmap=mymap,
                  sc=True,vmin=200,vmax=300):
        
            """ This class method will plot the input data 
            input       :    Tb_arr -- array object
            sat_name    : name of the satellite if the data is from gridsat
            t           : title string
            bbox        : to set the extent
            sc          : show colorbar (boolean)"""
            
            
            if sat_name==None:
                coords = [i for i in Tb_arr.coords]
                if 'longitude' and 'latitude' in coords:
                    lat = Tb_arr.latitude
                    lonn = Tb_arr.longitude
                elif 'lon' and 'lat' in coords:
                    lat = Tb_arr.lat
                    lonn = Tb_arr.lon
                else:
                    raise KeyError("Could not find any coordinate named lat/latitude or lon/longitude")
                #X,Y= np.meshgrid(lonn,lat)
            
            elif sat_name=='grds':
                lat=np.linspace(-70.0,69.93001, 2000)
                lonn=np.linspace(-179.93, 179.87  ,5143)
                
            X,Y= np.meshgrid(lonn,lat)   
            fig = plt.figure(dpi=500)
            ax = plt.subplot(projection=ccrs.PlateCarree())
            ax.add_feature(feature.NaturalEarthFeature(category='cultural',
                                                       name='admin_1_states_provinces_lines',
                                                       scale='50m',facecolor='none'))
            ax.coastlines()
            gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=2, color='grey', 
                              alpha=0.3, linestyle='-', draw_labels=True)
            fs=11
            gl.top_labels = False
            gl.right_labels = False
            gl.xformatter = LONGITUDE_FORMATTER
            gl.yformatter = LATITUDE_FORMATTER
            gl.xlabel_style = {'size': fs}
            gl.ylabel_style = {'size': fs}
            plt.title(t,size=15)
            iax = ax.pcolormesh(X,Y,Tb_arr, cmap=cmap,vmin=vmin,vmax=vmax)
            ax.set_extent(bbox,crs=ccrs.PlateCarree())
            if sc:
                cax,kw = clb.make_axes(ax,location='right',pad=0.02,shrink=0.5,
                                       fraction=0.10,aspect=12)
                fig.colorbar(iax, cax=cax)#extend = 'both',drawedges=True)
                #fig.subplots_adjust(wspace=0.2,top=0.8)
    
    
#%% 
          

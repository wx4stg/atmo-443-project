#!/usr/bin/env python3


import xarray as xr
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import numpy as np


targetLat = 42.493
targetLon = -79.272

snow1 = xr.open_dataset("sfav2_CONUS_6h_2022111818.nc").sel(lat=slice(targetLat+1, targetLat-1), lon=slice(targetLon-1, targetLon+1))
snow2 = xr.open_dataset("sfav2_CONUS_6h_2022111900.nc").sel(lat=slice(targetLat+1, targetLat-1), lon=slice(targetLon-1, targetLon+1))
snow3 = xr.open_dataset("sfav2_CONUS_6h_2022111906.nc").sel(lat=slice(targetLat+1, targetLat-1), lon=slice(targetLon-1, targetLon+1))
print(snow1.Data)
summedData = (snow1.Data.data + snow2.Data.data + snow3.Data.data)

vminTarget = np.nanmin(summedData)
vmaxTarget = np.nanmax(summedData)

print(vminTarget)
print(vmaxTarget)

fig = plt.figure()
gs = GridSpec(3, 2, figure=fig, height_ratios=[1, 1, 0.2])
ax1 = fig.add_subplot(gs[0, 0], projection=ccrs.PlateCarree())
pcm = ax1.pcolormesh(snow1.lon, snow1.lat, snow1.Data.data*1000, transform=ccrs.PlateCarree(), vmin=vminTarget, vmax=vmaxTarget, zorder=1)

ax2 = fig.add_subplot(gs[0, 1], projection=ccrs.PlateCarree())
ax2.pcolormesh(snow2.lon, snow2.lat, snow2.Data*1000, transform=ccrs.PlateCarree(), vmin=vminTarget, vmax=vmaxTarget,  zorder=1)

ax3 = fig.add_subplot(gs[1, 0], projection=ccrs.PlateCarree())
ax3.pcolormesh(snow3.lon, snow3.lat, snow3.Data*1000, transform=ccrs.PlateCarree(), vmin=vminTarget, vmax=vmaxTarget, zorder=1)

ax4 = fig.add_subplot(gs[1, 1], projection=ccrs.PlateCarree())
ax4.pcolormesh(snow3.lon, snow3.lat, summedData*1000, transform=ccrs.PlateCarree(), vmin=vminTarget, vmax=vmaxTarget, zorder=1)


for ax in [ax1, ax2, ax3, ax4]:
    ax.plot(targetLon, targetLat, marker='*', color='red', markersize=5, transform=ccrs.PlateCarree())
    ax.set_extent([targetLon-1, targetLon+1, targetLat-1, targetLat+1], crs=ccrs.PlateCarree())
    ax.set_box_aspect(1)
    ax.add_feature(cfeat.COASTLINE, zorder=2)
    ax.add_feature(cfeat.STATES, zorder=2)

cbax = fig.add_subplot(gs[2, :])
cb = plt.colorbar(pcm, cax=cbax, orientation='horizontal')

fig.savefig("snow.png")
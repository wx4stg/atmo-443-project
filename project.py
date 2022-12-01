#!/usr/bin/env python3

from os import path, listdir
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import colors as pltcolors
from cartopy import crs as ccrs
from cartopy import feature as cfeat
import pyart
import pandas as pd



targetLat = 42.493
targetLon = -79.272
# Load files
inputDir = path.join(path.dirname(__file__), "radar22")
radarDataFiles = sorted(listdir(inputDir))

px = 1/plt.rcParams["figure.dpi"]
lowfilt = plt.cm.bone_r(np.linspace(0, 0.6, 80))
specR = plt.cm.Spectral_r(np.linspace(0, 1, 200))
pink = plt.cm.PiYG(np.linspace(0, .25, 40))
purple = plt.cm.PRGn_r(np.linspace(0.75, 1, 40))
cArrLowFilt = np.vstack((lowfilt, specR, pink, purple))
cmapLowFilt = pltcolors.LinearSegmentedColormap.from_list("chaseSpectral_LowFilt", cArrLowFilt)

# Collect scan times and dBZ at the farm as we go
scanTimes = list()
ZsAtFarm = list()
ZDRSAtFarm = list()
CCSAtFarm = list()
for i in range(len(radarDataFiles)):
        radarPostageStampFig = plt.figure()
        gs = GridSpec(2, 1, figure=radarPostageStampFig, height_ratios=[50, 1])
        radarAx = radarPostageStampFig.add_subplot(gs[0, :], projection=ccrs.PlateCarree())
        rdaFile = path.join(inputDir, radarDataFiles[i])
        radarDat = pyart.io.read(rdaFile)
        # Gets the latitude and longitude of all gates
        ##lat = radarDat.get_gate_lat_lon_alt(0)[0]
        lon = radarDat.get_gate_lat_lon_alt(0)[1]
        # Gets the distance in degrees lat/lon of each gate to the farm mesonet
        ##latDisplacement = np.abs(radarDat.get_gate_lat_lon_alt(0)[0] - targetLat)
        ##lonDisplacement = np.abs(radarDat.get_gate_lat_lon_alt(0)[1] - targetLon)
        ##totalDisplacement = latDisplacement+lonDisplacement
        # Find which azimuth and range is the closest gate to the farm site
        ##closestAz, closestRange = np.where(totalDisplacement == np.min(totalDisplacement))
        ##closestAz, closestRange = closestAz[0], closestRange[0]
        # Get the value of reflectivity at the farm site
        ##ZatFarm = radarDat.fields["reflectivity"]["data"][closestAz, closestRange]
        ##ZDRatFarm = radarDat.fields["differential_reflectivity"]["data"][closestAz, closestRange]
        ##CCatFarm = radarDat.fields["cross_correlation_ratio"]["data"][closestAz, closestRange]
        # Create pyart RadarMapDisplay handle
        rmd = pyart.graph.RadarMapDisplay(radarDat)
        # Create ADRAD plan-position indicator
        rmd.plot_ppi_map("reflectivity", sweep=0, cmap=cmapLowFilt, vmin=-10, vmax=80, title_flag=False, colorbar_flag=False, ax=radarAx, fig=radarPostageStampFig, embellish=False)
        plotHandle = radarAx.get_children()[0]
        # Add a transparent square over the farm site to highlight the reflectivity value there
        radarAx.scatter(targetLon, targetLat, edgecolor="white", color="none", marker="s", transform=ccrs.PlateCarree())
        # Zoom in on the mesonet site +/- 0.4 degrees lat/lon
        radarAx.set_extent([targetLon-1, targetLon+1, targetLat-1, targetLat+1], crs=ccrs.PlateCarree())
        radarAx.add_feature(cfeat.STATES, edgecolor="white", linewidth=0.5)
        cbax = radarPostageStampFig.add_subplot(gs[-1, 0])
        xl = cbax.set_xlabel("Reflectivity (dBZ)")
        radarPostageStampFig.colorbar(plotHandle, cax=cbax, orientation="horizontal", extend="neither", label="Reflectivity (dBZ)")
        radarPostageStampFig.set_size_inches(1024*px, 1024*px)
        # Add the reflectivity value at the farm site to the title
        scanTime = pyart.util.datetime_from_radar(radarDat)
        radarAx.set_title(scanTime.strftime("KBUF 0.3° Reflectivity PPI\n%Y-%m-%d %H:%M:%S")+" UTC")
        radarAx.set_facecolor("black")
        radarPostageStampFig.set_facecolor("none")
        radarPostageStampFig.savefig(f"radarOut22/{i}.png")
        # Save scan time (for later)
        ##scanTimes.append(scanTime)
        # Save farm dBZ (for later)
        ##ZsAtFarm.append(ZatFarm)
        ##ZDRSAtFarm.append(ZDRatFarm)
        ##CCSAtFarm.append(CCatFarm)
        print("====")
        print(f"{i}/{len(radarDataFiles)}")
        ##print(scanTime)
        ##print(ZatFarm)
        ##print(ZDRatFarm)
        ##print(CCatFarm)
        print("====")
        plt.close("all")
##radarDataExtracted = pd.DataFrame({"Z": ZsAtFarm, "ZDR": ZDRSAtFarm, "CC": CCSAtFarm}, index=scanTimes)
##radarDataExtracted.to_csv("radar.csv")

#!/usr/bin/env python3

from os import path
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import dates as pltdates
from datetime import datetime as dt, timedelta


dirToRead = "20221118" # 20221118 20201224
asosFilePath = path.join(dirToRead, "asos.csv")
asosData = pd.read_csv(asosFilePath, comment="#").set_index("valid")
asosData.index = pd.to_datetime(asosData.index)  + timedelta(minutes=7)
asosData = asosData[asosData.index.minute == 0]
radarFilePath = path.join(dirToRead, "radar.csv")
reflData = pd.read_csv(radarFilePath).set_index("dates").replace("--", np.nan).astype(float)
reflData.index = pd.to_datetime(reflData.index)
scanTimes = reflData.index



hours = []
obs = []
MP = []
WSC = []
CSES = []
CSWS = []
RT = []
HS = []
thirtyMinInteg = []
shT = []
winners = []
brkVar = False
for i in range(1, len(asosData.index)):
    Zs = reflData[(reflData.index >= asosData.index[i] - timedelta(minutes=60)) & (reflData.index < asosData.index[i])]["Z"].values.flatten()
    zs = 10**(Zs/10)

    rainRateMP = ((200*zs**(-1))**(-0.625))
    rainRateWSC = ((300*zs**(-1))**(-1/1.4))
    rainRateCSES = ((130*zs**(-1))**(-1/2))
    rainRateCSWS = ((75*zs**(-1))**(-1/2))
    rainRateRT = ((250*zs**(-1))**(-1/1.2))
    snowRateHS = (zs/427)**(1/1.09)
    snowRate30mininteg = (zs/554)**(1/.88)
    snowShiina = (zs/158)**(1/1.1)

    selectedScanTimes = scanTimes[(scanTimes >= asosData.index[i] - timedelta(minutes=60)) & (scanTimes < asosData.index[i])]
    deltaSeconds = np.array([(selectedScanTimes[i+1]-selectedScanTimes[i]).total_seconds() for i in range(len(selectedScanTimes)-1)])
    adradTotalMP = np.sum(np.nan_to_num(rainRateMP[:-1])*(deltaSeconds/3600))
    adradTotalWSC = np.sum(np.nan_to_num(rainRateWSC[:-1])*(deltaSeconds/3600))
    adradTotalCSES = np.sum(np.nan_to_num(rainRateCSES[:-1])*(deltaSeconds/3600))
    adradTotalCSWS = np.sum(np.nan_to_num(rainRateCSWS[:-1])*(deltaSeconds/3600))
    adradTotalRT = np.sum(np.nan_to_num(rainRateRT[:-1])*(deltaSeconds/3600))
    adradTotalHS = np.sum(np.nan_to_num(snowRateHS[:-1])*(deltaSeconds/3600))
    adradTotal30mininteg = np.sum(np.nan_to_num(snowRate30mininteg[:-1])*(deltaSeconds/3600))
    adradTotalShiina = np.sum(np.nan_to_num(snowShiina[:-1])*(deltaSeconds/3600))

    print(f"Hour: {asosData.index[i].hour}")
    print(f"Observed Rainfall: {asosData.iloc[i]['p01m']} mm")
    print(f"ADRAD Rainfall (Marshall-Palmer): {adradTotalMP} mm")
    print(f"ADRAD Rainfall (Warm Season Convective): {adradTotalWSC} mm")
    print(f"ADRAD Rainfall (Cold Season East Stratiform): {adradTotalCSES} mm")
    print(f"ADRAD Rainfall (Cold Season West Stratiform): {adradTotalCSWS} mm")
    print(f"ADRAD Rainfall (Rosenfeld Tropical): {adradTotalRT} mm")
    print(f"ADRAD Snowfall (1-minute snowfall average method): {adradTotalHS} mm")
    print(f"ADRAD Snowfall (30-minute snowfall average method): {adradTotal30mininteg} mm")
    print(f"ADRAD Snowfall (Shiina snowfall method): {adradTotalShiina} mm")
    differenceDict = {
        "Marshall-Palmer": np.abs(adradTotalMP - asosData.iloc[i]['p01m']),
        "Warm Season Convective": np.abs(adradTotalWSC - asosData.iloc[i]['p01m']),
        "Cold Season East Stratiform": np.abs(adradTotalCSES - asosData.iloc[i]['p01m']),
        "Cold Season West Stratiform": np.abs(adradTotalCSWS - asosData.iloc[i]['p01m']),
        "Rosenfeld Tropical": np.abs(adradTotalRT - asosData.iloc[i]['p01m']),
        "1-minute snowfall average method": np.abs(adradTotalHS - asosData.iloc[i]['p01m']),
        "30-minute snowfall average method": np.abs(adradTotal30mininteg - asosData.iloc[i]['p01m']),
        "Shiina snowfall method": np.abs(adradTotalShiina - asosData.iloc[i]['p01m'])
    }
    winningVal = np.min(list(differenceDict.values()))
    winningIndices = [i for i, x in enumerate(differenceDict.values()) if x == winningVal]
    if len(winningIndices) == 8:
        winner = "Tie between all"
    elif len(winningIndices) > 1:
        winner = "Tie between " + ", ".join([list(differenceDict.keys())[i] for i in winningIndices])
    else:
        winner = list(differenceDict.keys())[list(differenceDict.values()).index(winningVal)]
    print(f"Winner: {winner}")
    hours.append(asosData.index[i])
    obs.append(asosData.iloc[i]['p01m'])
    MP.append(adradTotalMP)
    WSC.append(adradTotalWSC)
    CSES.append(adradTotalCSES)
    CSWS.append(adradTotalCSWS)
    RT.append(adradTotalRT)
    HS.append(adradTotalHS)
    thirtyMinInteg.append(adradTotal30mininteg)
    shT.append(adradTotalShiina)
    winners.append(winner)

outputData = pd.DataFrame({"Observed": obs, "Marshall-Palmer": MP, "Warm Season Convective": WSC, "Cold Season East Stratiform": CSES, "Cold Season West Stratiform": CSWS, "Rosenfeld Tropical": RT, "1-minute snowfall average method" : HS, "30-minute snowfall average method" : thirtyMinInteg, "Shiina snowfall method" : shT, "Winner" : winners}, index=hours)
outputData.to_csv(path.join(dirToRead, "output.csv"))


Zs = reflData["Z"].values.flatten()
zs = 10**(Zs/10)

rainRateMP = ((200*zs**(-1))**(-0.625))
rainRateWSC = ((300*zs**(-1))**(-1/1.4))
rainRateCSES = ((130*zs**(-1))**(-1/2))
rainRateCSWS = ((75*zs**(-1))**(-1/2))
rainRateRT = ((250*zs**(-1))**(-1/1.2))
snowRateHS = (zs/427)**(1/1.09)
snowRate30mininteg = (zs/554)**(1/.88)
snowShiina = (zs/158)**(1/1.1)

comparisonFig = plt.figure()
comparisonAx = comparisonFig.gca()
constantHrTimes = []
[constantHrTimes.append([datePoint-timedelta(minutes=59), datePoint]) for datePoint in asosData.index]
constantHrTimes = list(np.array(constantHrTimes).flatten())
constantHrVals = []
[constantHrVals.append([asosData.iloc[i]['p01m'], asosData.iloc[i]['p01m']]) for i in range(len(asosData))]
constantHrVals = list(np.array(constantHrVals).flatten())
comparisonAx.plot(constantHrTimes, constantHrVals, color="black", label="Observed Rainfall Rate", linewidth=1.5)
comparisonAx.plot(scanTimes, rainRateMP, color="maroon", label="KBUF (Marshall-Palmer)", linewidth=0.75)
comparisonAx.plot(scanTimes, rainRateWSC, color="red", label="KBUF (Warm Season Convective)", linewidth=0.75)
comparisonAx.plot(scanTimes, rainRateCSES, color="orange", label="KBUF (Cold Season East Stratiform)", linewidth=0.75)
comparisonAx.plot(scanTimes, rainRateCSWS, color="gold", label="KBUF (Cold Season West Stratiform)", linewidth=0.75)
comparisonAx.plot(scanTimes, rainRateRT, color="lime", label="KBUF (Rosenfeld Tropical)", linewidth=0.75)
comparisonAx.plot(scanTimes, snowRateHS, color="magenta", label="KBUF (1-minute snowfall average method)", linewidth=0.75)
comparisonAx.plot(scanTimes, snowRate30mininteg, color="cyan", label="KBUF (30-minute snowfall average method)", linewidth=0.75)
comparisonAx.plot(scanTimes, snowShiina, color="mediumslateblue", label="KBUF (Shiina snowfall method)", linewidth=0.75)

comparisonAx.legend()
comparisonAx.grid(True)
comparisonAx.xaxis.set_major_formatter(pltdates.DateFormatter("%H:%M"))
comparisonAx.set_xlabel("Time (UTC)")
comparisonAx.set_ylabel("Rainfall Rate (mm/hr)")
px = 1/plt.rcParams['figure.dpi']
comparisonFig.set_size_inches(1024*px, 1024*px)
comparisonFig.savefig(path.join(dirToRead, "rainfallRateComparison.png"), transparent=True)
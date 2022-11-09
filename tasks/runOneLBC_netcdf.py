from boundaryConditions import LBC_writer
import sys
import yaml
import bunch
import glob
import os
import xarray as xr
import datetime
from utils import getConfigurationByID,nextDay
from inputFieldProvider import selectInputFile
import datetime

ids = list(map(float, sys.argv[1].replace('[', '').replace(',', '').split()))
pointName = sys.argv[2].replace(',', '')
wd = sys.argv[3].replace(',', '')
confPath = sys.argv[4].replace(',', '')
ds = sys.argv[5:]

days = []
for d in ds:
    days.append(d.replace('[', '').replace(',', '').replace(']', ''))
days.append(nextDay(days[-1]))

print('id', ids)
print('pointName', pointName)
print('wd', wd)
print('conf', confPath)
print('days', days)

conf = getConfigurationByID(confPath, 'config')
yesterday = (datetime.datetime.strptime(days[0], '%Y%m%d') - datetime.timedelta(days=1)).strftime('%Y%m%d')
parent = []
wind = []

for day in days:
    print(f'selecting wind file for {day}')
    #windFile=selectInputFile(conf,'meteoData',day)

    wnds = selectInputFile(conf, 'meteoData', day)
    if isinstance(wnds, list):
        for wnd in wnds:
            if wnd not in wind:
                wind.append(wnd)
    else:
        wind.append(wnds)

    print(f'selecting wave file for {day}')

    prnts = selectInputFile(conf, 'parentWave', day)
    if isinstance(prnts, list):
        for prnt in prnts:
            if prnt not in parent:
                parent.append(prnt)
    else:
        parent.append(prnts)

LBC_writer(wd, confPath, wind, parent).onePointBC(ids, pointName)

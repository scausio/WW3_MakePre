import xarray as xr
import os
import numpy as np
import glob
import pandas as pd
import subprocess
from datetime import datetime,timedelta,date
from glob import glob
from utils import dayBefore,dateFormatter,checkFunctionCompleted,seaoverland,getDaysBetweenDates, kelvin_to_celsius,splitDate
from natsort import natsorted

def getUnstHydroData(path,date, conf, fieldType):
    """
    get data from UNSTructuded input
    Parameters
    ----------
    path
    conf
    fieldType

    Returns
    -------

    """
    #field type can be: waterVelocity,waterLevel,surfaceTemp

    fo = xr.open_mfdataset(path,combine='by_coords').sel({conf.copernicusFiles.parentHydro.time: splitDate(date)})

    fo_sub = fo.isel(**{ conf.copernicusFiles.parentHydro.depth:0})
    #fo[conf.copernicusFiles.parentHydro.time] = fo[conf.copernicusFiles.parentHydro.time].dt.floor('H')  # round hours

    if fieldType=='waterVelocity':
        print ('getting currents')
        try:
            fo_sub = fo_sub.compute()
        except:
            pass
        fo_sub.to_netcdf(f'CUR_{date}.nc')

    if fieldType == 'waterLevel':

        try:
            fo_sub = fo_sub.compute()
        except:
            pass
        fo_sub.to_netcdf(f'LEV_{date}.nc')

    if fieldType == 'surfaceTemp':

        try:
            fo_sub = fo_sub.compute()
        except:
            pass
    return fo_sub

def selectInputFile(conf,inputType,day):
    fi = glob(buildFilePath(conf, inputType, day))
    return fi

def cutField(conf, fieldName, lon, lat):
    tempFold = os.getcwd()
    cut = conf.cutArea

    if not os.path.exists(tempFold):
        os.makedirs(tempFold)
    latId = (lat >= cut.lat[0]) & (lat <= cut.lat[1])
    lonId = (lon >= cut.lon[0]) & (lon <= cut.lon[1])

    return lonId, latId

def buildFilePath(conf, dataset,day,yesterday=False):
    if dataset in ['waterVelocity','waterLevel','surfaceTemp']:
        fileconf = conf.copernicusFiles.parentHydro.fileConf
        filename = fillPathTemplates(os.path.join(fileconf.fileTemplate.base,fileconf.fileTemplate.name), fileconf, day, yesterday)
    else:
        fileconf = conf.copernicusFiles[dataset].fileConf
        filename = fillPathTemplates(os.path.join(fileconf.fileTemplate.base,fileconf.fileTemplate.name), fileconf, day, yesterday)
    return filename

def fillPathTemplates(template, fileConf, startdate,yesterday):
    d = {'date': str(startdate),
         'year': str(startdate)[:4],
         'month': str(startdate)[4:6],
         'day': str(startdate)[6:8],
         'yesterday': dayBefore(str(startdate)) if not yesterday else yesterday,
         'freq': fileConf.freq,
         'producer': fileConf.producer,
         'parameter': fileConf.parameter,
         'config': fileConf.config,
         'region': fileConf.region,
         'version': fileConf.version
         }
    print (template)
    return template.format(d=d)

def getHydroData(path,date, conf, fieldType):
    #field type can be: waterVelocity,waterLevel,surfaceTemp

    fo = xr.open_mfdataset(path,combine='by_coords').sel({conf.copernicusFiles.parentHydro.time:splitDate(date)})

    lat,lon = fo[conf.copernicusFiles.parentHydro.lat][:], fo[conf.copernicusFiles.parentHydro.lon][:]

    # get data inside box
    lonId, latId = cutField(conf, 'oceanMsk', lon, lat)
    fo_sub = fo.isel(**{conf.copernicusFiles.parentHydro.lat:np.where(latId)[0], conf.copernicusFiles.parentHydro.lon:np.where(lonId)[0], conf.copernicusFiles.parentHydro.depth:0})

    if fieldType=='waterVelocity':
        print ('getting currents')
        u = fo_sub[conf.copernicusFiles.parentHydro.waterVelocity.u]
        u.values[u.values==0]=np.nan
        u.values = np.array([seaoverland(np.ma.masked_array(f, mask=np.isnan(f)), 20) for f in u.values])
        print ('u current sol', u.values.shape)
        v = fo_sub[conf.copernicusFiles.parentHydro.waterVelocity.v]
        v.values[v.values==0]=np.nan
        v.values = np.array([seaoverland(np.ma.masked_array(f, mask=np.isnan(f)), 20) for f in v.values])
        print ('v current sol', v.values.shape)
        fo_sub[conf.copernicusFiles.parentHydro.waterVelocity.u].values = u
        fo_sub[conf.copernicusFiles.parentHydro.waterVelocity.v].values = v
        try:
            fo_sub = fo_sub.compute()
        except:
            pass
        fo_sub.to_netcdf('CUR.nc')

    if fieldType == 'waterLevel':
        lev= fo_sub[conf.copernicusFiles.parentHydro.waterLevel.ssh]
        print ('lev', lev.values.shape)
        lev.values[lev.values==0]=np.nan
        lev.values = np.array([seaoverland(np.ma.masked_array(f, mask=np.isnan(f)), 20) for f in lev.values])
        fo_sub = fo.isel(lat=latId, lon=lonId)
        fo_sub[conf.copernicusFiles.parentHydro.waterLevel.ssh].values = lev
        try:
            fo_sub = fo_sub.compute()
        except:
            pass
        fo_sub.to_netcdf('LEV.nc')

    if fieldType == 'surfaceTemp':
        sst = fo_sub[conf.copernicusFiles.parentHydro.surfaceTemp.sst]
        sst.values[sst.values==0]=np.nan
        sst.values = np.array([seaoverland(np.ma.masked_array(f, mask=np.isnan(f)), 20) for f in sst.values])
        print ('sst current sol', sst.values.shape)
        fo_sub[conf.copernicusFiles.parentHydro.surfaceTemp.sst].values = sst
        try:
            fo_sub = fo_sub.compute()
        except:
            pass

    return fo_sub


def getMeteoData(path, conf,date):
    print (path)

    fo = xr.open_mfdataset(path,combine='by_coords').sel({conf.copernicusFiles.meteoData.time:splitDate(date)})

    lon, lat = fo[conf.copernicusFiles.meteoData.lon][:], fo[conf.copernicusFiles.meteoData.lat][:]

    lonId, latId = cutField(conf, 'meteoMsk', lon, lat)
    #fo_sub = fo.isel(**{conf.copernicusFiles.parentHydro.lat:np.where(latId)[0], conf.copernicusFiles.parentHydro.lon:np.where(lonId)[0]})
    land=fo['LSM'][:, latId, lonId].compute()
    #land=fo_sub['LSM'][:]
    msk_land=land.data>0.5
    msk_land=np.logical_not(msk_land).astype(int)
    
    #fo[conf.copernicusFiles.meteoData.time] = fo[conf.copernicusFiles.meteoData.time].dt.floor('H')  # round hours
    u=fo[conf.copernicusFiles.meteoData.u][:, latId, lonId].compute()
    #u=fo_sub[conf.copernicusFiles.meteoData.u][:].compute()
    u_data=u.values
    u_data=u_data*msk_land
    u_data[u_data==0]=np.nan
    print ('u', u.shape)
    u_data=np.array([seaoverland(np.ma.masked_array(f, mask=np.isnan(f)),20)for f in u_data])
    print ('sol u completed')

    v=fo[conf.copernicusFiles.meteoData.v][:, latId, lonId].compute()
    #v=fo_sub[conf.copernicusFiles.meteoData.v][:].compute()
    v_data=v.values
    v_data=v_data*msk_land
    v_data[v_data==0]=np.nan
    v_data=np.array([seaoverland(np.ma.masked_array(f, mask=np.isnan(f)),20)for f in v_data])
    print ('sol v completed')

    t2m=fo[conf.copernicusFiles.meteoData.T2M][:, latId, lonId].compute()
    #t2m=fo_sub[conf.copernicusFiles.meteoData.T2M][:].compute()
    t2m_data=t2m.values
    t2m_data=t2m_data*msk_land
    t2m_data[t2m_data==0]=np.nan
    t2m_data=np.array([seaoverland(np.ma.masked_array(f, mask=np.isnan(f)),20)for f in t2m_data])
    print ('sol t2m completed')

    fo_sub=fo.isel(lat= latId,lon= lonId)
    fo_sub[conf.copernicusFiles.meteoData.v].values  = v_data #fix
    fo_sub[conf.copernicusFiles.meteoData.u].values  = u_data #fix
    fo_sub[conf.copernicusFiles.meteoData.T2M].values = t2m_data

    fo_sub=fo_sub.compute()
    fo_sub=fo_sub.sortby(conf.copernicusFiles.meteoData.lat, ascending=True)
    fo_sub.to_netcdf(f'WND_{date}.nc')

    return fo_sub

class WNSgetter():
    def __init__(self, conf, date, fillValue=-9999):
        self.workPath = os.getcwd()
        self.conf = conf
        self.date = date
        self.fillValue = fillValue
        self.duration=1
        print('getting WNS')

    def manageSST(self):

        sstFiles=[]
        ssts=selectInputFile(self.conf, 'surfaceTemp', self.date)
        if isinstance(ssts, list):
            for sst in ssts:
                if sst not in sstFiles:
                    sstFiles.append(sst)
        else:
            sstFiles.append(ssts)
        if self.conf.forcings.deltaT.type.upper()=='REG':
            SST= getHydroData(sstFiles,self.date, self.conf, 'surfaceTemp')
        elif self.conf.forcings.deltaT.type.upper()=='UNST':
            SST = getUnstHydroData(sstFiles,self.date, self.conf, 'surfaceTemp')
        elif self.conf.forcings.deltaT.type.upper() == 'F':
            pass
        else:
            exit ('please check.conf.forcings.wind.type')
        return SST


        return SST

    def manageWind(self):
        date = dateFormatter(self.date)
        day_before=(date - timedelta(days=1)).strftime('%Y%m%d')

        days=getDaysBetweenDates(day_before,self.duration+2)


        windFiles = []
        for day in days:
            wnds = selectInputFile(self.conf, 'meteoData', day)
            if isinstance(wnds, list):
                for wnd in wnds:
                    if wnd not in windFiles:
                        windFiles.append(wnd)
            else:
                windFiles.append(wnds)

        if self.conf.forcings.wind.type.upper()=='REG':
            WND= getMeteoData(windFiles, self.conf, 'meteoData')
        elif self.conf.forcings.wind.type.upper()=='UNST':
            exit ('.conf.forcings.wind.type == UNST not yet implemented')
        elif self.conf.forcings.wind.type.upper() == 'F':
            pass
        else:
            exit ('please check.conf.forcings.wind.type')
        return WND


    def getDeltaT(self, wind, hydro,conf):
        print ('computing DT')

        #wind['T2M'].values[wind['T2M'].values==self.fillValue]=np.nan
        
        sst=hydro[conf.copernicusFiles.parentHydro.surfaceTemp.sst] 
        sst.values[sst.values==0]=np.nan
        sst.values= np.array([seaoverland(np.ma.masked_array(f, mask=np.isnan(f)), 20) for f in sst])

        regridded_T=sst.interp(**{conf.copernicusFiles.parentHydro.lat:wind[conf.copernicusFiles.meteoData.lat].values,
                                                                          conf.copernicusFiles.parentHydro.lon:wind[conf.copernicusFiles.meteoData.lon].values})
        sea_T=regridded_T.interp(method='nearest',**{conf.copernicusFiles.parentHydro.time: wind[conf.copernicusFiles.meteoData.time]})
        #conf.copernicusFiles.parentHydro.time: conf.copernicusFiles.meteoData.time}
        #nearest_T.values[nearest_T.values==0]=np.nan
        sea_T= np.array([seaoverland(np.ma.masked_array(f, mask=np.isnan(f)), 20) for f in sea_T])
        air_T= wind[conf.copernicusFiles.meteoData.T2M] 
        air_T.values[air_T.values==self.fillValue]=np.nan
        if air_T.attrs['units'] in ['K', 'Kelvin', 'k', 'degK']:
            air_T=kelvin_to_celsius(air_T)
        dT = air_T.values- sea_T
        wind["DT"]=((conf.copernicusFiles.meteoData.time, conf.copernicusFiles.meteoData.lat, conf.copernicusFiles.meteoData.lon),dT)
        wind["DT"].attrs["units"]="Celsius"
        wind=wind.drop(conf.copernicusFiles.meteoData.T2M)
        wind.to_netcdf(F'WNS_{date}.nc')

    def getWNS(self):
        WND = self.manageWind()  # ,T2M,latW,lonW
        SST = self.manageSST()  # ,latT,lonT # u can force computation with force=True
        self.getDeltaT(WND, SST,self.conf)


class WNDgetter():
    def __init__(self, conf, date, fillValue=-9999):
        self.workPath = os.getcwd()
        self.conf = conf
        self.date = date
        self.fillValue = fillValue
        self.duration=1
        print('getting WND')

    def manageWind(self):

        windFiles = []
        wnds = natsorted(selectInputFile(self.conf, 'meteoData', self.date))

        if isinstance(wnds, list):
            for wnd in wnds:
                if wnd not in windFiles:
                    windFiles.append(wnd)
        else:
            windFiles.append(wnds)
        WND=getMeteoData(windFiles, self.conf,self.date)
        return WND

    def getWND(self):
        WND = self.manageWind()
        print('... done ...')


class CURgetter():

    def __init__(self, conf, date, fillValue=-9999):
        self.workPath = os.getcwd()
        self.conf = conf
        self.date = date
        self.fillValue = fillValue
        self.duration=1
        print('getting Currents')


    def manageCUR(self):
        curFiles=[]

        curs = selectInputFile(self.conf, 'waterVelocity', self.date)
        if isinstance(curs, list):
            for cur in curs:
                if cur not in curFiles:
                    curFiles.append(cur)
        else:
            curFiles.append(curs)
        if self.conf.forcings.currents.type.upper()=='REG':
            CUR= getHydroData(curFiles,self.date, self.conf, 'waterVelocity')
        elif self.conf.forcings.currents.type.upper()=='UNST':
            CUR= getUnstHydroData(curFiles,self.date, self.conf, 'waterVelocity')
        elif self.conf.forcings.currents.type.upper() == 'F':
            pass
        else:
            exit ('please check.conf.forcings.currents.type')
        return CUR

    def getCUR(self):
        CUR = self.manageCUR()  # ,latT,lonT # u can force computation with force=True
        print('... done ...')


class LEVgetter():
    def __init__(self,  conf, date, fillValue=-9999):
        self.workPath = os.getcwd()
        self.conf = conf
        self.date = date
        self.fillValue = fillValue
        self.duration=1
        print('getting Sea Level')

    def manageLEV(self):

        levFiles = []
        levs = selectInputFile(self.conf, 'waterLevel', self.date)
        if isinstance(levs, list):
            for lev in levs:
                if lev not in levFiles:
                    levFiles.append(lev)
        else:
            levFiles.append(levs)

        if self.conf.forcings.water_level.type.upper()=='REG':
            LEV= getHydroData(levFiles,self.date, self.conf, 'waterLevel')
        elif self.conf.forcings.water_level.type.upper()=='UNST':
            LEV= getUnstHydroData(levFiles,self.date, self.conf, 'waterLevel')
        elif self.conf.forcings.water_level.type.upper() == 'F':
            pass
        else:
            exit ('please check.conf.forcings.water_level.type')
        return LEV


    def getLEV(self):
        outname = 'LEV.data'
        LEV = self.manageLEV()
        print('... done ...')

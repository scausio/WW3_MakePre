import numpy as np
from datetime import datetime, timedelta
import xarray as xr
from munch import DefaultMunch,munchify
import yaml
import os
from glob import glob
from natsort import natsorted

def splitDate(date):
    return f"{date[:4]}-{date[4:6]}-{date[6:8]}"

def kelvin_to_celsius(k):
    c=k-273.15
    c.attrs['units']='Celsius'
    return  c


def myMFdataset(fileList,bbox=False):
    if bbox:
        try:
            out = [xr.open_dataset(f).isel(latitude=range(int(bbox[2]),int(bbox[3])+1),longitude=range(int(bbox[0]),int(bbox[1])+1)) for f in fileList]
        except:
            out = [xr.open_dataset(f).isel(lat=range(int(bbox[2]), int(bbox[3]) + 1), lon=range(int(bbox[0]), int(bbox[1]) + 1)) for
               f in fileList]
    else:
        out=[xr.open_dataset(f) for f in fileList]
    #check CI
    out_=[]
    for f in out:
        try:
            out_.append(f.drop(['CI','precip']))
        except:
            out_.append(f)

    return xr.concat(out_,dim='time')

def nautical2cart(degree):
    return (90-degree)%360

def getFreqs(minFreq, nFreq, ratio):
    freqs = []
    for f in range(nFreq):
        freqs.append(minFreq)
        minFreq = minFreq * ratio
    return np.array(freqs)


def getDirs(nDir):
    return np.arange(0, 360, (360 / nDir))
    #step=360 / nDir
    #return  np.arange(0, 360, step)+(step/2)


def degToRad(deg):
    return np.pi * deg / 180


def stringAllocation(n, string, d):
    string = str(string)

    dsize = int(n - len(string))
    if dsize < 0:
        string = string[:n]
    if d == 'back':
        return (' ' * dsize) + string
    else:
        return string + (' ' * dsize)


def cart2pol(x, y):
    return np.arctan2(y, x)

def bearingToStandAngle(deg):
    deg= 90-deg
    deg[deg<0]+=360
    deg[deg > 360] -= 360
    return deg


def reverseDir(deg):
    deg= deg-180
    deg[deg<0]+=360
    deg[deg > 360] -= 360
    return deg

def reverse90(deg):
    deg= deg-90
    deg[deg<0]+=360
    deg[deg > 360] -= 360
    return deg

def add90(deg):
    deg= deg+90
    deg[deg<0]+=360
    deg[deg > 360] -= 360
    return deg

def fixBCdir(deg, negativeVal=False):
    if negativeVal:
        deg=deg0_360(deg)
    #if isinstance(deg,np.ndarray):

    #    m_low=deg <= 180
    #    m_hi=deg > 180
    #
    #    deg[m_low]=np.abs(deg[m_low]- 180)
    #    deg[m_hi]=360 - (np.abs(deg[m_hi] - 180))

    #else:
    #    if deg <= 180:
    #        deg = np.abs(deg - 180)
    #    else:
    #        deg = 360 - (np.abs(deg - 180))
    #return deg
    return (180-deg)%360


def deg0_360(deg):
    if type(deg).__module__ == np.__name__:
        deg[deg<0]=360+deg[deg<0]
    else:
        if deg<0:
            deg+=360
    return deg


def antiCwiseSorter(coords):
    origin = np.mean(np.array(coords), axis=0)
    coordsCentered = coords - origin
    phi = cart2pol(coordsCentered[:, 0], coordsCentered[:, 1])
    return np.argsort(phi)


def duplicatesId(bc):
    #bc=np.array([bc])
    uniquesId = []
    uniques = []
    print (bc, 'BC')
    for n, i in enumerate(bc[:, :-1]):
        if tuple(i) not in uniquesId:
            uniquesId.append(tuple(i))
            uniques.append(n)
    return uniques


def seaoverland(data, iterations=1, copy=False):
    """
    Python implementation of G. Girardi's seaoverland.f90.

    Extends grids defined only at sea onto the land such that bilinear
    interpolation of the grid is also possible close to the coastline. The
    procedure determines the value of an undefined grid point by averaging the
    defined neighbour points among the 8 closest neighbours. Every iteration
    therefore extends the grid with one additional row of points towards the
    coastline.

    With copy set to True a copy of the data is returned leaving the original
    untouched. Otherwise the data is modified in place and returned.

    Parameters
    ----------
    data : numpy.ma.masked_array
        Grid (2 dimensions) to be extrapolated.
    iterations : int, optional, default: 1
        Number of iterations (i.e. extent of extrapolation).
    copy : boolean, optional, default: False
        Create a copy of data instead of modifying it in place.
    """
    if copy:
        data = np.ma.copy(data)

    if not isinstance(data, np.ma.masked_array) or not data.mask.any():
        return data

    for _ in range(iterations):
        shifted = []
        ni, nj = data.shape
        for i in range(-1, 2):
            for j in range(-1, 2):
                if i != 0 or j != 0:
                    # Shift the grid by i horizontally and j vertically and
                    # append it to the array. Shifted grids are 2 units smaller
                    # in both dimensions to accomodate the shift.
                    shifted.append(data[1 + i:ni - 1 + i, 1 + j:nj - 1 + j])

        # Calculate the mean value of the shifted grids to obtain the
        # approximated values. Only non-masked entries are taken into account.
        approx = np.ma.mean(shifted, axis=0)

        # Create a view without the outer points (so it is the same size as the
        # shifted grids), then copy the approximated values for the cells that
        # are masked.
        view = data[1:-1, 1:-1]
        np.copyto(view, approx, where=(view.mask & ~approx.mask))

        # Combine the two masks, unmasking values that were masked in view but
        # have been successfully approximated.
        view.mask &= approx.mask

    return data


def dateFormatter(date):
    try:
        return datetime.strptime(date, "%Y%m%d")
    except:
        return date


def getConfigurationByID(path,confId):
    globalConf = yaml.load(open(os.path.join(path,"conf.yaml")))
    return munchify(globalConf[confId])

def dayBefore(date):
    date=dateFormatter(date)
    return  (date - timedelta(days=1)).strftime('%Y%m%d')
def nextDay(date):
    date=dateFormatter(date)
    return  (date + timedelta(days=1)).strftime('%Y%m%d')


class Paths():
    def __init__(self,conf,args):
        self.workingDir=args.workingDir
        self.tasksDir=conf.run_specifics.base
        self.outDir=os.path.join(args.workingDir,'inputs')
        print ('path saved')


class Checks():
    def __init__(self,rundir):
        self.rundir=rundir

    def RUNcomplete(self):
        if os.path.isfile(os.path.join(self.rundir, 'RUN_complete')):
            exit('run complete')
        else:
            print('Run starting ...')
    def LBCcomplete(self):
         return os.path.isfile(os.path.join(self.rundir, 'LBC_complete'))

    def SBCcomplete(self):
        return os.path.isfile(os.path.join(self.rundir, 'SBC_complete'))

    def TMPLcomplete(self):
        return os.path.isfile(os.path.join(self.rundir, 'TMPL_complete'))

def checkFunctionCompleted(checks,funct,check=True):
    if check:
        try:
            checks()
        except:
            print ('%s not implemented'%checks)
        if checks():
            return
    funct()

def checkCompleted(paths, day):
    flag= os.path.exists(os.path.join(paths.outDir,'completed'))
    print (flag,'flag')
    if flag:
        print ('{} run already done'.format(day))
    return flag

def getDaysBetweenDates(startDate,duration):
    date=dateFormatter(startDate)
    return [(date + timedelta(days=n)).strftime('%Y%m%d') for n in range(0,duration+1)]

def getStations(path):
    print (os.path.join(path, '*.nc'))
    return np.unique([os.path.basename(f).split('_')[1] for f in natsorted(glob(os.path.join(path, '*.nc')))])





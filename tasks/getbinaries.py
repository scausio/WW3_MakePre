import os
import shutil
import subprocess
from glob import glob
from natsort import natsorted
from boundaryConditions import LBCindex_getter,LBC_writer,LBC_extractor
from inputFieldProvider import  WNSgetter, WNDgetter, \
    CURgetter, LEVgetter,selectInputFile
from templates import Templates
from utils import checkFunctionCompleted,getDaysBetweenDates, myMFdataset,getStations
from bjobs import Bjobs
import numpy as np
import xarray as xr


class GetBinaries():
    def __init__(self,SBC,duration, conf,rundir,args, checks,submitter):
        self.checks = checks
        self.conf = conf
        self.rundir = rundir
        self.startdate = args.startDate
        self.workingDir = args.workingDir
        self.submitter=submitter
        self.runDuration=duration
        self.fields=SBC

    def getSymLinks(self):
        if not os.path.islink(os.path.join(self.workingDir,'wd',os.path.basename(self.conf.model.grid ))):
            os.symlink(self.conf.model.grid,os.path.join(self.workingDir,'wd',os.path.basename(self.conf.model.grid )))
        if not os.path.islink(os.path.join(self.workingDir, 'wd','WND.nc')):
            os.symlink(os.path.join(self.workingDir, 'inputs','WND.nc'),
                       os.path.join(self.workingDir, 'wd','WND.nc'))


    def processTemplates(self):
        #try:
        if self.conf.lateralBoundaries.type.upper() in ['FROMPARENT', 'COMPUTED', "AVERAGED", "PARTITIONED"] :
            print (os.path.join(self.workingDir, 'inputs', 'id_*spec.nc'))
            spectra = natsorted(glob(os.path.join(self.workingDir, 'inputs', 'id_*spec.nc')))
            spectra_def=xr.open_dataset(spectra[0])
        else:
            spectra=False
            spectra_def=False
        Templates(self.conf, self.startdate,self.runDuration, os.path.join(self.workingDir,'wd'),spectra,spectra_def).generate()


    def processINP(self):
        print (f'I m here {os.getcwd()}')

        if not os.path.exists('mod_def.ww3'):
            self.submitter.systemCommand("{}".format(os.path.join(self.conf.executable, 'ww3_grid')))
        for field in self.fields:
            self.writePrncINP(field)
            self.submitter.systemCommand(os.path.join(self.conf.executable, 'ww3_prnc'))

        if self.conf.lateralBoundaries.type.upper() in ['FROMPARENT', 'COMPUTED', "AVERAGED", "PARTITIONED"]:
            if not os.path.exists('nest.ww3'):
                self.submitter.systemCommand("{}".format(os.path.join(self.conf.executable, 'ww3_bounc')))

    def writePrncINP(self, field):
        input_type = 'LL' if field in ['WND', 'WNS'] else 'AI'
        with open( 'ww3_prnc.inp', 'w') as out:
            out.write('$\n ')
            out.write(f"'{field}' '{input_type}' T T\n")
            out.write('$\n ')
            if field in ['CUR', 'LEV']:
                out.write("node \n")
            else:
                out.write("{lon} {lat}\n".format(lon=self.conf.copernicusFiles.meteoData.lon,
                                                 lat=self.conf.copernicusFiles.meteoData.lat))
            out.write('$\n ')
            if field == 'WND':
                out.write("{U} {V}\n".format(U=self.conf.copernicusFiles.meteoData.u,
                                             V=self.conf.copernicusFiles.meteoData.v))
                out.write(' WND.nc\n')
            elif field == 'WNS':
                out.write("{U} {V} DT\n".format(U=self.conf.copernicusFiles.meteoData.u,
                                                V=self.conf.copernicusFiles.meteoData.v))
                out.write(' WNS.nc\n')
            elif field == 'CUR':
                out.write("{U} {V} \n".format(U=self.conf.copernicusFiles.parentHydro.waterVelocity.u,
                                              V=self.conf.copernicusFiles.parentHydro.waterVelocity.v))
                out.write(' CUR.nc\n')
            elif field == 'LEV':
                out.write("{SSH} \n".format(SSH=self.conf.copernicusFiles.parentHydro.waterLevel.ssh))
                out.write(' LEV.nc\n')
            out.write('$\n ')

    def moveBins(self):
        os.chdir(os.path.join(self.workingDir,'wd'))
        print (glob(os.path.join(os.getcwd(),'*.ww3')))
        [shutil.move(f, os.path.join(self.workingDir,'inputs',os.path.basename(f))) for f in glob(os.path.join(os.getcwd(),'*.ww3'))]


    def main(self):
        print ('preprocessing starting')
        self.getSymLinks()
        os.chdir(os.path.join(self.workingDir,'wd'))
        self.processTemplates()
        self.processINP()
        self.moveBins()


class LateralBoundaryCondition():
    def __init__(self,conf,runningdir,workingdir,startdate):
        self.conf=conf
        self.workingdir=workingdir
        self.startdate=startdate
        self.runDuration=1
        self.rundir=runningdir

    def compute(self):
        print('Computing Indexes for LBC')
        day=self.startdate
        windFiles = []

        wnds = selectInputFile(self.conf, 'meteoData', day)
        if isinstance(wnds, list):
            for wnd in wnds:
                if wnd not in windFiles:
                    windFiles.append(wnd)
        else:
            windFiles.append(wnds)

        parentFiles = []
        prnts = selectInputFile(self.conf, 'parentWave', day)
        if isinstance(prnts, list):
            for prnt in prnts:
                if prnt not in parentFiles:
                    parentFiles.append(prnt)
        else:
            parentFiles.append(prnts)

        # wd, conf, bc, parent, wind
        lbc = LBCindex_getter(os.path.dirname(os.path.dirname(self.workingdir)),
                              self.conf.lateralBoundaries.path,
                              parentFiles,
                              windFiles)


        pointNames=[i for i in range(1,len(lbc.ids)+1)]

        LBC_writer(self.workingdir, self.conf, windFiles, parentFiles).allPointsBC(lbc.ids, pointNames,self.startdate)

        self.spectra= natsorted(glob(os.path.join(self.rundir, 'id_*.spc')))

    def getIDXSboxes(self,idxs):
        """

        :param idxs: indices from LBCindex_getter
        :return: parent, wind arrays with box indices to cut parent and wind files. Each of them is [minX,maxX,minY,maxY]. The size of the box is +2 for X and Y
        """

        latParent=idxs[:,1]
        lonParent = idxs[:, 0]
        latWind = idxs[:, 4]
        lonWind = idxs[:,3]

        boxParent=[np.nanmin(lonParent)-2,np.nanmax(lonParent)+2,np.nanmin(latParent)-2,np.nanmax(latParent)+2]
        boxWind = [np.nanmin(lonWind)-2, np.nanmax(lonWind)+2, np.nanmin(latWind)-2, np.nanmax(latWind)+2]

        return boxParent,boxWind

    def converterBoxIdxs(self,ids):
        ids.T[0]-=np.nanmin(ids.T[0])
        ids.T[1]-=np.nanmin(ids.T[1])
        ids.T[3]-=np.nanmin(ids.T[3])
        ids.T[4]-=np.nanmin(ids.T[4])
        
        ids.T[0]+=2
        ids.T[1]+=2
        ids.T[3]+=2
        ids.T[4]+=2
        return ids

    def extract(self):
        print('Extracting Indexes for LBC')
        days = getDaysBetweenDates(self.startdate, self.runDuration + 1)

        windFiles = []
        for day in days:
            wnds = selectInputFile(self.conf, 'meteoData', day)
            if isinstance(wnds, list):
                for wnd in wnds:
                    if wnd not in windFiles:
                        windFiles.append(wnd)
            else:
                windFiles.append(wnds)

        spectraFiles = []
        for day in days:
            spcs = selectInputFile(self.conf, 'parentWave', day)
            if isinstance(spcs, list):
                for spc in spcs:
                    if spc not in spectraFiles:
                        spectraFiles.append(spc)
            else:
                spectraFiles.append(spcs)

        # wd, conf, bc, parent, wind
        lbc = LBCindex_getter(self.rundir,
                              os.path.join(self.conf.base, self.conf.model.lateralBoundaries.path),
                              spectraFiles,
                              windFiles)

        boxParent, boxWind=self.getIDXSboxes( lbc.ids)

        winds=myMFdataset(windFiles,boxWind)
        spectra=myMFdataset(spectraFiles,boxParent)
        pointNames=[i for i in range(1,self.converterBoxIdxs(lbc.ids)+1)]
        writer=LBC_extractor(self.rundir, self.workingdir, winds, spectra)

        writer.getSpectraBC(self.converterBoxIdxs(lbc.ids), pointNames,self.startdate)

    def main(self):
        print ('LBC preprocessing starting')
        if self.conf.lateralBoundaries.type.upper() not in ["AVERAGED", "PARTITIONED", "FROMSPECTRA", 'F', 'FALSE']:
            print (self.conf.lateralBoundaries.type)
            exit('please check input at model.lateralBoundaries.type in conf.yaml')
        os.chdir(self.rundir)
        if self.conf.lateralBoundaries.type.upper() in ["AVERAGED", "PARTITIONED"]:
            LBC = self.compute()
            print ('computing boundary conditions')

        if self.conf.lateralBoundaries.type.upper() == 'FROMSPECTRA':
            print ('extracting boundary conditions')
            LBC = self.extract()


class SurfaceBoundaryCondition():
    def __init__(self,conf, startdate,workingDir):
        self.conf=conf
        self.startdate=startdate
        self.workingDir=workingDir
        os.chdir(self.workingDir)

    def main(self):
        print('starting with preprocessing Surface Boundaries')
        self.fields=[]
        if str(self.conf.forcings.deltaT.type).upper()!='F':
            self.fields.append('WNS')
            if not os.path.exists(os.path.join(self.workingDir ,f"WNS_{self.startdate}.nc")):
                WNSgetter(self.conf, self.startdate).getWNS()

                print('WNS processed')
        if str(self.conf.forcings.wind.type).upper()!='F':
            self.fields.append('WND')
            if not os.path.exists(os.path.join(self.workingDir,f"WND_{self.startdate}.nc")):
                WNDgetter(self.conf, self.startdate).getWND()

                print('WND processed')
        else:
            print('No WIND')
        if str(self.conf.forcings.currents.type).upper() !='F':
            self.fields.append('CUR')
            if not os.path.exists(os.path.join(self.workingDir, f"CUR_{self.startdate}.nc")):
                CURgetter(self.conf, self.startdate).getCUR()

            print('CUR processed')
        if str(self.conf.forcings.water_level.type).upper()!='F':
            if not os.path.exists(os.path.join(self.workingDir, f"LEV_{self.startdate}.nc")):
                LEVgetter(self.conf, self.startdate).getLEV()
                self.fields.append('LEV')
                print('LEV processed')
        os.chdir(self.workingDir)
        return self.fields

def rebuildSBC(wd):
    print ('concatenating SBC ')
    fields=['WND','WNS', 'LEV', 'CUR']
    for field in fields:
        if not os.path.exists(os.path.join(wd,'inputs', f'{field}.nc')):
            filelist=natsorted(glob (os.path.join(wd,'wd','SBC',f'{field}*nc')))
            if len(filelist)>0:
                buffer=[xr.open_dataset(f) for f in filelist]
                xr.concat(buffer,dim='time').to_netcdf(os.path.join(wd,'inputs', f'{field}.nc'))

def rebuildLBC(wd):
    print ('concatenating LBC ')
    stations=getStations(os.path.join(wd,'wd','LBC'))
    for station in stations:
        if not os.path.exists(os.path.join(wd,'inputs', f'id_{station}_spec.nc')):
            filelist=natsorted(glob (os.path.join(wd,'wd','LBC',f'id_{station}_*_spec.nc')))
            if len(filelist)>0:
                buffer=[xr.open_dataset(f) for f in filelist]
                xr.concat(buffer,dim='time').to_netcdf(os.path.join(wd,'inputs', f'id_{station}_spec.nc'))




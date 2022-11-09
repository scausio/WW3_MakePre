import os
from datetime import timedelta
import pybars
from utils import dateFormatter
import numpy as np


class Templates():
    def __init__(self,conf,startdate,duration,rundir,spectra,spectra_def):
        self.conf=conf
        self.startdate=startdate
        self.rundir=rundir
        self.spectra=spectra
        self.compiler = pybars.Compiler()
        self.runDuration=duration
        self.data=self.getDict(conf,spectra_def)

    def getDict(self,conf,spectra_def):

        startRunDate = dateFormatter(self.startdate).strftime('%Y%m%d %H%M%S')

        endRunDate = (dateFormatter(self.startdate) + timedelta(days=self.runDuration) + timedelta(
            seconds=conf.model.outputTime)).strftime('%Y%m%d %H%M%S')

        if self.conf.lateralBoundaries.type.upper() in ['AVERAGED','PARTITIONED', 'F', "FALSE"]:
            frequencyIncrement=conf.model.frequencyIncrement
            minFrequency=conf.model.minFrequency
            NoFrequencies=conf.model.NoFrequencies
            NoDirections=conf.model.NoDirections
        else:
            freqs=spectra_def[conf.copernicusFiles.parentWave.freq].values
            frequencyIncrement=np.round(freqs[1]/freqs[0],2)
            minFrequency=np.min(freqs)
            NoFrequencies=len(freqs)
            NoDirections=len(spectra_def[conf.copernicusFiles.parentWave.dir].values)

        data = {'model': conf.model.name,
                'namelist': [{'n': i} for i in conf.model.namelist],
                'frequencyIncrement': frequencyIncrement,
                'minFrequency': minFrequency,
                'NoFrequencies': NoFrequencies,
                'NoDirections': NoDirections,
                'gridType': conf.model.gridType,
                'depthLimit': conf.model.depthLimit,
                'minDepth': conf.model.minDepth,
                'itype': 3,
                'inputField': [{'field': 'CUR'},
                               {'field': conf.copernicusFiles.meteoData.ww3Name},
                               {'field': 'LEV'}],
                'global_TS': conf.model.timeStep.global_TS,
                'CFL_TS': conf.model.timeStep.CFL_TS,
                'refraction_TS': conf.model.timeStep.refraction_TS,
                'minimum_TS': conf.model.timeStep.minimum_TS,
                'timeInterval': conf.model.outputTime,
                'startDate': startRunDate,
                'endDate': endRunDate,
                'bottomFile': os.path.basename(conf.model.grid),
                'curFlag': conf.forcings.currents.type,
                'levFlag': conf.forcings.water_level.type,
                'windFlag': conf.forcings.wind.type,
                'bc_file': conf.lateralBoundaries.path,
                'scaleFactor': conf.model.depthScaleFactor,
                'restart': 86400* self.runDuration,
                'wd': self.rundir,
                }

        if conf.model.gridType == 'UNST':
            self.ext= '.unst'
            data['noPoints'] = conf.model.noNodes
            data['obstFile_local'] = ''
            data['obstFile_shadow'] = ''

        else:
            self.ext = '.reg'
            data['lonResolution'] = conf.grid.lonResolution
            data['latResolution'] = conf.grid.latResolution
            data['minLon'] = conf.grid.minLon
            data['minLat'] = conf.grid.minLat
            data['lonSize'] = conf.grid.lonSize
            data['latSize'] = conf.grid.latSize
            data['maskFile'] = os.path.basename(conf.grid.maskFile)
            data['obstFile_regular'] = os.path.basename(conf.grid.obstFile_regular)

        boundNodes = []
        try:
            fo = open(data['bc_file']).readlines()
            boundNodes = [{'nodes': '  %s 1 F' % str(i).replace('\n', '')} for i in fo]
        except Exception as e:
            print(e)

        data['bc'] = boundNodes
        data['spcs'] = []
        try:
            [data['spcs'].append({'spc': spc}) for spc in self.spectra]
        except:
            pass
        return data

    def generate(self):
        listTemplates = [f for f in os.listdir(self.conf.templates) if f.endswith(self.ext) or f.endswith('.tmpl')]
        [self.fillTemplate( tmpl)for tmpl in listTemplates]

    def fillTemplate(self, tmpl):
        print('compiling template %s' % tmpl)
        print (self.rundir)
        with open(os.path.join(self.conf.templates, tmpl), 'r') as s:
            source = s.read()
            template = self.compiler.compile(source)
            outTmpl = template(self.data)
            tmplName = tmpl.split('.')[0]
            if tmplName == 'submitModel':
                ext = ".sh"
            else:
                ext = ".inp"
            outputFile = open(os.path.join(self.rundir, tmplName + ext), "w")
            [outputFile.write(str(line)) for line in outTmpl]
            outputFile.close()

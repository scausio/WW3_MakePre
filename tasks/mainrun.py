#!/usr/bin/env python

import os
import subprocess
from getbinaries import GetBinaries,SurfaceBoundaryCondition, rebuildSBC,LateralBoundaryCondition,rebuildLBC
from utils import Paths, Checks
from datetime import datetime,timedelta
from utils import getConfigurationByID,checkFunctionCompleted,checkCompleted,dateFormatter,nextDay
import argparse
import shutil
from bjobs import Bjobs

parser = argparse.ArgumentParser()
parser.add_argument('--startDate', '-s', help='first run day')
parser.add_argument('--duration', '-d', help='run duration')
parser.add_argument('--workingDir', '-wd', help='run working dir')
args = parser.parse_args()


class NextRun():
    def __init__(self,conf,totalDuration,runDuration, paths,args,submitter):
        self.conf=conf
        self.submitter=submitter
        self.totalDuration=totalDuration
        self.runDuration=runDuration
        self.workingdir=paths.workingDir
        self.paths=paths

        if self.totalDuration>self.runDuration:
            self.buildCMD(args)
            self.lastRun=False
        else:
            self.lastRun=True
            self.cmd='run over'

    def buildCMD(self, args):
        nextDuration = self.totalDuration - self.runDuration
        nextStartDay = dateFormatter(args.startDate) + timedelta(days=self.runDuration)

        wd = args.workingDir
        cmd= setRunCMD(args)
        self.cmd = cmd.format(sla='' ,p=self.conf.model.project_queue,queue=self.conf.model.serial_queue, startdate=nextStartDay.strftime('%Y%m%d'), duration=nextDuration, wd=wd)
        return self.cmd
    def run(self):
        if self.totalDuration >self.runDuration:
            os.chdir(os.path.join(self.conf.base, 'tasks'))
            self.submitter.systemCommand(self.cmd)
        else:
            print (self.totalDuration, self.runDuration, self.workingdir)
            #os.remove(glob(os.path.join(self.workingdir,'restart*.ww3'))[-1])
            if self.conf.post.regrid.flag.upper() in ['T','TRUE','Y','YES']:
                regridder=Regridding(self.conf,self.paths)
                for area in self.conf.post.regrid.area:
                    print(area)
                    vars = self.conf.post.regrid.area[area].variables
                    if vars:
                       regridder.regriddingMain(area)
            # light-traffic file for formatter
            open(os.path.join(self.workingdir,'output','run_complete'), 'w').close() 
            print ('run over')


class Submit():
    def __init__(self,conf,paths):
        self.paths=paths
        self.conf=conf

    def systemCommand(self,cmd):
        print ('Submitting %s'%cmd)
        subprocess.run(str(cmd), check=True,shell=True)


def main():
    print (args)
    startDate=str(args.startDate)
    duration=int(args.duration)
    conf = getConfigurationByID(args.workingDir, 'config')

    runDir = os.path.dirname(args.workingDir)
    runDuration=duration

    paths= Paths(conf,args)
    submitter = Submit(conf,paths)
    nextrun = NextRun(conf,duration,runDuration,paths, args, submitter)
    #if not os.path.exists( os.path.join(runDir,'inputs','WND.nc')):
    checks = Checks(runDir)
    os.chdir(paths.outDir)

    runningDate=startDate
    for i in range(duration):
        print (runningDate)
        print (f'computing SBC for {runningDate}')
        SBC=SurfaceBoundaryCondition(conf, runningDate, os.path.join(args.workingDir,'wd','SBC')).main()
        print (f'computing LBC for {runningDate}')
        #LBC = LateralBoundaryCondition(conf,runDir, os.path.join(args.workingDir, 'wd', 'LBC'),runningDate).main()
        runningDate = nextDay(runningDate)
    rebuildSBC(args.workingDir)
    rebuildLBC(args.workingDir)
    GetBinaries(SBC,runDuration, conf, runDir, args, checks, submitter).main()
    # Getting binaries




    # getting binaries



    #nextrun.run()

#if __name__=='main':
main()


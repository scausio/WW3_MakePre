#!/usr/bin/env python
import argparse
from datetime import datetime
import subprocess
from utils import getConfigurationByID
import os
import shutil

parser = argparse.ArgumentParser()
parser.add_argument('--startDate', '-s', help='first run day')
parser.add_argument('--duration', '-d', help='run duration')
parser.add_argument('--endDate', '-e', help='last run day')

args = parser.parse_args()

def get_runRange(args):
    startdate = args.startDate

    try:
        duration = int(args.duration)
        pass
    except:
        try:
            endDate = str(args.endDate)
            dt = datetime.strptime(endDate, "%Y%m%d") - datetime.strptime(startdate, "%Y%m%d")
            duration=dt.days
        except Exception as e:
            print (e, 'please add the end-date or duration of the run')
    return startdate, duration

def setWorkingDir( startdate,duration):
    conf_path = os.path.dirname(os.getcwd())
    conf = getConfigurationByID(conf_path, 'config')
    wd = os.path.join(conf.run_specifics.wd,  f"{startdate}_{duration}d")
    os.makedirs(wd,exist_ok=True)

    shutil.copy(os.path.join(conf_path, 'conf.yaml'), wd)
    os.makedirs(os.path.join(wd, 'inputs'), exist_ok=True)
    os.makedirs(os.path.join(wd, 'wd'), exist_ok=True)
    os.makedirs(os.path.join(wd,  'wd', 'SBC'), exist_ok=True)
    os.makedirs(os.path.join(wd,  'wd', 'LBC'), exist_ok=True)
    return conf,wd

def main():

    startDate, duration = get_runRange(args)
    conf,wd = setWorkingDir(startDate,duration)
    subprocess.call(f'./mainrun.py -s {startDate} -d {duration} -wd {wd}' , shell=True)

main()

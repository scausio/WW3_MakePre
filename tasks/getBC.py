import os

import numpy as np

from toSpectra import JONSWAP

from datetime import datetime, timedelta
from utils import getFreqs,getDirs,degToRad


class JONSWAP():
    def __init__(self,hs,tp,dir,ndir,freqs,gamma=3.3):
        self.hs= hs
        self.tp=tp
        self.dir=dir
        self.ndir=ndir
        self.freqs=freqs
        self.gamma=gamma


    def jsEnergy(self):
        mo = []
        for f in self.freqs:
            freqPeak = 1. / self.tp  # fpk
            freqPeak4 = freqPeak ** 4  # fpk4
            alpha = ((self.hs ** 2) * freqPeak4) / ((0.06533 * ((self.gamma ** 0.8015) + 0.13467)) * 16)
            f4 = f ** 4
            f5 = f ** 5
            cpshap = 1.25 * freqPeak4 / f4
            if cpshap > 10:
                ra = 0
            else:
                ra = (alpha / f5) * np.exp(-cpshap)
            if f <= freqPeak:
                sigma = 0.07
            else:
                sigma = 0.09
            apshap = 0.5 * ((f - freqPeak) / (sigma * freqPeak)) ** 2
            if apshap > 10:
                syf = 1.
            else:
                ppshap = np.exp(-apshap)
                syf = self.hs ** ppshap
            specDens = syf * ra / (2 * np.pi * f)
            mo.append(specDens)
        return np.array(mo)

    def gammLn(self,val):
        coefs=[76.18009173,-86.50532033,24.01409822, -1.231739516,.120858003e-2,-.536382e-5]
        stp=2.50662827465
        fpf=4.5
        tmp=(val+0.5+fpf)*np.log(val+fpf)-(val+fpf)
        ser=1
        for coef in coefs:
            ser=ser+coef/(val)
            val=+1
        return tmp+(np.log(stp*ser))

    def gammaF(self,val):
        return np.exp(np.clip(self.gammLn(val),-30,30))


    def main(self):
        energySpectra = self.jsEnergy()
        self.directions = [np.pi * i / 180 for i in np.arange(0, 360, (360 / self.ndir))]
        dirSpread = 360 / self.ndir
        rad = np.pi * self.dir / 180

        if dirSpread < 12:
            ctot = (2. ** dirSpread) * (self.gammaF(0.5 * dirSpread + 1.)) ** 2 / (np.pi * self.gammaF(dirSpread + 1.))
        else:
            ctot = np.sqrt(0.5 * dirSpread / np.pi) / (1. - 0.25 / dirSpread)
        dirSpectra = []
        for d in self.directions:
            acos = np.cos(d - rad)

            if acos > 0:
                cdir = ctot * np.clip((acos ** dirSpread), 1.e-10, None)
            else:
                cdir = 0

            for spectrum in energySpectra:
                dirSpectra.append(cdir * spectrum)

        return energySpectra, np.array(dirSpectra).reshape(self.ndir, len(self.freqs))



def removeListduplicates(listSeq):
    unique=[]
    for i in listSeq:
        if i in unique:
            pass
        else:
            unique.append(i)
    return np.array( unique)

# def stringAllocation(n, string, reducingSize=0):
#     string=str(string)
#     size= n * 8
#     dsize=int(size-len(string)-1)-reducingSize
#     if dsize<0:
#         string=string[:7]
#     return string+ (' '* dsize)

# def stringAllocation(n, string):
#     string=str(string)
#     buf= n *' '
#     buf[:len(string)]=string
#     return buf

def stringAllocation(n, string,d):
    string = str(string)

    dsize = int(n - len(string) )
    if dsize < 0:
        string = string[:n]
    if d=='back':
        return (' ' * dsize)+string
    else:
        return string + (' ' * dsize)



class BC():

    def __init__(self,wd,lon,lat, hs,tp,dir,minFreq,ratio,nFreq,nDir,days):
        self.wd=wd
        self.hs=hs
        self.tp=tp
        self.dir=dir
        self.minFreq=minFreq
        self.ratio=ratio
        self.nFreq=nFreq
        self.nDir=nDir
        self.days=days
        self.lon=lon
        self.lat=lat

    def bc_header(self):
        outname = os.path.join(self.wd, 'bc_header.spc')
        self.freqs = getFreqs(self.minFreq, self.nFreq,
                              self.ratio)
        if os.path.exists(outname):
            return
        else:
            dirs = degToRad(getDirs(self.nDir))
            fbins = [np.format_float_scientific(f, precision=3, exp_digits=2, unique=False) for f in self.freqs]
            f_chunks = (' '.join(str(l) + '\n' * (n % 8 == 7) for n, l in enumerate(fbins)))

            dbins = [np.format_float_scientific(f, precision=3, exp_digits=2, unique=False) for f in dirs]
            d_chunks = (' '.join(str(l) + '\n' * (n % 8 == 7) for n, l in enumerate(dbins)))

            if d_chunks[-1] != '\n':
                d_chunks = d_chunks + '\n'

            with open(outname, 'w')as o:
                o.write("'WAVeWATCH III SPeCTRA'     %02d    %02d     1 'MULTIGRID GLOBAL05+OTHeRS     '\n " % (
                    self.nFreq,
                    self.nDir))

                o.write(f_chunks)
                o.write('\n ')
                o.write(d_chunks)
                o.close()


    def onePointBC(self):
        self.bc_header()
        pointName='%s_%s'%(self.lon,self.lat)
        spName = 'id_{}'.format(pointName)
        outFile = os.path.join(self.wd, '{}.spc'.format(spName))

        with open(outFile, 'w') as o:

            with open(os.path.join(self.wd, 'bc_header.spc')) as head:
                for line in head:
                    o.write(line)
                head.close()

            for t in self.days:

                wSpeed = 5
                wDir = 150
                try:
                    timestamp = datetime.strptime(str(t.data).split('.')[0], "%Y-%m-%dT%H:%M:%S").strftime(
                        '%Y%m%d %H%M%S')
                except:
                    timestamp = t.strftime(
                        '%Y%m%d %H%M%S')

                hs = self.hs
                tp = self.tp
                dir = self.dir
                nonDimSpec, dimSpec = JONSWAP(hs, tp, dir, nDir, self.freqs).main()

                o.write('%s\n' % timestamp)
                o.write("'{pointId}'{lat}{lon}{depth}{wSpeed}{wDir}   0.00 270.0\n  ".format(
                    pointId=stringAllocation(10, spName[:10], 'forth'),
                    lat='%7.2f' % self.lat,
                    lon='%7.2f' % self.lon,
                    depth='%9.1f' % float(500),
                    wSpeed='%7.2f' % np.round(wSpeed, 2),
                    wDir='%6.1f' % np.round(wDir, 1)))

                spc = [np.format_float_scientific(i, precision=3, exp_digits=2, unique=False) for i in
                       dimSpec.ravel()]

                chunkSpc = ('  '.join(str(l) + '\n' * (n % 7 == 6) for n, l in enumerate(spc)))
                if chunkSpc[-1] != '\n':
                    chunkSpc = chunkSpc + '\n'
                o.write(chunkSpc)
            o.write('\n')
        o.close()






hs=5
tp=10
dir=100

minFreq=0.0373
ratio=1.1
nFreq=32
nDir=24

wd='/Users/scausio/Desktop/hawaySpecsTest'

startday='20080225 000000'
h=3
endDay='20080227 000000'

days=[]

d=datetime.strptime(startday, '%Y%m%d %H%M%S')

while d!=datetime.strptime(endDay, '%Y%m%d %H%M%S'):
    days.append(d)
    d=d+ timedelta(hours=h)

lons=[-155,-154, -153]
lats=[15,16,17]

dirs=getDirs(nDir)
freqs=getFreqs(minFreq,nFreq,ratio)

if not os.path.exists(wd):
    os.makedirs(wd)
for lon,lat in zip(lons,lats):
    a=BC(wd,lon,lat, hs,tp,dir,minFreq,ratio,nFreq,nDir,days).onePointBC()
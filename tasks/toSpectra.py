import numpy as np
from utils import getDirs, getFreqs, fixBCdir,nautical2cart


class SPECTRA():
    def __init__(self,hs,tp,dir,ndir,freqs,spc_type):
        self.hs= hs
        self.tp=tp
        self.dir=dir
        self.ndir=ndir
        self.freqs=freqs
        self.spc_type=spc_type

    def gaussian_density(self,gaus_width=0.5): #CHECK GAUSWIDTH
        self.gaus_width=gaus_width
        freqPeak = 1. / self.tp  # fpk
        freqPeak4 = freqPeak ** 4  # fpk4
        aux1 = (self.hs ** 2) / (16 * np.sqrt(np.pi * 2) * self.gaus_width)
        aux3 = 2 * self.gaus_width ** 2
        mo = []
        for f in self.freqs:
            aux2= (f-((np.pi*2)/self.tp))**2
            ra= aux1* np.exp(-1 * aux2/aux3)/f
            mo.append(ra)
        return np.array(mo)

    def pm(self):

        mo = []
        freqPeak = 1. / self.tp  # fpk
        freqPeak4 = freqPeak ** 4  # fpk4
        alpha = ((self.hs ** 2) * freqPeak4)*5 / 16
        for f in self.freqs:
            sf=f/(2*np.pi)
            sf4 = sf ** 4
            sf5 = sf ** 5

            specDens = (alpha / sf5) * np.exp(-(5*sf4)/(4*sf4))/(f)

            mo.append(specDens)
        return np.array(mo)

    def jonswap_energy(self,gamma=3.3):
        self.gamma = gamma
        freqPeak = 1. / self.tp  # fpk
        freqPeak4 = freqPeak ** 4  # fpk4
        alpha = ((self.hs ** 2) * freqPeak4) / ((0.06533 * ((self.gamma ** 0.8015) + 0.13467)) * 16)
        mo = []
        for f in self.freqs:

            f4 = f ** 4
            f5 = f ** 5

            # f4 = f ** 4
            # f5 = f ** 5
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
            specDens = syf * ra / (f*2*np.pi)
            mo.append(specDens)
        return np.array(mo)

    def gammLn(self, val):
        coefs = [76.18009173, -86.50532033, 24.01409822, -1.231739516, .120858003e-2, -.536382e-5]
        stp = 2.50662827465
        fpf = 4.5
        x = val - 1
        tmp = x + fpf
        tmp = (x + 0.5) * np.log(tmp) - (tmp)
        ser = 1
        for coef in coefs:
            val += 1
            ser = ser + coef / (x)

        return tmp + (np.log(stp * ser))

    def gammaF(self,val):
        return np.exp(np.clip(self.gammLn(val),-30,30))

    def main(self, energy_coefficient=False):

        """

        :param energy_coefficient: this is to change gamma in jonsap or gaus_width in gaussian
        :return:
        """
        if self.spc_type=='jonswap':
            if energy_coefficient:
                energySpectra = self.jonswap_energy(energy_coefficient)
            else:
                energySpectra = self.jonswap_energy()
        elif self.spc_type=='gaussian':
            if energy_coefficient:
                energySpectra = self.gaussian_density(energy_coefficient)
            else:
                energySpectra = self.gaussian_density()
        elif self.spc_type == 'pm':
            energySpectra = self.pm()
        else:
            exit('not valid spectra type')
        #self.directions = np.deg2rad(nautical2cart(getDirs(self.ndir)))
        self.directions = np.deg2rad(getDirs(self.ndir))
        dirSpread = 360 / self.ndir
        adir = np.deg2rad((180-self.dir)%360)
        #adir = np.deg2rad(self.dir)
        #if DSHAPL==1:
        #    ms = np.nanmax([np.deg2rad(dirSpread)**(-2) - 2, 1])
        #else:
        ms= dirSpread
        if ms < 12:
            ctot = (2. ** ms) * (self.gammaF(0.5 * ms + 1.)) ** 2 / (np.pi * self.gammaF(ms + 1.))
        else:
            ctot = np.sqrt(0.5 * ms/ np.pi) / (1. - 0.25 / ms)
        dirSpectra = []

        for d in self.directions:
            acos = np.cos(d - adir)

            if acos > 0:
                cdir = ctot * np.clip((acos ** ms), 1.e-10, None)
            else:
                cdir = 1.e-10

            for spectrum in energySpectra:
                dirSpectra.append(cdir * spectrum)

        return energySpectra, np.array(dirSpectra).reshape(self.ndir, len(self.freqs))


#   IF YOU WANT TO SEE THE SPECTRA WITH PURE PYTHON SCRIPTING:
#   for d in self.directions: SHOULD BE
#   for d in self.directions[::-1]:
#   AND
#   for f in self.freqs[::-1]: SHOULD BE
#   for f in self.freqs:
#
#
# def mo(spec,freqs):
#     freqs=np.insert(freqs, 0,0)
#     bandWidth = freqs[1:] - freqs[:-1]
#     return np.sum((spec*bandWidth)/2)
# def hmo0(spec,freqs):
#     return 4*np.sqrt(mo(spec, freqs))*0.95
#
# def getFreqs(minFreq, nFreq, ratio):
#     freqs = []
#     for f in range(nFreq):
#         freqs.append(minFreq)
#         minFreq = minFreq * ratio
#     return np.array(freqs)
#
#
# def getDirs(nDir):
#     #return np.arange(0, 360, (360 / nDir))
#     step=360 / nDir
#     return  np.arange(0, 360, step)+(step/2)
#
# def test():
#     import matplotlib.pyplot as plt
#
#     #0.76500005 6.3030005 193.61
#     minFreq=0.05
#     ratio=1.1
#     nFreq=32
#
#     nDir=24
#
#
#     hs=0.76
#     tp=6.3
#     dir=90
#
#
#
#     dirs=getDirs(nDir)
#     freqs=getFreqs(minFreq,nFreq,ratio)
#
#     nonDimSpec, dimSpec = SPECTRA(hs, tp, dir, nDir, freqs,'jonswap').main(3.3)
#
#     plt.plot(nonDimSpec)
#     plt.show()
#     plt.imshow(dimSpec.T, origin='lower')
#     plt.colorbar()
#     plt.show()
#     #
#     # plt.plot(dimSpec.T[:,12])
#     # plt.show()
#     #js.header()
#     #plt.scatter(freqs,mo,c='r')
#
#     #a=hmo0(dirS,freqs)
#
#     #print mo(dirS,freqs)- mo(enerS,freqs)
#
#
# if __name__ == "__main__":
#     test()

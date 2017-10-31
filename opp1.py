"""Oppgave 1: Tidserie av havbolger"""
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy import stats

#Fra Read1D:
FichData = 'BayOfBiscay.mat'
file_ = h5py.File(FichData, 'r')
t1 = file_['t']
eta = file_['eta']
t1 = t1[0,:]
eta1 = eta[0,:]
file_.close()

FichData = 'NewYearWave.mat'
file_ = h5py.File(FichData, 'r')
t2 = file_['t']
eta = file_['eta']
t2 = t2[0,:]
eta2 = eta[0,:]
file_.close()

#a) Plott overflatehevningen som en tidsserie.

def plot_over(t,eta, dataname):
    plt.plot(t, eta)
    plt.title('Oppgave 1a - %s' % dataname)
    plt.xlabel('$t$[s]')
    plt.ylabel('$\eta_1$[m]')
    plt.show()

plot_over(t1,eta1, 'BayOfBiscay')
plot_over(t2,eta2, 'NewYearWave')

#b) Beregn signifikant bolgehoyde Hs fra standardavviket til overflatehevningen

def sign_wave1(eta, dataname):
    Hs = 4*np.std(eta, dtype=np.float64) #signifikant bolgehoyde
    print "Signifikant bolgehoyde til %s : %f " % (dataname, Hs)

sign_wave1(eta1, 'BayOfBiscay')
sign_wave1(eta2, 'NewYearWave')

#c) Beregn frekvensopplosningen df og Nyquistfrekvensen.

def freq(t, dataname):
    T = len(t)
    df = float(1./T)
    Nqf = df/2.
    print "%s: Frekvensopplosning=%f Nyquistfrekvens=%f" % (dataname, df, Nqf)

freq(t1, 'BayOfBiscay')
freq(t2, 'NewYearWave')

#d) Estimer varians/effekt spekteret S(f) fra tidsserien for overflatehevningen.

def spectrum(t, eta, dataname):
    T = len(t)
    X = np.fft.ifft(eta)
    S = np.abs(X)**2
    fs = T/(t[-1]-t[0])         #samplingfrekvensen
    f = np.linspace(0, fs, T)
    Sn = S/(f+fs)               #normaliseringsbetingelse
    plt.plot(f, Sn)
    plt.title('Oppgave 1d - %s' % dataname)
    plt.xlabel('Frekvens[Hz]')
    plt.ylabel('Spektrum')
    plt.show()
    return S, f

S1, f1 = spectrum(t1, eta1, 'BayOfBiscay')
S2, f2= spectrum(t2, eta2, 'NewYearWave')

#e) Beregn signifikant bolgehoyde Hs fra effekt spekteret S(f).

def sign_wave2(f, dataname, S):
    Hspec = 4*np.sqrt(np.sum(S)*(f[-1]-f[0])/len(f))
    print "Signifikant bolgehoyde til %s : %f " % (dataname, Hspec)

sign_wave2(f1, 'BayOfBiscay', S1)
sign_wave2(f2, 'NewYearWave', S2)

#f) Plott det en-sidige spekteret for f >= 0 for frekvenser opp til Nyquistfrekvensen.
# Bruk dimensjonelle akser, forste akse skal vare oppgitt i Hertz (Hz).
#Plottet skal ha dobbelt logaritmiske akser, dvs. lages med loglog.

def ensidig_spec(f,S,dataname):
    S_en = 2*(S[:int(len(S)+1)/2])
    f_en = f[:int(len(f)+1)/2]
    plt.loglog(f_en, S_en)
    plt.title("Oppgave 1f - Ensidig spektrum for %s" % dataname)
    plt.xlabel("Frekvens[Hz]")
    plt.ylabel("spektrum S(f)")
    plt.xlim(10E-4, 10)
    plt.ylim(10E-10, 10E-1)
    plt.show()

ensidig_spec(f1,S1,'BayOfBiscay')
ensidig_spec(f2,S2,'NewYearWave')

# g) Gjett en passende potens

def bay_potens(f,S,dataname):
    S_en = 2*(S[:int(len(S)+1)/2])
    f_en = f[:int(len(f)+1)/2]

    first_third = int(len(S_en)*3/4.)
    slope, intercept, r_value, p_value, std_err = stats.linregress(np.log(f_en[first_third:]),np.log(S_en[first_third:]))
    print slope,intercept

    plt.loglog(f_en, S_en)
    #plt.plot(np.log(f_en),slope*np.log(f_en) + intercept)
    plt.title("Oppgave 1g %s" % dataname)
    plt.xlabel("Frekvens[Hz]")
    plt.ylabel("spektrum S(f)")
    plt.xlim(0.05, 1)       #zoomer inn der det ser ut som spekteret folger en potenslov
    plt.ylim(10E-10, 10E-1)
    plt.hold('on')

    y_line1 = (f_en**(-2.89628725475))*1E-6 #hvor potensen er hentet fra resultatet fra slope
                                            #og ganget med 1E-6 for a legge linjen over dataen
    plt.loglog(f_en, y_line1)
    plt.legend(["S(f)", "n=-2.8962"])
    plt.show()

def new_potens(f,S,dataname):
    S_en = 2*(S[:int(len(S)+1)/2])
    f_en = f[:int(len(f)+1)/2]

    first_third = int(len(S_en)*1/4.)
    slope, intercept, r_value, p_value, std_err = stats.linregress(np.log(f_en[first_third:]),np.log(S_en[first_third:]))
    print slope,intercept

    plt.loglog(f_en, S_en)
    plt.title("Oppgave 1g- %s" % dataname)
    plt.xlabel("Frekvens[Hz]")
    plt.ylabel("spektrum S(f)")
    plt.xlim(0.05, 1)
    plt.ylim(10E-10, 1)
    plt.hold('on')

    y_line1 = (f_en**(-2.83802220014))*1E-6
    plt.loglog(f_en, y_line1)
    plt.legend(["S(f)", "n=-2.8380"])
    plt.show()

bay_potens(f1,S1,'BayOfBiscay')
new_potens(f2,S2,'NewYearWave')

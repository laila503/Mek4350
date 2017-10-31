import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import h5py
from mpl_toolkits.mplot3d import Axes3D

#Lager tomme lister
x_bimodal = []
y_bimodal = []
eta_bimodal = []
x_swell = []
y_swell = []
eta_swell = []
x_windsea = []
y_windsea = []
eta_windsea = []

#Leser data og legger inn i de tomme listene:
for i in range(0,30):
    FichData = 'ReWFLin00000_000%02d.mat' % i
    file_ = h5py.File(FichData, 'r')
    x = file_['x']
    y = file_['y']
    eta = file_['eta']
    x1 = x[0,:]
    y1 = y[0,:]
    eta1 = eta[:,:]
    if 0<i<10:
        x_bimodal.append(x1)
        y_bimodal.append(y1)
        eta_bimodal.append(eta1)
    if 10<i<20:
        x_swell.append(x1)
        y_swell.append(y1)
        eta_swell.append(eta1)
    if 20<i<30:
        x_windsea.append(x1)
        y_windsea.append(y1)
        eta_windsea.append(eta1)
    file_.close()

#a) beregn opplosningene i Fourier rom og Nyquist bolgetall
def opplosning(x_data, y_data):
    kxs = (2*np.pi)*len(x_data)/(x_data[-1]-x_data[0])
    kys = (2*np.pi)*len(y_data)/(y_data[-1]-y_data[0])

    dkx = kxs/len(x_data)
    dky = kys/len(y_data)

    Nqx = dkx/2.
    Nqy = dky/2.
    return dkx, dky, Nqx, Nqy

def pen_utskrift1(dataname, x_data, y_data):
    print "%s:" % dataname
    print "  dkx    |   dky    |   Nqx    |   Nqy"
    print "-----------------------------------------"
    for j in range(0,9):
        dkx, dky, Nqx, Nqy =  opplosning(x_data[j], y_data[j])
        print " %f | %f | %f | %f " % (dkx, dky, Nqx, Nqy)
    print " "

# pen_utskrift1("BINOMAL", x_bimodal, y_bimodal)
# pen_utskrift1("SWELL", x_swell, y_swell, eta_swell)
# pen_utskrift1("WINDSEA", x_windsea, y_windsea, eta_windsea)

#b) estimer det to-dimensjonale spektret
def spectrum(x_data, y_data, eta_data):
    X = np.fft.fft2(eta_data)#/(2*np.pi)**2
    Xk = np.fft.fftshift(X)

    kxs = (2*np.pi)*len(x_data)/(x_data[-1]-x_data[0])
    kys = (2*np.pi)*len(y_data)/(y_data[-1]-y_data[0])
    kx = np.linspace(-kxs/2.,kxs/2.*(1-2./(len(x_data))),len(x_data))
    ky = np.linspace(-kys/2.,kys/2.*(1-2./(len(y_data))),len(y_data))

    S = np.absolute(Xk)**2*kxs*kys/(len(x_data)*len(y_data))
    # k = [kx,ky]
    return kx,ky,S


def plot_3d(x_data, y_data, eta_data):
    eta_mean = np.mean(np.array(eta_data),axis=0)
    kx,ky,S = spectrum(x_data[-1], y_data[-1], eta_mean)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    X,Y = np.meshgrid(kx,ky)
    ax.plot_surface(X,Y,S,cmap=cm.coolwarm)
    plt.show()
    return eta_mean

#plot_3d(x_bimodal, y_bimodal, eta_bimodal)
#plot_3d(x_swell, y_swell, eta_swell)
#plot_3d(x_windsea, y_windsea, eta_windsea)

#e
def sign_wave(data, x_data, y_data,eta_data):
    eta_mean = np.mean(np.array(eta_data),axis=0)
    kx,ky,S = spectrum(x_data[-1], y_data[-1], eta_mean)

    Hs_stand = 4*np.std(eta_mean)
    Hs_spec = 4*np.sqrt(np.sum(S)*(kx[-1] - kx[0])/len(kx)*(ky[-1] - ky[0])/len(ky))
    print data
    print "Signifikant bolgehoyde fra standardaviket: %f" % Hs_stand
    print "Signifikant bolgehoyde fra spekteret: %f" % Hs_spec

# sign_wave("BIMODAL", x_bimodal, y_bimodal, eta_bimodal)
# sign_wave("SWELL", x_swell, y_swell, eta_swell)
# sign_wave("WINDSEA", x_windsea, y_windsea, eta_windsea)

#f
def plot_wave(x_data, y_data,eta_data):
    eta_mean = np.mean(np.array(eta_data),axis=0)
    eta = eta_mean
    x = x_data[0]
    y = y_data[0]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    X,Y = np.meshgrid(x,y)
    ax.plot_surface(X,Y,eta,cmap=cm.coolwarm)
    plt.show()

#plot_wave(x_bimodal, y_bimodal, eta_bimodal)
#plot_wave(x_swell, y_swell, eta_swell)
#plot_wave(x_windsea, y_windsea, eta_windsea)

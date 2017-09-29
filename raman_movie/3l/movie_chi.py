import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from yambopy import *
import argparse
from time import time
from math import pi

plt.rc('font', size=16)

def abs2(x):
    return x.real**2+x.imag**2

def normalize(a):
    return a/max(a)

def plot_reciprocal_bz(pathcolor):
    #plot reciprocal lattice
    x,y,z = lat.rlat
    x1,x2,_ = x
    y1,y2,_ = y
    G = [     0,      0, 0. ]
    M = [   0.0,    0.5, 0. ]
    M = [   0.5,    0.0, 0. ]
    K1 = [ 1./3.,  1./3., 0. ]
    K2 = [ 2./3., -1./3., 0. ]
    K3 = [ 1./3., -2./3., 0. ]
    K4 = [-1./3., -1./3., 0. ]
    K5 = [-2./3.,  1./3., 0. ]
    K6 = [-1./3.,  2./3., 0. ]
    K7 = [ 1./3.,  1./3., 0. ]
    ax = plt.gca()
    m = red_car([M], lat.rlat)[0]
    k = red_car([K1],lat.rlat)[0]
    ax.text(m[0], m[1], 'M',      color=pathcolor, fontsize=15)
    ax.text(k[0], k[1], 'K',      color=pathcolor, fontsize=15)
    ax.text(-0.02, 0, '$\Gamma$', color=pathcolor, fontsize=15)
    path = np.array(red_car([G,M,K1,K2,K3,K4,K5,K6,K7,G],lat.rlat))

    plt.plot(path[:,0], path[:,1], pathcolor)

#parse options
parser = argparse.ArgumentParser(description='Plot raman in the BZ.')
parser.add_argument('-e' ,'--energy', default=1.96, help='Energy to make the plot')
args = parser.parse_args()
e = float(args.energy)

modes = ['aprime1','eprime']
#modes = ['aprime3']

#read the lattice
stime = time()
lat = YamboLatticeDB()
print lat.rlat
print "took %3.2lfs"%(time()-stime)

pmode = 'phase'
#plot the figure
colora = {'aprime1':'blue',
          'eprime':'red'}
colori = {'aprime1':'lightblue',
          'eprime':'salmon'}
colorp = {'aprime1':'darkblue',
          'eprime':'darkred'}
label = {'aprime1':'$\\alpha^{\\rm xx}_\\mathbf{k}$ A$^\prime_1$',
         'eprime': '$\\alpha^{\\rm xx}_\\mathbf{k}$ E$^\prime$'}

fig = plt.figure(figsize=(15,8))
d = 0.08
f = 0.10
markersize = 30

#plot chi
filename = "undisplaced/chi.npz"
chi_file = np.load(filename) 
energies = chi_file['energies']
chi      = chi_file['chi']
kpts     = chi_file['kpt']
kfilter  = [nk for nk,k in enumerate(kpts) if np.linalg.norm(k) < 0.15]
kpts     = kpts[kfilter]
ychi     = np.sum(chi.imag,axis=0)

for ne,e in enumerate(energies):

    if not (100 <= ne <= 250) : 
        continue
    print ne, e

    #select energy    
    chie = chi[kfilter]
    chie = chie[:,ne]
    chie = chie.imag
    chie = chie/max(chie) 

    #plot bz
    ax = plt.axes([0.02,0.1,.4,.8])
    plt.title('$E_L$ = %5.2lf eV'%e)
    c2b = ax.scatter(kpts[:,0],kpts[:,1],marker='H',lw=0,
                     s=markersize,c=chie,vmin=0,cmap='viridis',
                     rasterized=True)

    ax.set_aspect('equal')
    dim = 0.15
    plt.xlim([-dim,dim])
    plt.ylim([-dim,dim])
    plt.axis('off')
    plot_reciprocal_bz('w')

    #plot chi
    ax = plt.axes([.55,0.1,.4,.8])
    plt.plot(energies,ychi)    
    plt.axvline(e)
    plt.xlim(0,2.5)
    plt.xlabel('$E_L$ eV')

    #colorbar viridis
    cax = plt.axes([.43, .05, .02, .9])
    plt.colorbar(cax=cax, mappable=c2b, ticks=[0.0, 0.5, 1.0])
    cax.set_ylabel('Im($\\chi^{\\rm xx}_{\\mathbf{k}}$)')

    plt.savefig('chibz_e%03d.pdf'%ne,dpi=200)
    plt.clf()



import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from yambopy import *
import argparse
from time import time
from math import pi

plt.rc('font', size=18)

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

modes = ['aprime1','aprime3']

#read the lattice
stime = time()
lat = YamboLatticeDB()
print lat.rlat
print "took %3.2lfs"%(time()-stime)

pmode = 'phase'
#plot the figure
colora = {'aprime1':'blue',
          'aprime3':'red'}
colori = {'aprime1':'lightblue',
          'aprime3':'salmon'}
colorp = {'aprime1':'darkblue',
          'aprime3':'darkred'}
label = {'aprime1':'$\\alpha^{\\rm xx}_\\mathbf{k}$ A$^\prime_1$',
         'aprime3':'$\\alpha^{\\rm xx}_\\mathbf{k}$ A$^\prime_3$'}

fig = plt.figure(figsize=(25,8))
d = 0.08
f = 0.10
markersize = 25

#plot chi
filename = "undisplaced/chi.npz"
chi_file = np.load(filename) 
energies = chi_file['energies']
chi      = chi_file['chi']
ychi     = np.sum(chi.imag,axis=0)

cbars  = 0.10
width  = 0.25
height = 0.8

for ne,e in enumerate(energies):

    if not (100 <= ne <= 250) : 
        continue
    print ne, e

    raman_mode = {}
    for n,mode in enumerate(modes):
        #load the two files
        filenamem = "%s/m.npz"%mode
        filenamep = "%s/p.npz"%mode
        chip_file = np.load(filenamep)
        chim_file = np.load(filenamem)
        energies = chip_file['energies']
        chip     = chim_file['chi']
        chim     = chip_file['chi']
        kpts     = chip_file['kpt']
        kfilter  = [nk for nk,k in enumerate(kpts) if np.linalg.norm(k) < 0.15]
        kpts     = kpts[kfilter]

        raman = chip-chim

        #select only points close to main gamma
        ramanf = raman[kfilter]

        #select energy
        e_index = min(range(len(energies)), key=lambda i: abs(energies[i]-e))

        ramane = ramanf[:,e_index]
        raman_mode[mode] = {}
        raman_mode[mode]['raman_mode'] = raman
        raman_mode[mode]['raman']     = ramane
        raman_mode[mode]['intensity'] = np.absolute(ramane)
        raman_mode[mode]['phase']     = np.angle(ramane)
        raman_mode[mode]['kpts']      = kpts

    #get max intensity of the plots
    max_intensity = 0
    max_phase = 0
    print "\ngetting max..."
    for n,mode in enumerate(modes):
        print 'mode', mode
        intensities = raman_mode[mode]['intensity']
        max_i = np.max(intensities)
        max_intensity = max(max_i,max_intensity)
        print "max_intensity", max_intensity

        phase = raman_mode[mode]['raman']
        max_p = np.max(phase)
        max_phase = max(max_p,max_phase)
        print "max_phase", max_phase/abs(max_phase)

    #normalize
    for n,mode in enumerate(modes):
        #raman_mode[mode]['intensity'] /= 1
        raman_mode[mode]['intensity'] /= max_intensity
        raman_mode[mode]['raman']     /= max_phase*1j


    #bz plot
    print "\nplotting..."
    for n,mode in enumerate(modes):

        print 'mode', mode
        x = 0.05+n*.3
        y = 0.05
        print x, y
        ax = plt.axes([x,y,width,height])
        plt.title(label[mode])

        intensities = raman_mode[mode]['intensity']
        #phases = raman_mode[mode]['phase']
        phases = np.angle(raman_mode[mode]['raman'])
        kpts   = raman_mode[mode]['kpts']

        smap = mpl.cm.ScalarMappable(cmap='hsv')
        print 'max_intensity', max(intensities)
        c = smap.to_rgba(phases)
        c[:,3] = np.sqrt(intensities)
        cax = ax.scatter(kpts[:,0],kpts[:,1],marker='H',lw=0,
                         s=markersize,c=c,cmap='hsv',rasterized=True)
        ax.set_aspect('equal')
        dim = 0.15
        plt.xlim([-dim,dim])
        plt.ylim([-dim,dim])
        plot_reciprocal_bz('k')
        plt.axis('off')

    #plot raman graph
    for n,mode in enumerate(modes):
        raman = raman_mode[mode]['raman_mode']
        raman = np.sum(raman,axis=0)
        raman = raman.real**2+raman.imag**2

        ax = plt.axes([.65,0.1,.3,.8])
        plt.plot(energies,raman)
        plt.axvline(e)
        plt.xlim(0,2.5)
        plt.xlabel('$E_L$ (eV)')

    plt.savefig('ramanbz_e%03d.pdf'%ne,dpi=200)
    plt.clf()

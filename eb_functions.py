import numpy as np
import healpy as hp

def tqu2qeueqbub(T, Q, U, lmax, getEfamily=True, getBfamily=True):
    nside = hp.npix2nside(len(T))
  
    alms = hp.map2alm([T, Q, U], lmax=lmax, pol=True)
    alms0 = np.zeros(len(alms[0]), dtype='complex')
  
    if getEfamily:
        Efamily = hp.alm2map([alms[0], alms[1], alms0], nside=nside, lmax=lmax, pol=True, verbose=False)
        QE = Efamily[1]
        UE = Efamily[2]
    if getBfamily:
        Bfamily = hp.alm2map([alms[0], alms0, alms[2]], nside=nside, lmax=lmax, pol=True, verbose=False)
        QB = Bfamily[1]
        UB = Bfamily[2]
      
    if getEfamily and getBfamily:
        return QE, UE, QB, UB
    elif getEfamily:
        return QE, UE
    elif getBfamily:
        return QB, UB

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

    
    
import pyssht as ssht # ssht package required for fast computation of spherical harmonics

def sYlm(s, l, m, L, rotation=None):
    flm = np.zeros(L*L, dtype=complex)
    ind = ssht.elm2ind(l, m)
    flm[ind] = 1.0
    if rotation:
        flm = ssht.rotate_flms(flm, rotation[0], rotation[1], rotation[2], L)
    f = ssht.inverse(flm, L, Spin=s)
    return f

def glm(pm, n, l, m, nside, L, Ylm, mYlm):
    # pm is +1 for g_+ and -1 for g_-
    # n is healpix pixel number
    th, ph = hp.pix2ang(nside, n)
    thin = ssht.theta_to_index(th, L)
    phin = ssht.phi_to_index(ph, L)
    return 0.25 * (Ylm[thin, phin] * np.conjugate(mYlm) + pm * mYlm[thin, phin] * np.conjugate(Ylm))
    

def g(pm, n, lmax, nside, L, verbose=True):
    ll = 0
    G = 0 
    while ll <= lmax:
        mm = -1 * ll
        while mm <= ll:
            Ylm = sYlm(2, ll, mm, L)
            mYlm = sYlm(-2, ll, mm, L)
            G = G + glm(pm, n, ll, mm, nside, L, Ylm, mYlm)
            mm = mm + 1
        if verbose:
            print '%d/%d\r' % (ll**2 / 2, lmax**2/2),
        ll = ll + 1
    return G

def g(pm, n, lmax, nside, L, Ylms, mYlms, verbose=True):
    ll = 0
    G = 0 
    while ll <= lmax:
        mm = -1 * ll
        while mm <= ll:
            #Ylm = sYlm(2, ll, mm, L)
            #mYlm = sYlm(-2, ll, mm, L)
            Ylm = Ylms[ll][mm + ll]
            mYlm = mYlms[ll][mm + ll]
            G = G + glm(pm, n, ll, mm, nside, L, Ylm, mYlm)
            mm = mm + 1
        if verbose:
            print '%d/%d\r' % (ll**2 / 2, lmax**2/2),
        ll = ll + 1
    return G

def getYlms(L, lmax):
    Ylms = []
    mYlms = []
    ll = 0
    while ll <= lmax:
        Ylmrow = []
        mYlmrow = []
        mm = -1 * ll
        while mm <= ll:
            Ylmrow.append(sYlm(2, ll, mm, L))
            mYlmrow.append(sYlm(-2, ll, mm, L))
            mm = mm + 1
        Ylms.append(Ylmrow)
        mYlms.append(mYlmrow)
        ll = ll + 1
    return Ylms, mYlms

def ssht_to_healpix(ssht_array, L, nside):
    out = np.zeros(hp.nside2npix(nside))
    for i in range(hp.nside2npix(nside)):
        th, ph = hp.pix2ang(nside, i)
        thin = ssht.theta_to_index(th, L)
        phin = ssht.phi_to_index(ph, L)
        out[i] = ssht_array[thin, phin]
    return out

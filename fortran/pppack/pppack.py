import numpy as np

import pypppack
__pyp__ = pypppack.pypppack

def interv(x, t):
	li_nk = len(t)
	li_left, li_mflag = __pyp__.pyinterv(x, t, li_nk)
	return li_left, li_mflag

def bvalue(x, p, t, bcoef, jderiv=1):
	li_nk = len(t)
	li_n = len(bcoef)
        lr_val = __pyp__.pybvalue(x, jderiv, p, t, bcoef, li_nk, li_n)
	return lr_val

def bsplvb(x, p, t, left=None):
	li_nk = len(t)
        if left is None:
            left, mflag = interv(x, t)
        biatx = __pyp__.pybsplvb(x, left, t[p:], li_nk-p-1, p+1)
	return biatx

def bsplvd(x, p, nderiv, t, left=None):
	li_nk = len(t)
        li_k = p + 1
        if left is None:
            left, mflag = interv(x, t)
        nderiv
        dbiatx = __pyp__.pybsplvd(x, left, nderiv, t, li_k, nderiv+1)
#        dbiatx = __pyp__.pybsplvd(x, left, t, li_nk, p+1)
	return dbiatx

def splint(k, t, tau, g):
	li_dim = len(g)
	li_n = len(g[0])
	li_k = k
	li_p = k-1
	li_nk = li_n + li_p + 1
	li_s = li_n - li_k
	li_npoint = li_n

	lpr_t = t

	lpr_tau = tau
	lpr_gtau = np.zeros((li_npoint, li_dim))
	for (i,Lg) in enumerate(g):
		lpr_gtau[:,i] = Lg

	lpr_P = np.zeros((li_npoint, li_dim))
	lpr_P = __pyp__.pysplint(li_n, li_p, lpr_t, lpr_tau, lpr_gtau, li_nk, li_npoint, li_dim)
	return zip(*lpr_P)

def splint2d(kx, taux, ky, tauy, g):

    nx = len(taux) ; ny = len(tauy)
    nkx = nx + kx ; nky = ny + ky
#    print nx, kx, taux
#    print ny, ky, tauy
#    print g
#    print nx, ny, nkx, nky
    return __pyp__.pysplint2d(kx, taux, ky, tauy, g, nx, ny, nkx, nky)





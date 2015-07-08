import numpy as np
from pigasus.geometry.pppack import interv, bvalue, bsplvb, splint, splint2d, bsplvd

def test_interv():
        print "*******  test_interv *******"
        li_s = 5
        li_k = 3
        t = np.asarray([0.]*li_k + list(np.linspace(0.,1.,li_s+2)[1:-1]) + [1.]*li_k)
        x = 0.3
        left, mflag = interv(x, t)
        print left, mflag

def test_bvalue():
        print "*******  test_bvalue *******"
        li_s = 3
        li_k = 3
        t = np.asarray([0.]*li_k + list(np.linspace(0.,1.,li_s+2)[1:-1]) + [1.]*li_k)
#        bcoef = np.random.random(li_k+li_s)
        bcoef = np.ones(li_k+li_s)
        x = 0.3
        val = bvalue(x, li_k-1, t, bcoef, jderiv=0)
        print val

def test_bsplvb():
        print "*******  test_bsplvb *******"
        li_s = 3
        li_k = 3
        t = np.asarray([0.]*li_k + list(np.linspace(0.,1.,li_s+2)[1:-1]) + [1.]*li_k)
        x = 0.3
        biatx = bsplvb(x, li_k-1, t)
        print biatx

def test_bsplvd():
        print "*******  test_bsplvd *******"
        li_s = 3
        li_k = 3
        p = li_k - 1
        nderiv = 2
        t = np.asarray([0.]*li_k + list(np.linspace(0.,1.,li_s+2)[1:-1]) + [1.]*li_k)
        x = 0.3
        dbiatx = bsplvd(x, p, nderiv, t)
        print dbiatx

def test_splint():
        print "*******  test_splint *******"
	from scipy.interpolate import splprep, splev, splrep
	k = 3
	n = 5

	x = np.linspace(0.,1.,n)
	y = np.sin(2*np.pi * x)

	# ...
	# using splrep
	# ...
	tck = splrep(x,y,s=0, k=k-1)
	t = tck[0]
	print "=== tck ==="
	print tck[1][:-k]
	# ...

	# ...
	# using splint
	# ...
	li_n = n
	li_k = k
	li_p = k-1
	li_nk = li_n + li_p + 1
	li_s = li_n - li_k
	li_npoint = n
	li_dim = 2
	#lpr_t = np.asarray([0.]*li_k + list(np.linspace(0.,1.,li_s+2)[1:-1]) + [1.]*li_k)
	lpr_t = t

	lpr_tau = x
	lpr_P = splint(k, lpr_t, lpr_tau, [y])
	print "=== splint ==="
	print lpr_P
	# ...

def test_splint2d():
    print "*******  test_splint2d *******"
    kx = 3 ; ky = 3
    nx = 5 ; ny = 5

    ux = np.linspace(0.,1.,nx) ; uy = np.linspace(0.,1.,ny)
    x,y = np.meshgrid(ux,uy)
    z = np.sin(2*np.pi * x) * np.sin(2*np.pi * y)

    # ...
    # using splrep
    # ...
    #	tck = splrep(x,y,s=0, k=k-1)
    #	t = tck[0]
    #	print "=== tck ==="
    #	print tck[1][:-k]
    # ...

    # ...
    # using splint
    # ...
    P, tx, ty = splint2d(kx, ux, ky, uy, z)
    #	print "=== splint2d ==="
    print "P  = ", P
    print "tx = ", tx
    print "ty = ", ty
    print P.shape

    c = np.zeros(25)
    c = np.zeros((5,5))
#    P = P.reshape(25)
    tck = [tx, ty, c, kx-1, ky-1]
    from scipy.interpolate import bisplev
    xnew, ynew = bisplev(ux,uy, tck)

    import pylab as pl
    pl.plot(xnew, ynew); pl.show()
    # ...



#test_interv()
#test_bvalue()
#test_bsplvd()
###test_bsplvb() # not working
#test_splint()
test_splint2d()

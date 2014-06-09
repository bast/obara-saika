
from math import sqrt, pi, exp, gamma
from mpmath import gammainc
import sys

#-------------------------------------------------------------------------------

class XF:
    def __init__(self,
                 scale=1,
                 prefactors=[],
                 q=[0,0,0,0,0,0,0,0,0,0,0,0],
                 order=0):
        self.scale = scale
        self.prefactors = prefactors
        self.q = q
        self.order = order

#-------------------------------------------------------------------------------

def list_is_flat(l):
    for x in l:
        if isinstance(x, list):
            return False
    return True

#-------------------------------------------------------------------------------

def flatten_sub(l):
    return [item for sublist in l for item in sublist]

#-------------------------------------------------------------------------------

def flatten(l):
    l_out = l[:]
    for i in range(50):
        if list_is_flat(l_out):
            return l_out
        else:
            l_out = flatten_sub(l_out)

    sys.stderr.write('ERROR: max depth reached in flatten\n')
    sys.exit(1)

#-------------------------------------------------------------------------------

def find_fun_to_lower(q):

    l = []
    for i in range(4):
        l.append(q[i*3] + q[i*3 + 1] + q[i*3 + 2])

    # find function to lower
    # start with lowest angular momentum above s
    fun = -1
    kmax = max(l) + 1
    for i in range(4):
        k = l[i]
        if k > 0:
            if k < kmax:
                kmax = k
                fun = i

    if fun == -1:
        sys.stderr.write('ERROR in find_fun_to_lower\n')
        sys.exit(1)

    return fun

#-------------------------------------------------------------------------------

def find_component_to_lower(fun):

    for i, c in enumerate(fun):
        if c > 0:
            return i

    sys.stderr.write('ERROR in find_component_to_lower\n')
    sys.exit(1)

#-------------------------------------------------------------------------------

def apply_os(xf):

    if sum(xf.q) == 0:
        return [xf]

    fun = find_fun_to_lower(xf.q)
    component = find_component_to_lower([xf.q[fun*3], xf.q[fun*3 + 1], xf.q[fun*3 + 2]])

    if component == 0:
        i1 = [ 0, 1, 2, 3][fun]
        i2 = [ 4, 4, 5, 5][fun]
    if component == 1:
        i1 = [ 6, 7, 8, 9][fun]
        i2 = [10,10,11,11][fun]
    if component == 2:
        i1 = [12,13,14,15][fun]
        i2 = [16,16,17,17][fun]

    bra = [0, 1]
    ket = [2, 3]

    pre = []
    pre.append(i1)
    pre.append(i2)
    if fun in bra:
        pre.append(18)
        pre.append(20)
        pre.append(18)
        pre.append(20)
    else:
        pre.append(19)
        pre.append(21)
        pre.append(19)
        pre.append(21)
    pre.append(22)
    pre.append(22)

    a = fun
    if a in bra:
       bra.remove(a)
       b = bra.pop()
       c = ket.pop()
       d = ket.pop()
    else:
       ket.remove(a)
       b = ket.pop()
       c = bra.pop()
       d = bra.pop()

    xf_copy = []
    scale = xf.scale
    order = xf.order
    for term in range(8):
        xf_new = XF(scale=scale,
                    prefactors=xf.prefactors[:],
                    q=xf.q[:],
                    order=order)
        xf_copy.append(xf_new)
        xf_copy[term].q[fun*3 + component] -= 1

    for term in [1, 3, 5, 6, 7]:
        xf_copy[term].order += 1

    xf_copy[2].q[a*3 + component] -= 1
    xf_copy[3].q[a*3 + component] -= 1
    xf_copy[4].q[b*3 + component] -= 1
    xf_copy[5].q[b*3 + component] -= 1
    xf_copy[6].q[c*3 + component] -= 1
    xf_copy[7].q[d*3 + component] -= 1

    n = []
    n.append(1)
    n.append(1)
    n.append(xf.q[a*3 + component] - 1)
    n.append(xf.q[a*3 + component] - 1)
    n.append(xf.q[b*3 + component])
    n.append(xf.q[b*3 + component])
    n.append(xf.q[c*3 + component])
    n.append(xf.q[d*3 + component])

    xf_list = []
    for term in range(8):
        if n[term] > 0:
            if all(i >= 0 for i in xf_copy[term].q):
                if n[term] > 1:
                    xf_copy[term].scale *= n[term]
                xf_copy[term].prefactors.append(pre[term])
                xf_list.append(xf_copy[term])

    xf_final = []
    for xf in xf_list:
        if all(i == 0 for i in xf.q):
            xf_final.append([xf])
        else:
            xf_final.append(apply_os(xf))

    return flatten(xf_final)

#-------------------------------------------------------------------------------

def get_r12_squared(r1, r2):
    return (r1[0] - r2[0])**2.0 + (r1[1] - r2[1])**2.0 + (r1[2] - r2[2])**2.0

#-------------------------------------------------------------------------------

def get_k(z1, z2, r1, r2):
    r12 = get_r12_squared(r1, r2)
    f0 = z1 + z2
    if r12 > 0.0:
        f1 = -z1*z2*r12/f0
        f2 = exp(f1)
    else:
        f2 = 1.0
    return sqrt(2.0)*f2*pi**(5.0/4.0)/f0

#-------------------------------------------------------------------------------

def get_rho(za, zb, zc, zd):
    z = za + zb
    n = zc + zd
    return z*n/(z + n)

#-------------------------------------------------------------------------------

def get_bi_center(z1, z2, r1, r2):
    z = z1 + z2
    rx = (z1*r1[0] + z2*r2[0])/z
    ry = (z1*r1[1] + z2*r2[1])/z
    rz = (z1*r1[2] + z2*r2[2])/z
    return [rx, ry, rz]

#-------------------------------------------------------------------------------

def boys(n, x):
    if x > 0.0:
        f = 2.0*x**(n + 0.5)
        g = gamma(n + 0.5)
        gi = 1.0 - gammainc(n + 0.5, x, regularized=True)
        return g*gi/f
    else:
        return 1.0/(n*2 + 1)

#-------------------------------------------------------------------------------

def get_aux(za, zb, zc, zd, ra, rb, rc, rd):
    k1 = get_k(za, zb, ra, rb)
    k2 = get_k(zc, zd, rc, rd)
    return k1*k2/sqrt(za + zb + zc + zd)

#-------------------------------------------------------------------------------

def get_integral(za, zb, zc, zd, ra, rb, rc, rd, c):

    rp = get_bi_center(za, zb, ra, rb)
    rq = get_bi_center(zc, zd, rc, rd)
    rw = get_bi_center((za + zb), (zc + zd), rp, rq)
    t  = get_rho(za, zb, zc, zd)*get_r12_squared(rp, rq)
    s  = get_aux(za, zb, zc, zd, ra, rb, rc, rd)

    z = za + zb
    n = zc + zd
    rho = z*n/(z + n)

    prefac = []
    prefac.append(rp[0] - ra[0])
    prefac.append(rp[0] - rb[0])
    prefac.append(rq[0] - rc[0])
    prefac.append(rq[0] - rd[0])
    prefac.append(rw[0] - rp[0])
    prefac.append(rw[0] - rq[0])
    prefac.append(rp[1] - ra[1])
    prefac.append(rp[1] - rb[1])
    prefac.append(rq[1] - rc[1])
    prefac.append(rq[1] - rd[1])
    prefac.append(rw[1] - rp[1])
    prefac.append(rw[1] - rq[1])
    prefac.append(rp[2] - ra[2])
    prefac.append(rp[2] - rb[2])
    prefac.append(rq[2] - rc[2])
    prefac.append(rq[2] - rd[2])
    prefac.append(rw[2] - rp[2])
    prefac.append(rw[2] - rq[2])
    prefac.append(0.5/z)
    prefac.append(0.5/n)
    prefac.append(-0.5*rho/(z*z))
    prefac.append(-0.5*rho/(n*n))
    prefac.append(0.5/(z + n))

    fun = XF(q=c)
    expansion = apply_os(fun)
    integral = 0.0
    for i in range(sum(c) + 1):
        b = boys(i, t)*s
        for f in expansion:
            if f.order == i:
                g = 1.0
                for k in f.prefactors:
                    g *= prefac[k]
                integral += float(f.scale)*b*g
    return integral

#-------------------------------------------------------------------------------

def test_get_integral():

    za = 1.1
    zb = 1.2
    zc = 1.3
    zd = 1.4

    ra = [1.0, 0.0, 1.0]
    rb = [0.0, 1.0, 2.0]
    rc = [0.0, 0.0, 3.0]
    rd = [0.0, 0.0, 4.0]

    ref = 1.71817807954e-05
    integral = get_integral(za, zb, zc, zd, ra, rb, rc, rd, [2, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0])

    assert abs(integral - ref) < 1.0e-16

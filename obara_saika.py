
from math import sqrt, pi, exp, gamma
from mpmath import gammainc
import sys

#-------------------------------------------------------------------------------

class X2:
    def __init__(self,
                 scale=1,
                 prefactors=[],
                 q=[0,0,0,0,0,0],
                 kind='S',
                 order=0):
        self.scale = scale
        self.prefactors = prefactors
        self.q = q
        self.kind = kind
        self.order = order

#-------------------------------------------------------------------------------

class X4:
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

def find_fun_to_lower(q, n):

    l = []
    for i in range(n):
        l.append(q[i*3] + q[i*3 + 1] + q[i*3 + 2])

    # find function to lower
    # start with lowest angular momentum above s
    fun = -1
    kmax = max(l) + 1
    for i in range(n):
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

def apply_os4(x4):

    if sum(x4.q) == 0:
        return [x4]

    fun = find_fun_to_lower(x4.q, 4)
    component = find_component_to_lower([x4.q[fun*3], x4.q[fun*3 + 1], x4.q[fun*3 + 2]])

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

    x4_copy = []
    scale = x4.scale
    order = x4.order
    for term in range(8):
        x4_new = X4(scale=scale,
                    prefactors=x4.prefactors[:],
                    q=x4.q[:],
                    order=order)
        x4_copy.append(x4_new)
        x4_copy[term].q[fun*3 + component] -= 1

    for term in [1, 3, 5, 6, 7]:
        x4_copy[term].order += 1

    x4_copy[2].q[a*3 + component] -= 1
    x4_copy[3].q[a*3 + component] -= 1
    x4_copy[4].q[b*3 + component] -= 1
    x4_copy[5].q[b*3 + component] -= 1
    x4_copy[6].q[c*3 + component] -= 1
    x4_copy[7].q[d*3 + component] -= 1

    n = []
    n.append(1)
    n.append(1)
    n.append(x4.q[a*3 + component] - 1)
    n.append(x4.q[a*3 + component] - 1)
    n.append(x4.q[b*3 + component])
    n.append(x4.q[b*3 + component])
    n.append(x4.q[c*3 + component])
    n.append(x4.q[d*3 + component])

    x4_list = []
    for term in range(8):
        if n[term] > 0:
            if all(i >= 0 for i in x4_copy[term].q):
                if n[term] > 1:
                    x4_copy[term].scale *= n[term]
                x4_copy[term].prefactors.append(pre[term])
                x4_list.append(x4_copy[term])

    x4_final = []
    for x4 in x4_list:
        if all(i == 0 for i in x4.q):
            x4_final.append([x4])
        else:
            x4_final.append(apply_os4(x4))

    return flatten(x4_final)

#-------------------------------------------------------------------------------

def apply_os2(x, kind='S'):

    if sum(x.q) == 0:
        x.kind = kind
        return [x]

    fun = find_fun_to_lower(x.q, 2)
    component = find_component_to_lower([x.q[fun*3], x.q[fun*3 + 1], x.q[fun*3 + 2]])

    if component == 0:
        i1 = [ 0, 1][fun]
    if component == 1:
        i1 = [ 2, 3][fun]
    if component == 2:
        i1 = [ 4, 5][fun]

    if kind == 'S':
        pre = []
        pre.append(i1)
        pre.append(6)
        pre.append(7)

    if kind == 'T':
        pre = []
        pre.append(i1)
        pre.append(6)
        pre.append(7)
        pre.append(8)
        i2 = [ 9, 10][fun]
        pre.append(i2)

    if kind == 'N':
        if component == 0:
            i2 = 6
        if component == 1:
            i2 = 7
        if component == 2:
            i2 = 8
        pre = []
        pre.append(i1)
        pre.append(i2)
        pre.append(9)
        pre.append(10)
        pre.append(9)
        pre.append(10)

    l = [0, 1]
    a = fun
    l.remove(a)
    b = l[0]

    if kind == 'S':
        num_terms = 3
    elif kind == 'T':
        num_terms = 5
    elif kind == 'N':
        num_terms = 6
    else:
        sys.stderr.write('ERROR: unexpected kind\n')

    x_copy = []
    scale = x.scale
    order = x.order # not used for S and T
    for term in range(num_terms):
        x_new = X2(scale=scale,
                   prefactors=x.prefactors[:],
                   q=x.q[:],
                   kind=kind,
                   order=order)
        x_copy.append(x_new)

    if kind == 'N':
        for term in [1, 3, 5]:
            x_copy[term].order += 1

    if kind == 'T':
        x_copy[3].kind = 'S'
        x_copy[4].kind = 'S'

    if kind == 'S':
        x_copy[0].q[fun*3 + component] -= 1
        x_copy[1].q[fun*3 + component] -= 1
        x_copy[2].q[fun*3 + component] -= 1
        x_copy[1].q[a*3 + component] -= 1
        x_copy[2].q[b*3 + component] -= 1

    if kind == 'T':
        x_copy[0].q[fun*3 + component] -= 1
        x_copy[1].q[fun*3 + component] -= 1
        x_copy[2].q[fun*3 + component] -= 1
        # term 4 (x_copy[3]) is just a copy
        x_copy[4].q[fun*3 + component] -= 2
        x_copy[1].q[a*3 + component] -= 1
        x_copy[2].q[b*3 + component] -= 1

    if kind == 'N':
        x_copy[0].q[fun*3 + component] -= 1
        x_copy[2].q[fun*3 + component] -= 1
        x_copy[4].q[fun*3 + component] -= 1
        x_copy[2].q[a*3 + component] -= 1
        x_copy[4].q[b*3 + component] -= 1
        x_copy[1].q = x_copy[0].q[:]
        x_copy[3].q = x_copy[2].q[:]
        x_copy[5].q = x_copy[4].q[:]

    n = []
    n.append(1)
    if kind == 'N':
        n.append(1)
    n.append(x.q[a*3 + component] - 1)
    if kind == 'N':
        n.append(x.q[a*3 + component] - 1)
    n.append(x.q[b*3 + component])
    if kind == 'N':
        n.append(x.q[b*3 + component])

    if kind == 'T':
        n.append(1)
        n.append(x.q[a*3 + component] - 1)

    x_list = []
    for term in range(num_terms):
        if n[term] > 0:
            if all(i >= 0 for i in x_copy[term].q):
                if n[term] > 1:
                    x_copy[term].scale *= n[term]
                x_copy[term].prefactors.append(pre[term])
                x_list.append(x_copy[term])

    x_final = []
    for y in x_list:
        if all(i == 0 for i in y.q):
            x_final.append([y])
        else:
            x_final.append(apply_os2(y, kind=y.kind))

    return flatten(x_final)

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

def get_coulomb(za, zb, zc, zd, ra, rb, rc, rd, c):

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

    fun = X4(q=c)
    expansion = apply_os4(fun)
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

def test_get_coulomb():

    za = 1.1
    zb = 1.2
    zc = 1.3
    zd = 1.4

    ra = [1.0, 0.0, 1.0]
    rb = [0.0, 1.0, 2.0]
    rc = [0.0, 0.0, 3.0]
    rd = [0.0, 0.0, 4.0]

    ref = 1.71817807954e-05
    integral = get_coulomb(za, zb, zc, zd, ra, rb, rc, rd, [2, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0])

    assert abs(integral - ref) < 1.0e-16

#-------------------------------------------------------------------------------

def get_overlap(za, zb, ra, rb, c):

    rp = get_bi_center(za, zb, ra, rb)

    z = za + zb
    e = za*zb/(za + zb)
    f = (ra[0] - rb[0])**2 + (ra[1] - rb[1])**2 + (ra[2] - rb[2])**2
    aux = exp(-e*f)*(pi/z)**1.5

    prefac = []
    prefac.append(rp[0] - ra[0])
    prefac.append(rp[0] - rb[0])
    prefac.append(rp[1] - ra[1])
    prefac.append(rp[1] - rb[1])
    prefac.append(rp[2] - ra[2])
    prefac.append(rp[2] - rb[2])
    prefac.append(0.5/z)
    prefac.append(0.5/z)

    fun = X2(q=c)
    expansion = apply_os2(fun)
    integral = 0.0
    for f in expansion:
        g = 1.0
        for k in f.prefactors:
            g *= prefac[k]
        integral += float(f.scale)*aux*g
    return integral

#-------------------------------------------------------------------------------

def test_get_overlap():

    za = 1.8
    zb = 2.8
    ra = [0.0, 0.0, 0.0]
    rb = [0.5, 0.8, -0.2]

    integral = get_overlap(za, zb, ra, rb, [0, 0, 0, 0, 0, 0])
    assert abs(integral - 0.20373275913014607) < 1.0e-16

    integral = get_overlap(za, zb, ra, rb, [1, 0, 0, 0, 0, 0])
    assert abs(integral - 0.062005622343957505) < 1.0e-16

    integral = get_overlap(za, zb, ra, rb, [1, 1, 0, 1, 1, 0])
    assert abs(integral - -0.00043801221837779696) < 1.0e-16

    integral = get_overlap(za, zb, ra, rb, [2, 1, 0, 1, 1, 0])
    assert abs(integral - -0.0002385994651113168) < 1.0e-16

#-------------------------------------------------------------------------------

def get_kinetic(za, zb, ra, rb, c):

    rp = get_bi_center(za, zb, ra, rb)

    z = za + zb
    e = za*zb/(za + zb)
    ab = (ra[0] - rb[0])**2 + (ra[1] - rb[1])**2 + (ra[2] - rb[2])**2
    aux = get_overlap(za, zb, ra, rb, [0, 0, 0, 0, 0, 0])

    prefac = []
    prefac.append(rp[0] - ra[0])
    prefac.append(rp[0] - rb[0])
    prefac.append(rp[1] - ra[1])
    prefac.append(rp[1] - rb[1])
    prefac.append(rp[2] - ra[2])
    prefac.append(rp[2] - rb[2])
    prefac.append(0.5/z)
    prefac.append(0.5/z)
    prefac.append(2.0*e)
    prefac.append(-e/za)
    prefac.append(-e/zb)

    fun = X2(q=c)
    expansion = apply_os2(fun, kind='T')
    integral = 0.0
    for f in expansion:
        if f.kind == 'T':
            g = e*(3.0 - 2.0*e*ab)
            for k in f.prefactors:
                g *= prefac[k]
            integral += float(f.scale)*aux*g
        if f.kind == 'S':
            g = 1.0
            for k in f.prefactors:
                g *= prefac[k]
            integral += float(f.scale)*aux*g
    return integral

#-------------------------------------------------------------------------------

def test_get_kinetic():

    za = 1.8
    zb = 2.0
    ra = [0.0, 0.0, 0.0]
    rb = [0.5, 0.8, -0.2]

    integral = get_kinetic(za, zb, ra, rb, [0, 0, 0, 0, 0, 0])
    assert abs(integral - 0.3652714583525358) < 1.0e-16

    integral = get_kinetic(za, zb, ra, rb, [1, 0, 0, 0, 0, 0])
    assert abs(integral - 0.2514265587836556) < 1.0e-16

    integral = get_kinetic(za, zb, ra, rb, [2, 2, 2, 2, 2, 2])
    assert abs(integral - -7.40057384314e-05) < 1.0e-16

#-------------------------------------------------------------------------------

def get_nuclear(za, zb, ra, rb, rc, c):

    z = za + zb
    rp = get_bi_center(za, zb, ra, rb)
    pc = (rp[0] - rc[0])**2 + (rp[1] - rc[1])**2 + (rp[2] - rc[2])**2

    u = z*pc
    aux = -2.0*(z/pi)**0.5*get_overlap(za, zb, ra, rb, [0, 0, 0, 0, 0, 0])

    prefac = []
    prefac.append(rp[0] - ra[0])
    prefac.append(rp[0] - rb[0])
    prefac.append(rp[1] - ra[1])
    prefac.append(rp[1] - rb[1])
    prefac.append(rp[2] - ra[2])
    prefac.append(rp[2] - rb[2])
    prefac.append(-rp[0] + rc[0])
    prefac.append(-rp[1] + rc[1])
    prefac.append(-rp[2] + rc[2])
    prefac.append(0.5/z)
    prefac.append(-0.5/z)

    fun = X2(q=c)
    expansion = apply_os2(fun, kind='N')
    integral = 0.0
    for f in expansion:
        g = 1.0
        for k in f.prefactors:
            g *= prefac[k]
        integral += float(f.scale)*aux*g*boys(f.order, u)
    return integral

#-------------------------------------------------------------------------------

def test_get_nuclear():

    za = 1.8
    zb = 2.0
    ra = [0.0, 0.0, 0.0]
    rb = [0.5, 0.8, -0.2]
    rc = [0.5, 0.8, 0.2]

    integral = get_nuclear(za, zb, ra, rb, rc, [0, 0, 0, 0, 0, 0])
    assert abs(integral - -0.49742209545104593) < 1.0e-16

    integral = get_nuclear(za, zb, ra, rb, rc, [1, 0, 0, 0, 0, 0])
    assert abs(integral - -0.15987439458254471) < 1.0e-16

    integral = get_nuclear(za, zb, ra, rb, rc, [2, 2, 2, 0, 0, 0])
    assert abs(integral - -0.003801373531942607) < 1.0e-16

    integral = get_nuclear(za, zb, ra, rb, rc, [1, 1, 1, 1, 1, 1])
    assert abs(integral - 8.8415484347060993e-5) < 1.0e-16


from itertools import product

from shell import get_shell
from obara_saika import get_coulomb

za = 1.1
zb = 1.2
zc = 1.3
zd = 1.4

ra = [1.0, 0.0, 1.0]
rb = [0.0, 1.0, 2.0]
rc = [0.0, 0.0, 3.0]
rd = [0.0, 0.0, 4.0]

# get all integrals up to pppp
for p in product('01', repeat=4):
    for c in get_shell(int(p[0]), int(p[1]), int(p[2]), int(p[3])):
        if sum(c) > 0:
            print c, get_coulomb(za, zb, zc, zd, ra, rb, rc, rd, c)

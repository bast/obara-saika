
def get_ijk_list(m):
    l = []
    for a in range(1, m + 2):
        for b in range(1, a + 1):
            i = m + 1 - a
            j = a - b
            k = b - 1
            l.append([i, j, k])
    return l

#-------------------------------------------------------------------------------

def get_shell(a, b, c, d):
    components = []
    for p in get_ijk_list(a):
        for q in get_ijk_list(b):
            for r in get_ijk_list(c):
                for s in get_ijk_list(d):
                    components.append(p + q + r + s)
    return components

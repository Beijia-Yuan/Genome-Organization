import numpy as np
import random
from tqdm import tqdm
import matplotlib.pylab as plt
import math
from scipy.integrate import quad

# Caculate Rg, Rrel

def in_loop_interaction(b,loop_n):
    return (b**2-loop_n)/(loop_n*(loop_n-1))
def A_ld(m,l,n,b):
    A = (l*n + l*m + m*n)/(2*m*l*n)
    def ld(b):
        return l/(l*n + l*m + m*n)*(m*b)
    return A, ld
def B2_fun_1D(m,l,n,b):
    A, ld = A_ld(m,l,n,b)
    def B2(b):
        return 1/2/A + ld(b)**2
    return B2
def pair_correlation_intersect(b,l_min,l_max,ep_min,ep_max):
    if l_min < ep_min:
        m = ep_min - l_min
        l = l_max - ep_min
        n = ep_max - l_max
    else:
        n = l_min - ep_min
        l = ep_max - l_min
        m = l_max - ep_max
    def in_loop(i):
        if l_min < i <= l_max:
            return True
        return False
    def in_ep(i):
        if ep_min < i <= ep_max:
            return True
        return False
    A, ld = A_ld(m,l,n,b)
    B2 = B2_fun_1D(m,l,n,b)
    def Imm(b):
        return (B2(b)-m)/(m*(m-1))
    def Ill(b):
        return (B2(b)-l)/(l*(l-1))
    def Inn(b):
        return (b*b-2*ld(b)*b+B2(b)-n)/(n*(n-1))
    def Imn(b):
        return (-ld(b)*b+B2(b))/(m*n)
    def Iml(b):
        return (-B2(b))/(m*l)
    def Iln(b):
        return (ld(b)*b-B2(b))/(n*l)
    def locate(i,j,b):
        if i==j:
            return 1
        if (not (in_loop(i) or in_ep(i))):
            return 0
        if (not (in_loop(j) or in_ep(j))):
            return 0
        if in_loop(i) and (not in_ep(i)):
            if in_loop(j) and (not in_ep(j)):
                return Imm(b)
            elif in_loop(j) and in_ep(j):
                return Iml(b)
            else:
                return Imn(b)
        elif in_loop(i) and in_ep(i):
            if in_loop(j) and (not in_ep(j)):
                return Iml(b)
            elif in_loop(j) and in_ep(j):
                return Ill(b)
            else:
                return Iln(b)
        else:
            if in_loop(j) and (not in_ep(j)):
                return Imn(b)
            elif in_loop(j) and in_ep(j):
                return Iln(b)
            else:
                return Inn(b)
    return locate, m*l/(m+l)+n
def pair_correlation_loop_in(b,l_min,l_max,ep_min,ep_max): 
    m = l_min - ep_min
    l = l_max - l_min
    n = ep_max - l_max
    def in_loop(i):
        if l_min < i <= l_max:
            return True
        return False
    def in_ep_only(i):
        if ep_min < i <= l_min:
            return True
        elif l_max < i <= ep_max:
            return True
        return False 

    def Imn(b):
        Imn = in_loop_interaction(b,m+n)
    Ill = in_loop_interaction(0,l)
    
    def locate(i,j,b):
        if i==j:
            return 1
        if (in_loop(i) and in_loop(j)):
            return Ill
        elif in_ep_only(i) and in_ep_only(j):
            return Imn(b)
        else:
            return 0
    return locate, m+n
def pair_correlation_loop_out(b,l_min,l_max,ep_min,ep_max): 
    m = ep_min - l_min
    l = ep_max - ep_min
    n = l_max - ep_max
    def in_ep(i):
        if ep_min < i <= ep_max:
            return True
        return False
    def in_loop_only(i):
        if l_min < i <= ep_min:
            return True
        elif ep_max < i <= l_max:
            return True
        return False 

    def Imn(b):
        return in_loop_interaction(b,m+n)
    def Ill(b):
        return in_loop_interaction(b,l)
    def Iml(b):
        return -b**2/l/(m+n)
    
    def locate(i,j,b):
        if i==j:
            return 1
        if in_loop_only(i):
            if in_loop_only(j):
                return Imn(b)
            elif in_ep(j):
                return Iml(b)
            else:
                return 0
        elif in_ep(i):
            if in_loop_only(j):
                return Iml(b)
            elif in_ep(j):
                return Ill(b)
            else:
                return 0
        else:
            return 0
    return locate, l*(m+n)/(l+m+n)
def pair_correlation_independent(b,l_min,l_max,ep_min,ep_max):
    l = l_max-l_min
    n = ep_max-ep_min
    def in_ep(i):
        if ep_min < i <= ep_max:
            return True
        return False
    def in_loop(i):
        if l_min < i <= l_max:
            return True
        return False 
    Ill = in_loop_interaction(0,l)
    def Inn(b):
        return in_loop_interaction(b,n)
    def locate(i,j,b):
        if i==j:
            return 1
        if in_loop(i) and in_loop(j):
            return Ill
        elif in_ep(i) and in_ep(j):
            return Inn(b)
        else:
            return 0
    return locate, n

def pair_correlation_no_loop(b,ep_min,ep_max):
    n = ep_max-ep_min
    def in_ep(i):
        if ep_min < i <= ep_max:
            return True
        return False
    def Inn(b):
        return in_loop_interaction(b,n)
    def locate(i,j,b):
        if i==j:
            return 1
        elif in_ep(i) and in_ep(j):
            return Inn(b)
        else:
            return 0
    return locate, n

def configure(b,l_min,l_max,ep_min,ep_max):
    if (l_max <= ep_min) or (l_min >= ep_max):
        return pair_correlation_independent(b,l_min,l_max,ep_min,ep_max)
    elif (l_max <= ep_max) and (l_min >= ep_min):
        return pair_correlation_loop_in(b,l_min,l_max,ep_min,ep_max)
    elif (l_max >= ep_max) and (l_min <= ep_min):
        return pair_correlation_loop_out(b,l_min,l_max,ep_min,ep_max)
    else:
        return pair_correlation_intersect(b,l_min,l_max,ep_min,ep_max)
def Rg(N,RiRj,b):
    NRg2 = 0
    for i in range(1,N+1):
        for j in range(1,N+1):
            NRg2 += (N-max(i-1,j-1)+1-(N-i)*(N-j)/N)*RiRj(i,j,b)
    return (NRg2/N)**0.5

# Calculate p

def Probden(b,ep_min,ep_max):
    RiRj, n_eff = configure(b,l_min,l_max,ep_min,ep_max)
    E = b**2/2/n_eff + v*(N**2)/(Rg(N,RiRj,b)**3)
    return np.exp(-E)
def Probden_no_loop(b,ep_min,ep_max):
    RiRj, n_eff = pair_correlation_no_loop(b,ep_min,ep_max)
    E = b**2/2/n_eff + v*(N**2)/(Rg(N,RiRj,b)**3)
    return np.exp(-E)
def Prob(Probden,ep_min,ep_max,N,v):
    N = quad(Probden, -N, N, args=(ep_min,ep_max))
    I = quad(Probden, -1, 1, args=(ep_min,ep_max))
    return I[0]/N[0]

# pairwise interactions

if __name__ == "__main__":
    N=558; v=0.0015; l_min=230; l_max=328
    w = open("HiC_simulation_no_loop.txt","w")
    for i in range(130,430,5):
        for j in range(i,430,5):
            if i == j:
                p = 1
            else:
                p = Prob(Probden_no_loop,i,j,N,v)
            w.write(str(i) + "\t" + str(j) + "\t" + str(p) +"\n")
    w.close()
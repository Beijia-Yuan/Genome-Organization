#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import random
from scipy.integrate import quad
import matplotlib.pylab as plt
from tqdm import tqdm 
import multiprocess as mp
from itertools import repeat

# In[250]:


def scalar_prod(c,x):
    try:
        iter(x)
        return(np.array(tuple(c*i for i in x)))
    
    except TypeError:
        return c*x


# In[199]:


def dot_prod(x,y):
    try:
        iter(x)
        if len(x) != len(y):
            print("error_dot_product")
            return
        return (sum(x[i]*y[i] for i in range(len(x))))
    except TypeError:
        return x*y


# In[234]:


def A_ld(m,l,n,x,y):
    A = (l*n + l*m + m*n)/(2*m*l*n)
    ld = scalar_prod(l/(l*n + l*m + m*n),scalar_prod(n,x)+scalar_prod(m,y))
    return A, ld


# In[206]:


def B2_fun(A,ld):
    return 3/2/A + dot_prod(ld,ld)


# In[256]:


def pair_correlation(m,l,n,x,y):
    A, ld = A_ld(m,l,n,x,y)
    B2 = B2_fun(A, ld)
    Imm = (dot_prod(x,x) - 2*dot_prod(ld,x) + B2 - m)/(m*(m-1))
    Ill = (B2-l)/(l*(l-1))
    Inn = (dot_prod(y,y) - 2*dot_prod(ld,y) + B2 - n)/(n*(n-1))
    Imn = (dot_prod(x,y) - dot_prod(ld,x) - dot_prod(ld,y) + B2)/(m*n)
    Iml = (dot_prod(ld,x)-B2)/(m*l)
    Inl = (dot_prod(ld,y)-B2)/(n*l)
    return Imm, Ill, Inn, Imn, Iml, Inl


# In[227]:


def generateData(distribution):
    import numpy as np
    def dot_prod(x,y):
        try:
            iter(x)
            if len(x) != len(y):
                print("error_dot_product")
                return
            return (sum(x[i]*y[i] for i in range(len(x))))
        except TypeError:
            return x*y
    def generate_step(steps):
        nxt = np.random.randint(6)
        if nxt == 0:
            steps.append(steps[-1]+np.array((1,0,0)))
        elif nxt == 1:
            steps.append(steps[-1]+np.array((-1,0,0)))
        elif nxt == 2:
            steps.append(steps[-1]+np.array((0,1,0)))
        elif nxt == 3:
            steps.append(steps[-1]+np.array((0,-1,0)))
        elif nxt == 4:
            steps.append(steps[-1]+np.array((0,0,1)))
        elif nxt == 5:
            steps.append(steps[-1]+np.array((0,0,-1)))
    steps = [np.array((0,0,0))]
    for j in range(1000):
        generate_step(steps)
    a = steps[300]
    b = steps[500]-steps[300]
    c = steps[1000]-steps[500]
    x = tuple(a+b)
    y = tuple(b+c)
    print(x,y)
    Imme = []
    for j in range(300):
        print("hi")
        sj = steps[j+1]-steps[j]
        for k in range(j+1,300):
            sk = steps[k+1]-steps[k]
            Imme.append(dot_prod(sj,sk))
        print("-", endl="")
    x = tuple(a+b)
    y = tuple(b+c)
    print(x,y)
    if (x,y) in distribution:
        distribution[(x,y)].append(np.mean(Imme))
    else:
        distribution[(x,y)] = [np.mean(Imme)]

def init_pool_processes():
    from numpy.random import seed
    seed()

def new_fun(d):
    import numpy as np
    a = np.random.randint(6)
    d[a] = a
    print(str(a)+"hi"+str(d))
# In[231]:

<<<<<<< HEAD
if __name__ == "__main__":
    with mp.Manager() as manager:
        distribution = manager.dict()
        with manager.Pool(initializer=init_pool_processes) as pool:
            pool.imap(generateData, repeat(distribution,10))
            print(distribution)
=======

with mp.Manager() as manager:
    import numpy as np
    distribution = manager.dict()
    with manager.Pool(16) as pool:
        tqdm(pool.imap_unordered(generateData, repeat(distribution, 1_000)))
>>>>>>> a7ee8fe64f3189e281dda8ec2ae9d796858e2137

# In[232]:



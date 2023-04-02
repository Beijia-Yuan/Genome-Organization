#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import random
from scipy.integrate import quad
import matplotlib.pylab as plt


# In[2]:


def scalar_prod(c,x):
    try:
        iter(x)
        return(tuple(c*i for i in x))
    
    except TypeError:
        return c*x

def dot_prod(x,y):
    try:
        iter(x)
        if len(x) != len(y):
            print("error_dot_product")
            return
        return (sum(x[i]*y[i] for i in range(len(x))))
    except TypeError:
        return x*y


# In[4]:


def A_ld(m,l,n,x,y):
    A = (l*n + l*m + m*n)/(2*m*l*n)
    ld = l * (scalar_prod(n,x)+scalar_prod(m,y))/(l*n + l*m + m*n)
    return A, ld

def B2_fun(A,ld):
    return 3/2/A + dot_prod(ld,ld)


# In[5]:


def pair_correlation(m,l,n,x,y):
    A, ld = A_ld(m,l,n,x,y)
    B2 = B2_fun(A, ld)
    Imm = (dot_prod(x-ld,x-ld) + 3/2/A - m)/(m*(m-1))
    Ill = (B2-l)/(l*(l-1))
    Inn = (dot_prod(y-ld,y-ld) + 3/2/A - n)/(n*(n-1))
    Imn = (dot_prod(x,y) - dot_prod(ld,(x+y)) + B2)/(m*n)
    Iml = (dot_prod(ld,x)-B2)/(m*l)
    Inl = (dot_prod(ld,y)-B2)/(n*l)
    return Imm, Ill, Inn, Imn, Iml, Inl


# In[6]:


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


# In[ ]:


m = 100
distribution = {}
for i in range(2000000):
    steps = [np.array((0,0,0))]
    for j in range(m):
        generate_step(steps)
    a = steps[30]
    b = steps[50]-steps[30]
    c = steps[100]-steps[50]
    Imme = []
    for j in range(30):
        sj = steps[j+1]-steps[j]
        for k in range(j+1,30):
            sk = steps[k+1]-steps[k]
            Imme.append(dot_prod(sj,sk))
    x = tuple(a+b)
    y = tuple(b+c)
    if (x,y) in distribution:
        distribution[(x,y)] += Imme
    else:
        distribution[(x,y)] = Imme


# In[ ]:


stat = {}
for eed in distribution:
    if len(distribution[eed]) > 10*(30+1)*30/2:
        stat[eed] = np.mean(distribution[eed])


# In[ ]:


distribution4 = []
for eed in stat:
    x = eed[0]; y = eed[1]
    mm = pair_correlation(30,20,50,x,y)[0]
    distribution4.append((stat[eed], mm))


# In[ ]:


plt.scatter(*zip(*distribution4),s=5)
plt.savefig("correlation_3D_mm.png")


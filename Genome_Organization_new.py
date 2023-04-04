#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import random
from scipy.integrate import quad
import matplotlib.pylab as plt
from tqdm import tqdm


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


# In[231]:


m = 100
distribution = {}
for i in tqdm(range(100_000_000)):
    steps = [np.array((0,0,0))]
    for j in range(m):
        generate_step(steps)
    a = steps[30]
    b = steps[50]-steps[30]
    c = steps[99]-steps[50]
    Imme = []
    for j in range(30):
        sj = steps[j+1]-steps[j]
        for k in range(j+1,30):
            sk = steps[k+1]-steps[k]
            Imme.append(dot_prod(sj,sk))
    x = tuple(a+b)
    y = tuple(b+c)
    if (x,y) in distribution:
        distribution[(x,y)].append(np.mean(Imme))
    else:
        distribution[(x,y)] = [np.mean(Imme)]


# In[232]:


stat = {}
for eed in distribution:
    if len(distribution[eed]) > 10:
        stat[eed] = np.mean(distribution[eed])


# In[258]:


distribution4 = []
with open("distribution4.csv", "w") as w:
    for eed in stat:
        x = eed[0]; y = eed[1]
        A, ld = A_ld(30,20,50,x,y)
        mm = pair_correlation(30,20,50,x,y)[0]
        distribution4.append((stat[eed], mm))
        w.write(str(x) + "," str(y) + "," + str(stat[eed]) + "," + str(mm) + "\n")


# In[259]:


try:
    x,y = zip(*distribution4)
except Exception:
    print("Not enough samples, subsitute with dummy samples")
    x, y = [1,2], [1,2]
plt.scatter(x=x, y=y,s=5)
# plt.show()
plt.savefig("correlation_3D_mm.png")


# In[ ]:


a,b = np.polyfit(x, y , 1)
print(a,b)


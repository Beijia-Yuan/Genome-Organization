import concurrent.futures
import random
import numpy as np
import random
from scipy.integrate import quad
import matplotlib.pylab as plt
from tqdm import tqdm
import time

def scalar_prod(c,x):
    try:
        iter(x)
        return(np.array(tuple(c*i for i in x)))
    
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

def A_ld(m,l,n,x,y):
    A = (l*n + l*m + m*n)/(2*m*l*n)
    ld = scalar_prod(l/(l*n + l*m + m*n),scalar_prod(n,x)+scalar_prod(m,y))
    return A, ld

def B2_fun(A,ld):
    return 3/2/A + dot_prod(ld,ld)

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


def generate_step(steps):
    nxt = random.randint(0,5)
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

def gen(_):
    m = 1000
    for _ in range(10):
        steps = [np.array((0,0,0))]
        for j in range(m):
            generate_step(steps)
        a = steps[300]
        b = steps[500]-steps[300]
        c = steps[1000]-steps[500]
        Imme = []
        for j in range(300):
            sj = steps[j+1]-steps[j]
            for k in range(j+1,300):
                sk = steps[k+1]-steps[k]
                Imme.append(dot_prod(sj,sk))
        x = tuple(a+b)
        y = tuple(b+c)
        if (x,y) in distribution:
            distribution[(x,y)].append(np.mean(Imme))
        else:
            distribution[(x,y)] = [np.mean(Imme)]

def process_distribution():
    for eed in distribution:
        if len(distribution[eed]) > 10:
            stat[eed] = np.mean(distribution[eed])
    distribution4 = []
    if len(stat) == 0:
        print("stat empty")
    with open("distribution4.csv", "w") as w:
        for eed in tqdm(stat):
            x = eed[0]; y = eed[1]
            A, ld = A_ld(300,200,500,x,y)
            mm = pair_correlation(300,200,500,x,y)[0]
            distribution4.append((stat[eed], mm))
            w.write(str(x) + "," + str(y) + "," + str(stat[eed]) + "," + str(mm) + "\n")
    try:
        x,y = zip(*distribution4)
    except Exception:
        print("Not enough samples, subsitute with dummy samples")
        x, y = [1,2,3], [1,2,3]
    return stat

def run(f, my_iter):
	l = len(my_iter)
	with tqdm(total=l) as pbar:
	    # let's give it some more threads:
	    with concurrent.futures.ThreadPoolExecutor(max_workers=32) as executor:
	        futures = {executor.submit(f, arg): arg for arg in my_iter}
	        results = {}
	        for future in concurrent.futures.as_completed(futures):
	            arg = futures[future]
	            results[arg] = future.result()
	            pbar.update(1)
	print("finished threading run")

if __name__ == "__main__":
    distribution = {}
    stat = {}
    run(gen, range(100000))
    process_distribution()
    print("Successfully executed")
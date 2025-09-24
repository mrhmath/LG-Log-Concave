import multiprocessing as mp
from sympy import * #load symbolic library
from collections import Counter
import pandas as pd
import numpy as np
from functools import lru_cache

t0 = Symbol('t0')
t1 = Symbol('t1')
z = Symbol('z')
p = Symbol('p')
q = Symbol('q')

#%%% Positivize polynomial

def evalneg(f):
    return expand(f.subs([(t0,-t0),(t1,-t1)]))

def is_positive(f):
    f=Poly(f)
    if not (f-f.abs()).is_zero:
        return 0
    return 1

def make_positive(f):
    f0 = f
    if is_positive(f) == 0:
        f = evalneg(f)
        if is_positive(f) == 0:
            f = -f
            if is_positive(f) == 0:
                f = evalneg(f)
                if is_positive(f) == 0:
                    print(f0,'is not positive')
                    return 0
    return f

#%%% Normalization

def homog(f):
    return Poly(Poly(f).homogenize(z),(t0,t1,z)) #(t0,t1,z) assigns the lexography

def multifactorial(v): #v is a vector
    mf = 1
    for d in v:
        mf *= factorial(d)
    return mf

def normalize(f):
    pf = homog(f)
    monodegs = pf.monoms()
    nf = 0
    for m in monodegs:
        nf+= pf.coeff_monomial(m)*multifactorial(m)**-1*t0**m[0]*t1**m[1]*z**m[2]
    return Poly(nf,(t0,t1,z))

#%%% Set of all length n unordered combinations of (t0,t1,z)

def partition(n, d, depth=0):
    if d == depth:
        return [[]]
    return [
        item + [i]
        for i in range(n+1)
        for item in partition(n-i, d, depth=depth+1)
        ]

@lru_cache(maxsize=None)
def d2cube(n):
    lst = [[n-sum(p)] + p for p in partition(n, 3-1)]
    total = []
    for item in lst:
        term = []
        for i in range(1,item[0]+1):
            term.append(t0)
        for i in range(1,item[1]+1):
            term.append(t1)
        for i in range(1,item[2]+1):
            term.append(z)
        total.append(term)
    return total

#%%% Count positive eigenvalues of Hessian

def countpos(a):
    ct = 0
    for x in a:
        if x>0:
            ct += 1
    return ct

def Hesscount(f):
    Hess = Matrix([[(f.as_expr()).diff(x).diff(y) for x in [t0,t1,z]] for y in [t0,t1,z]])
    npHess = matrix2numpy(Hess,dtype=float)
    eigens = np.linalg.eigvals(npHess)
    return np.sum(eigens > 0)

#%%% Lorentzian test

def shiftvars(poly): #f is a Laurant polynomial in t0 and t1, we shift f to polynomial in t0 and t1
    mpoly = Poly(poly.subs([(t0,p),(t1,q)]).as_expr())
    ppoly = ((mpoly * p**degree(mpoly,p**-1) * q**degree(mpoly,q**-1)).as_expr()).subs([(p,t0),(q,t1)])
    return ppoly

def Lorentz(f,name): #Test if polynomial f is Lorenztian, constraint on exactly one positive eigenvalue
    fpos = make_positive(f)
    if fpos == 0:
        print(name,'is not positive')
        return 0
    fnorm = normalize(fpos)
    deg = fnorm.total_degree()
    if deg <=1:
        print(name,'is low degree')
        return 0
    if deg == 2:
        if Hesscount(fnorm) > 1:
            print("[] is not a good derivative of", f)
            return 0
        return 1
    cube = d2cube(deg-2)
    for v in cube:
        fdeg2 = fnorm.diff(*v)
        if not fdeg2.is_quadratic:
            print("degree 2 error")
            return 0
        if Hesscount(fdeg2) > 1:
            print(v,"is not a good derivative of", name)
            return 0
    return 1

Lorentzlist = []
Notlist = []
for file in ['A03','A04','A05','A06','A07','A08','A09','A10','A11','A12','A13','A14','A15','A16']:
    T = pd.read_csv('data/A_whole/'+file+'.csv')
    for index, row in T.iterrows():
        name = "".join(T.loc[index,'key'].split()[2:5])
        poly = Poly(sympify(T.loc[index,'newExpr'])).as_expr()
        ppoly = shiftvars(poly)
        if Lorentz(ppoly,name):
            print(name,'is Lorentzian')
            Lorentzlist.append(name)
        else:
            Notlist.append(name)
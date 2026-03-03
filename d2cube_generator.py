import sympy as sp
import pickle
import os

d2cube_directory = "./d2cubes/"

def create_fold(fpath):
    if not os.path.exists(fpath):
        os.makedirs(fpath)
        print(f"Created folder: {fpath}")
    else:
        print(f"Folder already exists: {fpath}")

t0, t1 = sp.symbols("t0 t1")
z, p, q = sp.symbols("z p q")

def disk_load(fname):
    if os.path.exists(fname):
        with open(fname, "rb") as f:
            return pickle.load(f)
    return None

def disk_save(fname, obj):
    with open(fname, "wb") as f:
        pickle.dump(obj, f)

def partition(n, d, depth=0):
    if d == depth:
        return [[]]
    return [
        item + [i]
        for i in range(n+1)
        for item in partition(n-i, d, depth=depth+1)
        ]

def d2cube(n):
    fname = os.path.join(d2cube_directory, f"d2cube_{n}.pkl")
    output = disk_load(fname)
    if output is not None:
        return output
    lst = [[n-sum(p)] + p for p in partition(n, 3-1)]
    output = []
    for item in lst:
        term = []
        for i in range(1,item[0]+1):
            term.append(t0)
        for i in range(1,item[1]+1):
            term.append(t1)
        for i in range(1,item[2]+1):
            term.append(z)
        output.append(term)
    disk_save(fname, output)
    return output

if __name__ == "__main__":
    create_fold(d2cube_directory)
    for n in range(40):
        d2cube(n)

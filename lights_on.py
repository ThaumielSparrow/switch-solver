from operator import add
from itertools import chain, combinations
from functools import reduce

import math
import numpy as np
from scipy import ndimage
from tkinter import *


class GF2(object):  
    def __init__(self, a=0):
        self.value = int(a) & 1
    
    def __add__(self, rhs):
        return GF2(self.value + GF2(rhs).value)
    
    def __mul__(self, rhs):
        return GF2(self.value * GF2(rhs).value)
    
    def __sub__(self, rhs):
        return GF2(self.value - GF2(rhs).value)
    
    def __truediv__(self, rhs):
        return GF2(self.value / GF2(rhs).value)
    
    def __repr__(self):
        return str(self.value)
    
    def __eq__(self, rhs):
        if isinstance(rhs, GF2):
            return self.value == rhs.value
        return self.value == rhs
    
    def __le__(self, rhs):
        if isinstance(rhs, GF2):
            return self.value <= rhs.value
        return self.value <= rhs
    
    def __lt__(self, rhs):
        if isinstance(rhs, GF2):
            return self.value < rhs.value
        return self.value < rhs
    
    def __int__(self):
        return self.value
    
    def __long__(self):
        return self.value
    

GF2array = np.vectorize(GF2)


def gjel(A):
    nulldim = 0
    for i, row1 in enumerate(A):
        pivot = A[i:, i].argmax() + i
        if A[pivot, i] == 0:
            nulldim = len(A) - i
            break
        new_row = A[pivot] / A[pivot, i]
        A[pivot] = A[i]
        row1[:] = new_row
        
        for j, row2 in enumerate(A):
            if j == i:
                continue
            row2[:] -= new_row*A[j, i]
    return A, nulldim


def GF2inv(A):
    n = len(A)
    assert n == A.shape[1], "Matrix must be square"
    
    A = np.hstack([A, np.eye(n)])
    B, nulldim = gjel(GF2array(A))
    
    inverse = np.int_(B[-n:, -n:])
    E = B[:n, :n]
    null_vectors = []
    if nulldim > 0:
        null_vectors = E[:, -nulldim:]
        null_vectors[-nulldim:, :] = GF2array(np.eye(nulldim))
        null_vectors = np.int_(null_vectors.T)
    
    return inverse, null_vectors


def lightsoutbase(n):
    a = np.eye(n*n)
    a = np.reshape(a, (n*n, n, n))
    a = np.array(list(map(ndimage.binary_dilation, a)))
    return np.reshape(a, (n*n, n*n))


def powerset(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))


class LightsOut(object):
    def __init__(self, size=5):
        self.n = size
        self.base = lightsoutbase(self.n)
        self.invbase, self.null_vectors = GF2inv(self.base)
    
    def solve(self, b):
        b = np.asarray(b)
        assert b.shape[0] == b.shape[1] == self.n, "incompatible shape"
        
        if not self.issolvable(b):
            raise ValueError("The given setup is not solvable")
        
        first = np.dot(self.invbase, b.ravel()) & 1
        
        solutions = [(first + reduce(add, nvs, 0)) & 1 for nvs in powerset(self.null_vectors)]
        final = min(solutions, key=lambda x: x.sum())
        return np.reshape(final, (self.n, self.n))
    
    def issolvable(self, b):
        b = np.asarray(b)
        assert b.shape[0] == b.shape[1] == self.n, "incompatible shape"
        b = b.ravel()
        p = [np.dot(x, b) & 1 for x in self.null_vectors]
        return not any(p)


def text_to_mat(gridtxt, invert=True):
    gridlist = [int(s) for s in list(gridtxt)]

    shape = np.sqrt(len(gridlist))
    if shape%1 != 0:
        print("input matrix is not square.")
        return 1
    
    shape = int(shape)

    matlist = [gridlist[i: i+shape] for i in range(0, len(gridlist), shape)]
    mat = np.array(matlist)
    if invert: 
        mat = 1-mat
    
    return mat


def text_solver(gridtxt):
    global bsol
    mat_inv = text_to_mat(gridtxt, True)
    if type(mat_inv) == int:
        return 1
    lo = LightsOut(3)
    try:
        bsol = lo.solve(mat_inv)
    except:
        print("Error in determining solution")
        return 1
    
    return bsol


master = Tk()
master.title("it's really just simple linear algebra")
master.geometry("250x125")

check_size = 25
check_on = PhotoImage(width=check_size, height=check_size)
check_off = PhotoImage(width=check_size, height=check_size)
check_on.put(("green"), to=(0,0,check_size,check_size))
check_off.put(("red"), to=(0,0,check_size,check_size))

master_gridtxt = "000000000"
def update_gridtxt():
    master_gridtxt = ""
    for i in range(9):
        s = str(globals()[f"b_state{i}"].get())
        master_gridtxt += s


def reset_boxes():
    for i in range(9):
        globals()[f"b_state{i}"].set(0)


for i in range(9):
    j = i+1
    col = i%3
    row = math.ceil(j/3)
    globals()[f"b_state{i}"] = IntVar()
    globals()[f"b{i}"] = Checkbutton(master, variable=globals()[f"b_state{i}"], 
                                    image=check_off, selectimage=check_on, indicatoron=False,
                                    onvalue=1, offvalue=0, command=update_gridtxt)
    globals()[f"b{i}"].grid(row=row, column=col, padx=1, pady=1)


b_solve = Button(master, text="Solve", command=text_solver(master_gridtxt))
b_solve.grid(row=1, column=4, padx=1, pady=1)

b_reset = Button(master, text="Reset", command=reset_boxes)
b_reset.grid(row=2, column=4, padx=1, pady=1)

master.mainloop()

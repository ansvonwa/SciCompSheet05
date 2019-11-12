#!/usr/bin/env python3.6

import scipy
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def gen_mesh(n):
    nodes = range(0, (n+1)*(n+1))
    coordinates = list(map(lambda i: [i // (n+1), i % (n+1)], nodes))
    def corners(i):
        row = i // (2*n)
        lower_left = i % (2*n) // 2 + (n+1)*row
        if (i % 2 == 0): # even
            return [lower_left, lower_left+1, lower_left + n+2]
        else: # odd
            return [lower_left, lower_left + n+2, lower_left + n+1]
    elements = list(map(corners, range(0, 2*n*n)))
    dirichletboundary = set(range(0, n+1)) | set(range(n+1, n*(n+1), n+1)) | set(range(2*n+1, (n+1)*(n+1), n+1)) | set(range(n*(n+1), (n+1)*(n+1), 1))
    return (nodes, coordinates, elements, dirichletboundary)

#print(gen_mesh(3))

def assemble_stiffness_local(a0, a1, a2):
    42

def assemble_mass_local():
    42

def foo(coordinates, elements, dirichletboundary):
    # ...
    42



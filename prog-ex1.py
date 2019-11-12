#!/usr/bin/env python3.6

import scipy
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def gen_mesh(n):
    nodes = range(0, (n+1)*(n+1))
    coordinates = list(map(lambda i: [(i // (n+1))/n, (i % (n+1))/n], nodes))
    def corners(i):
        row = i // (2*n)
        lower_left = i % (2*n) // 2 + (n+1)*row
        if (i % 2 == 0): # even
            return [lower_left, lower_left+1, lower_left + n+2]
        else: # odd
            return [lower_left, lower_left + n+2, lower_left + n+1]
    elements = list(map(corners, range(0, 2*n*n)))
    dirichletboundary = set(range(0, n+1)) | set(range(n+1, n*(n+1), n+1)) | \
            set(range(2*n+1, (n+1)*(n+1), n+1)) | set(range(n*(n+1), (n+1)*(n+1), 1))
    return (nodes, coordinates, elements, dirichletboundary)

#print(gen_mesh(3))
nodes, coordinates, elements, dirichletboundary = gen_mesh(3)

def assemble_stiffness_local(triangle):
    g1 = np.array([[1] + coordinates[c] for c in elements[triangle]])
    g2 = np.array([[0, 0], [1, 0], [0, 1]])
    g = np.linalg.inv(np.transpose(g1)).dot(g2)
    x = g1[:,1]
    y = g1[:,2]
    a = 1/2 * abs((x[1]-x[0])*(y[2]-y[0]) - (y[1]-y[0])*(x[2]-x[0])) * \
            g.dot(np.transpose(g))
    return a

print(assemble_stiffness_local(0))

def assemble_mass_local(triangle):
    g1 = np.array([[1] + coordinates[c] for c in elements[triangle]])
    x = g1[:,1]
    y = g1[:,2]
    m = 1/24 * abs((x[1]-x[0])*(y[2]-y[0]) - (y[1]-y[0])*(x[2]-x[0])) * \
            (np.ones((3,3)) + np.diag(np.ones(3)))
    return m

print(assemble_mass_local(0))

def assemble_load_local(triangle):
    # ??
    43

def foo(coordinates, elements, dirichletboundary):
    # ...
    42



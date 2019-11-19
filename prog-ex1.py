#!/usr/bin/env python3.6

import scipy
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

def gen_mesh(n):
    nodes = range(0, (n+1)*(n+1))
    #This is written as compact as possible. Hopefully its still understandable
    coordinates = list(map(lambda i: [(i // (n+1))/n, (i % (n+1))/n], nodes))

    #Function to list all the corners for a given triangle i. Always start bottom left and count counterclockwise
    def corners(i):
        row = i // (2*n)
        lower_left = i % (2*n) // 2 + (n+1)*row
        if (i % 2 == 0): # even
            return [lower_left, lower_left+1, lower_left + n+2]
        else: # odd
            return [lower_left, lower_left + n+2, lower_left + n+1]

    elements = list(map(corners, range(0, 2*n*n)))

    #Dirichletboundary can be broken into 4 parts. Top, bottom, left and right boundary. Merge those together.
    dirichletboundary = set(range(0, n+1)) | set(range(n+1, n*(n+1), n+1)) | \
            set(range(2*n+1, (n+1)*(n+1), n+1)) | set(range(n*(n+1), (n+1)*(n+1), 1))
    return (nodes, coordinates, elements, dirichletboundary)




#Choose the n already here!!!
n = 25
#Choose the n already here!!!



#Those variables are needed for the functions below
nodes, coordinates, elements, dirichletboundary = gen_mesh(n)

#These functions follow the line of the exercise sheet. No magic happening here
def assemble_stiffness_local(triangle):
    g1 = np.array([[1] + coordinates[c] for c in elements[triangle]])
    g2 = np.array([[0, 0], [1, 0], [0, 1]])
    g = np.linalg.inv(np.transpose(g1)).dot(g2)
    x = g1[:,1]
    y = g1[:,2]
    a = abs((x[1]-x[0])*(y[2]-y[0]) - (y[1]-y[0])*(x[2]-x[0])) * \
            g.dot(np.transpose(g))
    return a


def assemble_mass_local(triangle):
    g1 = np.array([[1] + coordinates[c] for c in elements[triangle]])
    x = g1[:,1]
    y = g1[:,2]
    m = 1/24 * abs((x[1]-x[0])*(y[2]-y[0]) - (y[1]-y[0])*(x[2]-x[0])) * \
            (np.ones((3,3)) + np.diag(np.ones(3)))
    return m


def assemble_load_local(triangle, f):
    [xT, yT] = map(lambda v: v/3, map(sum, zip(*[coordinates[c] for c in elements[triangle]])))
    return f(xT, yT) * 1/(6*n*n) * np.ones(3)


#use lil_matrix for highest efficiency
def assemble_stiffness():
    num_nodes = (n+1)*(n+1)
    a = lil_matrix((num_nodes, num_nodes))
    for t in range(0, 2*n*n):
        loc = assemble_stiffness_local(t)
        for c1 in range(0,3):
            for c2 in range(0,3):
                a[elements[t][c1], elements[t][c2]] += loc[c1][c2]
    return a.tocsr()

def assemble_mass():
    num_nodes = (n+1)*(n+1)
    m = lil_matrix((num_nodes, num_nodes))
    for t in range(0, 2*n*n):
        loc = assemble_mass_local(t)
        for c1 in range(0,3):
            for c2 in range(0,3):
                m[elements[t][c1], elements[t][c2]] += loc[c1][c2]
    return m.tocsr()

def assemble_load(f):
    v = np.zeros((n+1)*(n+1))
    for t in range(0, 2*n*n):
        loc = assemble_load_local(t, f)
        for c1 in range(0,3):
            v[elements[t][c1]] += loc[c1]
    return v

#Solve the linear system
def solve(f):
    non_zero = list(filter(lambda n: not n in dirichletboundary, nodes))
    a = assemble_stiffness()
    v = assemble_load(f)
    res = spsolve(a[non_zero,:][:,non_zero], v[non_zero])
    full_res = np.zeros((n+1)*(n+1))
    full_res[non_zero] = res
    return full_res


def L2diff(f, coeff):
    #Use (f - v)^2 = f^2 + v^2 -2*f*v. Setting v = sum v_i * phi_i
    num_nodes = (n+1)*(n+1)
    mass = assemble_mass()
    v = assemble_load(f)

    #The L2 norm of the analytic solution u is 1/4
    temp = 0.25
    for i in range(num_nodes):
        temp -= 2*coeff[i]*v[i]
        for j in range(num_nodes):
            temp += coeff[i]*coeff[j]*mass[i, j]
    return temp


#Functions from the sheet
f = lambda x,y: 2*math.pi*math.pi*np.sin(math.pi*x)*np.sin(math.pi*y)
u = lambda x,y: np.sin(math.pi*x)*np.sin(math.pi*y)
u_n = lambda x,y: u(x/n, y/n)

#get the x and y coordinates of the nodes
coordinatesnp = np.array(coordinates)
x = coordinatesnp[:,0]
y = coordinatesnp[:,1]

#First the analytic solution:
z = np.fromfunction(u_n, (n+1, n+1), dtype = int)
z = z.flatten()

plot1 = plt.figure(1)
plt.triplot(x, y, elements)
plt.tricontourf(x, y, elements, z)
plot1.show()

#Second the numerical solution:
z = solve(f)

plot2 = plt.figure(2)
plt.triplot(x, y, elements)
plt.tricontourf(x, y, elements, z)
plot2.show()


print("L2 difference is: " + str(L2diff(u, z)))



"""
Due to the design of our program it is not simply possible to automatically repeat that process with different n.
So after setting n = 2, 4, 8, 16, 32, 64 manually I get the L2 distances:

0.12621363230682697
0.024064878810094007
0.004794883104977338
0.0011042694547959426
0.00026983720968328786
6.706458894656713e-05

Plotting these numbers on a log-log graph yields a fairly straight line. So our convergence rate is of the form N^k.
Looking at the ratio of the smallest and biggest numbers tells us that k = -4.5 approximately.

Use the code below to see the plot of the L2 errors:

index = np.array([2, 4, 8, 16, 32, 64])
norms = np.array([0.12621363230682697, 0.024064878810094007, 0.004794883104977338, 0.0011042694547959426, 0.00026983720968328786, 6.706458894656713e-05])

plt.plot(index, norms)
plt.xscale('log')
plt.yscale('log')
plt.show()
"""















# -*- coding: utf-8 -*-
# Tímadæmi í Varmaflutningsfræði

from numpy import *
from scipy import *
import math
from matplotlib import *
import matplotlib
from mpl_toolkits.mplot3d import Axes3D

# Breytur

T1 = 40
T2 = 10
Tinf = 0
h = 75
q_dot = 1000000
k = 20
B = 0.05
H = 0.1

# Reikninet
dx = 0.0025
dy = dx
M = B/dx + 1
N = H/dy + 1

# Skilgreining stærðar fylkja og vigra
A = zeros([M*N,M*N])
B = zeros([M+M*(N-1),1])

for i in range(1,int(M)+1):      # Kannski vandamál ??????????
    for j in range(1,int(N)+1):
        
        s=i+M*(j-1)-1
        s_iplus1=i+1+M*(j-1)-1
        s_iminus1=i-1+M*(j-1)-1
        s_jplus1=i+M*(j+1-1)-1
        s_jminus1=i+M*(j-1-1)-1
        
        # Vinstri hlið - Varmaburður
        if i==1 and j>1 and j<N:
            A[s,s] = 2+(h*dx)/k
            A[s,s_jplus1] = -.5
            A[s,s_jminus1] = -.5
            A[s,s_iplus1] = -1
            B[s] = (dx**2/2/k)*q_dot+h*Tinf*dx/k

        # Hægri hlið - Einangruð
        if i==M and j>1 and j<N:
            A[s,s] = 4
            A[s,s_jplus1] = -1
            A[s,s_jminus1] = -1
            A[s,s_iminus1] = -2
            B[s] = (dx**2/2/k)*q_dot
        
        # Innri punktar
        if i>1 and i<M and j>1 and j<N:
            A[s,s] = 4
            A[s,s_iplus1] = -1
            A[s,s_iminus1] = -1
            A[s,s_jplus1] = -1
            A[s,s_jminus1] = -1
            B[s] = (dx**2/k)*q_dot
            
        # Botn - T1 = fasti
        if j==1:
            A[s,s] = 1
            B[s] = T1
          
        # Toppur - T2 = fasti
        if j==N:
            A[s,s] = 1
            B[s] = T2
        
# Lausn fundin
T_list = linalg.solve(A,B)

# Röðun hitastigs á réttan stað
T_grid = zeros([N,M])

for i in range(1,int(M)+1):
    for j in range(1,int(N)+1):
        s = i+M*(j-1)-1
        T_grid[j-1,i-1] = T_list[int(s)]

# Teikna niðurstöður


x = linspace(0,dx*(M-1),M)
y = linspace(0,dy*(N-1),N)
X,Y = meshgrid(x,y)
Z = T_grid

# Teikna jafnhæðalínur fyrir hita
figure
pyplot.contour(X,Y,Z)
pyplot.colorbar()
pyplot.grid()
#pyplot.xlim(min(x),max(x))
#pyplot.ylim(min(y),max(y))
pyplot.xlabel('x [m]')
pyplot.ylabel('y [m]')
pyplot.axis('equal')
pyplot.show()

# Teiknar þrívíddaryfirborð
fig = pyplot.figure()
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X,Y,Z,rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0.1, antialiased=False)
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_zlabel('T [C]')
fig.colorbar(surf, shrink=0.5, aspect=5)
pyplot.show()
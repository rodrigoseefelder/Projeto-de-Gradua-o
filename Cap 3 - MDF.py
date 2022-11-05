# -*- coding: utf-8 -*-
"""
Created on Sat Jul 23 18:37:02 2022

@author: rodri
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

def gerador2d(nx,ny,e='tri'):
    x=[]
    for j in range(0,ny):
        for k in range(0,nx):
            x.append(k/(nx-1))
    y=[]
    for j in range(0,ny):
        for k in range(0,nx):
            y.append(j/(ny-1))
    if e=='quad':
        add=0
        IENQ=[]
        for j in range(0,ny-1):
            for k in range(0,nx-1):
                IENQ.append([add+k,add+k+1,add+k+nx+1,add+k+nx])
            add+=nx
            
        IENQ=np.array(IENQ)
        return[x,y,IENQ]
    elif e=='tri':
        Tri = matplotlib.tri.Triangulation(x,y)
        IENT = Tri.triangles  
        IENT=np.array(IENT)
        return[x,y,IENT]
    else:
        print('Parametro "e" deve ser uma string ("quad" ou "tri")' )
        

#PLOTANDO SOLUCOES POR MDF PARA MALHAS VARIADAS

npoints_min = 5
npoints_max = 5

for i in range(npoints_min,npoints_max+1):

    nx = i
    ny = i
    npoints = nx*ny
    ne = (nx-1)*(ny-1)
    dx = 1.0/ne
    dy = 1.0/ne
    
    # geracao de malha
    X = np.array(gerador2d(nx,ny)[0])
    Y = np.array(gerador2d(nx,ny)[1])
    
    
    cc2=[]
    cc4=[]
    for i in range(0,npoints):
        if X[i]==0:
            cc4.append(i)
        if X[i]==1:
            cc2.append(i)
            
    cc = cc2+cc4
    
    cc1=[]
    cc3=[]
    for i in range(0,npoints):
        if i not in cc:
            if Y[i]==0:
                cc1.append(i)
            if Y[i]==1:
                cc3.append(i)
                
    cc+= cc1+cc3
    
    inner=[]
                
    for i in range(0,npoints):
        if i not in cc:
            inner.append(i)
    
    
    #--------------------------------------------------
    # plt.plot(X,Y,'bo')
    # plt.plot(X[cc],Y[cc],'ko')
    # plt.show()
    #-------------------------------------------------- 
    
    bc = 0*np.ones( (npoints),dtype='float' )
    for i in cc1:
     bc[i] = 0
    for i in cc2:
     bc[i] = 4*Y[i] - 4*Y[i]*Y[i]
    for i in cc3:
     bc[i] = 0
    for i in cc4:
     bc[i] = 0
    
    #construindo matrizes
    
    A = np.zeros( (npoints,npoints),dtype='float' )
    b = np.zeros( (npoints),dtype='float' )
    
    # Eq. da cond. de contorno de Dirichlet
    for i in cc:
     A[i,i] = 1.0
     b[i]   = bc[i]
    
    # Eqs. pontos internos (inner):
    for i in inner:
     A[i,i-1] = 1/(dx*dx)
     A[i,i] = -2/(dx*dx) -2/(dy*dy)
     A[i,i+1] = 1/(dx*dx)
     A[i,i-nx] = 1/(dy*dy)
     A[i,i+nx] = 1/(dy*dy)
    
    
    # solucao do sistema linear Ax=b
    #Ainv = np.linalg.inv(A)
    #T = Ainv@b
    
    T = np.linalg.solve(A,b)
    
    
    # plot 2D cor (quadrilatero)
    Z = T.reshape(ny,nx)
    surf = plt.imshow(Z, interpolation='quadric', origin='lower',
                      cmap=matplotlib.cm.jet, extent=(X.min(),
                      X.max(), Y.min(), Y.max()))
    plt.colorbar(surf,shrink=1.0, aspect=20)
    plt.grid(color='black', linestyle='solid', linewidth=0.5)
    labx = np.linspace(X.min(),X.max(),nx)
    laby = np.linspace(Y.min(),Y.max(),ny)
    plt.xticks(labx,"")
    plt.yticks(laby,"")
    plt.title("T(x,y)")
    plt.show()
    
    
print(T)
    
    

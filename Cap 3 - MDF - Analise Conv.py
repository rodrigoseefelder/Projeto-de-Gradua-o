# -*- coding: utf-8 -*-
"""
Created on Sat Jul 23 18:37:02 2022

@author: rodri
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import time

print("INÍCIO DA ANÁLISE DE CONVERGÊNCIA PELO MDF")
print("")

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
        
        
listanpontos = []
listanelementos = []
listaerromax = []
listaduracao = []
listaTmed = []

#PLOTANDO SOLUCOES POR MDF

for i in range(3,40):
    
    start = time.time()

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
    
    
    end = time.time()
    duracao=end-start
    listaduracao.append(duracao)
    listaTmed.append(np.average(T))
    
    
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
    
    print(duracao)
    
    #definindo funcoes auxiliares da solucao analitica
    
    nmax=101 #numero maximo de iteracoes do somatorio dentro de cn
    
    def a1(n):
        return 2/(np.sinh(n*np.pi))
    
    def a2(n):
        return 8 - 4*np.pi*n*np.sin(np.pi*n) - 8*np.cos(np.pi*n)
    
    def a3(n):
        return np.pi*np.pi*np.pi*n*n*n
    
    def cn(n):
        return a1(n)*(a2(n)/a3(n))
    
    
    def termosomatorio(n,x,y):
        return cn(n)*np.sinh(n*np.pi*x)*np.sin(n*np.pi*y)
    
    #solucao analitica: retorna o valor da temperatura para x e y definidos baseado no somatorio de 1 ate n = nmax
    
    def T_x_y(nmax,x,y):
        n=1
        T=0
        for i in range(0,nmax):
            T+=termosomatorio(n,x,y)
            n+=1
        return T
    
    
            
    #definindo coordenadas dos pontos nos eixos X e Y
    
    X = np.linspace(0,1,nx)
    Y = np.linspace(0,1,ny)
    
    T_analit=[]
    
    for y in Y:
        for x in X:
            T_analit.append(T_x_y(nmax,x,y))
    T_analit = np.array(T_analit)
    
    erro = abs(T-T_analit)
    erromax = max(erro)
    
    listanpontos.append(nx)
    listanelementos.append(ne)
    listaerromax.append(erromax)
    
    
    print(erromax)
    
plt.rcParams["font.family"] = "serif"
plt.plot(listanpontos,listaerromax,'ko',markersize=2.5)
plt.title("Nº de pontos em cada eixo x Erro absoluto máximo")
plt.grid(True)
plt.savefig('analiseconvmdf.png', dpi=300)
plt.show()

plt.plot(listanpontos,listaduracao,'ko',markersize=2.5)
plt.title("Nº de pontos em cada eixo x Tempo de Execução")
plt.grid(True)
plt.savefig('tempoexecmdf.png', dpi=300)
plt.show()

print("")
print("Erro Máximo por Iteração: ")
print(listaerromax)
print("")
print("")
print("Duração por Iteração: ")
print(listaduracao)
print("Temperatura Média por Iteração: ")
print(listaTmed)
print("")

print("FIM DA ANÁLISE DE CONVERGÊNCIA PELO MDF")
    
    

# -*- coding: utf-8 -*-
"""
Created on Sat Jul 23 18:37:02 2022

@author: rodri
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import time

print("INÍCIO DA ANÁLISE DE CONVERGÊNCIA PELO MEF")
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
        
def printacolorido(nx,ny,IEN,T,vmin,vmax,titulo=None):
    plt.rcParams["font.family"] = "serif"
    Z = T.reshape(ny,nx)
    ax = plt.triplot(X,Y,IEN,color='k',linewidth=0.3)
    surf = plt.imshow(Z, interpolation='quadric', origin='lower',
                      cmap=matplotlib.cm.jet, extent=(X.min(),
                      X.max(), Y.min(), Y.max()))
    cb=plt.colorbar(surf,shrink=1.0, aspect=20)
    plt.clim(0,1)
    labx = np.linspace(X.min(),X.max(),nx)
    laby = np.linspace(Y.min(),Y.max(),ny)
    plt.title(titulo,size=12)
    plt.xticks([])
    plt.yticks([])
    plt.show() 
        
        
listanpontos = []
listanelementos = []
listaerromax = []
listaduracao = []
listaTmed = []

#PLOTANDO SOLUCOES POR MDF

for i in range(3,40):

#parametros malha

    start = time.time()

    nx = i
    ny = i
    npoints = nx*ny
    ne = 2*(nx-1)*(ny-1)
    dx = 1.0/ne
    dy = 1.0/ne
    
    # geracao de malha
    X = np.array(gerador2d(nx,ny)[0])
    Y = np.array(gerador2d(nx,ny)[1])
    IEN = np.array(gerador2d(nx,ny)[2])
    
    # gerando listas com nós de contorno
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
            
    #assembly da matriz K
    K=np.zeros((npoints,npoints))
    for e in range(0,ne):
        [v1,v2,v3] = IEN[e]
        xi=X[v1]
        xj=X[v2]
        xk=X[v3]
        yi=Y[v1]
        yj=Y[v2]
        yk=Y[v3]
        area=(1/2)*np.linalg.det([[1,xi,yi],[1,xj,yj],[1,xk,yk]])
        ai=xj*yk-xk*yj
        aj=xk*yi-xi*yk
        ak=xi*yj-xj*yj
        bi=yj-yk
        bj=yk-yi
        bk=yi-yj
        ci=xk-xj
        cj=xi-xk
        ck=xj-xi
        melem = (area/12)*np.array([ [2,1,1],
                                     [1,2,1],
                                     [1,1,2]])
        kx = (1/(4*area))*np.array([ [bi*bi,bi*bj,bi*bk],
                                     [bj*bi,bj*bj,bj*bk],
                                     [bk*bi,bk*bj,bk*bk]])
        ky = (1/(4*area))*np.array([ [ci*ci,ci*cj,ci*ck],
                                     [cj*ci,cj*cj,cj*ck],
                                     [ck*ci,ck*cj,ck*ck]])
    
        kelem = kx+ky
    
        for ilocal in range(0,3):
            iglobal=IEN[e,ilocal]
            for jlocal in range(0,3):
                jglobal=IEN[e,jlocal]
                K[iglobal,jglobal]+=kelem[ilocal,jlocal] 
    
    
    # distribuicao de fonte de calor
    Q0 = 0*np.ones( (npoints),dtype='float' )
    for i in range(0,npoints):
      if X[i]<0.5: 
          Q0[i] = 2
    
    
    #definindo parâmetros para solução transiente
    
    
    b= np.zeros(npoints)
    bcc=np.zeros(npoints)
    
    
    for i in cc1:
     bcc[i] = 0
    for i in cc2:
     bcc[i] = 4*Y[i] - 4*Y[i]*Y[i] 
    for i in cc3:
     bcc[i] = 0
    for i in cc4:
     bcc[i] = 0
     
    # imposicao das c.c.s de Dirichlet
    for i in cc:
     K[i,:] = 0.0
     K[i,i] = 1.0
     b[i] = bcc[i]
     
     T=np.linalg.solve(K,b)
     
    
    
    end = time.time()
    duracao = end-start
    listaduracao.append(duracao)
    listaTmed.append(np.average(T))
    
    
    # plot 2D cor
    titulo = "T(x,y)"
    printacolorido(nx,ny,IEN,T,0,2,titulo)
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
plt.savefig('analiseconvmef.png', dpi=300)
plt.show()

plt.plot(listanpontos,listaduracao,'ko',markersize=2.5)
plt.title("Nº de pontos em cada eixo x Tempo de Execução")
plt.grid(True)
plt.savefig('tempoexecmef.png', dpi=300)
plt.show()


print("")
print("Erro Máximo por Iteração: ")
print(listaerromax)
print("")
print("Duração por Iteração: ")
print(listaduracao)
print("")
print("Temperatura Média por Iteração: ")
print(listaTmed)
print("")

print("FIM DA ANÁLISE DE CONVERGÊNCIA PELO MEF")
    
    

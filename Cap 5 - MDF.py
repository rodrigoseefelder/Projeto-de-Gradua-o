# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 20:36:24 2022

@author: rodri
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import time


col1 = "T_mef"
col2 = "T_mdf"


#PROBLEMA TÉRMICO 2D DIRICHLET elementos triangulares

k=1
cv=1
rho=1
kappa_x = 1.0
kappa_y = 1.0
alpha_x = kappa_x/(rho*cv)
alpha_y = kappa_y/(rho*cv)
alpha=alpha_x

#definindo funcoes 
    

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

#parametros malha

malhas = [25,41]

for i in malhas:
    
    print("Início da resolução para a malha com nx = ny = "+str(i))
    print("") 
    
    nx = i
    ny = i
    npoints = nx*ny
    
            
    def printacolorido(nx,ny,T,titulo=None):
        plt.rcParams["font.family"] = "serif"
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
        plt.title(titulo)
        plt.show() 
    
    #DEFININDO PARÂMETROS DA MALHA
    ne = (nx-1)*(ny-1)
    dx = 1.0/ne
    dy = 1.0/ne
    
    start_general = time.time()
    
    #GERAÇÃO DE MALHA
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
    
    
    # condicoes de contorno espaciais
    bc = 0*np.ones( (npoints),dtype='float' )
    for i in cc1:
     bc[i] = 0
     #bc[i] = 0
    for i in cc2:
     bc[i] = 4*Y[i] - 4*Y[i]*Y[i] 
     #bc[i] = 1
    for i in cc3:
     bc[i] = 0
     #bc[i] = 2
    for i in cc4:
     bc[i] = 0
     #bc[i] = 1
    
    # condicao inicial (temporal)
    T = 0*np.ones( (npoints),dtype='float' )
    for i in cc:
     T[i] = bc[i]
     
          
    q=2000
    # distribuicao de fonte de calor
    #Qvec = q/(rho*cv)*np.ones( (npoints),dtype='float' )
    Qvec = np.zeros( (npoints),dtype='float' )
    for i in range(0,npoints):
       if X[i]<0.5: 
           Qvec[i] = q/rho*cv
    
    
    #Temperatura inicial:
    
    printacolorido(nx,ny,T,'Temperatura (Condição Inicial)')
    
    #Resolvendo o problema transiente:
        
    dt = 1e-8
    tempo=0
    
    lista_normas_MDF=[]
    lista_Tmedio_MDF=[]
    lista_Tmedio_MDF.append(np.average(T))
    
    duracao_plot = 0
    
    print("Soluções MDF com dt = "+str(dt)+"\n")
    
    n_it=0
    #while normaTplus - normaT > 1E-9:
    while n_it < 60000002:
     for i in inner:
      T[i] = T[i] + (alpha*dt/(dx*dx) ) *( T[i+1] - 2*T[i] + T[i-1] ) + (alpha*dt/(dy*dy) ) *( T[i+nx] - 2*T[i] + T[i-nx] ) + (dt * Qvec[i]/(rho*cv))
      tempo+=dt
      n_it+=1
      start_plot = time.time()
      if n_it in [2000000, 4000000, 8000000, 12000000, 15000000, 18000000, 25000000, 50000000, 60000000]:
          tempo=round(tempo,2)
          titulo2 = 'Temperatura (t = {})'.format(str(tempo))
          printacolorido(nx,ny,T,titulo2)
          tmediomdf=np.average(T)
          lista_Tmedio_MDF.append(tmediomdf)
      end_plot = time.time()
      duracao_plot+=end_plot-start_plot
      
    end_general = time.time()
    
    duracao_geral = end_general - start_general
    
    duracao_geral_sem_plot = duracao_geral - duracao_plot
    
    print("")
    print(duracao_geral_sem_plot)
    print("")
    print(duracao_plot)
    print("")
    
    print("Fim da resolução para a malha com nx = ny = "+str(nx))
    print("") 
    
     
      
      
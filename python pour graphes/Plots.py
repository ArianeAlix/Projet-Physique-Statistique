import matplotlib.pyplot as plt
import numpy as np


def read_file (fn ):
    with open (fn) as f:
        names =f. readline ()
        X=np. array ([[ float (x) for x in l. strip (). split (" ")] for l in f. readlines ()])
    return X[: ,0] ,X[: ,1],X[: ,2],names.split (" ")
    
    
    
plt.figure(figsize=(6,6))    
r,pot_prim,f,names=read_file('./Derivee.txt')
plt.plot(r,pot_prim,'blue',label='-dV/dr')
plt.plot(r,f,'red',label='f')
plt.axhline(0, color='black',linewidth=0.5)
plt.legend()
plt.show()

plt.figure(figsize=(8,6))
Etot,Ep,Ec,names=read_file('./Energie.txt')
plt.plot(np.arange(0,np.size(Etot),1),Etot,'blue',label='Etot',linewidth=0.5)
plt.plot(np.arange(0,np.size(Etot),1),Ep,'red',label='Ep',linewidth=0.5)
plt.plot(np.arange(0,np.size(Etot),1),Ec,'green',label='Ec',linewidth=0.5)
plt.legend()
plt.show()


def positions():
    with open ('./Positions.txt') as f:
        names =f. readline ()
        X=np. array ([[ float (x) for x in l. strip (). split (" ")] for l in f. readlines ()])
        
        plt.figure(figsize=(8,6))
        for i in range (0,np.shape(X)[1]):
            plt.plot(X[:,i],np.arange(0,np.size(X[:,i]),1),'blue',label='p'+str(i+1),linewidth=0.5)
        plt.legend()
        plt.show()
        
        
positions()
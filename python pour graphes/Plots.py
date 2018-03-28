import matplotlib.pyplot as plt
import numpy as np


def read_file (fn ):
    with open (fn) as f:
        names =f. readline ()
        X=np. array ([[ float (x) for x in l. strip (). split (" ")] for l in f. readlines ()])
    return X[: ,0] ,X[: ,1],X[: ,2],names.split (" ")
    
    
    
plt.figure(figsize=(5,5))    
r,pot_prim,f,names=read_file('Derivee.txt')
plt.plot(r,pot_prim,'blue',label='-dV/dr')
plt.plot(r,f,'red',label='f')
plt.axhline(0, color='black',linewidth=0.5)
plt.legend()
plt.show()

plt.figure(figsize=(10,8))
Etot,Ep,Ec,names=read_file('Energie.txt')
plt.plot(np.arange(0,np.size(Etot),1),Etot,'blue',label='Etot',linewidth=0.5)
plt.plot(np.arange(0,np.size(Etot),1),Ep,'red',label='Ep',linewidth=0.5)
plt.plot(np.arange(0,np.size(Etot),1),Ec,'green',label='Ec',linewidth=0.5)
plt.legend()
plt.show()

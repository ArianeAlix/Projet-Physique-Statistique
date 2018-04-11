import matplotlib.pyplot as plt
import numpy as np


def read_file (fn ):
    with open (fn) as f:
        names =f. readline ()
        X=np. array ([[ float (x) for x in l. strip (). split (" ")] for l in f. readlines ()])
    return X[: ,0] ,X[: ,1],X[: ,2],names.split (" ")
    
    
#Tracé de - la dérivée du potentiel et de f pour vérifier la correspondance
plt.figure(figsize=(6,6))    
r,pot_prim,f,names=read_file('./build/Derivee.txt')
plt.plot(r,pot_prim,'blue',label='-dV/dr')
plt.plot(r,f,'red',label='f')
plt.axhline(0, color='black',linewidth=0.5)
plt.title("Force et opposé de la dérivée du potentiel")
plt.legend()
plt.show()

#Tracé des énergies cinétique, potentielle et totale
plt.figure(figsize=(8,6))
Etot,Ep,Ec,names=read_file('./build/Energie.txt')
plt.plot(np.arange(0,np.size(Etot),1),Etot,'blue',label='Etot',linewidth=0.5)
plt.plot(np.arange(0,np.size(Etot),1),Ep,'red',label='Ep',linewidth=0.5)
plt.plot(np.arange(0,np.size(Etot),1),Ec,'green',label='Ec',linewidth=0.5)
plt.title("Evolution de l'énergie du système")
plt.legend()
plt.show()


def positions():
    with open ('./build/Positions.txt') as f:
        names =f. readline ()
        X=np. array ([[ float (x) for x in l. strip (). split (" ")] for l in f. readlines ()])
        
        plt.figure(figsize=(8,6))
        for i in range (0,np.shape(X)[1]):
            plt.plot(X[:,i],np.arange(0,np.size(X[:,i]),1),'blue',label='p'+str(i+1),linewidth=0.5)
        plt.legend()
        plt.show()
        
        
positions()


#Tracé de l'erreur en focntion de log_10(deltaT)

def erreur():
    with open ('./build/Erreur_energie.txt') as f:
        names =f. readline ()
        X=np. array ([[ float (x) for x in l. strip (). split (" ")] for l in f. readlines ()])
    plt.figure(figsize=(8,6))
    plt.plot(X[:,0],X[:,1],'.',label='Erreur totale')
    
    fit = np.polyfit(X[:,0],X[:,1],1)
    fit_fn = np.poly1d(fit)
    
    print(fit_fn);
    
    plt.plot(X[:,0],fit_fn(X[:,0]),'-',label='Regression linéaire')
    plt.title("Pourcentage de l'erreur en fonction du carré du pas de temps")
    plt.legend()
    plt.show()
   
erreur()













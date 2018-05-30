#coding=utf-8

import matplotlib.pyplot as plt
import numpy as np


def read_file (fn ):
    with open (fn) as f:
        names =f. readline ()
        X=np. array ([[ float (x) for x in l. strip (). split (" ")] for l in f. readlines ()])
    return X[: ,0] ,X[: ,1],X[: ,2],names.split (" ")
    
    
def energies():
    #Tracé des énergies cinétique, potentielle et totale
    plt.figure(figsize=(8,6))
    Etot,Ep,Ec,names=read_file('./Energie_2d.txt')
    times=10**(-3)*np.arange(0,np.size(Etot),1)#car pas de simulation de 10^-4 et de releve de valeurs 10^-3
    plt.plot(times,Etot,'blue',label='Etot',linewidth=0.5)
    plt.plot(times,Ep,'red',label='Ep',linewidth=0.5)
    plt.plot(times,Ec,'green',label='Ec',linewidth=0.5)
    plt.title("Evolution de l'énergie du système")
    plt.legend()
    plt.show()

energies()


def temperatures():
    with open ('./Temp_2d.txt') as f:
        names =f. readline ()
        X=np. array ([[ float (x) for x in l. strip (). split (" ")] for l in f. readlines ()])
        for i in range (0,np.shape(X)[1]):
            plt.figure(figsize=(8,6))
            plt.plot(10**(-3)*np.arange(0,np.size(X[:,i]),1),X[:,i],i*'red'+(i-1)*'blue',label='T'+(1-i)*'c'+i*'p',linewidth=0.5)
            plt.legend()
            plt.title('Température '+(1-i)*'cinétique'+i*'potentielle')
            plt.show()

temperatures()


def positions():
    with open ('./Positions_2d.txt') as f:
        names =f. readline ()
        X=np. array ([[ float (x) for x in l. strip (). split (" ")] for l in f. readlines ()])
        
        plt.figure(figsize=(8,6))
        for i in range (0,int(np.shape(X)[1]/2)):
            plt.plot(X[:,2*i+1],X[:,2*i],label='p'+str(i+1),linewidth=0.5)
        #plt.legend()
        plt.show()
        
positions()



#Tracé de l'erreur en fonction de log_10(deltaT)

def erreur():
    with open ('./Erreur_energie_2d.txt') as f:
        names =f. readline ()
        X=np. array ([[ float (x) for x in l. strip (). split (" ")] for l in f. readlines ()])
    plt.figure(figsize=(8,6))
    plt.plot(X[:,0],X[:,1],'.',label='Emax-Emin')
    
    fit= np.polyfit(X[:,0],X[:,1],1,full=True)
    fit_fn = np.poly1d(fit[0])
    
    print(fit_fn);
    print("Erreur résiduelle: ",fit[1])
    
    plt.plot(X[:,0],fit_fn(X[:,0]),'-',label='Regression linéaire')
    plt.title("Amplitude des variations d'énergie en fonction du carré du pas de temps")
    plt.legend()
    plt.show()
   
erreur()










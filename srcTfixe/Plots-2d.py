#coding=utf-8

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.mlab as mlab

def read_file (fn ):
    with open (fn) as f:
        names =f. readline ()
        X=np. array ([[ float (x) for x in l. strip (). split (" ")] for l in f. readlines ()])
    return X[: ,0] ,X[: ,1],X[: ,2],names.split (" ")
    
    
def energies():
    #Tracé des énergies cinétique, potentielle et totale
    plt.figure(figsize=(8,6))
    Etot,Ep,Ec,names=read_file('./build/Energie_2d.txt')
    times=10**(-3)*np.arange(0,np.size(Etot),1)#car pas de simulation de 10^-4 et de releve de valeurs 10^-3
    plt.plot(times[100:],Etot[100:],'blue',label='Etot',linewidth=0.5)
    plt.plot(times[100:],Ep[100:],'red',label='Ep',linewidth=0.5)
    plt.plot(times[100:],Ec[100:],'green',label='Ec',linewidth=0.5)
    plt.title("Evolution de l'énergie du système")
    plt.legend()
    plt.show()

# energies()


def temperatures():
    with open ('./build/Temp_2d.txt') as f:
        names =f. readline ()
        X=np. array ([[ float (x) for x in l. strip (). split (" ")] for l in f. readlines ()])
        for i in range (0,np.shape(X)[1]):
            plt.figure(figsize=(8,6))
            plt.plot(10**(-3)*np.arange(0,np.size(X[:,i]),1),X[:,i],i*'red'+(i-1)*'blue',label='T'+(1-i)*'c'+i*'p',linewidth=0.5)
            plt.legend()
            plt.title('Température '+(1-i)*'cinétique'+i*'potentielle')
            plt.show()

temperatures()


def pression():
    plt.figure(figsize=(8,6))
    plt.title("Pression en fonction de la densité \n à température fixée")
    
    with open ('./build/Press_2d.txt') as f:
        names =f. readline ()
        X=np. array ([[ float (x) for x in l. strip (). split (" ")] for l in f. readlines ()])
        times=10**(-3)*np.arange(0,np.size(X),1)#car pas de simulation de 10^-4 et de releve de valeurs 10^-3
        
        plt.plot(X[:,0],X[:,1],'blue',linewidth=0.5)
        
        #Selon la formule PV=NkT
        T=1.2
        k=1
        plt.plot(np.arange(0,1.1,0.1),k*T*np.arange(0,1.1,0.1),'green',linewidth=0.5,label="PV=NkT")
        
        
        '''
        #Selon la formule P= (rho kT)/(1 - b rho) - a rho^2
        T=1.2
        k=1
        Na=6.02*10**23
        Tc=150.7/120
        Pc=48.98 
        a=(27 * (k**2) * (Tc**2))/(64*Pc)
        b= (k*Tc) /(8*Pc)
        print(Pc)
        print(a,b)
        plt.plot(np.arange(0,1,0.01),(k*T*np.arange(0,1,0.01))/(1-b*np.arange(0,1,0.01))-a*np.arange(0,1,0.01)**2,'green',linewidth=0.5,label="van der Waals")
        '''
        
        plt.title("Pression en fonction de la densité \n à température fixée")
        
        
    with open ('./build/Press_2d_ecart_type.txt') as f:
        X=np. array ([[ float (x) for x in l. strip (). split (" ")] for l in f. readlines ()])
        moyenne=np.mean(X[:,1])
        ecart_type=np.std(X[:,1])
        plt.plot([X[0,0],X[0,0]],[moyenne+ecart_type,moyenne-ecart_type],'red')
        
        
    with open ('./build/Press_2d_ecart_type1.txt') as f:
        X=np. array ([[ float (x) for x in l. strip (). split (" ")] for l in f. readlines ()])
        moyenne=np.mean(X[:,1])
        ecart_type=np.std(X[:,1])
        plt.plot([X[0,0],X[0,0]],[moyenne+ecart_type,moyenne-ecart_type],'red')
        
    
    with open ('./build/Press_2d_ecart_type2.txt') as f:
        X=np. array ([[ float (x) for x in l. strip (). split (" ")] for l in f. readlines ()])
        moyenne=np.mean(X[:,1])
        ecart_type=np.std(X[:,1])
        plt.plot([X[0,0],X[0,0]],[moyenne+ecart_type,moyenne-ecart_type],'red')
        
        
    plt.xlabel('Rho')
    plt.ylabel('P')
    plt.show()

pression()



def EP():
    plt.figure(figsize=(8,6))
    plt.title("Energie potentielle en fonction de la densité \n à température fixée")
    
    with open ('./build/Ep_2d.txt') as f:
        names =f. readline ()
        X=np. array ([[ float (x) for x in l. strip (). split (" ")] for l in f. readlines ()])
        times=10**(-3)*np.arange(0,np.size(X),1)#car pas de simulation de 10^-4 et de releve de valeurs 10^-3
        
        plt.plot(X[:,0],X[:,1],'blue',linewidth=0.5)
        
        
    with open ('./build/Ep_2d_ecart_type.txt') as f:
        X=np. array ([[ float (x) for x in l. strip (). split (" ")] for l in f. readlines ()])
        moyenne=np.mean(X[:,1])
        ecart_type=np.std(X[:,1])
        plt.plot([X[0,0],X[0,0]],[moyenne+ecart_type,moyenne-ecart_type],'red')
        
        
    with open ('./build/Ep_2d_ecart_type1.txt') as f:
        X=np. array ([[ float (x) for x in l. strip (). split (" ")] for l in f. readlines ()])
        moyenne=np.mean(X[:,1])
        ecart_type=np.std(X[:,1])
        plt.plot([X[0,0],X[0,0]],[moyenne+ecart_type,moyenne-ecart_type],'red')
        
    with open ('./build/Ep_2d_ecart_type2.txt') as f:
        X=np. array ([[ float (x) for x in l. strip (). split (" ")] for l in f. readlines ()])
        moyenne=np.mean(X[:,1])
        ecart_type=np.std(X[:,1])
        plt.plot([X[0,0],X[0,0]],[moyenne+ecart_type,moyenne-ecart_type],'red')
        
        
    plt.xlabel('Rho')
    plt.ylabel('Ep')
    plt.show()

        
EP()



def positions():
    with open ('./build/Positions_2d.txt') as f:
        names =f. readline ()
        X=np. array ([[ float (x) for x in l. strip (). split (" ")] for l in f. readlines ()])
        
        plt.figure(figsize=(8,6))
        for i in range (0,int(np.shape(X)[1]/2)):
            plt.plot(X[:,2*i+1],X[:,2*i],label='p'+str(i+1),linewidth=0.5)
        #plt.legend()
        plt.show()
        
# positions()


'''
#Tracé de l'erreur en fonction de deltaT

def erreur():
    with open ('./build/Erreur_energie_2d.txt') as f:
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

'''

def histos():
    with open ('./build/Histo_2d.txt') as f:
        names =f. readline ()
        X=np. array ([[ float (x) for x in l. strip (). split (" ")] for l in f. readlines ()])
        for i in range (0,np.shape(X)[1]):
            plt.figure(figsize=(8,6))
            plt.hist(X[:,i],bins=np.arange(-4,4,0.2),normed=True)
            plt.title('Histogramme des vitesses en '+(1-i)*'x'+i*'y')
            
            mu = 0
            k=1
            m=1
            T=1.2
            sigma = k*T/m
            x = np.linspace(mu - 3*sigma, mu + 3*sigma, 100)
            plt.plot(x,mlab.normpdf(x, mu, sigma),'red',label="N(0,1.44)")
            
            plt.legend()
            plt.show()

histos()







import numpy as np
import matplotlib.pyplot as plt
from math import exp,atan,cos,sin
import random
import pandas as pd
from statistics import median,mean,stdev
from copy import copy
from random import randrange 
from functions import generate_random,uniform,cosRule,tempScaling
import time



def acceptChange(cord,energy,new_energy,matrix,n_matrix,which,Ts):
    delta = new_energy-energy
    if delta >0.0:
        accept_the_change = uniform(1) # generating a random number to decide if we accept the change
        if accept_the_change < exp(-delta/Ts):
            cord[0][which] = cord[0][-1]
            cord[1][which] = cord[1][-1]
            energy = new_energy
            matrix = n_matrix
     
    else:
        cord[0][which] = cord[0][-1]
        cord[1][which] = cord[1][-1]
        energy = new_energy
        matrix = n_matrix

        
    return cord,energy,matrix


def moveCharge(number,cord,step,which,r):
    n_theta=2*uniform(1)*np.pi
    n = randrange(1,3)
    delta_radius = (-1)**n *step
    n_radius =delta_radius+cord[0][which]
    if n_radius>r or n_radius<0.0:
        n_radius -=delta_radius
    n_cord= np.array([[n_radius],[n_theta]])
    if cord.shape[1] == number+1:
        cord[0][-1] = n_radius
        cord[1][-1]= n_theta
    else:
        cord = np.append(cord,n_cord,1)
    return cord
   

def partial_energy(no,radii,thetas,enMatrix):
    """
    no- ordinary number of the charge that you moved and calculate the change in energy as a result of displacement" 
    """
    radiusA = radii[no]
    thetaA  = thetas[no]
    for key,theta in enumerate(thetas):
        if not key==no:
            radiusB = radii[key]
            thetaB = theta
            energy = cosRule(radiusB,radiusA,thetaB,thetaA)
            enMatrix[key][no]= 0.5*energy
            enMatrix[no][key] = enMatrix[key][no]
    return enMatrix

def recPartial(no,cord,matrix):
    cord = copy(cord)
    cord[0][no] = cord[0][-1]
    cord[1][no] =cord[1][-1]
    cord = np.delete(cord,-1,1)
    radiusA = cord[0][no]
    thetaA = cord[1][no]
    for key,theta in enumerate(cord[1]):
        if key!=no:
            radiusB = cord[0][key]
            thetaB = theta
            energy = cosRule(radiusB,radiusA,thetaB,thetaA)
            matrix[key][no]= 0.5*energy
            matrix[no][key] = matrix[key][no]
    return matrix
    
        
def total_energy(n,cord,enMatrix=None,which=None):
    if enMatrix is None:
        enMatrix = np.zeros([n,n])
        for i in range(n):
            enMatrixNew= partial_energy(i,cord[0],cord[1],enMatrix=enMatrix)
            energy = sum(enMatrixNew).sum()
        return energy,enMatrix
    else:
        enMatrix = copy(enMatrix)
        enMatrixNew=recPartial(which,cord,enMatrix)
        energy = sum(enMatrixNew).sum()
        return energy,enMatrixNew
       

def initialPos(number):
    radius,theta = generate_random(number)
    cord= np.array([radius,theta])
    energy,matrix = total_energy(number,cord)
    return cord,energy,matrix
 

Ts,Tf = tempScaling(0.2)
#### HERE THE SIMULATION BEGINS, WE LOOK FOR GLOBAL MINIMUM AS A FUNCTION OF ALL CHARGES POSITION ON THE DISC
def energy_find(number,Ts=Ts):
    r=10.0
    cord,energy,matrix= initialPos(number)
    step = 1e-1
    k = int(((2*r)**2/step**2)*number)
    divider =int(k/100)
    for i in range(k):
        which= randrange(number)
        cord = moveCharge(number,cord,step,which,r)
        new_energy,n_matrix= total_energy(number,cord,matrix,which)
        cord,energy,matrix= acceptChange(cord,energy,new_energy,matrix,n_matrix,which,Ts)
        if i%divider==0:
            Ts = Ts/1.3
    cord = np.delete(cord,-1,1)
    return energy,cord



for i,j in enumerate([23]):
    print(j)
    number = j
    energies=[]
    cords = []
    for i in range(30):
        energy, cord  = energy_find(number)  
        energies.append(energy)
        cords.append(cord) 
        r = 10.0
        c = "orangered"
        ### PLOTING CHARGES' POSITION ON THE POLAR PLOT AND CHECKING THE FINAL ENERGY OF THE SYSTEM
        fig, ax = plt.subplots(1, 1,figsize=(10,10), subplot_kw=dict(projection='polar'))
        ax.scatter(cord[1],cord[0],marker='o',s=100,c=c)
        ax.set_ylim([0,10.3])
        ax.set_title("Energy of this system is "+str(energy)+" eV",y=1.08)
        plt.show()
    index = energies.index(min(energies))
    final = cords[index]
    final_e = energies[index]   
    
    r = 10.0
    c = "orangered"
    ### PLOTING CHARGES' POSITION ON THE POLAR PLOT AND CHECKING THE FINAL ENERGY OF THE SYSTEM
    fig, ax = plt.subplots(1, 1,figsize=(10,10), subplot_kw=dict(projection='polar'))
    ax.scatter(final[1],final[0],marker='o',s=100,c=c)
    ax.set_ylim([0,10.3])
    ax.set_title("Energy of this system is "+str(final_e)+" eV",y=1.08)
    #plt.savefig("correcting/{}b.png".format(number))

    






import numpy as np
from math import exp,cos,floor,sqrt,log
import random
from copy import copy
from random import randrange 


r = 10
def uniform(n):
    global _first_random
    if n==0: random.seed(1234); _first_random=0
    if _first_random==1: random.seed(None); _first_random=0
    if n==1: return random.random()
    else: return floor(n*random.random()+1)
   
_first_random=1


def generate_random(number):
    """
    Function for creating random position for n charges
    
    Parameters
    ---------
    number -- number of charges in the system. Takes integer values.
    
    Outputs 
    --------
    charges_radius,
    charges_theta
    
    """
    charges_radius = []
    charges_theta = []
    for i in range(number):
        radius = randrange(r)
        theta= np.random.random() * 2.0 * np.pi
        while theta in charges_theta and radius in charges_radius or radius==0 and radius in charges_radius:
            radius = uniform(10)
            theta= np.random.random() * 2.0 * np.pi
        charges_radius.append(radius)
        charges_theta.append(theta)
    return charges_radius,charges_theta
    


def cosRule(rad1,rad2,ang1,ang2):
    """
    Function for calculating energy of a pair of charges using cosine rule 
    
    Parameters
    ---------
    ang1,ang2 --the angular coordinates of two charges in radians.
    
    rad1,rad2 -- radial coordinates of two charges.
    
    Outputs 
    --------
    resulting energy. 
    
    """
    net= ang2-ang1
    try:
        net_distance = sqrt(rad1**2+rad2**2-2*rad1*rad2*cos(net))
    except ValueError:
        print(net)
        print(rad1**2+rad2**2-2*rad1*rad2*cos(net))
        raise ValueError("math domain error")
    try:
        energy =(1.0/net_distance)  
    except ZeroDivisionError:
        energy = 1e12
    return energy




def tempScaling(a):
    """
    Function for estimating start and final temperature for cooling schedule
    
    Parameters
    ---------
    a - average change in energy of the system
    
    Outputs 
    --------
    Ts - Start energy
    Tf - Final Energy
    
    """
    T_s = -a/log(0.7) 
    T_f = -a/log(0.01) 
    return T_s, T_f
    



def acceptChange(cord,energy,new_energy,matrix,n_matrix,which,Ts):
    """
    Function for deciding whether to keep or reject new change in coordinates and resulting energy of the system
    
    Parameters
    ---------
    cord --array with radial and angular coordinates of all charges.
    
    energy -- previously found energy of the system.
    
    new_energy -- newly found energy of the system.
    
    matrix -- matrix of previously found energy of the system.
    
    n_matrix -- matrix of newly found energy of the system.
    
    which -- index of the charge with changed position.
    
    Ts -- temperature of the system.
    
    Outputs 
    --------
    cord --  array with new radial and angular coordinates of all charges.
    
    energy -- accepted energy of the system 
    
    matrix -- matrix of accepted energy of the system
    
    """
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
    """
    Function for moving a charge.
    
    Parameters
    ---------
    number -- number of charges in the system.
    
    cord --array with radial and angular coordinates of all charges.
    
    step -- radial step size
    
    which -- index of the charge with changed position.
    
    r -- radius of the disc 
    
    Outputs 
    --------
    cord --new array with radial and angular coordinates of all charges.
    
    """
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
    Function for calculating partial energy of the system from one charge.
    
    Parameters
    ---------
    no -- index of charge we are calculating partial energies of with relation to the system.
    
    radii -- radial coordinates of the charge
    
    thetas -- angular coordinates of the charge
    
    enMatrix -- energy matrix of the system
    
    Outputs 
    --------
    enMatrix -- updated energy matrix of the system
    
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
    """
    Function for recalculating partial energy of the system from one charge.
    
    Parameters
    ---------
    no -- index of charge we are calculating partial energies of with relation to the system.

    cord --array with radial and angular coordinates of all charges.
    
    matrix -- energy matrix of the system
    
    Outputs 
    --------
    matrix -- updated energy matrix of the system
    
    """
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
    """
    Function for calculating total energy of the system.
    
    Parameters
    ---------
    n -- number of the charges in the system

    cord --array with radial and angular coordinates of all charges.
    
    enMatrix -- previous energy matrix of the system
    
    which -- index of last changed charge's position.
    
    Outputs 
    --------
    energy -- new total energy of the system.
    
    matrix -- updated energy matrix of the system.
    
    """
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
    """
    Function for creating initial positions of charges
    
    Parameters
    ---------
    number -- number of charges in the system. Takes integer values.
    
    Outputs 
    --------
    energy -- total energy of the system.
    
    matrix -- energy matrix of the system.
    
    cord --array with radial and angular coordinates of all charges.
    
    """
    radius,theta = generate_random(number)
    cord= np.array([radius,theta])
    energy,matrix = total_energy(number,cord)
    return cord,energy,matrix



Ts,Tf = tempScaling(0.2)
#### HERE THE SIMULATION BEGINS, WE LOOK FOR GLOBAL MINIMUM AS A FUNCTION OF ALL CHARGES POSITION ON THE DISC
def energy_find(number,Ts=Ts):
    """
    Function for implementing cooling schedule procedure in order to find the global energy minimum.
    
    Parameters
    ---------
    number -- number of charges in the system. Takes integer values.
    
    Ts -- starting temperature
    
    Outputs 
    --------
    energy -- total energy of the system.
    
    cord --array with radial and angular coordinates of all charges.
    
    """
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














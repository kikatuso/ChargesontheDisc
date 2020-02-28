import numpy as np
import matplotlib.pyplot as plt
from random import randrange
from math import cos,sqrt,exp,floor,radians
import random
import pandas as pd
import matplotlib.animation as animation
from statistics import mean

def uniform(n):
    global _first_random
    if n==0: random.seed(1234); _first_random=0
    if _first_random==1: random.seed(None); _first_random=0
    if n==1: return random.random()
    else: return floor(n*random.random()+1)
   
_first_random=1
r = 10

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
    

def calculate_energy(n,charges_radius,charges_theta):
    """Function for calculating the net distances between charges and the consequent energy of the system. 
    
    Parameters
    ---------
    n -- number of charges in the system. Takes integer values (default 3)
    
    Outputs 
    --------
    energy -- The resultant energy of the system. Expressed as a float. 
    
    """
    q = 1
    net_distances=[]
    energy=0.0
    for i in range(n-1):
        radius_1=charges_radius[i]
        theta_1=charges_theta[i]
        radii = list(charges_radius[i+1:]) # all radii with index higher than the current one
        thetas = list(charges_theta[i+1:]) # all thetas with index higher than the current one
        for i,val in enumerate(thetas):
            net = theta_1-val
            if net == np.inf or net==-np.inf:
                print(theta_1)
                print(val)
                raise ValueError('infinity!')
            net_distance = sqrt(radius_1**2+radii[i]**2-2*radii[i]*radius_1*cos(net))
            net_distances.append(net_distance)
            energy += q*q*(1.0/net_distance)
            
    return energy,net_distances


def scale(c):
    scalling=1.0
    if c==int(0.8*336001):
        scalling = 5
    elif c==int(0.9*336001):
        scalling = 10.0
    return scalling

 

#### HERE THE SIMULATION BEGINS, WE LOOK FOR GLOBAL MINIMUM AS A FUNCTION OF ALL CHARGES POSITION ON THE DISC
def energy_find(number):
    val_dic = {}
    delta = []
    r = 10
    c = 0 
    if c==0:
        radius,theta = generate_random(number)
        energy,net_distances = calculate_energy(number,radius,theta)
        val_dic[c] = []
        val_dic[c] = {"radius":radius,"theta":theta,"energy":energy}

    T_0 = 1 # our initial temperature expressed in Kelvins 
    theta_incr = 0.0174533
    radial_incr = 0.10
    m = 1000 # number of repetitions per given temperature
    while T_0 > 1e-7:
        for i in range(m):
            c +=1
            old_theta = val_dic[c-1]["theta"]
            old_radius = val_dic[c-1]["radius"]
            energy= val_dic[c-1]["energy"]
            which_charge = randrange(number) # creating a random number from range 0-number to find a random charge by index 
            which_coordinate = randrange(2) # choosing a coordinate to be changed at random
            if which_coordinate==0:
                new_theta = list(old_theta)
                new_radius = list(old_radius)
                n = randrange(1,3)
                delta_theta = (-1)**n *theta_incr
                new_theta[which_charge] +=delta_theta
                if new_theta[which_charge]>2*np.pi:
                    new_theta[which_charge] -= 2*np.pi
                elif new_theta[which_charge]<0.0:
                    new_theta[which_charge] += 2*np.pi
            if which_coordinate==1:
                new_theta = list(old_theta)
                new_radius =list(old_radius)
                n = randrange(1,3)
                delta_radius = (-1)**n *radial_incr
                new_radius[which_charge] += delta_radius
                if new_radius[which_charge]>r or new_radius[which_charge]<0.0:
                    new_radius[which_charge]=old_radius[which_charge]
            new_energy,net_distances = calculate_energy(number,new_radius,new_theta)
            delta_energy = new_energy-energy
            delta.append(delta_energy)
            if delta_energy >=0.0:
                accept_the_change = uniform(1) # generating a random number to decide if we accept the change
                if accept_the_change < exp(-delta_energy/T_0):
                    val_dic[c]=[]
                    val_dic[c]={"radius":new_radius,"theta":new_theta,"energy":new_energy}

                else:
                    val_dic[c]=[]
                    val_dic[c]={"radius":old_radius,"theta":old_theta,"energy":energy}

            else:
                val_dic[c]=[]
                val_dic[c] = {"radius":new_radius,"theta":new_theta,"energy":new_energy}

        T_0= 0.95*T_0
        #ax = plt.subplot(111, projection='polar')
        #ax.plot(val_dic[c]["theta"],val_dic[c]["radius"],marker='o',markersize=7)
        #plt.show()
        
        
    df = pd.DataFrame(val_dic).T
    energy = df.energy.min()
    index =  pd.to_numeric(df.energy).idxmin()
    theta = df.loc[index,"theta"]
    radius = df.loc[index,"radius"]
    
    return df,energy,radius,theta,delta


df,energy,radius,theta,delta= energy_find(11)

### PLOTING CHARGES' POSITION ON THE POLAR PLOT AND CHECKING THE FINAL ENERGY OF THE SYSTEM
ax = plt.subplot(111, projection='polar')
ax.plot(theta,radius,marker='o',markersize=7)
ax.set_yticks(np.linspace(1.0,r,r//2))
ax.set_title("Energy of this system is "+str(energy)+" eV",y=1.08)
ax.set_xlabel("The energy difference is  " + str(energy-0.6881909602355869))
plt.show()




"""



fig = plt.figure(figsize=(4,4))
ax = fig.add_subplot(111, projection='polar')
ax.set_ylim(0,10)


l,  = ax.plot([],[],'*',markersize=12)

def update(i):
    global df
    radius= df.radius[df.shape[0]-1-i]
    theta = df.theta[df.shape[0]-1-i]
    #df.loc['radius',-1] = df.iloc['radius',0]
    #df.loc['theta',-1] = df.loc['theta',0]
    l.set_data(theta, radius )
    return l, 

ani = animation.FuncAnimation(fig, update, frames=1000, interval=200, blit=True)
ani.save('simulation.gif', writer = 'imagemagick')

"""








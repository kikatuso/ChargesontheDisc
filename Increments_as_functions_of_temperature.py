import numpy as np
import matplotlib.pyplot as plt
from random import randrange
from math import cos,sqrt,exp,floor,radians,e
import random
import pandas as pd

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
            net_distance = sqrt(radius_1**2+radii[i]**2-2*radii[i]*radius_1*cos(net))
            net_distances.append(net_distance)
            energy += q*q*(1.0/net_distance)
            
    return energy,net_distances


def scale(c):
    scalling=1.0
    if c>200:
        scalling = 1.5
    elif c>300:
        scalling = 2.0
    return scalling

 
    

#### HERE THE SIMULATION BEGINS, WE LOOK FOR GLOBAL MINIMUM AS A FUNCTION OF ALL CHARGES POSITION ON THE DISC
def energy_find(number):
    val_dic = {}
    r = 10
    c = 0 
    if c==0:
        radius,theta = generate_random(number)
        energy,net_distances = calculate_energy(number,radius,theta)
        val_dic[c] = []
        val_dic[c] = {"radius":radius,"theta":theta,"energy":energy}
        #entry=pd.DataFrame({"radius":[radius],"theta":[theta],"energy":energy})
        #values=values.append(entry,ignore_index=True)
    T_0 = 300 # our initial temperature expressed in Kelvins 
    theta_incr = 0.0174533
    radial_incr = 0.5
    m = 50000 # number of repetitions per given temperature
    while T_0 > 1e-2:
        for i in range(m):
            c +=1
            theta_incr = theta_incr*scale(c)
            radial_incr = radial_incr*scale(c)
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
                if new_theta[which_charge]<0 or new_theta[which_charge]>2*np.pi:
                    new_theta[which_charge]=np.random.random()*2*np.pi
            if which_coordinate==1:
                new_theta = list(old_theta)
                new_radius =list(old_radius)
                n = randrange(1,3)
                delta_radius = (-1)**n *radial_incr
                new_radius[which_charge] += delta_radius
                if new_radius[which_charge]>r or new_radius[which_charge]<0.0:
                    new_radius[which_charge]=np.random.random()*10.0
            new_energy,net_distances = calculate_energy(number,new_radius,new_theta)
            delta_energy = new_energy-energy
            if delta_energy >=0.0:
                accept_the_change = uniform(1) # generating a random number to decide if we accept the change
                if accept_the_change < exp(-delta_energy/T_0):
                    val_dic[c]=[]
                    val_dic[c]={"radius":new_radius,"theta":new_theta,"energy":new_energy}

                else:
                    val_dic[c]=[]
                    val_dic[c]={"radius":old_radius,"theta":old_theta,"energy":energy}
#

            else:
                val_dic[c]=[]
                val_dic[c] = {"radius":new_radius,"theta":new_theta,"energy":new_energy}

        T_0= (1/e)*T_0
        
    val_df = pd.DataFrame(val_dic).T
    energy = val_df.energy.min()
    index =  val_df.energy.idxmin()
    theta = val_df.loc[index,"theta"]
    radius = val_df.loc[index,"radius"]
    return val_df,energy,theta,radius


val_df,energy,theta,radius = energy_find(11)


### PLOTING CHARGES' POSITION ON THE POLAR PLOT AND CHECKING THE FINAL ENERGY OF THE SYSTEM
ax = plt.subplot(111, projection='polar')
ax.plot(theta,radius,marker='o',markersize=7)
ax.set_yticks(np.linspace(1.0,r,r//2))
ax.set_title("Energy of this system is "+str(energy)+" eV",y=1.08)
ax.set_xlabel("The energy difference is  " + str(energy-0.6881909602355869))
plt.show()









import numpy as np
import matplotlib.pyplot as plt
from random import randrange
from math import cos,sqrt,exp,floor
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

     
#### HERE THE SIMULATION BEGINS, WE LOOK FOR GLOBAL MINIMUM AS A FUNCTION OF ALL CHARGES POSITION ON THE DISC
def findingbest(z,l,number):
    r = 10
    
    charges_radius,charges_theta = generate_random(number)
    energy,net_distances = calculate_energy(number,charges_radius,charges_theta)
    T_0 = 300 # our initial temperature expressed in Kelvins 
    theta_incr = 0.0174533*z # expressed in radians, equivalent of 1 degree
    radial_incr = 0.05*l
    m = 1000 # number of repetitions per given temperature
    while T_0 > 1e-9:
        for i in range(m):
            kept_charges_radius=list(charges_radius)
            kept_charges_theta = list(charges_theta)
            which_charge = randrange(number) # creating a random number from range 0-number to find a random charge by index 
            which_coordinate = randrange(2) # choosing a coordinate to be changed at random
            if which_coordinate==0: 
                n = randrange(1,3)
                delta_theta = (-1)**n *theta_incr
                charges_theta[which_charge] +=delta_theta
                if charges_theta[which_charge]<0.0 or charges_theta[which_charge]>2*np.pi:
                    charges_theta[which_charge]=np.random.random() * 2.0 * np.pi
            if which_coordinate==1:
                n = randrange(1,3)
                delta_radius = (-1)**n *radial_incr
                charges_radius[which_charge] += delta_radius
                if charges_radius[which_charge]>r or charges_radius[which_charge]<0.0:
                    charges_radius[which_charge]=np.random.random()*10.0
    
            changed_radius = charges_radius[which_charge]
            changed_theta =  charges_theta[which_charge]
            except_one_radius = list(charges_radius)
            except_one_theta = list(charges_theta)
            del except_one_radius[which_charge]
            del except_one_theta[which_charge]
            if not changed_radius in except_one_radius and not changed_theta in except_one_theta and charges_radius.count(0.0)<=1:
                new_energy,net_distances = calculate_energy(number,charges_radius,charges_theta)
                delta_energy = new_energy-energy
                if delta_energy >=0.0:
                    accept_the_change = uniform(1) # generating a random number to decide if we accept the change
                    if accept_the_change < exp(-delta_energy/T_0):
                        energy += delta_energy
                        
                    else:
                        energy = energy
                        charges_theta = kept_charges_theta
                        charges_radius = kept_charges_radius
                else:
                    energy +=delta_energy
            
        T_0= 0.9*T_0
    return z,l,energy,charges_radius,charges_theta



"""
find_best_matrix = pd.DataFrame(columns=["first_para","second_para","energy"])
for i in np.arange(0.1,1.0,0.05):
    for y in np.arange(0.25,5.0,0.25):
        z,l,energy,charges_radius,charges_theta= findingbest(y,i)
        find_best_matrix = find_best_matrix.append({"first_para":z,"second_para":l,"energy":energy},ignore_index=True)


z,l,energy,charges_radius,charges_theta= findingbest(1.0,1.0,7)
### PLOTING CHARGES' POSITION ON THE POLAR PLOT AND CHECKING THE FINAL ENERGY OF THE SYSTEM
ax = plt.subplot(111, projection='polar')
ax.plot(charges_theta,charges_radius,marker='o',markersize=7)
ax.set_yticks(np.linspace(1.0,r,r//2))
ax.set_title("Energy of this system is "+str(energy)+" eV",y=1.08)
ax.set_xlabel("The energy difference is  " + str(energy-sqrt(3)/10.0))
plt.show()
"""

def update_table():
    energy_df = pd.read_csv("Energy_min_table.csv")
    energy_obtained = []
    energy_difference= []
    
    for i,val in enumerate(np.arange(2,7)):
        print(val)
        z,l,energy,charges_radius,charges_theta= findingbest(1.0,0.3,val)
        energy_obtained.append(energy)
        difference = energy-energy_df["min energy"][i] 
        energy_difference.append(difference)
    
    energy_df= energy_df.assign(energy_obtained=energy_obtained)
    energy_df= energy_df.assign(energy_difference=energy_difference)     
    energy_df.to_csv("Energy_min_comparision.csv",header=True,index=False, encoding='utf-8')
    return energy_df


    
    

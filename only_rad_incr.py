from copy import copy
from statistics import median
import pandas as pd
from random import randrange
from functions import generate_random,total_energy,tempScaling,uniform,moveCharge
import numpy as np 
from math import exp,sin,cos
from matplotlib import pyplot as plt

T_s, T_f = tempScaling(0.2)
 
r = 10



"""
def cartes(radius,theta):
    x=[]
    y=[]
    for i,val in enumerate(radius):
       x.append(val*cos(theta[i])) 
       y.append(val*sin(theta[i]))
    return x,y

def polar(x,y):
    radius=[]
    theta = []
    for i, val in enumerate(x):
        radius.append(np.sqrt(val**2+y[i]**2))
        theta.append(np.arctan2(y[i],val))
    return radius,theta


r = 10.0

def moveCharge(number,theta,radius,step,which):
    r=10.0
    x_old,y_old = cartes(radius,theta)
    angle = 2*uniform(1)*np.pi
    n = randrange(1,3)
    delta_radius = (-1)**n *step
    x_new = copy(x_old)
    y_new = copy(y_old)
    x_new[which]=x_old[which]+cos(angle)*delta_radius
    y_new[which] = y_old[which] +sin(angle)*delta_radius
    n_radius,n_theta = polar(x_new,y_new)
    if n_radius[which]>r or n_radius[which]<0.0:
        n_radius[which] =radius[which]
    return n_theta,n_radius

"""



#### HERE THE SIMULATION BEGINS, WE LOOK FOR GLOBAL MINIMUM AS A FUNCTION OF ALL CHARGES POSITION ON THE DISC
def energy_find(number,Ts=T_s,Tf=T_f):
    val_dic = {}
    delta = []
    c = 0 
    if c==0:
        radius,theta = generate_random(number)
        energy,enMatrix = total_energy(number,radius,theta)
        val_dic[c] = []
        val_dic[c] = {"radius":radius,"theta":theta,"energy":energy,"temp":Ts,"energies":enMatrix}
    radial_incr = 0.5
    m = 500 # number of repetitions per given temperature
    while Ts > Tf:
        for i in range(m):
            c +=1
            old_theta = copy(val_dic[c-1]["theta"])
            old_radius = copy(val_dic[c-1]["radius"])
            energy= copy(val_dic[c-1]["energy"])
            old_energies = copy(val_dic[c-1]["energies"])
            which_charge = randrange(number) # creating a random number from range 0-number to find a random charge by index 
            w=0
            while True:
                new_theta,new_radius=moveCharge(number,old_theta,old_radius,radial_incr,which_charge)
                new_energy,enMatrix= total_energy(number,new_radius,new_theta,old_energies,which_charge)
                delta_energy = new_energy-energy
                if delta_energy >0.0:
                    delta.append(delta_energy)
                    accept_the_change = uniform(1) # generating a random number to decide if we accept the change
                    if accept_the_change < exp(-delta_energy/Ts):
                        val_dic[c]=[]
                        val_dic[c]={"radius":new_radius,"theta":new_theta,"energy":new_energy,"temp":Ts,"energies":enMatrix}
                        w +=1
                     
                            
                    else:
                        val_dic[c]=[]
                        val_dic[c]={"radius":old_radius,"theta":old_theta,"energy":energy,"temp":Ts,"energies":old_energies}
                        w +=1
                        
                    
                else:
                    val_dic[c]=[]
                    val_dic[c] = {"radius":new_radius,"theta":new_theta,"energy":new_energy,"temp":Ts,"energies":enMatrix}
                    break
                if w%360==0:
                    radial_incr = np.random.normal(0.25,0.10)
                if w == 10*360:
                    break

        Ts=Ts*0.95
        
    df = pd.DataFrame(val_dic).T
    energy = df.energy.min()
    index =  pd.to_numeric(df.energy).idxmin()
    theta = df.loc[index,"theta"]
    radius = df.loc[index,"radius"]
    
    return df,energy,radius,theta,delta


df,energy,radius,theta,delta= energy_find(11)

### PLOTING CHARGES' POSITION ON THE POLAR PLOT AND CHECKING THE FINAL ENERGY OF THE SYSTEM
ax = plt.subplot(111, projection='polar')
ax.plot(theta,radius,marker='o',markersize=7,c="fuchsia")
ax.set_yticks(np.linspace(1.0,r,r//2))
ax.set_title("Energy of this system is "+str(energy)+" eV",y=1.08)
ax.set_xlabel("The energy difference is  " + str(energy-1.0964101615137758))
plt.show()














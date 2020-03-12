import numpy as np
import matplotlib.pyplot as plt
from math import exp,atan,cos,sin
import random
import pandas as pd
from statistics import median,mean
import copy
from random import randrange 
from functions import generate_random,total_energy,tempScaling,uniform,moveCharge



def initialPos(val_dic,number):
    radius,theta = generate_random(number)
    energy,matrix = total_energy(number,radius,theta)
    val_dic.setdefault(0,{"radius":radius,"theta":theta,"energy":energy,"energies":matrix})
    return val_dic

def acceptChange(newset,delta,val_dic,c,Ts,acceptanceR):
    if delta >0.0:
        accept_the_change = uniform(1) # generating a random number to decide if we accept the change
        if accept_the_change < exp(-delta/Ts):
            val_dic.setdefault(c,{"radius":newset[0],"theta":newset[1],"energy":newset[2],"energies":newset[3]})
            acceptanceR.append(1)
  
        else:
            old_set = copy.copy(val_dic[c-1])
            val_dic.setdefault(c,old_set)
            acceptanceR.append(0)
    else:
        val_dic.setdefault(c,{"radius":newset[0],"theta":newset[1],"energy":newset[2],"energies":newset[3]})
        acceptanceR.append(None)
    
    return val_dic,acceptanceR


 
def stepFind(energy,energy2,step):
    if energy2 is None:
        step = 0.5        
    else:
        grad = abs(np.gradient([energy2,energy])[0])
        decInt = uniform(1)*abs(np.gradient(df.energy)).max()
        if decInt < atan(grad):
            step=0.9*step
            if step<1e-4:
                step = 1e-4
        else:
            step = 1.1*step
            if step>0.5:
                step=0.5
        
    return step
        
    

def acceptR(accR,c):
    yes = 0 
    for i,val in enumerate(accR[:c]):
        if val==1:
            yes +=1
    try:
        divider = [x for x in accR[:c] if x is not None]
        ratio= yes/len(divider)
    except ZeroDivisionError:
        ratio= yes/1.0
    return ratio
    


Ts, Tf = tempScaling(0.2)
#### HERE THE SIMULATION BEGINS, WE LOOK FOR GLOBAL MINIMUM AS A FUNCTION OF ALL CHARGES POSITION ON THE DISC
def energy_find(number,Ts=Ts,Tf=Tf):
    val_dic = {}
    delta_list = []
    val_dic = initialPos(val_dic,number)
    k = 3000 # number of repetitions per given temperature
    acceptanceR = []
    ratio = []
    while Ts > Tf:
        for i in range(k):
            step = random.uniform(5*1e-3,0.5)
            c =i+1
            old_theta = copy.copy(val_dic[c-1]["theta"])
            old_radius = copy.copy(val_dic[c-1]["radius"])
            energy= copy.copy(val_dic[c-1]["energy"])
            old_energies = copy.copy(val_dic[c-1]["energies"])
            try:
                energy2= copy.copy(val_dic[c-2]["energy"])
            except KeyError:
                energy2 = None
            which= randrange(number)
            new_theta,new_radius = moveCharge(number,old_theta,old_radius,step,which)
            new_energy,enMatrix= total_energy(number,new_radius,new_theta,old_energies,which)
            delta_energy = new_energy-energy
            if delta_energy>=0.0:
                delta_list.append(delta_energy)
            newset = [new_radius,new_theta,new_energy,enMatrix]
            val_dic,acceptanceR = acceptChange(newset,delta_energy,val_dic,c,Ts,acceptanceR)
            ratio.append(acceptR(acceptanceR,i-1))

     
        Ts = Ts*0.90
    
    df = pd.DataFrame(val_dic).T
    energy = df.energy.min()
    index =  pd.to_numeric(df.energy).idxmin()
    theta = df.loc[index,"theta"]
    radius = df.loc[index,"radius"]
    return df,energy,radius,theta,delta_list,ratio,acceptanceR


df,energy,radius,theta,delta,ratio,acceptanceR= energy_find(3)






### PLOTING CHARGES' POSITION ON THE POLAR PLOT AND CHECKING THE FINAL ENERGY OF THE SYSTEM
ax = plt.subplot(111, projection='polar')
ax.plot(theta,radius,marker='o',markersize=7,c="fuchsia")
ax.set_yticks(np.linspace(1.0,r,r//2))
ax.set_title("Energy of this system is "+str(energy)+" eV",y=1.08)
ax.set_xlabel("The energy difference is  " + str(energy-1.0964101615137758))
plt.show()










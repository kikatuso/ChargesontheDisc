import numpy as np
import matplotlib.pyplot as plt
from random import randrange
from math import cos,sqrt,exp,floor,radians,log
import random
import pandas as pd
import matplotlib.animation as animation
from statistics import median

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
    


def cosRule(rad1,rad2,ang1,ang2):
    q = 1.0
    net= ang2-ang1
    net_distance = sqrt(rad1**2+rad2**2-2*rad1*rad2*cos(net))
    try:
        energy = q*q*(1.0/net_distance)
        
    except ZeroDivisionError:
        energy = 1e12
    return energy

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

def total_energy(n,radius,thetas,enMatrix=None,which=None):
    if enMatrix is None:
        enMatrix = np.zeros([n,n])
        for i in range(n):
            enMatrixNew= partial_energy(i,radius,thetas,enMatrix=enMatrix)
            energy = sum(enMatrixNew).sum()
        return energy,enMatrix
    else:
        enMatrixNew=partial_energy(which,radius,thetas,enMatrix)
        energy = sum(enMatrixNew).sum()
        return energy,enMatrixNew
       


def tempScaling(a):
    T_s = -a/log(0.7) 
    T_f = -a/log(0.01) 
    return T_s, T_f
    
 
T_s, T_f = tempScaling(0.2)
 


#### HERE THE SIMULATION BEGINS, WE LOOK FOR GLOBAL MINIMUM AS A FUNCTION OF ALL CHARGES POSITION ON THE DISC
def energy_find(number,Ts=T_s,Tf=T_f):
    val_dic = {}
    delta = []
    r = 10
    c = 0 
    if c==0:
        radius,theta = generate_random(number)
        energy,enMatrix = total_energy(number,radius,theta)
        val_dic[c] = []
        val_dic[c] = {"radius":radius,"theta":theta,"energy":energy,"temp":Ts,"energies":enMatrix}
    radial_incr = 0.5
    m = 5 # number of repetitions per given temperature
    while Ts > Tf:
        for i in range(m):
            radial_incr = 0.5
            c +=1
            old_theta = val_dic[c-1]["theta"]
            old_radius = val_dic[c-1]["radius"]
            energy= val_dic[c-1]["energy"]
            old_energies = val_dic[c-1]["energies"]
            which_charge = randrange(number) # creating a random number from range 0-number to find a random charge by index 
            new_theta = list(old_theta)
            new_radius = list(old_radius)
            w=0
            while True:
                new_theta[which_charge] =2*uniform(1)*np.pi
                n = randrange(1,3)
                delta_radius = (-1)**n *radial_incr
                new_radius[which_charge] += delta_radius
                if new_radius[which_charge]>r or new_radius[which_charge]<0.0:
                    new_radius[which_charge]=old_radius[which_charge]
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
                    radial_incr = 0.9*radial_incr
                if radial_incr< 5*1e-3:
                    break
                    


        Ts=Ts*0.95
        
    df = pd.DataFrame(val_dic).T
    energy = df.energy.min()
    index =  pd.to_numeric(df.energy).idxmin()
    theta = df.loc[index,"theta"]
    radius = df.loc[index,"radius"]
    
    return df,energy,radius,theta,delta


df,energy,radius,theta,delta= energy_find(5)

### PLOTING CHARGES' POSITION ON THE POLAR PLOT AND CHECKING THE FINAL ENERGY OF THE SYSTEM
ax = plt.subplot(111, projection='polar')
ax.plot(theta,radius,marker='o',markersize=7,c="fuchsia")
ax.set_yticks(np.linspace(1.0,r,r//2))
ax.set_title("Energy of this system is "+str(energy)+" eV",y=1.08)
ax.set_xlabel("The energy difference is  " + str(energy-1.0964101615137758))
plt.show()




def createGif(df):
    fig = plt.figure(figsize=(4,4))
    ax = fig.add_subplot(111, projection='polar')
    ax.set_ylim(0,10)
    l,  = ax.plot([],[],'*',markersize=12,c="fuchsia")
    
    def update(i):
        global df
        radius= df.radius[100*i]
        theta = df.theta[100*i]
        l.set_data(theta, radius )
        
        return l
    ani = animation.FuncAnimation(fig, update, frames=int(df.shape[0]/100), interval=100, blit=True)
    ani.save('simulation.gif')
    print("Movie saved succesfully")


def createMovie(df):
    fig = plt.figure(figsize=(4,4))
    ax = fig.add_subplot(111, projection='polar')
    ax.set_ylim(0,10)
    l,  = ax.plot([],[],'*',markersize=12)
    ax.legend()
    def update(i):
        radius= df.radius[i]
        theta = df.theta[i]
        temp = df.temp[i]
        l.set_data(theta, radius )
        l.set_label(temp)
        plt.legend(bbox_to_anchor=(1.20, 0.05))
        return l,  
    plt.rcParams['animation.ffmpeg_path'] = r'C:/Users/zuzan/ffmpeg-20200227-9b22254-win64-static/ffmpeg-20200227-9b22254-win64-static/bin/ffmpeg'
    FFwriter=animation.FFMpegWriter(fps=30, extra_args=['-vcodec', 'libx264'])
    ani = animation.FuncAnimation(fig, update, frames=int(df.shape[0]), interval=100, blit=True)
    ani.save('simulation.mp4', writer = FFwriter)














import numpy as np
import matplotlib.pyplot as plt
from random import randrange
from math import cos,sqrt,floor,log
import random
import matplotlib.animation as animation
from copy import copy


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
        enMatrix = copy(enMatrix)
        enMatrixNew=partial_energy(which,radius,thetas,enMatrix)
        energy = sum(enMatrixNew).sum()
        return energy,enMatrixNew
       


def tempScaling(a):
    T_s = -a/log(0.7) 
    T_f = -a/log(0.01) 
    return T_s, T_f
    


def moveCharge(number,thetas,radius,radial_incr,which):
    n_thetas = copy(thetas)
    n_radius = copy(radius)
    n_thetas[which] =2*uniform(1)*np.pi
    n = randrange(1,3)
    delta_radius = (-1)**n *radial_incr
    n_radius[which] +=delta_radius
    if n_radius[which]>r or n_radius[which]<0.0:
        n_radius[which] -=delta_radius
    return n_thetas,n_radius


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








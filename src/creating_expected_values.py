

from math import sqrt,cos,radians,sin,atan,acos
import matplotlib.pyplot as plt
import numpy as np

from random import randrange
from copy import copy

"""
def moveCharge(number,thetas,radius):
    step = 2.05
    which = randrange(number)
    thetas[which] =2*uniform(1)*np.pi
    n = randrange(1,3)
    delta_radius = (-1)**n *step
    radius[which] +=delta_radius
    if radius[which]>r or radius[which]<0.0:
        radius[which] -=delta_radius
    return thetas,radius,which
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

def moveCharge(number,theta,radius,which):
    r=10.0
    step = 5.5
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


from functions import total_energy,uniform
from math import radians

number =28
inners = 5 
r1 = 10.0
energies = []
radi2 = []
for i in np.arange(4.0,6.0,0.001):
    r2 =4.524000000000175
    angle= 360.0/float(number-inners)
    theta = [radians(angle*i) for i in range(number-inners)]
    radius=[r1 for i in range(number-inners)]
    angle2= 360.0/float(inners)
    theta2 = [radians(angle2*i+22) for i in range(inners)]
    radius2 = [r2 for i in range(inners)]
    radius.extend(radius2)
    theta.extend(theta2)
    energy, matrix = total_energy(number,radius,theta)
    energies.append(energy)
    radi2.append(i)
index = energies.index(min(energies))
r2 = radi2[index]





r = 10.0
c = "orangered"
### PLOTING CHARGES' POSITION ON THE POLAR PLOT AND CHECKING THE FINAL ENERGY OF THE SYSTEM
fig, ax = plt.subplots(1, 1,figsize=(10,10), subplot_kw=dict(projection='polar'))
ax.scatter(theta,radius,marker='o',s=100,c=c)
ax.set_ylim([0,10.3])
ax.set_title("Energy of this system is "+str(energy)+" eV",y=1.08)
plt.savefig("correcting/{}.png".format(number))

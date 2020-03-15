
from functions import total_energy,uniform
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




number =3
angle= 360.0/float(number)
theta = [radians(angle*i) for i in range(number)]
radius=[r for i in range(number)]
energy, matrix = total_energy(number,radius,theta)


ax = plt.subplot(111, projection='polar')
ax.plot(theta,radius,marker='o',markersize=7)
ax.set_yticks(np.linspace(1.0,r,r//2))
plt.show()
energy, matrix = total_energy(number,radius,theta)
print(energy)


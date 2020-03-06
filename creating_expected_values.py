

from math import sqrt,cos,radians
import matplotlib.pyplot as plt
import numpy as np
r=10.0
from random import randrange

def uniform(n):
    global _first_random
    if n==0: random.seed(1234); _first_random=0
    if _first_random==1: random.seed(None); _first_random=0
    if n==1: return random.random()
    else: return floor(n*random.random()+1)
   
_first_random=1

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
            enMatrix[key][no] =enMatrix[no][key]= 0.5*energy
    return enMatrix

def total_energy(n,radius,thetas,enMatrix=None,which=None):
    if enMatrix is None:
        enMatrix = np.zeros([n,n])
        for i in range(n):
            enMatrix= partial_energy(i,radius,thetas,enMatrix=enMatrix)

            energy = sum(enMatrix).sum()
        return energy,enMatrix
    else:
        enMatrixNew=partial_energy(which,radius,thetas,enMatrix)
        energy = sum(enMatrixNew).sum()
        return energy,enMatrixNew
    
number =3
angle= 120.0
theta1 = [radians(angle*i) for i in range(number)]
radius1=[r for i in range(number)]


#energy,matrix = total_energy(number,radius1,theta1)
"""
print(energy)
print(matrix)
print("================================================================================")
"""
def moveCharge(number,thetas,radius):
    step = 2.05
    which = randrange(number)
    thetas[which] =2*uniform(1)*np.pi
    n = randrange(1,3)
    delta_radius = (-1)**n *step
    radius[which] +=delta_radius
    if radius[which]>r or radius[which]<0.0:
        radius[which] +=(-1)**(n+1) * step
    return thetas,radius,which

#theta2,radius2,which = moveCharge(number,theta1,radius1)
#print(str(which)+" charge that was moved")
#energy1, matrixNew = total_energy(number,radius2,theta2,matrix,which)
#energy2, matrix = total_energy(number,radius2,theta2)
"""
print(energy1)
print(matrixNew)
print("------------------------------------------------------------------------------")
print(energy2)
print(matrix)
"""
energy2, matrix = total_energy(5,radius,theta)
ax = plt.subplot(111, projection='polar')
ax.plot(theta,radius,marker='o',markersize=7)
ax.set_yticks(np.linspace(1.0,r,r//2))
print(energy)




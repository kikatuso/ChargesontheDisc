

import math as m
import matplotlib.pyplot as plt
import numpy as np
r=10.0


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
            net_distance = m.sqrt(radius_1**2+radii[i]**2-2*radii[i]*radius_1*m.cos(net))
            net_distances.append(net_distance)
            energy += q*q*(1.0/net_distance)
            
    return energy,net_distances






n =8
angle= 360.0/8.0

theta = [m.radians(angle*i) for i in range(n)]
radius=[r for i in range(n)]
energy,net_distances  = calculate_energy(n,radius,theta)
ax = plt.subplot(111, projection='polar')
ax.plot(theta,radius,marker='o',markersize=7)
ax.set_yticks(np.linspace(1.0,r,r//2))
plt.show()
print(energy)

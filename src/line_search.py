





import numpy as np 
from scipy.optimize import line_search

x = [1,25,100,1300] # list of energies


def objective(p):
    return x[p]



def gradient(p):
    if not p==0:
        low_x  = x[p-1]
        current_x = x[p]
        grad = np.gradient([current_x,low_x])
    else:
        grad = [0.2, 0.2]
        
    return grad


p = 3
sd = -gradient(p)
a = line_search(objective, gradient, x, p)[0]
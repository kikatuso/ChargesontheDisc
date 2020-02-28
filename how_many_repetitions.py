from random import randrange



def repetitions(radial_incr,disc_rad,rang):
    radius =0 
    for j in range(int(rang)):
        for i in range(0,j):
            which_coordinate = randrange(2) # choosing a coordinate to be changed at random               
            if which_coordinate==1:
                x = randrange(1,3)
                radius += (-1)**x *radial_incr
            if radius >= 2*disc_rad:
                n_found = j 
                return n_found
    return "Increase the range"
                


n=repetitions(0.05,10.0,10000)
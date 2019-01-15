##Code to take in arbitrary shaped regions designated by 1s and 0s and calculates spatial statistics on them

import numpy as np
import numpy.random as rnd
import matplotlib.pyplot as plt

def is_wall(grid,current,step):
    """checks if grid + step exceeds grid.
    If so, return True
    """
    dim = np.shape(grid)
    stepped = current+step
    if np.any(stepped < 0) or stepped[0] > dim[0]-1 or stepped[1] > dim[1]-1:
        return True
    else:
        return False

def num_to_step(num):
    """
    converts number 1 to 4 into a step in either x or y direction
           2
        1     3
           4
    """
    if step == 1:
        return np.array([-1,0])
    elif step == 2:
        return np.array([0,1])
    elif step == 3:
        return np.array([1,0])
    elif step == 4:
        return np.array([0,-1])
    
def make_grid(n_areas,n_length,n_height):
    """
    generates an arbitrary shaped region through random walking starting in the centre
    """

    def random_walk(loc):
        """
        random walker direction
           2
        1     3
           4
        """
        step = num_to_step(rnd.randint(1,5))
        if is_wall(grid,loc,step):
            step = 
        

    def check_stay(grid,w):
        """
        check if walker walks or stays
        walker stays if current square is 0 and one neighbouring square is 1
        otherwise return false. 
        """
        if np.any([num_to_step(step) for i in range(1,5)
    grid = np.zeros((n_length,n_height))
    mid = np.array([(n_length-)/2, (n_height-1)/2])
    grid[mid[0],mid[1]] = 1

    while np.sum(grid < n_areas):
        hopper = np.array([mid[0],mid[1]])
        
    
    

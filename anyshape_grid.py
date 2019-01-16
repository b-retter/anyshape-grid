##Code to take in arbitrary shaped regions designated by 1s and 0s and calculates spatial statistics on them

import numpy as np
import numpy.random as rnd
import matplotlib.pyplot as plt


def num_to_step(num):
    """
    converts number 1 to 4 into a step in either x or y direction
           2
        1     3
           4
    """
    if num == 1:
        return np.array([-1,0])
    elif num == 2:
        return np.array([0,1])
    elif num == 3:
        return np.array([1,0])
    elif num == 4:
        return np.array([0,-1])

def is_wall(current,step,grid):
    """checks if grid + step exceeds grid.
    If so, return True
    """
    dim = np.shape(grid)
    stepped = current+step
    if np.any(stepped < 0) or stepped[0] > dim[0]-1 or stepped[1] > dim[1]-1:
        return True
    else:
        return False
        
def make_grid(n_areas,n_length,n_height):
    """
    generates an arbitrary shaped region through random walking starting in the centre
    """

    def b_condition(loc,step):
        if is_wall(loc,step,grid):
            # if it is a wall, add the length or height of the grid
            # making a torus
            return loc+step-np.array([n_length,n_height])*step
        else:
            return loc+step
        
    def random_walk(loc):
        """
        random walker direction
           2
        1     3
           4
        """
        step = num_to_step(rnd.randint(1,5))
        return b_condition(loc,step)

    def check_stay(w):
        """
        check if walker walks or stays
        walker stays if current square is 0 and one neighbouring square is 1
        otherwise return false. 
        """
        if np.sum([grid[b_condition(w,num_to_step(num))[0],b_condition(w,num_to_step(num))[1]] \
                   for num in range(1,5)]) > 0 and grid[w[0],w[1]] == 0:
            return True
        else:
            return False
        
    grid = np.zeros((n_length,n_height))
    mid = np.array([(n_length-1)/2, (n_height-1)/2])
    grid[mid[0],mid[1]] = 1

    while np.sum(grid) < n_areas:
        walker = np.array([mid[0],mid[1]])
        while check_stay(walker) == False:
            walker = random_walk(walker)
        grid[walker[0],walker[1]] = 1

    return grid

def circle(x0,y0,R,gx,gy,grid):
    shape = np.shape(grid)
    dx = gx[1]-gx[0]
    dy = gy[1]-gy[0]
    Rx = np.ceil(R/float(dx))
    Ry = np.ceil(R/float(dy))
    x0_i = (x0-gx[0])/dx-0.5
    y0_i = (y0-gy[0])/dy-0.5
        
    for i in range(shape[0]-Rx,shape[0]+R[y]):
        for j in range(np.shape(grid)[1]):
            grid[x,y] = np.sqrt((gx[i]-x0)**2 + (gy[j]-y0)**2) <= R
            

    return None



x_side = 10
y_side = 10
XMIN,XMAX = 0,10
YMIN,YMAX = 0,10

#central coordinates of the grid squares
gx = np.linspace(XMIN,XMAX,x_side,endpoint=False) + (XMAX-XMIN)/(2.0*x_side)
gy = np.linspace(YMIN,YMAX,y_side,endpoint=False) + (YMAX-YMIN)/(2.0*y_side)

grid = make_grid(50,x_side,y_side)
plt.pcolormesh(grid)
plt.show()

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

def index2x(index,axis,grid=None):
    """
    Returns the x (or y) value equal to the midpoint of
    the grid square indicated by index.
    """
    #allow a default value of grid to be the coverage map
    if grid == None:
        grid = coverage

    #check if index is within coverage map
    if axis == 'x':
        if index >= np.shape(grid)[0]:
            print('Index out of range for x')
            return None
        dx = x[1]-x[0]
        return index*dx + dx/2.0 + x[0]
    elif axis == 'y':
        if index >= np.shape(grid)[1]:
            print('Index out of range for y')
            return None
        dy = y[1]-y[0]
        return index*dy + dy/2.0 + y[0]
    else:
        print("axis designation required: 'x' or 'y'")
        return None

def x2index(v,axis):
    """
    Finds the grid index for a given point on x
    or y axis.
    """
    
    if axis == 'x':
        if v < x[0] or v > max(x):
            print('x value is outside of range')
            return None
        elif v == x[0]:
            return 0
        else:
            dx = x[1]-x[0]
            return int(round((v-x[0])/float(dx)-0.5))
    elif axis == 'y':
        if v < y[0] or v > max(y):
            print('y value is outside of range')
            return None
        if v == y[0]:
            return 0
        else:
            dy = y[1]-y[0]
            return int(round((v-y[0])/float(dy)-0.5))
    else:
        print('axis designation required')
        return None

def xy2grid(v1,v2):
    return x2index(v1,'x'),x2index(v2,'y')

def delDist2Grid(v2,v1=0,axis=None):
    if axis == 'x':
        return int((v2-v1)/float(dx))
    elif axis == 'y':
        return int((v2-v1)/float(dy))
    
def delGrid2Dist(i2,i1=0,axis=None):
    if axis == 'x':
        return (i2-i1)*dx
    if axis == 'y':
        return (i2-i1)*dy
    
def circle(xp,yp,R,grid=None):
    """
    Finds all the grid squares that are in a circle around
    xp,yp with radius R. 
    """
    #allow a default value of grid to be the coverage map
    if grid == None:
        grid = coverage

    if xp < x[0] or xp > max(x):
        print('xp outside of range')
    elif yp < y[0] or yp > max(y):
        print('yp outside of range')
        
    Rx = delDist2Grid(R,axis='x')
    Ry = delDist2Grid(R,axis='y')
    xg,yg = xy2grid(xp,yp)

    shape = np.shape(grid)
    
    x_min,x_max = max(xg-Rx,0), min(xg+Rx+1,shape[0])
    y_min,y_max = max(yg-Ry,0), min(yg+Ry+1,shape[1])

    coords = [[],[]]
    for i in range(x_min,x_max):
        for j in range(y_min,y_max):
            if grid[i,j] == 1 and np.sqrt((gx[i]-xp)**2 + (gy[j]-yp)**2) <= R:
                coords[0].append(i)
                coords[1].append(j)
            

    return coords

def yso_to_grid(yso,grid=None):
    """
    Make a new grid to place YSOs into using grid as a mask
    and basis of next grid.
    yso should by a 2xN array containing x and y values.
    """
    #allow a default value of grid to be the coverage map
    if grid == None:
        grid = coverage

    yso_map = np.zeros(np.shape(grid))
        
    N = np.shape(yso)[1]    
    for i in range(N):
        yso_map[xy2grid(yso[0,i],yso[1,i])] += 1

    return yso_map


    
x_side = 50
y_side = 50
XMIN,XMAX = 0,10
YMIN,YMAX = 0,10
dx = (XMAX-XMIN)/float(x_side)
dy = (YMAX-YMIN)/float(y_side)


#central coordinates of the grid squares
x = np.arange(XMIN,XMAX+dx,dx)
y = np.arange(YMIN,YMAX+dy,dy)
gx = np.linspace(XMIN,XMAX,x_side,endpoint=False) + (XMAX-XMIN)/(2.0*x_side)
gy = np.linspace(YMIN,YMAX,y_side,endpoint=False) + (YMAX-YMIN)/(2.0*y_side)

Nyso = 500
yso = np.vstack(((rnd.rand(Nyso)*(XMAX-XMIN)+XMIN),(rnd.rand(Nyso)*(YMAX-YMIN)+YMIN)))

#grid = make_grid(50,x_side,y_side)
coverage = np.ones((x_side,y_side))
coords = circle(5,5,4)

yso_map = yso_to_grid(yso)

plt.pcolormesh(yso_map)
plt.show()

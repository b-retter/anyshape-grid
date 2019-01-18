##Code to take in arbitrary shaped regions designated by 1s and 0s and calculates spatial statistics on them

import numpy as np
import numpy.random as rnd
import matplotlib.pyplot as plt
import sys
sys.path.append('/Users/bretter/Documents/StarFormation/RandomDistribution/spatialStats/Functions')
import allstats as alls

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
            

    return np.array(coords)

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

def random_ysos(val,mode='binomial',grid=None):
    """
    Function to populate a grid with random YSOs. YSOs can be placed anywhere
    with grid == 1.
    Two modes:
    If mode is 'binomial' randomly distribute val YSOs around the region.
    If mode is 'csr', place Poisson((val/study area)*pixel area) ysos in each pixel.

    Both return the coordinates of the ysos and the completed grid. 
    """

    if grid == None:
        grid = coverage

    shape = np.shape(grid)
    inside_pixels = np.array(np.where(grid == 1))
    n_pixels = np.shape(inside_pixels)[1]
    
    yso_map = np.zeros(shape)
    yso = [[],[]]
    if mode == 'csr':
        lmda = val/float((XMAX-XMIN)*(YMAX-YMIN))
        pixel_area = dx*dy
        for pixel in range(n_pixels):
            i,j = inside_pixels[0,pixel], inside_pixels[1,pixel]
            yso_map[i,j] = rnd.poisson(lmda*pixel_area)
        return yso, yso_map
    elif mode == 'binomial':
        while np.sum(yso_map) < val:
            rand_pixel = rnd.randint(0,n_pixels)
            i,j = inside_pixels[0,rand_pixel], inside_pixels[1,rand_pixel]
            yso_map[i,j] += 1
        return yso, yso_map

def kfunc(x,y,t,yso_map=None,grid=None):
    """
    Calculates K function for points with coords x,y.
    Most likely x and y are the positions of yso.
    """

    if yso_map == None:
        yso_map = yso_to_grid(np.array([x,y]),grid)
        
    if grid == None:
        grid = coverage

    #yso will each count themselves once over the course of the 
    #algorithm. This accounts for the self-counting.
    yso_sum = -np.sum(yso_map)
    area_sum = 0
    for i in range(len(x)):
        coords = circle(x[i],y[i],t,grid)
        n_coords = np.shape(coords)[1]
        area_sum += n_coords
        
        for j in range(n_coords):
            yso_sum+=yso_map[coords[0,j],coords[1,j]]
    
    area = dx*dy*area_sum
    lmda = np.sum(yso_map)/float(np.sum(grid)*dx*dy)
    K = np.pi*t**2/lmda*yso_sum/float(area)
    L = np.sqrt(K/np.pi) - t
    return K, L

def Oring(x,y,t,w,yso_map=None,grid=None):
    """
    Calculates Oring function for points with coords x,y.
    Most likely x and y are the positions of yso.
    """

    if yso_map == None:
        yso_map = yso_to_grid(np.array([x,y]),grid)
        
    if grid == None:
        grid = coverage
        
    yso_sum = 0
    area_sum = 0
    for i in range(len(x)):
        coords = ring(x[i],y[i],t,w,grid)
        n_coords = np.shape(coords)[1]
        area_sum += n_coords

        xg,yg = xy2grid(x[i],y[i])
        for j in range(n_coords):
            if coords[0,j]==xg and coords[1,j]==yg:
                continue
            yso_sum+=yso_map[coords[0,j],coords[1,j]]
    
    area = dx*dy*area_sum
    lmda = np.sum(yso_map)/float(np.sum(grid)*dx*dy)
    O = yso_sum/float(area)
    return O, O/lmda

def ring(xp,yp,R,w,grid=None):
    """
    Finds all the grid squares that are in an annulus around
    xp,yp with radius R and width w. 
    """
    #allow a default value of grid to be the coverage map
    if grid == None:
        grid = coverage

    Rout = R+w/2.0
    Rin = R-w/2.0
    if xp < x[0] or xp > max(x):
        print('xp outside of range')
    elif yp < y[0] or yp > max(y):
        print('yp outside of range')
        
    Rx = delDist2Grid(Rout,axis='x')
    Ry = delDist2Grid(Rout,axis='y')
    xg,yg = xy2grid(xp,yp)

    shape = np.shape(grid)
    
    x_min,x_max = max(xg-Rx,0), min(xg+Rx+1,shape[0])
    y_min,y_max = max(yg-Ry,0), min(yg+Ry+1,shape[1])

    coords = [[],[]]
    for i in range(x_min,x_max):
        for j in range(y_min,y_max):
            d = np.sqrt((gx[i]-xp)**2 + (gy[j]-yp)**2)
            if grid[i,j] == 1 and d <= Rout and d >= Rin:
                coords[0].append(i)
                coords[1].append(j)            

    return np.array(coords)

x_side = 100
y_side = 100
XMIN,XMAX = 0,30
YMIN,YMAX = 0,30
AREA = (XMAX-XMIN)*(YMAX-YMIN)
dx = (XMAX-XMIN)/float(x_side)
dy = (YMAX-YMIN)/float(y_side)

bounds = np.array([[XMIN,XMAX],[YMIN,YMAX]])

#central coordinates of the grid squares
x = np.arange(XMIN,XMAX+dx,dx)
y = np.arange(YMIN,YMAX+dy,dy)
gx = np.linspace(XMIN,XMAX,x_side,endpoint=False) + (XMAX-XMIN)/(2.0*x_side)
gy = np.linspace(YMIN,YMAX,y_side,endpoint=False) + (YMAX-YMIN)/(2.0*y_side)

Nyso = 200
coverage = np.ones((x_side,y_side))

x0,y0 = 15,15
R = 10
circ = circle(x0,y0,R)
coverage = np.zeros((x_side,y_side))
for i in range(np.shape(circ)[1]):
    coverage[circ[0,i],circ[1,i]] = 1

N = 0
yso = [[],[]]
while N < Nyso:
    xx = (rnd.rand(Nyso-N)-0.5)*2*R
    yy = (rnd.rand(Nyso-N)-0.5)*2*R
    for i in range(Nyso-N):
        d = np.sqrt(xx[i]**2 + yy[i]**2)
        if d <= R:
            yso[0].append(xx[i]+x0)
            yso[1].append(yy[i]+y0)
            N+=1

yso = np.array(yso)
yso_map = yso_to_grid(yso)

step = 20
r = np.linspace(1,15,step)
h = 0.5

O1,L1 = [], []
O2,L2 = [], []
for i,t in enumerate(r):
    w = h
    o,oo = Oring(yso[0,:],yso[1,:],t,w,yso_map=None,grid=None)
    O1.append(oo)
    k,kk = kfunc(yso[0,:],yso[1,:],t,yso_map=None,grid=None)
    L1.append(kk)

    o,oo = alls.Oring(yso[0,:],yso[1,:],t,w,AREA,bounds)
    O2.append(oo)
    k,kk = alls.kfunc(yso[0,:],yso[1,:],t,AREA,bounds)
    L2.append(kk)


plt.figure()
plt.plot(r,O1,'r')
plt.plot(r,O2,'b')

plt.figure()
plt.plot(r,L1,'r')
plt.plot(r,L2,'b')

plt.show()
## yso_map = random_ysos(Nyso,mode='binomial'):

## yso = np.vstack(((rnd.rand(Nyso)*(XMAX-XMIN)+XMIN),(rnd.rand(Nyso)*(YMAX-YMIN)+YMIN)))
## #grid = make_grid(50,x_side,y_side)
## coverage = np.ones((x_side,y_side))


## t = 5
## w = 1
## K1,L1 = Oring(yso[0,:],yso[1,:],t,w,yso_map=None,grid=None)
## K2,L2 = alls.Oring(yso[0,:],yso[1,:],t,w,AREA,bounds)

## print('O1 = {:.4f}, O/lambda = {:.4f}'.format(K1,L1))
## print('O1 = {:.4f}, O/lambda = {:.4f}'.format(K2,L2))

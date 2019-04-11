##Code to take in arbitrary shaped regions designated by 1s and 0s and calculates spatial statistics on them

import numpy as np
import numpy.random as rnd
import matplotlib.pyplot as plt
import sys
import os
from astropy.io import fits
from astropy import wcs

#/Users/bretter/Documents/StarFormation/RandomDistribution/spatialStats/Functions
sys.path.append('/Users/bretter/Documents/StarFormation/RandomDistribution/spatialStats/Functions')
import allstats as alls
from timeit import default_timer as timer

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
    current = array([x,y])
    step = array([x,y])
    """
    dim = np.shape(grid)
    stepped = current+step
    if np.any(stepped < 0) or stepped[0] > dim[0]-1 or stepped[1] > dim[1]-1:
        return True
    else:
        return False
def box_check(xg,yg,Lx,Ly=None,grid=None):
    """
    Checks if there are any cells not covered by coverage map, or
    exceed the bounds of the coverage map within box of side L,
    centered on xg,yg.
    Returns False if not entirely within map.
    """
    if grid is None:
        grid = coverage
    if Ly is None:
        Ly = Lx
        
    #Check if bounds nearby
    steps = map(num_to_step,range(1,5))
    xy = np.array([xg,yg])
    if is_wall(xy,steps[0]*(Lx/2),grid) or is_wall(xy,steps[2]*(Lx/2),grid):
        return False
    if is_wall(xy,steps[1]*(Ly/2),grid) or is_wall(xy,steps[3]*(Ly/2),grid):
        return False

    #Check if any cells not in coverage map
    x_min,x_max = xg-Lx/2, xg+Lx/2+1
    y_min,y_max = yg-Ly/2, yg+Ly/2+1
    if np.any(grid[x_min:x_max,y_min:y_max] == 0):
        return False
    else:
        return True

def circle_check(xp,yp,R,grid=None):
    """
    Checks if there are any cells not covered by coverage map, or
    exceed the bounds of the coverage map within circle of radius R,
    centered on xp,yp.
    Returns False if not entirely within map.
    """
    if grid is None:
        grid = coverage

    xg,yg = xy2grid(xp,yp)
    Rx = delDist2Grid(R,axis='x')
    Ry = delDist2Grid(R,axis='y')

    #Check if bounds nearby
    steps = map(num_to_step,range(1,5))
    xy = np.array([xg,yg])
    if is_wall(xy,steps[0]*(Rx),grid) or is_wall(xy,steps[2]*(Rx),grid):
        return False
    if is_wall(xy,steps[1]*(Ry),grid) or is_wall(xy,steps[3]*(Ry),grid):
        return False
    

    GX, GY = np.meshgrid(gx,gy,indexing='ij')
    dists = np.sqrt((GX-xp)**2 + (GY-yp)**2)
    co_x,co_y = np.where((dists <= R) & (grid == 0))
    if len(co_x) == 0:
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

def xy2grid(v1,v2,wcs_obj=None):
    """
    convert world coordinates (v1,v2) to grid coordinates
    (i,j).
    """
    if wcs_obj is None:
        wcs_obj = w_obj

    if inverted == True:
        j,i = wcs_obj.all_world2pix(v2,v1,0)
    else:
        i,j = wcs_obj.all_world2pix(v1,v2,0)
        
    return int(np.round(i)),int(np.round(j))

def ij2xy(i,j,wcs_obj=None):
    """
    convert grid coordinates (i,j) to world coordinates.
    """
    
    if wcs_obj is None:
        wcs_obj = w_obj

    if inverted == True:
        y,x = wcs_obj.all_pix2world(j,i,0)
    else:
        x,y = wcs_obj.all_pix2world(i,j,0)
    return x,y

def delDist2Grid(v2,v1=0,axis=None):
    if axis == 'x':
        return int((v2-v1)/float(dx))
    elif axis == 'y':
        return int((v2-v1)/float(dy))

def inside_check(v1,v2,wcs_obj=None):
    """
    Checks if coordinates (RA,Dec) are inside of coverage
    map given by wcs_obj.
    Returns False if not inside coverage map.
    """
    if wcs_obj is None:
        wcs_obj = w_obj
    if inverted == True:
        cdec,cra = wcs_obj.all_world2pix(v2,v1,0)
    else:
        cra,cdec = wcs_obj.all_world2pix(v1,v2,0)
    if  0 <= cra < ra_axis and 0 <= cdec < dec_axis:
        return True
    else:
        return False
    
def circle(xp,yp,R,grid=None,relative=False):
    """
    Finds all the grid squares that are in a circle around
    xp,yp with radius R.
    If relative is true, provide relative differences in
    grid coords between xp,yp and circle cells.
    Otherwise provide the absolute references.
    """
    #allow a default value of grid to be the coverage map
    if grid is None:
        grid = coverage

    if not inside_check(xp,yp):
        print('world coordinate outside coverage map')
        
    #reduce distance search to more immediate values
    il,ir,jl,jr = angle2box(xp,yp,R)
    dists = gcircle((gx[il:ir,jl:jr],gy[il:ir,jl:jr]),(xp,yp))
    co_x,co_y = np.where((dists <= R) & (grid[il:ir,jl:jr] == 1))
    if relative == False:
        return np.array([co_x+il,co_y+jl])
    elif relative == True:
        return np.array([co_x,co_y])-np.array([xg,yg]).reshape(2,1)
    
def yso_to_grid(yso,grid=None,yso_return=False):
    """
    Make a new grid to place YSOs into using grid as a mask
    and basis of next grid.
    yso should by a 2xN array containing x and y values.
    Optional yso_return function. Returns yso coordinates
    that were inside coverage map.
    """
    #Check yso_return is boolean
    if not type(yso_return) == bool:
        print('yso_return must be a boolean')
        return None
    
    #allow a default value of grid to be the coverage map
    if grid is None:
        grid = coverage

    yso_map = np.zeros(np.shape(grid))
    filtered_ysos = [[],[]]
    N = np.shape(yso)[1]
    fail_count = 0
    for i in range(N):
        x,y = xy2grid(yso[0,i],yso[1,i])
        if x is None or y is None:
             fail_count += 1
             continue

        if yso_return == True:
            filtered_ysos[0].append(yso[0,i])
            filtered_ysos[1].append(yso[1,i])
        yso_map[xy2grid(yso[0,i],yso[1,i])] += 1

    if fail_count > 0:
        print('{:d} YSOs failed to position'.format(fail_count))

    if yso_return == True:
        return yso_map*grid, np.array(filtered_ysos)
    elif yso_return == False:
        return yso_map*grid
    
def random_ysos(val,mode='binomial',grid=None):
    """
    Function to populate a grid with random YSOs. YSOs can be placed anywhere
    with grid == 1.
    Two modes:
    If mode is 'binomial' randomly distribute val YSOs around the region.
    If mode is 'csr', place Poisson((val/study area)*pixel area) ysos in each pixel.

    Both return the coordinates of the ysos and the completed grid.
    Assumes a uniform probability across the coverage map.
    """

    if grid is None:
        grid = coverage

    shape = np.shape(grid)
    inside_pixels = np.array(np.where(grid == 1))
    n_pixels = np.shape(inside_pixels)[1]
    
    yso_map = np.zeros(shape)
    yso_x = []
    yso_y = []
    if mode == 'csr':
        total_area = get_area()
        lmda = val/total_area
        for pixel in range(n_pixels):
            i,j = inside_pixels[0,pixel], inside_pixels[1,pixel]
            Nyso = rnd.poisson(lmda*area_array[i,j])
            yso_map[i,j] = Nyso

            if Nyso > 0:
                if inverted:
                    y,x = w_obj.all_pix2world(rnd.rand(Nyso)+j,rnd.rand(Nyso)+i,0)
                else:
                    x,y = w_obj.all_pix2world(rnd.rand(Nyso)+i,rnd.rand(Nyso)+j,0)

                yso_x = np.append(yso_x,x)
                yso_y = np.append(yso_y,y)
            
        return np.vstack([yso_x,yso_y]), yso_map
    
    elif mode == 'binomial':
        while np.sum(yso_map) < val:
            rand_pixel = rnd.randint(0,n_pixels)
            i,j = inside_pixels[0,rand_pixel], inside_pixels[1,rand_pixel]
            yso_map[i,j] += 1
            if inverted:
                y,x = w_obj.all_pix2world(rnd.rand()+j,rnd.rand()+i,0)
            else:
                x,y = w_obj.all_pix2world(rnd.rand()+i,rnd.rand()+j,0)
            yso[0].append(x)
            yso[1].append(y)
            
        return np.array(yso), yso_map

def kfunc(x,y,t,yso_map=None,grid=None,opti=False,diag=False):
    """
    Calculates K function for points with coords x,y.
    Most likely x and y are the positions of yso.
    """

    if yso_map is None:
        yso_map = yso_to_grid(np.array([x,y]),grid)
        
    if grid is None:
        grid = coverage

    #yso will each count themselves once over the course of the 
    #algorithm. This accounts for the self-counting.
    yso_sum = -np.sum(yso_map)
    area = 0

    #Generate relative circle coords for approx centre of map
    shape = np.shape(yso_map)
    x_mid,y_mid = ij2xy(shape[0]/2,shape[1]/2)
    
    for i in range(len(x)):

        coords = circle(x[i],y[i],t,grid)
                
        n_coords = np.shape(coords)[1]
        
        for j in range(n_coords):
            yso_sum+=yso_map[coords[0,j],coords[1,j]]
            area += area_array[coords[0,j],coords[1,j]]

    total_area = get_area()
    lmda = np.sum(yso_map)/total_area
    K = (np.pi*t**2/lmda)*yso_sum/float(area)
    L = np.sqrt(K/np.pi) - t
    if diag == True:
        return K,L,yso_sum,float(area)
    else:
        return K, L

def Oring(x,y,t,w,yso_map=None,grid=None,opti=False,diag=False):
    """
    Calculates Oring function for points with coords x,y.
    Most likely x and y are the positions of ysos.
    """

    if yso_map is None:
        yso_map = yso_to_grid(np.array([x,y]),grid)
        
    if grid is None:
        grid = coverage

    #Generate relative circle coords for approx centre of map
    shape = np.shape(yso_map)
    x_mid,y_mid = ij2xy(shape[0]/2,shape[1]/2)
    
    yso_sum = 0
    area = 0
    for i in range(len(x)):
        self_count = False

        coords = ring(x[i],y[i],t,w,grid)
        n_coords = np.shape(coords)[1]

        xg,yg = xy2grid(x[i],y[i])
        for j in range(n_coords):

            #Allow o-ring to skip one point within its own grid square.
            if coords[0,j]==xg and coords[1,j]==yg and self_count == False:
                yso_sum-=1

            area += area_array[coords[0,j],coords[1,j]]
            yso_sum+=yso_map[coords[0,j],coords[1,j]]

    total_area = get_area()
    lmda = np.sum(yso_map)/total_area
    O = yso_sum/float(area)
    if diag == True:
        return O, O/lmda, yso_sum, float(area)
    else:
        return O, O/lmda
    
def ring(xp,yp,R,w,grid=None,relative=False):
    """
    Finds all the grid squares that are in an annulus around
    xp,yp with radius R and width w.
    If relative is True, provide relative differences in
    grid coords between xp,yp and circle cells.
    Otherwise provide the absolute references.
    """
    #allow a default value of grid to be the coverage map
    if grid is None:
        grid = coverage

    Rout = R+w/2.0
    Rin = R-w/2.0
    
    if not inside_check(xp,yp):
        print('world coordinate outside coverage map')

    #reduce distance search to more immediate values
    il,ir,jl,jr = angle2box(xp,yp,Rout)
    dists = gcircle((gx[il:ir,jl:jr],gy[il:ir,jl:jr]),(xp,yp))
    co_x,co_y = np.where((dists <= Rout) & (dists >= Rin) & (grid[il:ir,jl:jr] == 1))
    if relative == False:
        return np.array([co_x+il,co_y+jl])
    elif relative == True:
        return np.array([co_x,co_y])-np.array([xg,yg]).reshape(2,1)

def gfunc(x,y,t,yso_map=None,grid=None):
    """
    Calculates Diggle's G function by calculating how many events
    have a nearest neighbour within distance t.
    Applying the border method of edge correction.
    """
    print('gfunc no longer functions in this version. Requires rewrite.')
    return None
    
    if grid is None:
        grid = coverage
        
    if yso_map is None:
        yso_map,xy = yso_to_grid(np.array([x,y]),grid,yso_return = True)
        x,y = xy[0,:],xy[1,:]
        
    x = np.copy(x)
    y = np.copy(y)

    N = np.size(x)

    X,Y = np.meshgrid(x,y)
    x.resize((len(x),1))

    #collect nearest neighbour distances
    dists = np.sqrt((X-x)**2 + (Y-y)**2)
    dists[dists == 0] = np.max(dists) 
    nearest = np.min(dists,axis=0)

    #count number of points inside boundary with nn <= t
    n_inside = 0
    n_nearest = 0
    for i in range(N):
        if circle_check(x[i],y[i],t,grid=None):
            n_inside += 1
            if nearest[i] <= t:
                n_nearest += 1

    #if all points are excluded use backup
    if n_inside == 0:
        gsum = np.sum(nearest <= t)
        G = gsum/float(N)
    else:
        gsum = n_nearest
        G = gsum/float(n_inside)

    area = get_area(grid)
    lmda = N/area
    E = 1 - np.exp(-lmda*np.pi*t**2)
    return G,E

def ffunc(x,y,t,a=None,b=None,yso_map=None,grid=None):
    """
    Calculates free space function by calculating how many events
    have a point-event nearest neighbour within distance t.
    Applying the border method of edge correction.
    a and b are optional arguments for the positions used to
    probe the space.
    If not provided a and b will be randomly distributed over the
    map. 
    """
    print('ffunc no longer functions in this version. Requires rewrite.')
    return None
    if grid is None:
        grid = coverage
    
    if yso_map is None:
        yso_map, xy = yso_to_grid(np.array([x,y]),grid,yso_return=True)
        x,y = xy[0,:],xy[1,:]

    #deal with random positions
    if bool(a is None) != bool(b is None):
        print('Provide either a and b, or neither')
        return None
    elif a is None and b is None:
        #if no values given for a and b, generate N randomly
        #spaced ysos across effective area of coverage map.
        val = len(x)
        ab, ab_map = random_ysos(val,'binomial',grid)
        a,b = ab[0,:],ab[1,:]
    else:
        #len(a) must be equal to len(x)
        if len(x) != len(a):
            print('Currently x, y, a and b must all be of equal length.')
            return None
        else:
            #copy a and b, filter by coverage map and return
            a,b = np.copy(a),np.copy(b)
            ab_map, ab = yso_to_grid(np.array([a,b]),grid,yso_return=True)
            a,b = ab[0,:],ab[1,:]
        
    x = np.copy(x)
    y = np.copy(y)

    N = np.size(x)
    
    A,B = np.meshgrid(a,b)
    B = B.transpose()

    x.resize((np.size(x),1))
    y.resize((np.size(y),1))

    dists = np.sqrt((A-x)**2 + (B-y)**2)
    nearest = np.min(dists,axis=0)

    a.resize((len(a),1))
    b.resize((len(b),1))

    #count number of points inside boundary with nn <= t
    n_inside = 0
    n_nearest = 0
    for i in range(N):
        if circle_check(a[i],b[i],t,grid=None):
            n_inside += 1
            if nearest[i] <= t:
                n_nearest += 1

    #if all points are excluded use backup
    if n_inside == 0:
        fsum = np.sum(nearest <= t)
        F = fsum/float(N)
    else:
        fsum = n_nearest
        F = fsum/float(n_inside)

    area = get_area(grid)
    lmda = N/area
    E = 1 - np.exp(-lmda*np.pi*t**2)
    return F,E

def get_area(grid = None):
    """
    Return the effective area of the coverage map
    as float.
    """
    if grid is None:
        grid = coverage
        
    return float(np.sum(area_array*grid))

def gcircle(p1,p2):
    """
    Returns the central angle between two angular
    coordinates.
    p1 and p2 must be tuples of form (ra,dec).
    """
    #conversion from degrees to radians.
    d2r = lambda x: x*np.pi/180.0
    x0,y0 = d2r(p1[0]),d2r(p1[1])
    x1,y1 = d2r(p2[0]),d2r(p2[1])

    #Using Vincenty Great Circle Formula
    dra = np.abs(x0-x1)
    s1 = np.sqrt((np.cos(y1)*np.sin(dra))**2+(np.cos(y0)*np.sin(y1)-np.sin(y0)*np.cos(y1)*np.cos(dra))**2)
    s2 = np.sin(y0)*np.sin(y1)+np.cos(y0)*np.cos(y1)*np.cos(dra)
    sep = np.arctan(s1/s2)
    return sep*180/np.pi

def get_area_array(tan=True):
    """
    Return the celestial pixel areas for each pixel
    in the FITS file.
    If tan projection, calculates celestial area by calculating
    great circle angle to each pixel and multiplying reference 
    pixel angle by plate scale cos**3(theta). See wikipedia
    gnomonic projection. 

    Else, assumes pixels are small and can be approximated 
    by rectangles rotated with respect to lines of const ra.
    """
    if tan:
        #If tan projection use da_sphere = cos**3(theta)*da_plane
        angles = gcircle((gx,gy),(ra_ref,dec_ref))
        ref_area = wcs.utils.proj_plane_pixel_area(w_obj)
        return np.cos(d2r(angles))**3*ref_area
    else:
        #Find RA and Dec at each grid coordinate
        gx = np.arange(ra_axis+1)
        gy = np.arange(dec_axis+1)
        GX,GY = np.meshgrid(gx,gy,indexing='ij')
        GX,GY = GX.flatten(), GY.flatten()
        gy,gx = w_obj.all_pix2world(GY,GX,0)
        gx, gy = gx.reshape(ra_axis+1,dec_axis+1), gy.reshape(ra_axis+1,dec_axis+1)
        
        #Find dRA_ij = RA_{i+1},j - RA_ij and likewise for dec
        dRA = gx[1:ra_axis+1,:dec_axis]-gx[:ra_axis,:dec_axis]
        dDec = gy[:ra_axis,1:dec_axis+1]-gy[:ra_axis,:dec_axis]
        
        #Due to rotation a step in grid changes both ra and dec.
        #Calculate angle between lines of constant Ra and the
        #grid at constant i.
        dely = gy[1:ra_axis+1,:dec_axis]-gy[:ra_axis,:dec_axis]
        theta = np.arctan(dely/dRA)
        
        angle_part = np.sin(np.pi/2-gy[:ra_axis,:dec_axis]*np.pi/180.0)
        #celestial steradians for all pixels
        return angle_part*(dRA*dDec)/(np.cos(theta)**2)

def angle2box(xp,yp,t):
    """
    Calculate a VERY conservative estimate of the grid squares
    that cover an angular distance t.
    Returns the i and j values of bottom left and top right
    corners.
    """
    t = np.sqrt(2)*t
    if inverted:
        j,i = w_obj.all_world2pix(yp,xp,0)
        
        ##translate angle t in both ra and dec
        jl,il = w_obj.all_world2pix(yp-t,xp-t,0)
        jr,ir = w_obj.all_world2pix(yp+t,xp+t,0)
        
    else:
        i,j = w_obj.all_world2pix(xp,yp,0)

        ##translate angle t in both ra and dec
        il,jl = w_obj.all_world2pix(xp-t,yp-t,0)
        ir,jr = w_obj.all_world2pix(xp+t,yp+t,0)

    if jl < 0:
        jl = 0
    if il < 0:
        il = 0
    if ir > ra_axis:
        ir = ra_axis
    if jr > dec_axis:
        jr = dec_axis

    return int(il),int(ir),int(jl),int(jr)

"""
Extract the relevant data from the fits file. 
Array size.
Make a wcs object.
Array of grid square coordinates.
Array of grid centre coordinates.
Coverage map.
"""

fits_path = '/Users/bretter/Documents/StarFormation/SFR_data'
#fits_path = '../SFR_data'
fits_name = 'SERAQU_IRAC1234M1_cov_sm.fits'
coverage,header = fits.getdata(os.path.join(fits_path,fits_name), header=True)
w_obj = wcs.WCS(header)
##Find which axis is RA and which is Dec.
##inverted kayword defines if axes are RA and Dec or
##Dec and RA.
if 'DEC' in header['CTYPE1']:
    inverted = True
    dec_axis = header['NAXIS1']
    ra_axis = header['NAXIS2']
else:
    inverted = False
    coverage = coverage.T
    ra_axis = header['NAXIS1']
    dec_axis = header['NAXIS2']

##Getting celestial coordinates of pixel centres
gx = np.arange(0,ra_axis)
gy = np.arange(0,dec_axis)
GX,GY = np.meshgrid(gx,gy,indexing='ij')
GX,GY = GX.flatten(), GY.flatten()
gy,gx = w_obj.all_pix2world(GY,GX,0)
gx, gy = gx.reshape(ra_axis,dec_axis), gy.reshape(ra_axis,dec_axis)

##Clear GX and GY from memory
GX, GY = None, None

cov2 = np.zeros(np.shape(coverage))
cov2 += coverage == 1

coverage = cov2

##Getting pixel scales
area_array = get_area_array()
total_area = np.sum(area_array)
val = 200
print(np.shape(coverage))
yso, yso_map = random_ysos(val,mode='csr',grid=coverage)

steps = 20
r = np.linspace(0.1,2,steps)
w = 0.3*r

results = np.empty((2,steps))
for i,t in enumerate(r):
    o1, o2 = Oring(yso[0,:],yso[1,:],t,w[i],yso_map=None,grid=None,opti=False,diag=False)
    k1, k2 = kfunc(yso[0,:],yso[1,:],t)
    results[0,i] = o2
    results[1,i] = k2


fpath = '/Users/bretter/Documents/StarFormation/Meetings/04-04-2019/'

plt.pcolormesh(gx,gy,coverage)
plt.plot(yso[0,:],yso[1,:],'*')
plt.xlabel('RA')
plt.ylabel('Dec')
plt.title('200 YSOs randomly distibuted in coverage map')

fname = 'yso_coverage.png'
plt.savefig(fpath+fname)

plt.plot(r,results[0,:])
plt.xlabel('r (angle)')
plt.ylabel(r'$O/\lambda$')
plt.title('O-ring for YSOs randomly distibuted in coverage map')

fname = 'oring.png'
plt.savefig(fpath+fname)

plt.plot(r,results[1,:])
plt.xlabel('r (angle)')
plt.ylabel(r'$L$')
plt.title("Ripley's K for YSOs randomly distibuted in coverage map")

fname = 'oring.png'
plt.savefig(fpath+fname)

fname='results'
np.save(fpath+fname)


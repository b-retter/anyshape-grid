##Code to take in arbitrary shaped regions designated by 1s and 0s and calculates spatial statistics on them

import numpy as np
import numpy.random as rnd
import matplotlib.pyplot as plt
import multiprocessing as mp
import time
import sys
import os
from astropy.io import fits
from astropy import wcs
from scipy.special import factorial

#/Users/bretter/Documents/StarFormation/RandomDistribution/spatialStats/Functions
#sys.path.append('/Users/bretter/Documents/StarFormation/RandomDistribution/spatialStats/Functions')
#import allstats as alls
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
    if  0 <= round(cra) < ra_axis and 0 <= round(cdec) < dec_axis:
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
                    y,x = w_obj.all_pix2world(rnd.rand(Nyso)+j-0.5,rnd.rand(Nyso)+i-0.5,0)
                else:
                    x,y = w_obj.all_pix2world(rnd.rand(Nyso)+i-0.5,rnd.rand(Nyso)+j-0.5,0)

                yso_x = np.append(yso_x,x)
                yso_y = np.append(yso_y,y)
            
        return np.vstack([yso_x,yso_y]), yso_map
    
    if mode == 'binomial':
        while np.sum(yso_map) < val:
            rand_pixel = rnd.randint(0,n_pixels)
            i,j = inside_pixels[0,rand_pixel], inside_pixels[1,rand_pixel]
            yso_map[i,j] += 1
            if inverted:
                y,x = w_obj.all_pix2world(rnd.rand()+j-0.5,rnd.rand()+i-0.5,0)
            else:
                x,y = w_obj.all_pix2world(rnd.rand()+i-0.5,rnd.rand()+j-0.5,0)
            yso_x = np.append(yso_x,x)
            yso_y = np.append(yso_y,y)
            
        return np.vstack([yso_x,yso_y]), yso_map

    elif mode == 'sphere_binomial':    
        #Use sphere point picking to place ysos into coverage map using angles
        while np.sum(yso_map) < val:
            #Repeat until desired number of ysos have been added to the map
            d2r = lambda x: x*np.pi/180
            n_yso = val-int(np.sum(yso_map))
            x = (np.max(gx)-np.min(gx))*rnd.rand(val)+np.min(gx)
            
            r_theta = rnd.rand(val-int(np.sum(yso_map)))
            vmin = 0.5*(np.cos(np.pi/2-d2r(np.min(gy)))+1)
            vmax = 0.5*(np.cos(np.pi/2-d2r(np.max(gy)))+1)
            del_v = vmax-vmin
            V = del_v*r_theta+vmin
            y = 90-180/np.pi*np.arccos(2*V-1)

            for i in range(n_yso):
                if inside_check(x[i],y[i]):
                    if grid[xy2grid(x[i],y[i])] == 1:
                        yso_map[xy2grid(x[i],y[i])] += 1
                        yso_x = np.append(yso_x,x[i])
                        yso_y = np.append(yso_y,y[i])

        return np.vstack([yso_x,yso_y]), yso_map

def one_kfunc(x1,y1,t,yso_map,grid):
    """
    Returns the number of ysos and area in circle around 
    a given yso at position x1,y1.
    """
    coords = circle(x1,y1,t,grid)
                
    n_coords = np.shape(coords)[1]

    #subtract one from yso_sum to stop self-counting
    yso_sum = -1
    area = 0
    for j in range(n_coords):
        yso_sum+=yso_map[coords[0,j],coords[1,j]]
        area += area_array[coords[0,j],coords[1,j]]
    
    return np.array([area,yso_sum])

def kfunc(x,y,t,yso_map=None,grid=None,noP=None,diag=False):
    """
    Calculates K function for points with coords x,y.
    Most likely x and y are the positions of yso.
    noP determines the number of workers assigned to function for 
    multiprocessing. A value of None or 1 runs the function 
    without multiprocessing.
    """
    #if values not specified take global values.
    if yso_map is None:
        yso_map = yso_to_grid(np.array([x,y]),grid=grid)
        
    if grid is None:
        grid = coverage
        
    if noP > 1:
        ##Initialise pool of workers
        pool = mp.Pool(noP)

        ##Perform one_kfunc for each yso using workers
        results = []
        for i in range(len(x)):
            results.append(pool.apply_async(one_kfunc,(x[i],y[i],t,yso_map,grid)))

        #Close down workers to finish calculating results
        pool.close()
        pool.join()

        #collate results
        finished_results = np.empty([len(x),2])
        for i in range(len(x)):
            finished_results[i,:] = results[i].get()
    elif noP == 1 or noP == None:
        #If there are fewer than 2 workers specified, then ignore
        #multiprocessing.
        
        ##Perform one_kfunc for each yso
        finished_results = np.empty([len(x),2])
        for i in range(len(x)):
            finished_results[i,:] = one_kfunc(x[i],y[i],t,yso_map,grid)
            
    area = np.float(np.sum(finished_results[:,0]))
    yso_sum = np.sum(finished_results[:,1])
    
    total_area = get_area()
    lmda = np.sum(yso_map)/total_area
    K = (np.pi*t**2/lmda)*yso_sum/float(area)
    L = np.sqrt(K/np.pi) - t
    if diag == True:
        return K,L,yso_sum,float(area)
    else:
        return K, L

def one_oring(x1,y1,t,w,yso_map,grid):
    """
    Returns the number of ysos and area in annulus around 
    a given yso at position x1,y1.
    """
    self_count = False

    coords = ring(x1,y1,t,w,grid)
    n_coords = np.shape(coords)[1]

    xg,yg = xy2grid(x1,y1)
    area = 0
    yso_sum = 0
    for j in range(n_coords):

        #Allow o-ring to skip one point within its own grid square.
        if coords[0,j]==xg and coords[1,j]==yg and self_count == False:
            yso_sum-=1

        area += area_array[coords[0,j],coords[1,j]]
        yso_sum += yso_map[coords[0,j],coords[1,j]]
    return np.array([area,yso_sum])
    
def Oring(x,y,t,w,yso_map=None,grid=None,noP=None,diag=False):
    """
    Calculates Oring function for points with coords x,y.
    Most likely x and y are the positions of ysos.
    noP determines the number of workers assigned to function for 
    multiprocessing. A value of None or 1 runs the function 
    without multiprocessing.
    """

    #if values not specified take global values
    if yso_map is None:
        yso_map = yso_to_grid(np.array([x,y]),grid)
        
    if grid is None:
        grid = coverage
    
    if noP > 1:
        ##Initialise pool of workers
        pool = mp.Pool(noP)

        ##Perform one_oring for each yso using workers
        results = []
        for i in range(len(x)):
            results.append(pool.apply_async(one_oring,(x[i],y[i],t,w,yso_map,grid)))

        #Close down workers to finish calculating results
        pool.close()
        pool.join()

        #collate results
        finished_results = np.empty([len(x),2])
        for i in range(len(x)):
            finished_results[i,:] = results[i].get()
    elif noP == 1 or noP == None:
        #If there are fewer than 2 workers specified, then ignore
        #multiprocessing.
        
        ##Perform one_oring for each yso
        finished_results = np.empty([len(x),2])
        for i in range(len(x)):
            finished_results[i,:] = one_oring(x[i],y[i],t,w,yso_map,grid)
            
    area = np.float(np.sum(finished_results[:,0]))
    yso_sum = np.sum(finished_results[:,1])

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
        yso_map,xy = yso_to_grid(np.array([x,y]),grid=grid,yso_return = True)
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
        yso_map, xy = yso_to_grid(np.array([x,y]),grid=grid,yso_return=True)
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
            ab_map, ab = yso_to_grid(np.array([a,b]),grid=grid,yso_return=True)
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

def get_area_array(tan=True,grid=None):
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

    if grid is None:
        grid = coverage
        
    if tan:
        #If tan projection use da_sphere = cos**3(theta)*da_plane
        if inverted:
            ra_ref = header['CRVAL2']
            dec_ref = header['CRVAL1']
        else:
            dec_ref = header['CRVAL2']
            ra_ref = header['CRVAL1']

        d2r = lambda x: x*np.pi/180.0
        
        angles = gcircle((gx,gy),(ra_ref,dec_ref))
        ref_area = wcs.utils.proj_plane_pixel_area(w_obj)
        return np.cos(d2r(angles))**3*ref_area*coverage
    else:
        #Find RA and Dec at each grid coordinate
        grx = np.arange(ra_axis+1)
        gry = np.arange(dec_axis+1)
        GX,GY = np.meshgrid(grx,gry,indexing='ij')
        GX,GY = GX.flatten(), GY.flatten()
        gry,grx = w_obj.all_pix2world(GY,GX,0)
        grx, gry = grx.reshape(ra_axis+1,dec_axis+1), gry.reshape(ra_axis+1,dec_axis+1)
        
        #Find dRA_ij = RA_{i+1},j - RA_ij and likewise for dec
        dRA = grx[1:ra_axis+1,:dec_axis]-grx[:ra_axis,:dec_axis]
        dDec = gry[:ra_axis,1:dec_axis+1]-gry[:ra_axis,:dec_axis]
        
        #Due to rotation a step in grid changes both ra and dec.
        #Calculate angle between lines of constant Ra and the
        #grid at constant i.
        dely = gry[1:ra_axis+1,:dec_axis]-gry[:ra_axis,:dec_axis]
        theta = np.arctan(dely/dRA)
        
        angle_part = np.sin(np.pi/2-gry[:ra_axis,:dec_axis]*np.pi/180.0)
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

    #CDELT may be positive or negative
    #This is to ensure that slicing of the arrays works correctly
    ras = np.sort([il,ir])
    decs = np.sort([jl,jr])
    
    if ras[0] < 0:
        ras[0] = 0
    if decs[0] < 0:
        decs[0] = 0
    if ras[1] > ra_axis:
        ras[1] = ra_axis
    if decs[1] > dec_axis:
        decs[1] = dec_axis

    return int(round(ras[0])),int(round(ras[1])),int(round(decs[0])),int(round(decs[1]))

def collate(results):
    array = np.empty([4,LOOPS,len(r)])
    for i in range(LOOPS):
        array[:,i,:] = np.array(results[i].get())
    return array

def run_csr(val,r,w,mode='sphere_binomial',noP=None,grid=None):
    """
    Perform desired number of runs get Oring and Kfunc data for CSR envelopes.
    r is array of radial distances
    w is either single value or array of bin widths
    noP is number of desired workers to parallelise Oring and K seperately.
    (not recommended)
    """
    #if values not specified take global values
    if grid is None:
        grid = coverage

    #randomise seed
    rnd.rand()

    yso,yso_map = random_ysos(val,mode,grid)
    steps = len(r)
    results = np.empty((2,steps))
    for i,t in enumerate(r):
        ##Oring
        if type(w) == np.ndarray:
            w_i = w[i]
        else:
            w_i = w
        o,oo = Oring(yso[0,:],yso[1,:],t,w_i,yso_map,grid,noP)
        results[0,i] = oo
        k,kk = kfunc(yso[0,:],yso[1,:],t,yso_map,grid,noP)
        results[1,i] = kk
    
    return results

def allenv(val,r,w,LOOPS,mode='sphere_binomial',noP=None,grid=None):
    """
    perform the LOOPS number of realisations using multiprocessing 
    for increased speed producing results which can then be processed 
    into confidence envelopes for csr.
    """
    if grid is None:
        grid = coverage
        
    if noP > 1:
        #initialise multiprocessing pool with desired number of processes
        pool = mp.Pool(noP)
        results = []
        for loop in range(LOOPS):
            results.append(pool.apply_async(run_csr,(val,r,w,mode,None,grid),callback=callbackTimer))

        final_results = np.empty((2,LOOPS,len(r)))
        for loop in range(LOOPS):
            final_results[:,loop,:] = results[loop].get()
    elif noP == 1 or noP == None:
        start = time.time()
        final_results = np.empty((2,LOOPS,len(r)))
        for loop in range(LOOPS):
            final_results[:,loop,:] = run_csr(val,r,w,mode,None,grid)

            #estimate the time left for code to run
            completed = loop/float(LOOPS)*100
            est = (time.time()-start)/(60.0*loop) * (LOOPS-loop) 
            print('%f%% complete: ~ %f more minutes' %(completed,est))
    return final_results

def get_coords(wcs_obj,grid,centre=True):
    """
    Collect coordinates of pixel centres from a pre-sliced wcs object and 
    equivalently sliced coverage map.
    """
    ra_axis,dec_axis = np.shape(grid)
    if centre:
        gx = np.arange(ra_axis)
        gy = np.arange(dec_axis)
    elif centre==False:
        gx = np.arange(0.5,ra_axis)-1
        gy = np.arange(0.5,dec_axis)-1
        
    GX,GY = np.meshgrid(gx,gy,indexing='ij')
    GX,GY = GX.flatten(), GY.flatten()
    if inverted:
        gy,gx = wcs_obj.all_pix2world(GY,GX,0)
    else:
        gx,gy = wcs_obj.all_pix2world(GX,GY,0)
        
    gx, gy = gx.reshape(ra_axis,dec_axis), gy.reshape(ra_axis,dec_axis)
    
    return gx,gy

def extract_region(bounds,wcs_obj,grid):
    """
    Extract rectangular section of coverage map designated by ra and dec coordinates.
    bounds = (2x2) array containing the boundaries of the desired section. [[Ra_0,Ra_1],[Dec_0,Dec_1]]
    grid = array to be sliced using pixel coordinates given by wcs_obj.
    """
    #coordinates of bottom-left and top-right corners of RA, Dec box.
    bl = (bounds[0,0],bounds[1,0])
    tr = (bounds[0,1],bounds[1,1])
    
    tl,br = (bl[0],tr[1]),(tr[0],bl[1])
    if inverted:
        dec_box,ra_box = w_obj.all_world2pix(np.array([bl[1],br[1],tl[1],tr[1]]),np.array([bl[0],br[0],tl[0],tr[0]]),0)
    else:
        ra_box,dec_box = w_obj.all_world2pix(np.array([bl[0],br[0],tl[0],tr[0]]),np.array([bl[1],br[1],tl[1],tr[1]]),0)
        
    dec_box = np.round(dec_box)
    ra_box = np.round(ra_box)

    #slightly widen box to ensure complete extraction of region
    ra_lims,dec_lims = [np.min(ra_box)-1,np.max(ra_box)+1],[np.min(dec_box)-1,np.max(dec_box)+1]

    axes_lims = np.shape(coverage)
    #make sure no impossible slices are attempted
    ra_lims[0] = 0 if ra_lims[0] < 0 else ra_lims[0]
    ra_lims[1] = axes_lims[0] if ra_lims[1] > axes_lims[0] else ra_lims[1]

    dec_lims[0] = 0 if dec_lims[0] < 0 else dec_lims[0]
    dec_lims[1] = axes_lims[1] if dec_lims[1] > axes_lims[1] else dec_lims[1]
    
    #slice map and wcs_object
    grid = grid[int(ra_lims[0]):int(ra_lims[1]),int(dec_lims[0]):int(dec_lims[1])]

    #according to astropy documentation the order of wcs slices "should be reversed (as for the data)
    #compared to the natural WCS order." Which I have interpreted to mean the second axes of the WCS
    #object is sliced using the first slice and vice versa.
    if inverted:
        wcs_obj = wcs_obj[int(ra_lims[0]):int(ra_lims[1]),int(dec_lims[0]):int(dec_lims[1])]
    else:
        wcs_obj = wcs_obj[int(dec_lims[0]):int(dec_lims[1]),int(ra_lims[0]):int(ra_lims[1])]

    ##Getting celestial coordinates of pixel centres for extraction
    gx,gy = get_coords(wcs_obj,grid)
    
    remove_extra_coverage = np.where( (gx < bl[0]) | (gx > tr[0]) | (gy < bl[1]) | (gy > tr[1]))
    grid[remove_extra_coverage] = False
    return wcs_obj,grid

#time length of project and estimate completion time
A = []
def callbackTimer(x):
    toc = time.time() - tic
    A.append(toc)
    completed = len(A)/float(LOOPS)*100
    est = toc/(60.0*len(A)) * (LOOPS-len(A)) #estimates the time left for code to run
    print('%f%% complete: ~ %f more minutes' %(completed,est))

def reduce_map(grid,yso_map,L,p0,tol=1,yso_return=False):
    """
    uses a square window with side length L to remove sections of grid that 
    contain too few YSOs.
    The minimum number of ysos is specified by p0 using the method by 
    Wiegand et al 2004.
    Optional yso_return argument for when ysos have been excluded due to 
    coverage map being reduced.
    """
    def get_k(lmda,w,p0):
        """
        Find kmin assuming Poisson distribution of yso counts in areas.
        """
        Pk = lambda x,y: x**y*np.exp(-x)/factorial(y)
        Pk2 = lambda x,y: np.exp(y*np.log(x)-x-(0.5*np.log(2*np.pi*y)+y*np.log(y/np.exp(1))))
        k=-1
        p=0
        while p < p0:
            k+=1
            if k < 20:
                p+=Pk(lmda*w,k)
            else:
                p+=Pk2(lmda*w,k)
        return k

    shape = np.shape(grid)
    #estimate first-order intensity
    lmda = np.sum(yso_map)/get_area()
    half_width = int(round(L/2.0))

    #check preliminary k:
    w = float(L**2)/np.sum(grid)*get_area()
    k_min=get_k(lmda,w,p0)

    for i in range(shape[0]):
        i0 = i-half_width
        i1 = i+half_width
        if i0 < 0:
            i0 = 0
        if i1 >= shape[0]:
            i1 = shape[0]
            
        for j in range(shape[1]):
            ##skip point if already excluded from coverage map
            if grid[i,j] == 0:
                continue
            
            j0 = j-half_width
            j1 = j+half_width       
            if j0 < 0:
                j0 = 0
            if j1 >= shape[1]:
                j1 = shape[1]
            w = np.sum(area_array[i0:i1,j0:j1])
            k=get_k(lmda,w,p0)
            while k < tol:
                i0 -= 1
                i1 += 1
                j0 -= 1
                j1 += 1
                if i0 <0:
                    i0 = 0
                if i1 > shape[0]:
                    i1 = shape[0]
                if j0 < 0:
                    j0 = 0
                if j1 > shape[1]:
                    j1 = shape[1]
                
                w = np.sum(area_array[i0:i1,j0:j1])
                k=get_k(lmda,w,p0)
                
            if np.sum(yso_map[i0:i1,j0:j1]) < k:
                grid[i,j] = 0
    if yso_return == False:
        return grid
    else:
        return grid, yso_map*grid
        
def clean_map(grid,L,p0):
    """
    uses a square window with side length L to remove sections of grid that 
    contain too few YSOs.
    The minimum number of ysos is specified by p0 using the method by 
    Wiegand et al 2004.
    Optional yso_return argument for when ysos have been excluded due to 
    coverage map being reduced.
    """
    def get_k(lmda,w,p0,k0=0):
        """
        Find kmin assuming Poisson distribution of yso counts in areas.
        """
        Pk = lambda x,y: x**y*np.exp(-x)/factorial(y)
        Pk2 = lambda x,y: np.exp(y*np.log(x)-x-(0.5*np.log(2*np.pi*y)+y*np.log(y/np.exp(1))))
        k= k0-1
        p=0
        while p < p0:
            k+=1
            if k < 20:
                p+=Pk(lmda*w,k)
            else:
                p+=Pk2(lmda*w,k)
            
        return k

    end_grid = np.copy(grid)
    shape = np.shape(grid)
    
    #estimate first-order intensity
    lmda = float(np.sum(grid))/np.size(grid)
    half_width = int(round(L/2.0))

    #check preliminary k:
    w0 = float(L**2)
    k_min=get_k(lmda,w0,p0)
    k=k_min

    for i in range(shape[0]):
        i0 = i-half_width
        i1 = i+half_width
        if i0 < 0:
            i0 = 0
        if i1 >= shape[0]:
            i1 = shape[0]
            
        for j in range(shape[1]):
            ##skip point if already excluded from coverage map
            if grid[i,j] == 0:
                continue
            
            j0 = j-half_width
            j1 = j+half_width       
            if j0 < 0:
                j0 = 0
            if j1 >= shape[1]:
                j1 = shape[1]

            #check if w has changed, if it has recalculate w and k
            w = np.size(grid[i0:i1,j0:j1])
            if w != w0:
                w0=w
                k = get_k(lmda,w,p0,k0=k_min/4.0)
            
            if np.sum(grid[i0:i1,j0:j1]) < k:
                end_grid[i,j] = 0
                
    return end_grid


"""
Extract the relevant data from the fits file. 
Array size.
Make a wcs object.
Array of grid square coordinates.
Array of grid centre coordinates.
Coverage map.
"""

#fits_path = '/Users/bretter/Documents/StarFormation/SFR_data'
#fits_path = '../SFR_data'
fits_path = '.'
fits_name = 'PER_IRAC1234M1_cov.fits'
coverage,header = fits.getdata(os.path.join(fits_path,fits_name), header=True)
w_obj = wcs.WCS(header)

##Find which axis is RA and which is Dec.
##inverted kayword defines if axes are RA and Dec or
##Dec and RA.
if 'DEC' in header['CTYPE1']:
    inverted = True
    dec_axis = header['NAXIS1']
    ra_axis = header['NAXIS2']
    dec_ref = header['CRVAL1']
    ra_ref = header['CRVAL2']
else:
    inverted = False
    coverage = coverage.T
    ra_axis = header['NAXIS1']
    dec_axis = header['NAXIS2']

##Extracting sections of map
#desired sector of sky bottom-left and top-right.

bounds = np.array([[52.0,53.0],[31,32]])

#number of processes
noProcess = 1
w_obj,coverage = extract_region(bounds,w_obj,coverage)

#Remove non-binary values from coverage map
cov = np.zeros(np.shape(coverage))
cov += coverage == 1

coverage = cov.astype(bool)
#cov = None

ra_axis,dec_axis = np.shape(coverage)

##Getting celestial coordinates of pixel centres
gx,gy = get_coords(w_obj,coverage)

##Getting pixel scales
area_array = get_area_array()
total_area = np.sum(area_array)

##Getting ysos
dfile = 'dunham15_tab_ysos_coords.txt'
data = np.loadtxt(dfile,skiprows=1,usecols=(2,3,8))

#Class limits with corrected alpha 
#Class 0/I  alpha >= 0.3
#Flat-spectrum -0.3 <= alpha < 0.3
#Class II -1.6 <= alpha < - 0.3
#Class III alpha < -1.6

#combine class masks into one 3d array
class01_mask = data[:,2] >= 0.3
flat_mask = (-0.3 <= data[:,2]) & (data[:,2] < 0.3)
class2_mask = (-1.6 <= data[:,2]) & (data[:,2] < -0.3)
class3_mask = data[:,2] < -1.6

alpha_mask = np.array([class01_mask,flat_mask,class2_mask,class3_mask])
pos_data = data[:,:2]
pos_mask = (pos_data[:,0] > bounds[0,0]) & (pos_data[:,0] < bounds[0,1]) & (pos_data[:,1] > bounds[1,0]) & (pos_data[:,1] < bounds[1,1])

region = 'ngc1333'
fpath = '{:s}/'.format(region)
class_list = ['classI0','flat','classII','classIII']
#loop over each yso class

tic = time.time()
for a in range(4):
    #reset coverage map for each yso class
    coverage = cov.astype(bool)
    area_array = get_area_array()
    
    total_mask = pos_mask & alpha_mask[a]
    yso = data[total_mask,:2]
    yso = yso.T
    yso_map = yso_to_grid(yso)
    
    # #save image of each coverage map
    # plt.figure()
    # plt.pcolormesh(gx,gy,coverage)
    # plt.plot(yso[0,:],yso[1,:],'*')
    # plt.xlabel('RA')
    # plt.ylabel('Dec')
    # plt.title('Non-reduced coverage map for {:s} {:s} YSOs'.format(region,class_list[a]))
    # plt.axis('equal')
    # plt.savefig('{:s}{:s}_{:s}_map.png'.format(fpath,region,class_list[a]))

    # #uncomment to perform map reduction
    # #First pass
    # p0=0.01
    # lmda = np.sum(yso_map)/get_area()
    # wmin = -np.log(p0)/lmda
    # L_sqrd = wmin/total_area*np.size(coverage)
    # L = np.sqrt(L_sqrd)
    # coverage,yso_map = reduce_map(coverage,yso_map,int(1.5*L),p0,0,True)
    # area_array = get_area_array()

    # #Second pass
    # p0=0.01
    # lmda = np.sum(yso_map)/get_area()
    # wmin = -np.log(p0)/lmda
    # L_sqrd = wmin/total_area*np.size(coverage)
    # L = np.sqrt(L_sqrd)
    # coverage,yso_map = reduce_map(coverage,yso_map,L,p0,1,True)
    # area_array = get_area_array()

    # #save image of each reduced coverage map
    # plt.figure()
    # plt.pcolormesh(gx,gy,coverage)
    # plt.plot(yso[0,:],yso[1,:],'*')
    # plt.xlabel('RA')
    # plt.ylabel('Dec')
    # plt.title('Reduced coverage map for {:s} {:s} YSOs'.format(region,class_list[a]))
    # plt.axis('equal')
    # plt.savefig('{:s}{:s}_{:s}_reduced_map.png'.format(fpath,region,class_list[a]))
    
    #Get stats
    steps = 20
    r = np.linspace(0.01,0.5,steps)
    w = r*0.6

    results = np.empty((2,steps))
    for i,v in enumerate(r):
        w_i = w[i]
        o,oo = Oring(yso[0,:],yso[1,:],v,w_i,yso_map,coverage,noP=noProcess)
        k,kk = kfunc(yso[0,:],yso[1,:],v,yso_map,coverage,noP=noProcess)
        results[0,i] = oo
        results[1,i] = kk

        #estimate time remaining
        toc = time.time()
        cur = 1+i+a*steps
        completed = cur/float(4*steps)*100
        est = (toc-tic)/(60.0*cur) * (4*steps-cur) #estimates the time left for code to run
        print('%f%% complete: ~ %f more minutes' %(completed,est))

    np.save('{:s}{:s}_{:s}_stats'.format(fpath,region,class_list[a]),results)
    #np.save('{:s}{:s}_{:s}_map'.format(fpath,region,class_list[a]),coverage)
    
#coverage = clean_map(coverage,60,p0)
#area_array = get_area_array()
#plt.figure()
#plt.pcolormesh(gx,gy,coverage)
#plt.plot(yso[0,:],yso[1,:],'*')
#plt.axis('equal')

#p0=0.01
#lmda = np.sum(yso_map)/get_area()
#wmin = -np.log(p0)/lmda
#L_sqrd = wmin/get_area()*np.size(coverage)
#coverage,yso_map = reduce_map(coverage,yso_map,np.sqrt(L_sqrd),p0,True)
#area_array = get_area_array()

#plt.figure()
#plt.pcolormesh(gx,gy,coverage)
#plt.plot(yso[0,:],yso[1,:],'*')
#plt.axis('equal')

#plt.show()

import numpy as np
import numpy.random as rnd
from math import factorial as fct
import matplotlib.pyplot as plt

def cellfind(x,y,gx,gy): #takes in x,y coords of point and array of x and y coords for grid to return coordinate of parent cell
    #xcoord = sum((gx - x) <= 0) - 1
    #ycoord = sum((gy - y) <= 0) - 1
    xcoord = np.argmin(np.abs(gx-x))
    ycoord = np.argmin(np.abs(gy-y))
    return (xcoord,ycoord)

def populate(x,y,gx,gy,grid): #takes in x,y coords of point and array of x and y coords for grid to populate grid with points in cells
    for p in range(len(x)):
        grid[cellfind(x[p],y[p],gx,gy)] += 1
    return grid

def poisson(k,lmda,area):
    return ((lmda*area)**k)/float(fct(k)) * np.exp(-lmda*area)

def kmin(p0,lmda,area):
    k = 0
    p = 0

    while True:
        q = poisson(k,lmda,area)
        if p+q > p0:
            return k
        else:
            p = p+q
            k += 1

def circle(x0,y0,R,gx,gy,grid):
    for x in range(np.shape(grid)[0]):
        for y in range(np.shape(grid)[1]):
            grid[x,y] = np.sqrt((gx[x]-x0)**2 + (gy[y]-y0)**2) <= R

def ring(x0,y0,R,h,gx,gy,grid):
    for x in range(np.shape(grid)[0]):
        for y in range(np.shape(grid)[1]):
            grid[x,y] = R-h < np.sqrt((gx[x]-x0)**2 + (gy[y]-y0)**2) <= R+h
    
#Grid-based system of identifying homogeneous regions of point patterns

#need to program in the basic geometry of the plot.
XMIN = 0
XMAX = 50
YMIN = 0 
YMAX = 50

p0 = 0.01

sideL = float(100)
grid = np.zeros((sideL,sideL)) #grid of where points are located within the overall grid

plotArea = np.ones((sideL,sideL)) #grid of where can be considered as an area
cellarea = ((XMAX-XMIN)*(YMAX-YMIN))/float((sideL**2))

gx = np.linspace(XMIN,XMAX,sideL,endpoint=False) + (XMAX-XMIN)/(2.0*sideL)
gy = np.linspace(YMIN,YMAX,sideL,endpoint=False) + (YMAX-YMIN)/(2.0*sideL)



##
N = 1000
RADIUS = 30
## Generating N randomly distributed points with an r^-2 density power law
r = rnd.rand(N/2)*RADIUS
angle = rnd.rand(N/2)*2*np.pi
x = r*np.cos(angle) + XMAX/2
y = r*np.sin(angle) + YMAX/2

x = np.hstack((x,rnd.rand(N-len(x))*XMAX))
y = np.hstack((y,rnd.rand(N-len(y))*YMAX))



## Window
RADIUS = 25
#number of expected points in a window
sidex = np.ceil(RADIUS/((XMAX-XMIN)/sideL)-0.5)
sidey = np.ceil(RADIUS/((YMAX-YMIN)/sideL)-0.5)
expecWin = float((2*sidex)*(2*sidey))


area = np.sum(plotArea)*cellarea
lmda = np.sum(grid)/area

#populate(x,y,gx,gy,grid)
ring(25,25,10,2.5,gx,gy,grid)
Kmin = kmin(p0,lmda,expecWin*cellarea)

windowCount = np.zeros((sideL,sideL)) #grid of where points are located within the overall grid
#looks over with a rough window to find counts of points within a given radius of grid spaces
 
for i in range(len(gx)):
    xwin0 = np.argmin(np.abs(gx[i]-RADIUS-gx))
    xwin1 = np.argmin(np.abs(gx[i]+RADIUS-gx))
    for j in range(len(gy)):
        ywin0 = np.argmin(np.abs(gy[j]-RADIUS-gy))
        ywin1 = np.argmin(np.abs(gy[j]+RADIUS-gy))
        windowCount[i,j] = np.sum(grid[xwin0:xwin1+1,ywin0:ywin1+1])*(expecWin/np.sum(plotArea[xwin0:xwin1+1,ywin0:ywin1+1]))

plt.figure()
plt.imshow(windowCount)
plt.figure()
plt.imshow(grid)

plotArea = plotArea * (windowCount >= Kmin)
grid = grid * plotArea #remove isolated points

area = np.sum(plotArea)*cellarea
lmda = np.sum(grid)/area

Kmin = kmin(p0,lmda,expecWin*cellarea)

## Window
RADIUS = 9

#number of expected points in a window
sidex = np.ceil(RADIUS/((XMAX-XMIN)/sideL)-0.5)
sidey = np.ceil(RADIUS/((YMAX-YMIN)/sideL)-0.5)
expecWin = float((2*sidex)*(2*sidey))


windowCount = np.zeros((sideL,sideL)) #reset window count

for i in range(len(gx)):
    xwin0 = np.argmin(np.abs(gx[i]-RADIUS-gx))
    xwin1 = np.argmin(np.abs(gx[i]+RADIUS-gx))
    for j in range(len(gy)):
        ywin0 = np.argmin(np.abs(gy[j]-RADIUS-gy))
        ywin1 = np.argmin(np.abs(gy[j]+RADIUS-gy))
        windowCount[i,j] = np.sum(grid[xwin0:xwin1+1,ywin0:ywin1+1])*(expecWin/np.sum(plotArea[xwin0:xwin1+1,ywin0:ywin1+1]))


Kmin = kmin(p0,lmda,expecWin*cellarea)
plotArea = plotArea * (windowCount >= Kmin)
grid = grid * plotArea



plt.figure()
plt.imshow(windowCount*plotArea)
plt.figure()
plt.imshow(grid)
plt.show()

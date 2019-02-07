##Code to analyse grid centre results
import numpy as np
import numpy.random as rnd
import matplotlib.pyplot as plt

def LtoK(L,r):
    return None
    
##Read in data

results = np.load('grid_centre_all_stats_50_1.npy')
#(# of resolutions, {O-ring, K} , {grid, analytical}, {Stat, Count, area}, radial steps)
shape = np.shape(results)

fpath = '/Users/bretter/Documents/StarFormation/RandomDistribution/spatialStats/GridBased/figs/50'
for i in range(shape[0]):
    O = results[i,0,:,:,:]
    K = results[i,1,:,:,:]

    test = ['stat','count','area']
    res = [ 27,  81, 243, 729]
    for k in range(shape[3]):
        plt.figure()
        plt.plot(r,O[0,k,:],'r',lw=3)
        plt.plot(r,O[1,k,:],'b')
        plt.title('O-ring. Grid based (r), analytical (blue)')
        plt.savefig('{:s}/{:s}/O_{:d}'.format(fpath,test[k],res[i]))
    
        plt.figure()
        plt.plot(r,O[1,k,:]-O[0,k,:])
        plt.title('O-ring. Analytical-Grid based')
        plt.savefig('{:s}/{:s}/Odiff_{:d}'.format(fpath,test[k],res[i]))
    
        plt.figure()
        plt.plot(r,K[0,k,:],'r',lw=3)
        plt.plot(r,K[1,k,:],'b')
        plt.title('L. Grid based (r), analytical (blue)')
        plt.savefig('{:s}/{:s}/L_{:d}'.format(fpath,test[k],res[i]))
    
        plt.figure()
        plt.plot(r,K[1,k,:]-K[0,k,:])
        plt.title('L. Analytical-Grid based')
        plt.savefig('{:s}/{:s}/Ldiff_{:d}'.format(fpath,test[k],res[i]))
    
        plt.close('all')
        

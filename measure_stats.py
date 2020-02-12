#Code to test if anyshape_grid can count as a module

from anyshape_grid import *

def initialise_obj(astro_obj,bounds):
    
"""
Steps:
Make astro_box object using fits_file
"""

region = 'serpens_south'
bounds = np.array([[277.2, 277.7],[-2.25,-1.75]])
fits_name = 'SERAQU_IRAC1234M1_cov_sm.fits'
ext_name='HGBS_aquilaM2_hires_column_density_map.fits'
distance_to = 484
steps = 20
r = np.linspace(0.01,0.25,steps)

#number of processes
noProcess = 4

fits_path = '/Users/bretter/Documents/StarFormation/SFR_data'
#fits_path = '../../SFR_data'
#fits_path = '.'
fits_file = os.path.join(fits_path,fits_name)
w_obj = astro_box(fits_file,bounds)

#Remove non-binary values from coverage map and store as boolean.
w_obj.booleanise()

##Getting pixel scales
w_obj.get_area_array(dist=distance_to)

class_list = ['class0I','flat','classII','classIII','all']
#loop over each yso class

tic = time.time()
for a,cl in enumerate(class_list):
    if a > 0:
        break
    #extract ysos
    yso = get_yso_locs(bounds,cl,dpath='..')
    yso_map = yso_to_grid(yso,w_obj,yso_return=False)
    
    #Get stats
    w = r*0.6
    results = np.empty((2,steps))
    for i,v in enumerate(r):
        w_i = w[i]
        o,oo = Oring(yso[0,:],yso[1,:],v,w_i,w_obj,yso_map,noP=noProcess)
        k,kk = kfunc(yso[0,:],yso[1,:],v,w_obj,yso_map,noP=noProcess)
        
        results[0,i] = oo
        results[1,i] = kk

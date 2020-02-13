#code to make envelopes

from anyshape_grid import *
import copy 
"""
Extract the relevant data from the fits file. 
Array size.
Make a wcs object.
Array of grid square coordinates.
Array of grid centre coordinates.
Coverage map.
"""

region = 'serpens_south'
bounds = np.array([[277.2, 277.7],[-2.25,-1.75]])
fits_name = 'SERAQU_IRAC1234M1_cov_sm.fits'
ext_name='HGBS_aquilaM2_hires_column_density_map.fits'
distance_to = 484
steps = 10
r = np.linspace(0.01,0.25,steps)

fits_path = '/Users/bretter/Documents/StarFormation/SFR_data'
#fits_path = '../../SFR_data'
#fits_path = '.'
fits_file = os.path.join(fits_path,fits_name)
w_obj = astro_box(fits_file,bounds)

#Remove non-binary values from coverage map and store as boolean.
w_obj.booleanise()

##Getting pixel scales
w_obj.get_area_array(dist=distance_to)

##Read in extinction map
extinction_fits = os.path.join(fits_path,ext_name)
ext_obj = astro_box(extinction_fits,get_coords=False)

ext_obj = resample_fits(ext_obj,w_obj)
ext_obj.extract_region(bounds)

##Initialise envelope
noProcesses = 4

#loop over each yso class
tic = time.time()
fpath = '{:s}/'.format(region)
class_list = ['class0I','flat','classII','classIII','all']
for a,cl in enumerate(class_list):

    if a > 0:
        break
    #extract ysos
    yso = get_yso_locs(bounds,cl,dpath='..')
    
    steps = 20
    w = 0.6*r
    val = int(np.shape(yso)[1])
    prob_map = ext_obj.grid**2.05
    
    #produce bin edges
    binss = [np.log10(1e20),np.log10(np.max(prob_map)),11]
    binss = np.logspace(*binss)

    rnd.seed(1)

    
    LOOPS = 20
    results = np.empty((LOOPS,2,val))
    for i in range(LOOPS):
        yso,yso_map = random_ysos(val,w_obj,'nhpp',prob_map)
        results[i,:,:] = yso
        
    #results = allenv(val,r,w,LOOPS,w_obj,ext_obj.grid,mode='nhpp',noP=noProcesses,timer=False)
np.save('../new_random_yso',results)

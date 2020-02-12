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
steps = 20
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

ext_obj = resample_fits2(ext_obj,w_obj)

# ##Initialise envelope
# noProcesses = 30

# #loop over each yso class
# tic = time.time()
# fpath = '{:s}/'.format(region)
# class_list = ['classI0','flat','classII','classIII','all']
# for a,cl in enumerate(class_list):

#     #extract ysos
#     yso, yso_map = get_yso_locs(bounds,cl,dpath=None)
    
#     steps = 20
#     w = 0.6*r
#     val = int(np.shape(yso)[1])
    
#     #produce bin edges
#     binss = [np.log10(1e20),np.log10(np.max(lmda_map)),11]
#     binss = np.logspace(*binss)
#     map2 = extinction_prob(yso,binss,area_array,lmda_map,w_obj)
    
#     LOOPS = 99
#     results = allenv(val,r,w,LOOPS,mode='nhpp',noP=noProcesses,grid=coverage,density=map2,timer=False)

#     #estimate time remaining
#     toc = time.time()
#     cur = 1+a
#     completed = cur/float(len(class_list))*100
#     est = ((toc-tic)/60.0)*(float(len(class_list))/cur -1) #estimates the time left for code to run
#     print('{:s} complete: ~ {:.3f} more minutes'.format(cl,est))

#     np.save('{:s}{:s}_{:s}_herschel_probmap_extinction_nhpp.npy'.format(fpath,region,cl),results)



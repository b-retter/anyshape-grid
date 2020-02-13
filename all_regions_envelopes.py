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

noProcesses = 4
steps = 20
LOOPS = 2
fits_path = '/Users/bretter/Documents/StarFormation/SFR_data'
#fits_path = '../../SFR_data'
#fits_path = '.'

allregions = ['serpens_south','serpens_core','ophiuchus','ngc1333','ic348']
alldistances = [484,484,144,293,321]
allbounds = np.stack([np.array([[277.2, 277.7],[-2.25,-1.75]]),
                          np.array([[277.4, 277.6],[1.18,2.28]]),
                          np.array([[246.0, 248.5],[-25.2,-23.8]]),
                          np.array([[52.0, 52.8],[31,31.8]]),
                          np.array([[55.8, 56.4],[31.9,32.5]])],axis=0)

allfits_names = ['SERAQU_IRAC1234M1_cov_sm.fits',
                     'SER_IRAC1234M1_cov.fits',
                     'OPH_ALL_IRAC1234M1_cov_sm.fits',
                     'PER_IRAC1234M1_cov_sm.fits',
                     'PER_IRAC1234M1_cov_sm.fits']

allherschel = ['HGBS_aquilaM2_hires_column_density_map.fits',
                   'HGBS_serpens_hires_column_density_map.fits',
                   'HGBS_oph_l1688_hires_column_density_map.fits',
                   'HGBS_perseus_hires_column_density_map.fits',
                   'HGBS_perseus_hires_column_density_map.fits']

all_r = np.stack([np.linspace(0.01,0.25,steps),np.linspace(0.001,0.05,steps),
                      np.linspace(0.01,0.7,steps),np.linspace(0.01,0.5,steps),
                      np.linspace(0.01,0.3,steps)],axis=0)

for j,region in enumerate(allregions):
    print('{:s}'.format(region))
    fits_name = allfits_names[j]
    bounds = allbounds[j]
    fits_file = os.path.join(fits_path,fits_name)
    w_obj = astro_box(fits_file,bounds)

    #Remove non-binary values from coverage map and store as boolean.
    w_obj.booleanise()

    ##Getting pixel scales
    w_obj.get_area_array(dist=alldistances[j])

    ##Read in extinction map
    herschel_name = allherschel[j]
    extinction_fits = os.path.join(fits_path,herschel_name)
    ext_obj = astro_box(extinction_fits,get_coords=False)

    ext_obj = resample_fits(ext_obj,w_obj)
    ext_obj.extract_region(bounds)
    
    r = all_r[j]
    distance_to = alldistances[j]
    
    ##Initialise envelope
    #loop over each yso class
    tic = time.time()
    fpath = '{:s}/'.format(region)
    class_list = ['class0I','flat','classII','classIII','all']
    for a,cl in enumerate(class_list):
        #extract ysos
        yso = get_yso_locs(bounds,cl,dpath='..')
        w = 0.6*r
        val = int(np.shape(yso)[0])

        prob_map = ext_obj.grid**2.05
        
        results = allenv(val,r,w,LOOPS,w_obj,density=prob_map,mode='nhpp',noP=noProcesses,timer=False)

        #estimate time remaining
        toc = time.time()
        cur = 1+a
        completed = cur/float(len(class_list))*100
        est = ((toc-tic)/60.0)*(float(len(class_list))/cur -1) #estimates the time left for code to run
        print('{:s} complete: ~ {:.3f} more minutes'.format(class_list[a],est))

        #np.save('{:s}{:s}_{:s}_global_area_probmap_nhpp.npy'.format(fpath,region,cl),results)

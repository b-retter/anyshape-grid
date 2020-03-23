#Code to test if anyshape_grid can count as a module

from anyshape_grid import *
    
"""
Steps:
Make astro_box object using fits_file
"""
noProcess = 10
steps = 20

#directory for saving stats
#fpath = '/Users/bretter/Documents/StarFormation/RandomDistribution/Grid Based Spatial Stats/Data/SFR_stats'
fpath = '.'

#directory containing fits files
#fits_path = '/Users/bretter/Documents/StarFormation/SFR_data'
#fits_path = '../../SFR_data'
fits_path = 'HERSCHEL_SPITZER'

allregions = ['serpens_south','serpens_core','ophiuchus','ngc1333','ic348']
alldistances = [484,484,144,293,321]
allbounds = np.stack([np.array([[277.2, 277.7],[-2.25,-1.75]]),
                          np.array([[277.4, 277.6],[1.18,2.28]]),
                          np.array([[246.0, 248.5],[-25.2,-23.8]]),
                          np.array([[52.0, 52.8],[31,31.8]]),
                          np.array([[55.8, 56.4],[31.9,32.5]])],axis=0)

allfits_names = ['SERAQU_IRAC1234M1_cov.fits',
                     'SER_IRAC1234M1_cov.fits',
                     'OPH_ALL_IRAC1234M1_cov.fits',
                     'PER_IRAC1234M1_cov.fits',
                     'PER_IRAC1234M1_cov.fits']

all_r = np.stack([np.linspace(0.25/steps,0.25,steps),
                      np.linspace(0.05/steps,0.05,steps),
                      np.linspace(0.7/steps,0.7,steps),
                      np.linspace(0.5/steps,0.5,steps),
                      np.linspace(0.3/steps,0.3,steps)],axis=0)

#Prep for finding r steps in terms of parsecs
#define maximum distances in degrees and size of parsec steps.
max_r_values = np.array([0.25,0.05,0.7,0.5,0.3])
pc_step = 0.03

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
    
    distance_to = alldistances[j]

    #Recalculate radial distances by converting to parsecs, finding the steps in parsecs
    #and converting them back into angles.
    conversion_factor = distance_to*(np.pi/180)
    max_r_pc = max_r_values[j]*conversion_factor
    r_pc = np.arange(pc_step,max_r_pc,pc_step)
    r = r_pc/conversion_factor
    
    ##Initialise envelope
    #loop over each yso class
    tic = time.time()
    class_list = ['class0I','flat','classII','classIII','all']
    for a,cl in enumerate(class_list):
        if a > 0:
            break
        #extract ysos
        yso = get_yso_locs(bounds,cl,dpath='.')
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

        statsfile='{:s}_{:s}_0.03pc_statistics'.format(region,cl)
        statsdir = os.path.join(fpath,region,statsfile)
        np.save(statsdir,results)
        
        #estimate time remaining
        toc = time.time()
        cur = 1+a
        completed = cur/float(len(class_list))*100
        est = ((toc-tic)/60.0)*(float(len(class_list))/cur -1) #estimates the time left for code to run
        print('{:s} complete: ~ {:.3f} more minutes'.format(class_list[a],est))

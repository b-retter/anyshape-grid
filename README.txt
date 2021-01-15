Measuring O-ring and K for YSOs:
Make astro_box object using fits_file
	Optional: provide regions boundaries in array of form np.array([[Ra_0,Ra_1],[Dec_0,Dec_1]]) to use just a subset of region.

Tell astro_box to make an area_array using .get_area_array(distance_to_region)

Read in YSO locations.
	anyshape_grid has a function get_yso_locs which can extract YSOs from Dunham 2015 YSO catalogue.
	Otherwise just a 2XN numpy array of coords with [[RA_0,...,RA_n],[Dec_0,...,Dec_n]].

Make yso_map
	yso_map is an array the same shape as the FITS data.
	elements of a yso_map contain the number of YSOs within the equivalent cell in the FITS array. 

Define radial steps and annulus widths (for O-ring only).
	Note: code expects these in degrees.	

Call Oring() and kfunc().

Making envelopes:
Follow steps above except for calling Oring() and kfunc().

allenv() will make envelopes for O-ring and K. 
see allenv() docstring.

allenv() can accept a probability map for non-homogeneous models.

For Herschel data you can make an astro_box from Herschel FITS files and then astro_box.grid for the column density data.

